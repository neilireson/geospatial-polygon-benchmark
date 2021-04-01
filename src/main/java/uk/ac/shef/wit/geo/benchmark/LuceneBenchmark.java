package uk.ac.shef.wit.geo.benchmark;

import me.tongfei.progressbar.ProgressBar;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.LatLonPoint;
import org.apache.lucene.document.LatLonShape;
import org.apache.lucene.document.ShapeField;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.geo.LatLonGeometry;
import org.apache.lucene.geo.Point;
import org.apache.lucene.geo.Polygon;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreMode;
import org.apache.lucene.search.SimpleCollector;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.MMapDirectory;
import org.apache.lucene.store.NIOFSDirectory;
import org.apache.lucene.store.RAMDirectory;
import org.geotools.data.DataStore;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LinearRing;
import org.opengis.feature.simple.SimpleFeature;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Warmup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

import static java.lang.Math.abs;
import static org.apache.lucene.search.ScoreMode.COMPLETE_NO_SCORES;

@State(Scope.Thread)
public class LuceneBenchmark
        extends AbstractBenchmark {

    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    private final String fieldName = "location";
    private IndexSearcher indexSearcher = null;

    private enum LuceneType {
        ram, mmap, directory;

        public Path getIndexPath(File outputDirectory, String indexSource) {
            switch (this) {
                case ram:
                    return null;
                case mmap:
                case directory:
                    return Paths.get(outputDirectory.getAbsolutePath(), "lucene-polygons-index-" + indexSource);
                default:
                    throw new UnsupportedOperationException("LuceneType: " + this);
            }
        }

        public Directory getDirectory(File outputDirectory, String indexSource)
                throws IOException {
            switch (this) {
                case ram:
                    return new RAMDirectory();
                case mmap:
                    //noinspection ConstantConditions
                    return new MMapDirectory(getIndexPath(outputDirectory, indexSource));
                case directory:
                    //noinspection ConstantConditions
                    return NIOFSDirectory.open(getIndexPath(outputDirectory, indexSource));
                default:
                    throw new UnsupportedOperationException("LuceneType: " + this);
            }
        }
    }

    private final LuceneType luceneType = LuceneType.directory;

    @Setup
    public void setup() throws IOException {

        logger.info("Setting up Lucene {}", luceneType);

        Map.Entry<DataStore, SimpleFeatureCollection> dataStoreCollection = getIndexPolygons();
        DataStore dataStore = dataStoreCollection.getKey();
        SimpleFeatureCollection polygons = dataStoreCollection.getValue();

        try {
            boolean createIndex = false;
            Path indexPath = luceneType.getIndexPath(getOutputDirectory(), configName);
            if (indexPath != null && Files.exists(indexPath)) {
                logger.info("Reading Lucene index...");
                Directory directory = luceneType.getDirectory(getOutputDirectory(), configName);
                int count;
                try (final IndexReader indexReader = DirectoryReader.open(directory)) {
                    count = indexReader.numDocs();
                }

                if (count == 0) {
                    logger.error("Index is empty");
                    createIndex = true;
                } else if (polygons.size() != count) {
                    logger.error("Index contains incorrect number of documents. Expected {}, found {}",
                            polygons.size(), count);
                    if (abs(polygons.size() - count) <= config.getMissingDataThreshold()) {
                        logger.info("Missing number of documents is within acceptable limit, {} <= {}",
                                abs(polygons.size() - count), config.getMissingDataThreshold());
                    } else {
                        createIndex = true;
                    }
                } else {
                    logger.info("Index {} contains {} documents", indexPath, count);
                }
            }

            IndexWriterConfig indexWriterConfig = new IndexWriterConfig();
            if (createIndex || indexPath == null || !Files.exists(indexPath)) {
                if (createIndex) {
                    // this will delete all documents currently indexed
                    indexWriterConfig.setOpenMode(IndexWriterConfig.OpenMode.CREATE);
                }
                try (Directory directory = luceneType.getDirectory(getOutputDirectory(), configName);
                     IndexWriter indexWriter = new IndexWriter(directory, indexWriterConfig)) {

                    logger.info("Indexing features...");
                    try (ProgressBar progressBar = new ProgressBar("Features:", polygons.size());
                         SimpleFeatureIterator it = polygons.features()) {
                        Function<SimpleFeature, SimpleFeature> simplificationFunction =
                                getSimplificationFunction(polygons);
                        int id = 0;
                        while (it.hasNext()) {
                            SimpleFeature feature = it.next();
                            if (simplificationFunction != null) {
                                feature = simplificationFunction.apply(feature);
                            }
                            Geometry featureGeometry = (Geometry) feature.getDefaultGeometry();
                            try {
                                Document doc = new Document();
                                doc.add(new StoredField("id", ++id));
                                addGeometry(featureGeometry, doc);
                                indexWriter.addDocument(doc);
                            } catch (Exception e) {
                                logger.error("Failed to add polygon:{}", featureGeometry, e);
                            } finally {
                                progressBar.step();
                            }
                        }
                        indexWriter.commit();
                    }
                } catch (Exception e) {
                    logger.error("Failed to create Lucene index", e);
                    if (indexPath != null) {
                        try {
                            Files.walk(indexPath)
                                    .sorted(Comparator.reverseOrder())
                                    .forEach(path -> {
                                        try {
                                            Files.delete(path);
                                        } catch (IOException e2) {
                                            throw new RuntimeException(e2);
                                        }
                                    });
                            logger.info("Deleted Lucene index: {}", indexPath);
                        } catch (Exception e2) {
                            logger.error("Failed to delete Lucene index", e2);
                        }
                    }
                    throw (e instanceof IOException) ? (IOException) e : new IOException(e);
                }
            }
        } finally {
            if (dataStore != null) {
                dataStore.dispose();
            }
        }

        logger.info("Reading Lucene index...");
        Directory directory = luceneType.getDirectory(getOutputDirectory(), configName);
        final IndexReader indexReader = DirectoryReader.open(directory);
        indexSearcher = new IndexSearcher(indexReader);
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 1)
    public void pointIntersectsQuery() throws IOException {

        long candidateCount = 0;
        long nearestCount = 0;
        for (double[] latlon : getQueryPoints()) {
            Point point = new Point(latlon[0], latlon[1]);
            long totalHits = query(point);
            if (totalHits > 0) {
                candidateCount += totalHits;
                nearestCount++;
            }
        }
        candidateCounts.add(nearestCount == 0 ? 0 : candidateCount / nearestCount);
        nearestCounts.add(nearestCount);
    }


    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 1)
    public void polygonIntersectsQuery() throws IOException {

        long candidateCount = 0;
        long nearestCount = 0;
        for (double[][] latlons : getQueryPolygons()) {
            Polygon polygon = getPolygon(latlons);
            long totalHits = query(polygon);
            if (totalHits > 0) {
                candidateCount += totalHits;
                nearestCount++;
            }
        }
        candidateCounts.add(nearestCount == 0 ? 0 : candidateCount / nearestCount);
        nearestCounts.add(nearestCount);
    }

    private void addGeometry(Geometry geometry, Document doc) throws java.text.ParseException {
        org.locationtech.jts.geom.Point point = geometry.getCentroid();
        double lon = point.getX();
        double lat = point.getY();
        doc.add(new StoredField("latitude", lat));
        doc.add(new StoredField("longitude", lon));
        doc.add(new LatLonPoint("latlon", lat, lon));
        String geoJSON = geometryJSON.get().toString(geometry);
        doc.add(new StoredField("geoJson", geoJSON));
        for (Polygon poly : Polygon.fromGeoJSON(geoJSON)) {
            for (Field f : LatLonShape.createIndexableFields(fieldName, poly)) {
                doc.add(f);
            }
        }
    }

    private Polygon getPolygon(double[][] latlons) {
        double[] lats = new double[latlons.length];
        double[] lons = new double[latlons.length];
        for (int i = 0; i < latlons.length; i++) {
            lats[i] = latlons[i][0];
            lons[i] = latlons[i][1];
        }
        return new Polygon(lats, lons);
    }

    // FIXME does not work for complex geometries
    private Polygon getPolygon(SimpleFeature feature) {
        return getPolygon((Geometry) feature.getDefaultGeometry());
    }

    private Polygon getPolygon(Geometry geometry) {
        if (geometry.getNumGeometries() == 0) return null;

        switch (geometry.getGeometryType()) {
            case Geometry.TYPENAME_POLYGON: {
                org.locationtech.jts.geom.Polygon jtsPolygon = (org.locationtech.jts.geom.Polygon) geometry;
                double[][] latlons = getLatlons(jtsPolygon.getExteriorRing().getCoordinates());
                int numInteriorRings = jtsPolygon.getNumInteriorRing();
                if (numInteriorRings == 0)
                    return new Polygon(latlons[0], latlons[1]);
                else {
                    Polygon[] holes = new Polygon[numInteriorRings];
                    for (int i = 0; i < numInteriorRings; i++) {
                        holes[i] = getPolygon(jtsPolygon.getInteriorRingN(i));
                    }
                    return new Polygon(latlons[0], latlons[1], holes);
                }
            }
            case Geometry.TYPENAME_MULTIPOLYGON: {
                // First polygon is outer rest are holes
                if (geometry.getNumGeometries() == 1) {
                    return getPolygon(geometry.getGeometryN(0));
                }
                double[][] latlons = getLatlons(geometry.getGeometryN(0).getCoordinates());
                Polygon[] holes = new Polygon[geometry.getNumGeometries() - 1];
                for (int i = 1; i < geometry.getNumGeometries(); i++) {
                    holes[i - 1] = getPolygon(geometry.getGeometryN(i));
                }
                return new Polygon(latlons[0], latlons[1], holes);
            }
            case Geometry.TYPENAME_LINEARRING: {
                LinearRing linearRing = (LinearRing) geometry;
                double[][] latlons = getLatlons(linearRing.getCoordinates());
                return new Polygon(latlons[0], latlons[1]);
            }
            default:
                throw new IllegalArgumentException("Cannot create polygon from geometry: " + geometry.getClass());

        }
    }

    double[][] getLatlons(Coordinate[] coordinates) {
        double[] lats = new double[coordinates.length];
        double[] lons = new double[coordinates.length];
        for (int i = 0; i < coordinates.length; i++) {
            lats[i] = coordinates[i].y;
            lons[i] = coordinates[i].x;
        }
        return new double[][]{lats, lons};
    }

    private long query(LatLonGeometry geometry)
            throws IOException {
        String id = "0";
        float distance = -1;
        Query query = LatLonShape.newGeometryQuery(fieldName, ShapeField.QueryRelation.INTERSECTS, geometry);
        TotalHitCollector collector = new TotalHitCollector();
        indexSearcher.search(query, collector);
        if (collector.getTotalHits() != 0) {
            Document doc = indexSearcher.doc(collector.getDoc(0));
            id = doc.get("id");
        }
        results.add(new AbstractMap.SimpleImmutableEntry<>(id, (double) distance));
        return collector.getTotalHits();
    }

    @TearDown
    public void teardown() {
        super.teardown();
    }

    @Override
    public void close() throws IOException {
        indexSearcher.getIndexReader().close();
    }

    public static class TotalHitCollector extends SimpleCollector {

        private int base;
        private final List<Integer> docs = new ArrayList<>();

        public int getTotalHits() {
            return docs.size();
        }

        public int getDoc(int i) {
            return docs.get(i);
        }

//        public int[] getDocs() {
//            return docs.stream().mapToInt(i -> i).toArray();
//        }

        @Override
        public void collect(int doc) {
            doc += this.base;
            docs.add(doc);
        }

        @Override
        protected void doSetNextReader(LeafReaderContext context) {
            this.base = context.docBase;
        }

        @Override
        public ScoreMode scoreMode() {
            return ScoreMode.COMPLETE_NO_SCORES;
        }
    }

    

    public static void main(String[] args) {
        try (LuceneBenchmark benchmark = new LuceneBenchmark()) {
            benchmark.setup();
            benchmark.pointIntersectsQuery();
            benchmark.teardown();
            benchmark.polygonIntersectsQuery();
            benchmark.teardown();
        } catch (IOException e) {
            logger.error("", e);
        }
    }
}

package uk.ac.shef.wit.geo.benchmark;

import me.tongfei.progressbar.ProgressBar;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
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
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreMode;
import org.apache.lucene.search.SimpleCollector;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.MMapDirectory;
import org.apache.lucene.store.NIOFSDirectory;
import org.apache.lucene.store.RAMDirectory;
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

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import static org.apache.lucene.search.ScoreMode.COMPLETE_NO_SCORES;

@State(Scope.Thread)
public class LuceneBenchmark
        extends AbstractBenchmark {

    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    private final String fieldName = "location";
    private IndexSearcher indexSearcher = null;

    private enum LuceneType {
        ram, mmap, directory;

        public Path getIndexPath(String outputDirectoryName, int numberOfIndexPolygons) {
            switch (this) {
                case ram:
                    return null;
                case mmap:
                    return Paths.get(outputDirectoryName, "lucene-polygons-mm-index" + numberOfIndexPolygons);
                case directory:
                    return Paths.get(outputDirectoryName, "lucene-polygons-index" + numberOfIndexPolygons);
                default:
                    throw new UnsupportedOperationException("LuceneType: " + this);
            }
        }

        public Directory getDirectory(String outputDirectoryName, int numberOfIndexPolygons)
                throws IOException {
            switch (this) {
                case ram:
                    return new RAMDirectory();
                case mmap:
                    //noinspection ConstantConditions
                    return new MMapDirectory(getIndexPath(outputDirectoryName, numberOfIndexPolygons));
                case directory:
                    //noinspection ConstantConditions
                    return NIOFSDirectory.open(getIndexPath(outputDirectoryName, numberOfIndexPolygons));
                default:
                    throw new UnsupportedOperationException("LuceneType: " + this);
            }
        }
    }

    private final LuceneType luceneType = LuceneType.mmap;


    @Setup
    public void setup()
            throws IOException {

        logger.info("Setting up Lucene {}", luceneType);

        createDirectory(outputDirectoryName);
        Path indexPath = luceneType.getIndexPath(outputDirectoryName, numberOfIndexPolygons);
        if (indexPath != null && Files.exists(indexPath)) {
            logger.info("Reading Lucene index: {}...", indexPath);
            final IndexReader indexReader =
                    DirectoryReader.open(luceneType.getDirectory(outputDirectoryName, numberOfIndexPolygons));
            indexSearcher = new IndexSearcher(indexReader);
            return;
        }

        try (Directory directory = luceneType.getDirectory(outputDirectoryName, numberOfIndexPolygons);
             IndexWriter indexWriter = new IndexWriter(directory, new IndexWriterConfig())) {

            int id = 0;
            List<double[][]> polygons = getIndexPolygons();
            try (ProgressBar progressBar = new ProgressBar("Features:", numberOfIndexPolygons)) {
                for (double[][] latlons : polygons) {
                    try {
                        Polygon polygon = getPolygon(latlons);
                        Document doc = new Document();
                        doc.add(new StoredField("id", ++id));
                        for (Field f : LatLonShape.createIndexableFields(fieldName, polygon)) {
                            doc.add(f);
                        }
                        indexWriter.addDocument(doc);
                    } catch (Exception e) {
                        StringBuilder sb = new StringBuilder();
                        for (double[] latlon : latlons) {
                            sb.append(" ").append(latlon[1]).append(",").append(latlon[0]);
                        }
                        logger.error("Failed to add polygon:{}", sb, e);
                        throw e;
                    } finally {
                        progressBar.step();
                    }
                }
                indexWriter.commit();
            }
            final IndexReader indexReader = DirectoryReader.open(directory);
            indexSearcher = new IndexSearcher(indexReader);
        }
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
        results.clear();
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
        results.clear();
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

    private Polygon getPolygon(double[][] latlons) {
        double[] lats = new double[latlons.length];
        double[] lons = new double[latlons.length];
        for (int i = 0; i < latlons.length; i++) {
            lats[i] = latlons[i][0];
            lons[i] = latlons[i][1];
        }
        return new Polygon(lats, lons);
    }

    private long query(LatLonGeometry geometry)
            throws IOException {
        int id = 0;
        float distance = -1;
        Query query = LatLonShape.newGeometryQuery(fieldName, ShapeField.QueryRelation.INTERSECTS, geometry);
        TotalHitCollector collector = new TotalHitCollector();
        indexSearcher.search(query, collector);
        if (collector.getTotalHits() != 0) {
            Document doc = indexSearcher.doc(collector.getDoc(0));
            id = Integer.parseInt(doc.get("id"));
        }
        results.add(new AbstractMap.SimpleImmutableEntry<>(id, (double) distance));
        return collector.getTotalHits();
    }

    @TearDown
    public void teardown()
            throws IOException {
        super.teardown();
//        indexSearcher.getIndexReader().close();
    }


    public static void main(String[] args) throws IOException {
        LuceneBenchmark benchmark = new LuceneBenchmark();
        benchmark.setup();
        benchmark.pointIntersectsQuery();
        benchmark.teardown();
        benchmark.polygonIntersectsQuery();
        benchmark.teardown();
    }

    private static class TotalHitCollector extends SimpleCollector {

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
            docs.add(doc);
        }

        @Override
        public ScoreMode scoreMode() {
            return COMPLETE_NO_SCORES;
        }
    }
}

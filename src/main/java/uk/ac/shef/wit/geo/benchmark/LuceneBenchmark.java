package uk.ac.shef.wit.geo.benchmark;

import me.tongfei.progressbar.ProgressBar;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.LatLonShape;
import org.apache.lucene.document.ShapeField;
import org.apache.lucene.document.StoredField;
import org.apache.lucene.geo.Point;
import org.apache.lucene.geo.Polygon;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.search.FieldDoc;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.search.TopScoreDocCollector;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.NIOFSDirectory;
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
import java.util.List;
import java.util.concurrent.TimeUnit;

@State(Scope.Thread)
public class LuceneBenchmark
        extends AbstractBenchmark {

    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    private final String fieldName = "location";
    private IndexSearcher indexSearcher = null;

    @Setup
    public void setup()
            throws IOException {

        // MMapDirectory makes no difference to performance
//            Path tempDirectory = Files.createTempDirectory(this.getClass().getName());
        //            try (Directory directory = new MMapDirectory(tempDirectory);
        //            try (Directory directory = new RAMDirectory();
//                 IndexWriter indexWriter = new IndexWriter(directory, new IndexWriterConfig())) {

        createDirectory(outputDirectoryName);
        Path indexPath = Paths.get(outputDirectoryName, "lucene-polygons-index" + numberOfIndexPolygons);
        if (Files.exists(indexPath)) {
            final IndexReader indexReader = DirectoryReader.open(NIOFSDirectory.open(indexPath));
            indexSearcher = new IndexSearcher(indexReader);
            return;
        }

        try (Directory directory = NIOFSDirectory.open(indexPath);
             IndexWriter indexWriter = new IndexWriter(directory, new IndexWriterConfig())) {

            int id = 0;
            List<double[][]> polygons = getIndexPolygons();
            try (ProgressBar progressBar = new ProgressBar("Features:", numberOfIndexPolygons)) {
                for (double[][] latlons : polygons) {
                    try {
                        double[] lats = new double[latlons.length];
                        double[] lons = new double[latlons.length];
                        for (int i = 0; i < latlons.length; i++) {
                            lats[i] = latlons[i][0];
                            lons[i] = latlons[i][1];
                        }
                        Document doc = new Document();
                        doc.add(new StoredField("id", ++id));
                        Polygon polygon = new Polygon(lats, lons);
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
    @Measurement(iterations = 5)
    public void pointIntersectsQuery() {
        long candidateCount = 0;
        long nearestCount = 0;
        results.clear();
        for (double[] latlon : getQueryPoints()) {
            int id = 0;
            float distance = -1;
            try {
                Query query = LatLonShape.newGeometryQuery(fieldName,
                        ShapeField.QueryRelation.INTERSECTS, new Point(latlon[0], latlon[1]));

                TopScoreDocCollector collector = TopScoreDocCollector.create(1000, 1000);
                indexSearcher.search(query, collector);
                TopDocs topDocs = collector.topDocs();

                candidateCount += topDocs.totalHits.value;
                if (topDocs.totalHits.value != 0) {
                    ScoreDoc scoreDoc = topDocs.scoreDocs[0];
                    if (scoreDoc instanceof FieldDoc)
                        distance = ((Double) ((FieldDoc) scoreDoc).fields[0]).floatValue();
                    else
                        distance = scoreDoc.score;
                    Document doc = indexSearcher.doc(scoreDoc.doc);
                    id = Integer.parseInt(doc.get("id"));
                    nearestCount++;
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            results.add(new AbstractMap.SimpleImmutableEntry<>(id, (double) distance));
        }

        candidateCounts.add(nearestCount == 0 ? 0 : candidateCount / nearestCount);
        nearestCounts.add(nearestCount);
    }


    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 5)
    public void polygonIntersectsQuery() {

        long candidateCount = 0;
        long nearestCount = 0;
        results.clear();
        for (double[][] latlons : getQueryPolygons()) {
            int id = 0;
            float distance = -1;
            try {
                double[] lats = new double[latlons.length];
                double[] lons = new double[latlons.length];
                for (int i = 0; i < latlons.length; i++) {
                    lats[i] = latlons[i][0];
                    lons[i] = latlons[i][1];
                }
                Polygon polygon = new Polygon(lats, lons);
                Query query = LatLonShape.newGeometryQuery(fieldName,
                        ShapeField.QueryRelation.INTERSECTS, polygon);

                TopScoreDocCollector collector = TopScoreDocCollector.create(1000, 1000);
                indexSearcher.search(query, collector);
                TopDocs topDocs = collector.topDocs();


                candidateCount += topDocs.totalHits.value;
                if (topDocs.totalHits.value != 0) {
                    ScoreDoc scoreDoc = topDocs.scoreDocs[0];
                    if (scoreDoc instanceof FieldDoc)
                        distance = ((Double) ((FieldDoc) scoreDoc).fields[0]).floatValue();
                    else
                        distance = scoreDoc.score;
                    Document doc = indexSearcher.doc(scoreDoc.doc);
                    id = Integer.parseInt(doc.get("id"));
                    nearestCount++;
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            results.add(new AbstractMap.SimpleImmutableEntry<>(id, (double) distance));
        }

        candidateCounts.add(nearestCount == 0 ? 0 : candidateCount / nearestCount);
        nearestCounts.add(nearestCount);
    }

    @TearDown
    public void teardown()
            throws IOException {
        super.teardown();
        indexSearcher.getIndexReader().close();
    }


    public static void main(String[] args) throws IOException {
        LuceneBenchmark benchmark = new LuceneBenchmark();
        benchmark.setup();
        benchmark.pointIntersectsQuery();
        benchmark.polygonIntersectsQuery();
        benchmark.teardown();
    }
}

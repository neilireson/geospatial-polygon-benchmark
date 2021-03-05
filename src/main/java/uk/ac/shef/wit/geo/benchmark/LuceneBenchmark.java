package uk.ac.shef.wit.geo.benchmark;

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

import java.io.IOException;
import java.util.AbstractMap;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

@State(Scope.Thread)
public class LuceneBenchmark
        extends AbstractBenchmark {

    private final String fieldName = "location";
    private IndexSearcher indexSearcher = null;

    @Setup
    public void setup()
            throws IOException {

        // MMapDirectory makes no difference to performance
//            Path tempDirectory = Files.createTempDirectory(this.getClass().getName());
//            final Directory directory = new MMapDirectory(tempDirectory);
        final Directory directory = new RAMDirectory();
        IndexWriterConfig iwConfig = new IndexWriterConfig();
        IndexWriter indexWriter = new IndexWriter(directory, iwConfig);

        int id = 0;
        List<double[][]> polygons = getIndexPolygons();
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
                for (double[] latlon : latlons) {
                    System.out.println(" " + latlon[1] + "," + latlon[0]);
                }
                System.out.println();
                throw e;
            }
        }
        indexWriter.commit();
        indexWriter.close();
        final IndexReader indexReader = DirectoryReader.open(directory);
        indexSearcher = new IndexSearcher(indexReader);
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 1)
    public void pointIntersectsQuery() {
        pointQuery(latlon -> {
            Query query = LatLonShape.newGeometryQuery(fieldName,
                    ShapeField.QueryRelation.INTERSECTS, new Point(latlon[0], latlon[1]));

            TopScoreDocCollector collector = TopScoreDocCollector.create(1000, 1000);
            try {
                indexSearcher.search(query, collector);
                return collector.topDocs();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
    }

    public void pointQuery(Function<double[], TopDocs> getTopDocsFunction) {

        long candidateCount = 0;
        long nearestCount = 0;
        results.clear();
        for (double[] latlon : getQueryPoints()) {
            int id = 0;
            float distance = -1;
            try {
                TopDocs topDocs = getTopDocsFunction.apply(latlon);

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
    @Measurement(iterations = 1)
    public void polygonIntersectsQuery() {
        polygonQuery(polygon -> {
            Query query = LatLonShape.newGeometryQuery(fieldName,
                    ShapeField.QueryRelation.INTERSECTS, polygon);

            TopScoreDocCollector collector = TopScoreDocCollector.create(1000, 1000);
            try {
                indexSearcher.search(query, collector);
                return collector.topDocs();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
    }

    public void polygonQuery(Function<Polygon, TopDocs> getTopDocsFunction) {

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
                TopDocs topDocs = getTopDocsFunction.apply(polygon);

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
    public void teardown() {
        super.teardown();
    }
}

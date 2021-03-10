package uk.ac.shef.wit.geo.benchmark;

import com.mongodb.client.FindIterable;
import com.mongodb.client.model.Filters;
import com.mongodb.client.model.geojson.Polygon;
import com.mongodb.client.model.geojson.Position;
import com.mongodb.MongoClient;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import com.mongodb.client.model.IndexOptions;
import com.mongodb.client.model.Indexes;
import me.tongfei.progressbar.ProgressBar;
import org.bson.Document;
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
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

@State(Scope.Thread)
public class MongoDbBenchmark
        extends AbstractBenchmark {

    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    public static final String placesDbName = "BenchmarkPlacesDb";
    public static final String placesCollectionName = "places";
    // According to https://docs.mongodb.com/manual/tutorial/build-a-2d-index/
    // The (default) value 26 bots is equivalent tp a precision of ~60cm
    // which presumably is 40000/(2^26)
    public static final int geoHashBits = 26;
    private final String fieldName = "location";
    private final MongoClient mongoClient;

    public MongoDbBenchmark() {
        this.mongoClient = new MongoClient();
    }

    @Setup
    public void setup() {
        MongoDatabase db = mongoClient.getDatabase(placesDbName);
        for (String collectionName : db.listCollectionNames()) {
            logger.info("Collection: {}", collectionName);
        }
//        db.drop();
        MongoCollection<Document> collection = db.getCollection(placesCollectionName + numberOfIndexPolygons);
        long count = collection.countDocuments();
        if (count == numberOfIndexPolygons) return;
        else if (count > 0) {
            logger.error("Unexpected number of features in collection:{} expected={}, read={}",
                    placesCollectionName, numberOfIndexPolygons, count);
            collection.deleteMany(new Document());
        }

        IndexOptions indexOptions = new IndexOptions();
        indexOptions.bits(geoHashBits);
        collection.createIndex(Indexes.geo2dsphere(fieldName), indexOptions);

        try (ProgressBar progressBar = new ProgressBar("Features:", numberOfIndexPolygons)) {

            int id = 0;
            List<double[][]> polygons = getIndexPolygons();
            for (double[][] latlons : polygons) {
                try {
                    List<Position> positions = new ArrayList<>();
                    for (double[] latlon : latlons) {
                        positions.add(new Position(latlon[1], latlon[0]));
                    }
                    @SuppressWarnings("unchecked")
                    Polygon polygon = new Polygon(positions);

                    Document doc = new Document();
                    doc.append("id", ++id);
                    doc.append(fieldName, polygon);
                    collection.insertOne(doc);

                    progressBar.step();
                } catch (Exception e) {
                    for (double[] latlon : latlons) {
                        System.out.println(" " + latlon[1] + "," + latlon[0]);
                    }
                    System.out.println();
                    throw e;
                }
            }
        }
    }


    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 5)
    public void polygonQuery() {
        MongoDatabase db = mongoClient.getDatabase(placesDbName);
        MongoCollection<Document> collection = db.getCollection(placesCollectionName + numberOfIndexPolygons);
        long foundCount = 0;
        long intersectingCount = 0;
        results.clear();
        for (double[][] latlons : getQueryPolygons()) {
            int id = 0;
            float distance = -1;
            List<Position> positions = new ArrayList<>();
            for (double[] latlon : latlons) {
                positions.add(new Position(latlon[1], latlon[0]));
            }
            @SuppressWarnings("unchecked")
            Polygon polygon = new Polygon(positions);
            FindIterable<Document> found = collection.find(Filters.geoIntersects(fieldName, polygon));
            if (found.iterator().hasNext()) {
                foundCount++;
                for (Document doc : found) {
                    intersectingCount++;
                    id = doc.getInteger("id");
                }
            }

            results.add(new AbstractMap.SimpleImmutableEntry<>(id, (double) distance));
        }

        candidateCounts.add(foundCount == 0 ? 0 : intersectingCount);
        nearestCounts.add(foundCount);
    }

    @TearDown
    public void teardown() throws IOException {
        super.teardown();
        mongoClient.close();
    }


    public static void main(String[] args) throws IOException {
        MongoDbBenchmark benchmark = new MongoDbBenchmark();
        benchmark.setup();
        benchmark.polygonQuery();
        benchmark.teardown();
    }
}

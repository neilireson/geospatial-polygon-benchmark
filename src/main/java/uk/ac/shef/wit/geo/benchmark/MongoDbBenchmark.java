package uk.ac.shef.wit.geo.benchmark;

import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoIterable;
import com.mongodb.client.model.Filters;
import com.mongodb.client.model.geojson.Point;
import com.mongodb.client.model.geojson.Polygon;
import com.mongodb.client.model.geojson.Position;
import com.mongodb.MongoClient;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import com.mongodb.client.model.IndexOptions;
import com.mongodb.client.model.Indexes;
import me.tongfei.progressbar.ProgressBar;
import org.bson.Document;
import org.geotools.data.DataStore;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.json.simple.JSONValue;
import org.locationtech.jts.geom.Geometry;
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

import java.lang.invoke.MethodHandles;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
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
    MongoCollection<Document> collection;

    String collectionName = placesCollectionName + "_" + configName;

    public MongoDbBenchmark() {
        this.mongoClient = new MongoClient();
    }

    @Setup
    public void setup() {

        MongoDatabase db = mongoClient.getDatabase(placesDbName);
        //        db.drop();

        logger.info("Checking if database contains collection");
        MongoIterable<String> collections = db.listCollectionNames();
        boolean collectionExists = false;
        if (collections.iterator().hasNext()) {
            StringBuilder sb = new StringBuilder();
            for (String name : collections) {
                sb.append('\n').append("Collection: ").append(name);
                if (collectionName.equals(name)) {
                    collectionExists = true;
                    sb.append(" *");
                }
            }
            logger.info("Available collections:{}", sb);
        } else {
            logger.info("Database currently has no collections");
        }
        collection = db.getCollection(collectionName);

        Map.Entry<DataStore, SimpleFeatureCollection> dataStoreCollection = getIndexPolygons();
        DataStore dataStore = dataStoreCollection.getKey();
        SimpleFeatureCollection polygons = dataStoreCollection.getValue();

        if (collectionExists) {
            long count = collection.countDocuments();
            if (count == 0) {
                logger.error("Collection {} is empty", collectionName);
            } else if (config.getNumberOfIndexPoints() != null && config.getNumberOfIndexPoints() != count) {
                logger.error("Collection contains incorrect number of points. Expected {}, found {}",
                        config.getNumberOfIndexPoints(), collection.countDocuments());
                collection.deleteMany(new Document());
            } else {
                logger.info("Collection {} contains {} documents", collectionName, collection.countDocuments());
                return;
            }
        }

        IndexOptions indexOptions = new IndexOptions();
        indexOptions.bits(geoHashBits);
        collection.createIndex(Indexes.geo2dsphere(fieldName), indexOptions);

        logger.info("Indexing polygons...");
        try (ProgressBar progressBar = new ProgressBar("Features:", polygons.size());
             SimpleFeatureIterator it = polygons.features()) {
            int id = 0;
            while (it.hasNext()) {
                SimpleFeature feature = it.next();
                Geometry geometry = (Geometry) feature.getDefaultGeometry();
                try {
                    Document doc = new Document();
                    doc.append("id", String.valueOf(++id));
                    String geoJSON = geometryJSON.get().toString(geometry);
                    doc.append(fieldName, JSONValue.parse(geoJSON));
                    collection.insertOne(doc);

                    progressBar.step();
                } catch (Exception e) {
                    logger.error("Failed to add feature: " + feature);
                    throw e;
                }
            }
        } finally {
            if (dataStore != null) {
                dataStore.dispose();
            }
        }
    }


    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 1)
    public void pointQuery() {
        long foundCount = 0;
        long intersectingCount = 0;
        for (double[] latlon : getQueryPoints()) {
            String id = "0";
            float distance = -1;
            Point point = new Point(new Position(latlon[1], latlon[0]));
            FindIterable<Document> found = collection.find(Filters.geoIntersects(fieldName, point));
            if (found.iterator().hasNext()) {
                foundCount++;
                for (Document doc : found) {
                    intersectingCount++;
                    id = doc.getString("id");
                }
            }

            results.add(new AbstractMap.SimpleImmutableEntry<>(id, (double) distance));
        }

        candidateCounts.add(foundCount == 0 ? 0 : intersectingCount / foundCount);
        nearestCounts.add(foundCount);
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 1)
    public void polygonQuery() {
        long foundCount = 0;
        long intersectingCount = 0;
        for (double[][] latlons : getQueryPolygons()) {
            String id = "0";
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
                    id = doc.getString("id");
                }
            }

            results.add(new AbstractMap.SimpleImmutableEntry<>(id, (double) distance));
        }

        candidateCounts.add(foundCount == 0 ? 0 : intersectingCount / foundCount);
        nearestCounts.add(foundCount);
    }

    @TearDown
    public void teardown() {
        super.teardown();
    }


    public static void main(String[] args) {

        try (MongoDbBenchmark benchmark = new MongoDbBenchmark()) {
            benchmark.setup();
            benchmark.pointQuery();
            benchmark.teardown();
            benchmark.polygonQuery();
            benchmark.teardown();
        } catch (Exception e) {
            logger.error("", e);
        }
    }

    @Override
    public void close() {
        mongoClient.close();
    }
}

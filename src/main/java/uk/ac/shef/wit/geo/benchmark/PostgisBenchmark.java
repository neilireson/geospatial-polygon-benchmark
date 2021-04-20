package uk.ac.shef.wit.geo.benchmark;

import me.tongfei.progressbar.ProgressBar;
import org.geotools.data.DataStore;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
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
import org.postgresql.PGConnection;
import org.postgresql.util.PGobject;
import org.postgresql.util.PSQLException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.lang.invoke.MethodHandles;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.AbstractMap;
import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.stream.Collectors;

import static java.lang.Math.abs;


@State(Scope.Thread)
public class PostgisBenchmark
        extends AbstractBenchmark {

    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    public static final String placesDbName = "benchmarkplacesdb";
    public static final String placesTableName = "places";
    private final boolean overwrite = false;

    @Setup
    public void setup() {

        Statement stmt = null;
        try {
            Connection conn = getConnection("");
            logger.info("Creating database...");
            stmt = conn.createStatement();

            String sql = "CREATE DATABASE " + placesDbName;
            try {
                stmt.executeUpdate(sql);
                logger.info("Database created successfully...");
            } catch (PSQLException e) {
                if (("ERROR: database \"" + placesDbName + "\" already exists").equals(e.getLocalizedMessage()))
                    logger.warn("{}", e.getLocalizedMessage());
                else throw e;
            }

            logger.info("Closing connection and reopening on database...");
            conn.close();

            logger.info("Connecting to database...");
            conn = getConnection(placesDbName);
            stmt = conn.createStatement();

            logger.info("Adding EXTENSION postgis...");

            sql = "CREATE EXTENSION postgis";
            try {
                stmt.executeUpdate(sql);
                logger.info("Successfully added postgis");
            } catch (PSQLException e) {
                if ("ERROR: extension \"postgis\" already exists".equals(e.getLocalizedMessage()))
                    logger.warn("{}", e.getLocalizedMessage());
                else throw e;
            }

            if (overwrite) {
                logger.info("Dropping table...");
                try {
                    stmt.executeUpdate("DROP TABLE " + placesTableName);
                    logger.info("Table created successfully");
                } catch (PSQLException e) {
                    if (("ERROR: table \"" + placesTableName + "\" does not exist").equals(e.getLocalizedMessage()))
                        logger.warn("{}", e.getLocalizedMessage());
                    else throw e;
                }
            }
            //STEP 4: Execute create database table

            boolean tableExists = false;
            logger.info("Creating table...");
            sql = "CREATE TABLE " + placesTableName + " (id text, geometry geometry(MULTIPOLYGON,4326), " +
                    "CONSTRAINT uid UNIQUE (id))";
            try {
                stmt.executeUpdate(sql);
                logger.info("Table created successfully");

                logger.info("Creating GIST index...");
                sql = "CREATE INDEX geometry_idx ON " + placesTableName + " USING GIST (geometry)";
                stmt.executeUpdate(sql);
            } catch (PSQLException e) {
                if (("ERROR: relation \"" + placesTableName + "\" already exists").equals(e.getLocalizedMessage())) {
                    logger.warn("{}", e.getLocalizedMessage());
                    tableExists = true;
                } else throw e;
            }

            Map.Entry<DataStore, SimpleFeatureCollection> dataStoreCollection = getIndexPolygons();
            DataStore dataStore = dataStoreCollection.getKey();
            SimpleFeatureCollection polygons = dataStoreCollection.getValue();

            if (tableExists) {
                int count;
                ResultSet rs = stmt.executeQuery("SELECT COUNT(*) FROM " + placesTableName);
                rs.next();
                count = rs.getInt(1);

                if (count == 0) {
                    logger.error("Table {} is empty", placesTableName);
                } else if (polygons.size() != count) {
                    logger.error("Table contains incorrect number of rows. Expected {}, found {}",
                            polygons.size(), count);
                    if (abs(polygons.size() - count) <= config.getMissingDataThreshold()) {
                        logger.info("Missing number of rows is within acceptable limit, {} <= {}",
                                abs(polygons.size() - count), config.getMissingDataThreshold());
                        return;
                    } else {
                        logger.info("Deleting all documents in table {}", placesTableName);
                        stmt.executeUpdate("TRUNCATE " + placesTableName);
                    }
                } else {
                    logger.info("Table {} contains {} documents", placesTableName, count);
                    return;
                }
            }

            logger.info("Indexing polygons...");
            try (ProgressBar progressBar = new ProgressBar("Features:", polygons.size());
                 PreparedStatement statement = conn.prepareStatement(
                         "INSERT INTO " + placesTableName + "(id, geometry) " +
                                 "VALUES (?, ST_GeomFromText(?, 4326))");
                 SimpleFeatureIterator it = polygons.features()) {

                Function<SimpleFeature, SimpleFeature> simplificationFunction =
                        getSimplificationFunction(polygons);

                logger.info("Turn off table indexing...");
                sql = "ALTER TABLE " + placesTableName + " SET UNLOGGED";
                stmt.executeUpdate(sql);
                logger.info("Table altered successfully.");

                int id = 0;
                while (it.hasNext()) {
                    SimpleFeature feature = it.next();
                    if (simplificationFunction != null) {
                        feature = simplificationFunction.apply(feature);
                    }
                    Geometry geometry = (Geometry) feature.getDefaultGeometry();
                    try {
                        statement.setString(1, String.valueOf(++id));
                        statement.setString(2, geometry.toString());
                        statement.executeUpdate();

                        progressBar.step();
                    } catch (Exception e) {
                        logger.error("Failed to add feature: " + feature, e);
                    }
                }

                logger.info("Turn on table indexing");
                sql = "ALTER TABLE " + placesTableName + " SET LOGGED";
                stmt.executeUpdate(sql);
                logger.info("Indexes created successfully.");

            } catch (SQLException e) {
                logger.error("Failed to index", e);
            } finally {
                if (dataStore != null) {
                    dataStore.dispose();
                }
            }
        } catch (SQLException | ClassNotFoundException e) {
            throw new RuntimeException("Failed to setup PostGIS database", e);
        } finally {
            if (stmt != null) {
                try {
                    stmt.close();
                } catch (SQLException e) {
                    logger.warn("Failed to close statement", e);
                }
            }
        }
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 3)
    public void pointQuery() {
        long foundCount = 0;
        long intersectingCount = 0;
        try (Connection conn = getConnection(placesDbName);
             Statement s = conn.createStatement()) {

            for (double[] latlon : getQueryPoints()) {
                String id = "0";
                float distance = -1;
                ResultSet r = s.executeQuery("SELECT id FROM " + placesTableName + " WHERE " +
                        "ST_Intersects(geometry, 'SRID=4326;POINT(" + latlon[1] + " " + latlon[0] + " )');");
                while (r.next()) {
                    id = r.getString(1);
                    intersectingCount++;
                }
                if (!"0".equals(id)) {
                    foundCount++;
                }
                results.add(new AbstractMap.SimpleImmutableEntry<>(id, (double) distance));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        candidateCounts.add(foundCount == 0 ? 0 : intersectingCount / foundCount);
        nearestCounts.add(foundCount);
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 3)
    public void polygonQuery() {
        long foundCount = 0;
        long intersectingCount = 0;
        String query = "";
        try (Connection conn = getConnection(placesDbName);
             Statement s = conn.createStatement()) {

            for (double[][] latlons : getQueryPolygons()) {
                String id = "0";
                float distance = -1;

                query = "SELECT id FROM " + placesTableName + " WHERE " +
                        "ST_Intersects(geometry, 'SRID=4326;POLYGON((" +
                        Arrays.stream(latlons)
                                .map(latlon -> latlon[1] + " " + latlon[0])
                                .collect(Collectors.joining(",")) +
                        "))')";
                ResultSet r = s.executeQuery(query);
                while (r.next()) {
                    id = r.getString(1);
                    intersectingCount++;
                }
                if (!"0".equals(id)) {
                    foundCount++;
                }
                results.add(new AbstractMap.SimpleImmutableEntry<>(id, (double) distance));
            }
        } catch (Exception e) {
            logger.error("Query failed: {}", query, e);
        }

        candidateCounts.add(foundCount == 0 ? 0 : intersectingCount / foundCount);
        nearestCounts.add(foundCount);
    }

    @TearDown
    public void teardown() {
        super.teardown();
    }

    private static Connection getConnection(String database)
            throws ClassNotFoundException, SQLException {
        /*
         * Load the JDBC driver and establish a connection.
         */
        logger.info("Connecting to database...");
        Class.forName("org.postgresql.Driver");
        Connection conn = DriverManager.getConnection("jdbc:postgresql://localhost:6432/" + database, "postgres", "Q]b[+2=7E34SUZ2?");
        /*
         * Add the geometry types to the connection. Note that you
         * must cast the connection to the pgsql-specific connection
         * implementation before calling the addDataType() method.
         */
        //noinspection unchecked
        ((PGConnection) conn).addDataType("geometry",
                (Class<? extends PGobject>) Class.forName("org.postgis.PGgeometry"));

        return conn;
    }

    @Override
    public void close() {
    }


    public static void main(String[] args) {
        try (PostgisBenchmark benchmark = new PostgisBenchmark()) {
            benchmark.setup();
            benchmark.pointQuery();
            benchmark.teardown();
            benchmark.polygonQuery();
            benchmark.teardown();
        } catch (Exception e) {
            logger.error("", e);
        }
    }
}

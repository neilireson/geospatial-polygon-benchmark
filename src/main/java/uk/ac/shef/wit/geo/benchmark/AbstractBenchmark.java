package uk.ac.shef.wit.geo.benchmark;

import me.tongfei.progressbar.ProgressBar;
import org.geotools.data.DataStore;
import org.geotools.data.DataUtilities;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.Transaction;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geojson.geom.GeometryJSON;
import org.geotools.geometry.jts.JTS;
import org.geotools.map.FeatureLayer;
import org.geotools.map.Layer;
import org.geotools.map.MapContent;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.geotools.styling.SLD;
import org.geotools.styling.Style;
import org.geotools.swing.JMapFrame;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.hull.ConcaveHull;
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.results.RunResult;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import static java.lang.Math.*;

@State(Scope.Benchmark)
public abstract class AbstractBenchmark implements AutoCloseable {

    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    static final GeometryFactory gf = new GeometryFactory();
    private static final String outputDirectoryName = "out";
    private static final Random random = new Random();

    final int decimalPrecision = 8;
    // todo an alternative would be to use Jackson, which may be faster
    // see https://mvnrepository.com/artifact/org.n52.jackson/jackson-datatype-jts
    final ThreadLocal<GeometryJSON> geometryJSON = ThreadLocal.withInitial(() -> new GeometryJSON(decimalPrecision));

    //        @Param({"10000", "100000"})
    final String configName = "england-landuse";
    TrialConfiguration config;

    final List<Long> candidateCounts = new ArrayList<>();
    final List<Long> nearestCounts = new ArrayList<>();
    final List<Map.Entry<String, Double>> results = new ArrayList<>();

    final SimpleFeatureType polygonFeature;

    {
        SimpleFeatureTypeBuilder polyTypeBuilder = new SimpleFeatureTypeBuilder();
        polyTypeBuilder.setName("Polygon");
        polyTypeBuilder.setNamespaceURI("Polygon");
        polyTypeBuilder.setCRS(DefaultGeographicCRS.WGS84);
        // Note "the_geom" seems to be necessary name for shapefile read/write
        polyTypeBuilder.add("the_geom", Polygon.class);
        polyTypeBuilder.setDefaultGeometry("the_geom");
        polygonFeature = polyTypeBuilder.buildFeatureType();
    }

    public AbstractBenchmark() {
        File outputDirectory = createDirectory(outputDirectoryName);

        Pattern configFilePattern = Pattern.compile("^" + configName + "(\\.(cnf|conf|config|xml))?$");
        File[] configFiles = outputDirectory.listFiles(file -> {
            if (!file.isFile()) return false;
            return configFilePattern.matcher(file.getName()).matches();
        });

        if (configFiles == null || configFiles.length == 0)
            throw new IllegalArgumentException("Configuration file not found: " + configName);
        else if (configFiles.length > 1)
            throw new IllegalArgumentException("Multiple configuration files found: " + configName + ": " +
                    Arrays.toString(configFiles));

        logger.info("Reading trial configuration file: {}", configFiles[0].getPath());
        config = TrialConfiguration.create(configFiles[0].getPath());
    }

    synchronized Map.Entry<DataStore, SimpleFeatureCollection> getIndexPolygons() {

        DataStore dataStore = null;
        final SimpleFeatureCollection polygons;

        String fileNamePrefix = "benchmark-index-polygons-" + configName;
        Path generatedShapefile = Paths.get(fileNamePrefix + ".shp");

        if (config.getShapeFile() != null) {
            Path shapefile = Paths.get(config.getShapeFile());
            if (Files.exists(shapefile)) {
                logger.info("Reading from shapefile: {}", shapefile);
                try {
                    dataStore = FileDataStoreFinder.getDataStore(shapefile.toFile());
                    SimpleFeatureSource featureSource = dataStore.getFeatureSource(config.getTypeName());
                    polygons = featureSource.getFeatures();
                } catch (IOException e) {
                    if (dataStore != null) {
                        dataStore.dispose();
                    }
                    throw new RuntimeException("Failed to read features from shapefile: " + shapefile, e);
                }
            } else {
                throw new RuntimeException("Index source Shapefile does not exist: " + shapefile);
            }
        } else if (Files.exists(generatedShapefile)) {
            try {
                logger.info("Reading from generated shapefile: {}", generatedShapefile);
                dataStore = FileDataStoreFinder.getDataStore(generatedShapefile.toFile());
                SimpleFeatureSource featureSource = dataStore.getFeatureSource(fileNamePrefix);
                polygons = featureSource.getFeatures();
            } catch (IOException e) {
                if (dataStore != null) {
                    dataStore.dispose();
                }
                throw new RuntimeException("Failed to read features from shapefile: " + generatedShapefile, e);
            }
        } else {
            if (config.getNumberOfIndexPolygons() == null) {
                throw new NullPointerException("Configuration file specifies neither shapefile nor number of indexed polygons");
            }
            List<double[][]> indexPolygons = new ArrayList<>();
            for (int i = 0; i < config.getNumberOfIndexPolygons(); i++) {
                indexPolygons.add(createPolygon(config.getBoundingBox(), null));
            }
            ArrayList<SimpleFeature> features = new ArrayList<>();
            try (ProgressBar progressBar = new ProgressBar("Features:", config.getNumberOfIndexPolygons())) {
                for (double[][] latlons : indexPolygons) {
                    progressBar.step();

                    Coordinate[] coords = new Coordinate[latlons.length];
                    for (int i = 0; i < coords.length; i++) {
                        coords[i] = new Coordinate(latlons[i][1], latlons[i][0]);
                    }
                    Polygon polygon = gf.createPolygon(coords);
                    SimpleFeature feature = createSimpleFeature(polygonFeature, polygon);

                    if (feature != null) {
                        features.add(feature);
                    } else {
                        logger.error("Not a valid feature");
                    }
                }
            }
            try {
                writeToShapefile(generatedShapefile.toFile(), polygonFeature, features);
            } catch (IOException e) {
                logger.error("Failed to write shapefile", e);
            }

            logger.info("Wrapping features in a collection");
            polygons = DataUtilities.collection(features);
        }

        checkBoundingBox(polygons);

        return new AbstractMap.SimpleImmutableEntry<>(dataStore, polygons);
    }

    void checkBoundingBox(SimpleFeatureCollection indexFeatureCollection) {

        Path queryPointsFilePath = getQueryPointsFilePath();
        Path queryPolygonsFilePath = getQueryPolygonsFilePath();
        if (!Files.exists(queryPointsFilePath) || !Files.exists(queryPolygonsFilePath)) {
            logger.info("Creating points/polygon query files");
            Geometry convexHullGeometry = convexHull(indexFeatureCollection);
            Polygon convexHull = convexHullGeometry instanceof Polygon ? (Polygon) convexHullGeometry : null;
            Envelope boundingBox = convexHull == null ? config.getBoundingBox() : convexHullGeometry.getEnvelopeInternal();
            if (boundingBox == null) {
                throw new NullPointerException("Boundary box not available from configuration file or index feature collection");
            }

            // create query points and polygons that are within the index polygons
            getQueryPoints(boundingBox, convexHull);
            getQueryPolygons(boundingBox, convexHull);
        }
    }

    public static SimpleFeature createSimpleFeature(SimpleFeatureType schema, Geometry geometry) {
        if (geometry != null) {
            if (geometry.isValid()) {
                SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(schema);
                featureBuilder.add(geometry);
                return featureBuilder.buildFeature(null);
            } else {
                logger.error("Feature geometry is not valid: " + geometry);
                DrawShapes.draw(geometry);
                throw new IllegalArgumentException("Feature geometry is not valid");
            }
        }
        return null;
    }

    private void writeToShapefile(File shapefile, SimpleFeatureType type, ArrayList<SimpleFeature> features)
            throws IOException {

        ShapefileDataStore dataStore = null;
        try {
            ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();

            Map<String, Serializable> params = new HashMap<>();
            params.put("url", shapefile.toURI().toURL());
            params.put("create spatial index", Boolean.TRUE);

            dataStore = (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);

            /*
             * TYPE is used as a template to describe the file contents
             */
            dataStore.createSchema(type);

            /*
             * Write the features to the shapefile
             */
            Transaction transaction = new DefaultTransaction("create");

            String[] typeNames = dataStore.getTypeNames();
            String typeName = typeNames[0];
            SimpleFeatureSource featureSource = dataStore.getFeatureSource(typeName);
            SimpleFeatureType SHAPE_TYPE = featureSource.getSchema();
            /*
             * The Shapefile format has a couple limitations:
             * - "the_geom" is always first, and used for the geometry attribute name
             * - "the_geom" must be of type Point, MultiPoint, MuiltiLineString, MultiPolygon
             * - Attribute names are limited in length
             * - Not all data types are supported (example Timestamp represented as Date)
             *
             * Each data store has different limitations so check the resulting SimpleFeatureType.
             */
            logger.info("SHAPE:{}", SHAPE_TYPE);

            if (featureSource instanceof SimpleFeatureStore) {
                SimpleFeatureStore featureStore = (SimpleFeatureStore) featureSource;
                /*
                 * SimpleFeatureStore has a method to add features from a
                 * SimpleFeatureCollection object, so we use the ListFeatureCollection
                 * class to wrap our list of features.
                 */
                SimpleFeatureCollection collection = new ListFeatureCollection(type, features);
                featureStore.setTransaction(transaction);
                try {
                    featureStore.addFeatures(collection);
                    transaction.commit();
                } catch (Exception e) {
                    logger.error("Failed to add features to shapefile", e);
                    transaction.rollback();
                } finally {
                    transaction.close();
                }
            } else {
                logger.error("{} does not support read/write access", typeName);
            }
        } finally {
            if (dataStore != null)
                dataStore.dispose();
        }
    }

    Function<SimpleFeature, SimpleFeature> getSimplificationFunction(SimpleFeatureCollection collection) {
        switch (config.getSimplificationType()) {
            case DouglasPeucker:
            case TopologyPreserving:
                if (config.getSimplificationThreshold() <= 0) {
                    logger.error("Simplification type set to {} but threshold is {}, which is effectively a NoOp. " +
                                    "Setting type to {}",
                            config.getSimplificationType(), config.getSimplificationThreshold(),
                            TrialConfiguration.SimplificationType.None);
                    config.setSimplificationType(TrialConfiguration.SimplificationType.None);
                }
        }

        if (!config.getRemoveHoles() &&
                (config.getSimplificationType() == null ||
                        config.getSimplificationType() == TrialConfiguration.SimplificationType.None)) {
            return null;
        }

        final SimpleFeatureBuilder fb = new SimpleFeatureBuilder(collection.getSchema());
        return feature -> {
            for (Object attribute : feature.getAttributes()) {
                if (attribute instanceof Geometry) {
                    Geometry geometry = (Geometry) attribute;
                    if (config.getRemoveHoles()) {
                        geometry = geometry.getBoundary();
                    }
                    switch (config.getSimplificationType()) {
                        case BoundingBox:
                            geometry = geometry.getEnvelope();
                            break;
                        case ConvexHull:
                            geometry = geometry.convexHull();
                            break;
                        case ConcaveHull:
                            geometry = ConcaveHull.compute(geometry);
                            break;
                        case DouglasPeucker:
                            geometry = DouglasPeuckerSimplifier.simplify(geometry,
                                    config.getSimplificationThreshold());
                            break;
                        case TopologyPreserving:
                            geometry = TopologyPreservingSimplifier.simplify(geometry,
                                    config.getSimplificationThreshold());
                            break;
                    }
                    attribute = geometry;
                }
                fb.add(attribute);
            }
            return fb.buildFeature(feature.getID());
        };
    }


    List<double[]> getQueryPoints() {
        return getQueryPoints(config.getBoundingBox(), null);
    }

    List<double[]> getQueryPoints(Envelope boundingBox, Polygon polygon) {
        return getPoints(getQueryPointsFilePath(), config.getNumberOfQueryPoints(), boundingBox, polygon);
    }

    List<double[][]> getQueryPolygons() {
        return getQueryPolygons(config.getBoundingBox(), null);
    }

    List<double[][]> getQueryPolygons(Envelope boundingBox, Polygon polygon) {
        return getPolygons(getQueryPolygonsFilePath(), config.getNumberOfQueryPolygons(),
                boundingBox, polygon);
    }

    private Path getQueryPointsFilePath() {
        createDirectory(outputDirectoryName);
        String filename = "benchmark-query-points-" + configName + "-" + config.getNumberOfQueryPoints() + ".csv.gz";
        return Paths.get(outputDirectoryName, filename);
    }

    private Path getQueryPolygonsFilePath() {
        createDirectory(outputDirectoryName);
        String filename = "benchmark-query-polygons-" + configName + "-" + config.getNumberOfQueryPolygons() + ".csv.gz";
        return Paths.get(outputDirectoryName, filename);
    }

    private synchronized List<double[]> getPoints(Path pointFilePath, int numberOfPoints,
                                                  Envelope boundingBox, Polygon polygon) {
        List<double[]> indexPoints = new ArrayList<>();
        if (Files.exists(pointFilePath)) {
            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(new GZIPInputStream(new FileInputStream(pointFilePath.toFile()))))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    String[] latlon = line.split(",");
                    indexPoints.add(new double[]{Double.parseDouble(latlon[0]), Double.parseDouble(latlon[1])});
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            if (indexPoints.size() != numberOfPoints) {
                throw new RuntimeException("File contains incorrect number of points. Expected " +
                        numberOfPoints + " found " + indexPoints.size());
            }
        } else {
            logger.info("Creating {} points", numberOfPoints);
            try (BufferedWriter writer = new BufferedWriter(
                    new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(pointFilePath.toFile()))))) {
                for (int i = 0; i < numberOfPoints; i++) {
                    double[] latlon = createRandomLatLon(boundingBox, polygon);
                    indexPoints.add(latlon);
                    writer.write(latlon[0] + "," + latlon[1]);
                    writer.newLine();
                }
            } catch (IOException e) {
                try {
                    Files.delete(pointFilePath);
                } catch (IOException e2) {
                    logger.error("Failed to delete points file: {}", pointFilePath);
                }
                throw new RuntimeException("Failed to create points", e);
            }
        }
        return indexPoints;
    }

    private synchronized List<double[][]> getPolygons(Path polygonFilePath, int numberOfPolygons,
                                                      Envelope boundingBox,
                                                      Polygon polygon) {

        List<double[][]> indexPolygons = new ArrayList<>();
        if (Files.exists(polygonFilePath)) {
            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(new GZIPInputStream(new FileInputStream(polygonFilePath.toFile()))))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    String[] xys = line.split(" ");
                    double[][] latlons = new double[xys.length][];
                    for (int i = 0; i < latlons.length; i++) {
                        String[] xy = xys[i].split(",");
                        latlons[i] = new double[]{Double.parseDouble(xy[1]), Double.parseDouble(xy[0])};
                    }
                    indexPolygons.add(latlons);
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            if (indexPolygons.size() != numberOfPolygons) {
                throw new RuntimeException("File contains incorrect number of points. Expected " +
                        numberOfPolygons + " found " + indexPolygons.size());
            }
        } else {
            logger.info("Creating {} polygons", numberOfPolygons);
            try (BufferedWriter writer = new BufferedWriter(
                    new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(polygonFilePath.toFile()))))) {
                for (int i = 0; i < numberOfPolygons; i++) {
                    double[][] latlons = createPolygon(boundingBox, polygon);
                    indexPolygons.add(latlons);
                    writer.write(latlons[0][1] + "," + latlons[0][0]);
                    for (int j = 1; j < latlons.length; j++) {
                        writer.write(" " + latlons[j][1] + "," + latlons[j][0]);
                    }
                    writer.newLine();
                }
            } catch (IOException e) {
                try {
                    Files.delete(polygonFilePath);
                } catch (IOException e2) {
                    logger.error("Failed to delete polygon file: {}", polygonFilePath);
                }
                throw new RuntimeException("Failed to create points", e);
            }
        }
        return indexPolygons;
    }

    File getOutputDirectory() {
        return createDirectory(outputDirectoryName);
    }

    private File createDirectory(String directoryPath) {
        Path path = Paths.get(directoryPath);
        if (!Files.exists(path)) {
            if (!path.toFile().mkdirs()) {
                throw new RuntimeException("Failed to create output directory: " + path.toAbsolutePath());
            }
        } else if (!Files.isDirectory(path)) {
            throw new RuntimeException("Output directory is not a directory: " + path.toAbsolutePath());
        }
        return path.toFile();
    }

    private double[] createRandomLatLon(Envelope boundingBox) {
        return createRandomLatLon(boundingBox, null);
    }

    private double[] createRandomLatLon(Envelope boundingBox, Polygon polygon) {
        double latitude, longitude;
        do {
            latitude = ThreadLocalRandom.current().nextDouble(boundingBox.getMinY(), boundingBox.getMaxY());
            longitude = ThreadLocalRandom.current().nextDouble(boundingBox.getMinX(), boundingBox.getMaxY());
        } while (polygon != null && !polygon.contains(gf.createPoint(new Coordinate(longitude, latitude))));
        return new double[]{latitude, longitude};
    }

    double[][] createPolygon(Envelope boundingBox, Polygon polygon) {
        final double[] latlon = createRandomLatLon(boundingBox, polygon);
        return createPolygon(latlon[0], latlon[1]);
    }

    private double[][] createPolygon(double lat, double lon) {
        return createPolygon((float) lat, (float) lon,
                config.getMinVertices() + random.nextInt(config.getMaxVertices() - config.getMinVertices()),
                config.getIrregularity(), config.getSpikiness(), config.getMinRadius(), config.getMaxRadius(), config.getLambda());
    }

    /**
     * Start with the centre of the polygon at lon, lat,
     * then creates the polygon by sampling points on a circle around the centre.
     * Randon noise is added by varying the angular spacing between sequential points,
     * and by varying the radial distance of each point from the centre.
     * Adapted from: https://stackoverflow.com/questions/8997099/algorithm-to-generate-random-2d-polygon
     *
     * @param latitude     latitude of the "centre" of the polygon
     * @param longitude    longitude of the "centre" of the polygon
     * @param numVertices  number of vertices
     * @param irregularity [0,1] indicating how much variance there is in the angular spacing of vertices. [0,1] will map to [0, 2pi/numberOfVerts]
     * @param spikiness    [0,1] indicating how much variance there is in each vertex from the circle of radius averageRadius. [0,1] will map to [0, averageRadius]
     * @param minRadius    the minimum radius of this polygon, in metres
     * @param maxRadius    the maximum radius of this polygon, in metres
     * @param lambda       parameter used in bounded exponential function to determine polygon radius
     * @return an array of vertices forming a polygon.
     */
    public static double[][] createPolygon(float latitude, float longitude,
                                           int numVertices,
                                           float irregularity, float spikiness,
                                           float minRadius, float maxRadius, double lambda) {

        double radiusExp = -log(1 - (1 - exp(-lambda)) * random.nextDouble()) / lambda;
        float averageRadius = (float) (minRadius + (radiusExp * (maxRadius - minRadius)));

        if (minRadius <= 0) {
            throw new RuntimeException("Cannot create random polygons if the average radius is zero");
        }

        double lat = toRadians(latitude);
        double lon = toRadians(longitude);

        // Radius of Earth at given latitude
        double latRadius = WGS84EarthRadius(lat);

        irregularity = (float) (clip(irregularity, 0, 1) * 2 * PI / numVertices);
        spikiness = (float) (clip(spikiness, 0, 1) * averageRadius);

        // generate n angle steps
        double[] angleSteps = new double[numVertices];
        double sum = 0;
        for (int i = 0; i < numVertices; i++) {
            angleSteps[i] = 1 + (random.nextDouble() * irregularity);
            sum += angleSteps[i];
        }

        // normalize the steps so that point 0 and point n+1 are the same
        double sumRadians = sum / (2 * PI);
        for (int i = 0; i < numVertices; i++) {
            angleSteps[i] = angleSteps[i] / sumRadians;
        }
        // now generate the points
        double[][] points = new double[numVertices + 1][];
        double angle = 0;
        for (int i = 0; i < numVertices; i++) {
            angle += angleSteps[i];
            float vertLat, vertLon;
            do {
                float d = (float) (clip(random.nextGaussian() * spikiness + averageRadius, 0, 2 * averageRadius)
                        / latRadius);
                vertLat = (float) toDegrees(lat + (sin(angle) * d));
                vertLon = (float) toDegrees(lon + (cos(angle) * d));
            } while (latitude == vertLat && longitude == vertLon);
            points[i] = new double[]{vertLat, vertLon};
        }
        points[numVertices] = points[0];

        return points;
    }

    Geometry convexHull(SimpleFeatureCollection fc) {
        logger.info("Creating convex hull for feature collection...");
        ArrayList<Geometry> geoms = new ArrayList<>();
        try (ProgressBar pb = new ProgressBar("Features:", fc.size());
             SimpleFeatureIterator it = fc.features()) {
            if (it.hasNext()) {
                pb.step();
                geoms.add(((Geometry) it.next().getDefaultGeometry()).convexHull());
                while (it.hasNext()) {
                    pb.step();
                    geoms.add(JTS.toGeometry(it.next().getBounds()));
                    GeometryCollection geometryCollection = (GeometryCollection) gf.buildGeometry(geoms);
                    geoms.clear();
                    geoms.add(geometryCollection.union().convexHull());
                }
            }
        }
        return geoms.isEmpty() ? null : geoms.get(0);
    }

    // Semi-axes of WGS-84 geoidal reference
    private static final double WGS84_a = 6378137.0; // Major semiaxis [metres]
    private static final double WGS84_b = 6356752.3; // Minor semiaxis [metres]

    // Earth radius at a given latitude, according to the WGS-84 ellipsoid [metres]
    public static double WGS84EarthRadius(double lat) {
        // http://en.wikipedia.org/wiki/Earth_radius
        double An = WGS84_a * WGS84_a * cos(lat);
        double Bn = WGS84_b * WGS84_b * sin(lat);
        double Ad = WGS84_a * cos(lat);
        double Bd = WGS84_b * sin(lat);
        return sqrt((An * An + Bn * Bn) / (Ad * Ad + Bd * Bd));
    }

    private static double clip(double x, double min, double max) {
        return (x < min) ? min : min(x, max);
    }

    protected void teardown() {
        System.out.format("%n%s: average number of candidates per successful query = %.0f, " +
                        "number of intersecting polygons found = %.0f/%d%n",
                getClass().getSimpleName(),
                candidateCounts.stream().mapToLong(x -> x).average().orElse(0),
                nearestCounts.stream().mapToLong(x -> x).average().orElse(0),
                results.size());

        // write results
        String filename = "results" +
                "-" + this.getClass().getSimpleName() +
                "-" + configName +
                ".csv";
        Path path = Paths.get(outputDirectoryName, filename);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(path.toFile()))) {
            for (Map.Entry<String, Double> result : results) {
                writer.write(result.getKey());
                writer.write('\t');
                writer.write(String.valueOf(result.getValue()));
                writer.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        candidateCounts.clear();
        nearestCounts.clear();
        results.clear();
    }

    private void showFeatures(SimpleFeatureSource featureSource) {
        // Create a map content and add our shapefile to it
        MapContent map = new MapContent();
        map.setTitle("Quickstart");

        Style style = SLD.createSimpleStyle(featureSource.getSchema());
        Layer layer = new FeatureLayer(featureSource, style);
        map.addLayer(layer);

        // Now display the map
        JMapFrame.showMap(map);
    }


    public static void main(String[] args) {
//        LuceneBenchmark benchmark = new LuceneBenchmark();
//        DrawPolygons drawPolygons = new DrawPolygons(benchmark.minLon, benchmark.maxLon, benchmark.minLat, benchmark.maxLat);
//        List<double[][]> indexPolygons = benchmark.getIndexPolygons();
//        List<double[][]> queryPolygons = benchmark.getQueryPolygons();
//        drawPolygons.addPolygons(indexPolygons, Color.GRAY);
//        drawPolygons.addPolygons(queryPolygons, Color.RED);
//
//        benchmark.setup();
//        benchmark.polygonIntersectsQuery();

        Options opt = new OptionsBuilder()
//                .addProfiler(GCProfiler.class)
                .addProfiler(MaxMemoryProfiler.class)
//                .include(GeotoolsBenchmark.class.getSimpleName())
                .include(LuceneBenchmark.class.getSimpleName())
                .include(MongoDbBenchmark.class.getSimpleName())
                .build();

        try {
            Collection<RunResult> runResults = new Runner(opt).run();

            //        for (RunResult runResult : runResults) {
//            for (BenchmarkResult benchmarkResult : runResult.getBenchmarkResults()) {
//                Result primaryResult = benchmarkResult.getPrimaryResult();
//                BenchmarkParams params = benchmarkResult.getParams();
//                System.out.format("%s\t%s\t%.3f\t%s%n",
//                        params.getBenchmark(),
//                        params.getMode().shortLabel(),
//                        primaryResult.getScore(),
//                        primaryResult.getScoreUnit());
//            }
//        }
        } catch (RunnerException e) {
            logger.error("", e);
        }
    }
}
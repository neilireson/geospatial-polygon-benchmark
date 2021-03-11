package uk.ac.shef.wit.geo.benchmark;

import me.tongfei.progressbar.ProgressBar;
import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.DataUtilities;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.Transaction;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.data.collection.SpatialIndexFeatureCollection;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.factory.CommonFactoryFinder;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geometry.jts.JTS;
import org.geotools.map.FeatureLayer;
import org.geotools.map.Layer;
import org.geotools.map.MapContent;
import org.geotools.referencing.GeodeticCalculator;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.geotools.styling.SLD;
import org.geotools.styling.Style;
import org.geotools.swing.JMapFrame;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.distance.DistanceOp;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.filter.Filter;
import org.opengis.filter.FilterFactory2;
import org.opengis.referencing.operation.TransformException;
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
import java.io.Serializable;
import java.lang.invoke.MethodHandles;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

@State(Scope.Thread)
public class GeotoolsBenchmark
        extends AbstractBenchmark {

    private enum DatastoreType {shapefile, h2}

    private final DatastoreType datastoreType = DatastoreType.shapefile;
    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    private static final FilterFactory2 ff = CommonFactoryFinder.getFilterFactory2();
    private SpatialIndexFeatureCollection index;

    @Setup
    public void setup() {
        SimpleFeatureTypeBuilder polyTypeBuilder = new SimpleFeatureTypeBuilder();
        polyTypeBuilder.setName("Polygon");
        polyTypeBuilder.setNamespaceURI("Polygon");
        polyTypeBuilder.setCRS(DefaultGeographicCRS.WGS84);
        // Note "the_geom" seems to be necessary name for shapefile read/write
        polyTypeBuilder.add("the_geom", Polygon.class);
        polyTypeBuilder.setDefaultGeometry("the_geom");
        polyTypeBuilder.add("id", Integer.class);
        SimpleFeatureType polygonFeature = polyTypeBuilder.buildFeatureType();

        createDirectory(outputDirectoryName);
        String fileNamePrefix= "benchmark-polygons-" + numberOfIndexPolygons;
        File shapefile = new File(outputDirectoryName, fileNamePrefix + ".shp");
        File h2DbFile = new File(outputDirectoryName, fileNamePrefix + "-h2-index");

        SimpleFeatureCollection featureCollection = null;
        switch (datastoreType) {
            case shapefile:
                if (shapefile.exists()) {
                    try {
                        featureCollection = getFeatureCollection(datastoreType,
                                FileDataStoreFinder.getDataStore(shapefile),
                                fileNamePrefix);
                    } catch (IOException e) {
                        throw new RuntimeException("Failed to read features from shapefile: " + shapefile, e);
                    }
                }
                break;
            case h2:
                if (h2DbFile.exists()) {
                    Map<String, Object> params = new HashMap<>();
                    params.put("dbtype", "h2");
                    params.put("database", h2DbFile.getAbsolutePath());

                    try {
                        DataStore datastore = DataStoreFinder.getDataStore(params);
                        featureCollection = getFeatureCollection(datastoreType, datastore,"h2");
                    } catch (IOException e) {
                        throw new RuntimeException("Failed to read features from h2: " + h2DbFile, e);
                    }
                }
                break;

        }

        if (featureCollection == null) {
            logger.info("Creating {} features", numberOfIndexPolygons);

            ArrayList<SimpleFeature> features = new ArrayList<>();
            int id = 0;
            List<double[][]> polygons = getIndexPolygons();
            try (ProgressBar progressBar = new ProgressBar("Features:", numberOfIndexPolygons)) {
                for (double[][] latlons : polygons) {
                    Coordinate[] coords = new Coordinate[latlons.length];
                    for (int i = 0; i < coords.length; i++) {
                        coords[i] = new Coordinate(latlons[i][1], latlons[i][0]);
                    }
                    Polygon polygon = gf.createPolygon(coords);
                    SimpleFeature feature = createSimpleFeature(polygonFeature, polygon);

                    if (feature != null) {
                        feature.setAttribute("id", ++id);
                        features.add(feature);
                    } else {
                        logger.error("Not a valid feature");
                    }
                    progressBar.step();
                }
            }

            logger.info("Wrapping features in a collection");
            featureCollection = DataUtilities.collection(features);

            logger.info("Writing features to shapefile: {}", datastoreType);
            switch (datastoreType) {
                case shapefile:
                    try {
                        writeToShapefile(shapefile, polygonFeature, features);
                    } catch (IOException e) {
                        logger.error("Failed to write shapefile", e);
                    }
                    break;
                case h2:
                default:
                    throw new UnsupportedOperationException();
            }
        }
        index = new SpatialIndexFeatureCollection(featureCollection.getSchema());
        index.addAll(featureCollection);
    }

    private SimpleFeatureCollection getFeatureCollection(DatastoreType datastoreType, DataStore store, String typeName)
            throws IOException {

        logger.info("Reading features from {}", datastoreType);
        SimpleFeatureSource featureSource = store.getFeatureSource(typeName);
        SimpleFeatureCollection collection = featureSource.getFeatures();
        int size = collection.size();
        if (size != numberOfIndexPolygons) {
            throw new RuntimeException("Unexpected number of features in " +
                    datastoreType + ": expected=" + numberOfIndexPolygons + ", read=" + size);
        }
//                showFeatures(featureSource);

        logger.info("Wrapping features in a collection");
        SimpleFeatureCollection featureCollection = DataUtilities.collection(featureSource.getFeatures());
        store.dispose();

        logger.info("Finished reading {} features from {}", size, datastoreType);
        return featureCollection;
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

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 1)
    public void pointIntersectsQuery() throws TransformException {

        GeodeticCalculator gc = new GeodeticCalculator(DefaultGeographicCRS.WGS84);

        long candidateCount = 0;
        long nearestCount = 0;
        for (double[] latlon : getQueryPoints()) {
            SimpleFeature nearestFeature = null;
            double nearestDistance = Double.POSITIVE_INFINITY;
            Coordinate coord = new Coordinate(latlon[1], latlon[0]);
            Point point = gf.createPoint(coord);
            // get all features that are within maxSearchDistance of the query
            SimpleFeatureCollection candidates = getCandidateFeatures(point);
            if (!candidates.isEmpty()) {
                candidateCount += candidates.size();
                nearestCount++;
                // iterate through the candidates
                try (SimpleFeatureIterator itr = candidates.features()) {
                    while (itr.hasNext()) {
                        SimpleFeature feature = itr.next();
                        Geometry featureGeometry = (Geometry) feature.getDefaultGeometry();
                        double distance = DistanceOp.distance(point, featureGeometry);
                        if (nearestDistance > distance) {
                            nearestDistance = distance;
                            nearestFeature = feature;
                        }
                    }
                }
            }
            if (nearestFeature != null) {
                gc.setStartingPosition(JTS.toDirectPosition(coord, DefaultGeographicCRS.WGS84));
                gc.setDestinationPosition(JTS.toDirectPosition(((Geometry) nearestFeature.getDefaultGeometry()).getCoordinate(), DefaultGeographicCRS.WGS84));
                double distance = gc.getOrthodromicDistance();
                results.add(new AbstractMap.SimpleImmutableEntry<>((int) nearestFeature.getAttribute("id"), distance));
            } else {
                results.add(new AbstractMap.SimpleImmutableEntry<>(0, -1.0));
            }
        }

        candidateCounts.add(nearestCount == 0 ? 0 : (candidateCount / nearestCount));
        nearestCounts.add(nearestCount);
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 1)
    public void polygonIntersectsQuery() throws TransformException {

        GeodeticCalculator gc = new GeodeticCalculator(DefaultGeographicCRS.WGS84);

        long candidateCount = 0;
        long nearestCount = 0;
        for (double[][] latlons : getQueryPolygons()) {
            SimpleFeature nearestFeature = null;
            double nearestDistance = Double.POSITIVE_INFINITY;

            Coordinate[] coords = new Coordinate[latlons.length];
            for (int i = 0; i < coords.length; i++) {
                coords[i] = new Coordinate(latlons[i][1], latlons[i][0]);
            }
            Polygon polygon = gf.createPolygon(coords);

            // get all features that are within maxSearchDistance of the query
            SimpleFeatureCollection candidates = getCandidateFeatures(polygon);
            if (!candidates.isEmpty()) {
                candidateCount += candidates.size();
                nearestCount++;
                // iterate through the candidates
                try (SimpleFeatureIterator itr = candidates.features()) {
                    while (itr.hasNext()) {
                        SimpleFeature feature = itr.next();
                        Geometry featureGeometry = (Geometry) feature.getDefaultGeometry();
                        double distance = DistanceOp.distance(polygon, featureGeometry);
                        if (nearestDistance > distance) {
                            nearestDistance = distance;
                            nearestFeature = feature;
                        }
                    }
                }
            }
            if (nearestFeature != null) {
                gc.setStartingPosition(JTS.toDirectPosition(polygon.getCoordinate(), DefaultGeographicCRS.WGS84));
                gc.setDestinationPosition(JTS.toDirectPosition(((Geometry) nearestFeature.getDefaultGeometry()).getCoordinate(), DefaultGeographicCRS.WGS84));
                double distance = gc.getOrthodromicDistance();
                results.add(new AbstractMap.SimpleImmutableEntry<>((int) nearestFeature.getAttribute("id"), distance));
            } else {
                results.add(new AbstractMap.SimpleImmutableEntry<>(0, -1.0));
            }
        }

        candidateCounts.add(nearestCount == 0 ? 0 : (candidateCount / nearestCount));
        nearestCounts.add(nearestCount);
    }

    @TearDown
    public void teardown() throws IOException {
        super.teardown();
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

    private SimpleFeatureCollection getCandidateFeatures(Geometry geometry) {
        SimpleFeatureType schema = index.getSchema();
        Filter filter = ff.intersects(ff.property(schema.getGeometryDescriptor().getName()), ff.literal(geometry));
        return index.subCollection(filter);
    }

    private void writeToShapefile(File shapefile, SimpleFeatureType type, ArrayList<SimpleFeature> features) throws IOException {
        ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();

        Map<String, Serializable> params = new HashMap<>();
        params.put("url", shapefile.toURI().toURL());
        params.put("create spatial index", Boolean.TRUE);

        ShapefileDataStore newDataStore =
                (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);

        /*
         * TYPE is used as a template to describe the file contents
         */
        newDataStore.createSchema(type);

        /*
         * Write the features to the shapefile
         */
        Transaction transaction = new DefaultTransaction("create");

        String[] typeNames = newDataStore.getTypeNames();
        String typeName = typeNames[0];
        SimpleFeatureSource featureSource = newDataStore.getFeatureSource(typeName);
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
    }

    public static void main(String[] args) throws IOException {
//        for (Iterator<DataStoreFactorySpi> i = DataStoreFinder.getAvailableDataStores(); i.hasNext(); ) {
//            DataStoreFactorySpi factory = i.next();
//            System.out.println(factory.getDisplayName());
//        }

        GeotoolsBenchmark benchmark = new GeotoolsBenchmark();
        benchmark.setup();
        try {
            benchmark.pointIntersectsQuery();
            benchmark.teardown();
            benchmark.polygonIntersectsQuery();
            benchmark.teardown();
        } catch (TransformException e) {
            throw new IOException(e);
        }
    }
}

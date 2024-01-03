package uk.ac.shef.wit.geo.benchmark;

import org.geotools.api.data.DataStore;
import org.geotools.api.feature.simple.SimpleFeature;
import org.geotools.api.feature.simple.SimpleFeatureType;
import org.geotools.api.filter.Filter;
import org.geotools.api.filter.FilterFactory;
import org.geotools.api.referencing.operation.TransformException;
import org.geotools.data.collection.SpatialIndexFeatureCollection;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.factory.CommonFactoryFinder;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.GeodeticCalculator;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.distance.DistanceOp;
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
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

@State(Scope.Thread)
public class GeotoolsBenchmark
        extends AbstractBenchmark {

    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    private SpatialIndexFeatureCollection index;

    @Setup
    public void setup() {
        Map.Entry<DataStore, SimpleFeatureCollection> dataStoreCollection = getIndexPolygons();
        DataStore dataStore = dataStoreCollection.getKey();
        SimpleFeatureCollection polygons = dataStoreCollection.getValue();

        logger.info("Creating spatial index");
        index = new SpatialIndexFeatureCollection(polygons.getSchema());
        Function<SimpleFeature, SimpleFeature> simplificationFunction =
                getSimplificationFunction(polygons);
        if (simplificationFunction == null) {
            index.addAll(polygons);
        } else {
            SimpleFeatureIterator it = polygons.features();
            while (it.hasNext()) {
                SimpleFeature feature = simplificationFunction.apply(it.next());
                index.add(feature);
            }
        }

        dataStore.dispose();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Fork(value = 1)
    @Warmup(iterations = 0)
    @Measurement(iterations = 3)
    public void pointIntersectsQuery() {
        long candidateCount = 0;
        long nearestCount = 0;
        List<double[]> queryPoints = getQueryPoints();
        for (double[] latlon : queryPoints) {
            Coordinate coord = new Coordinate(latlon[1], latlon[0]);
            Point point = gf.createPoint(coord);
            // get all features that are within maxSearchDistance of the query
            SimpleFeatureCollection candidates = getIntersectingFeatures(point);
            if (!candidates.isEmpty()) {
                candidateCount += candidates.size();
                nearestCount++;
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
    @Measurement(iterations = 3)
    public void polygonIntersectsQuery() throws TransformException {

        GeodeticCalculator gc = new GeodeticCalculator(DefaultGeographicCRS.WGS84);

        long candidateCount = 0;
        long nearestCount = 0;
        List<double[][]> queryPolygons = getQueryPolygons();
        for (double[][] latlons : queryPolygons) {
            SimpleFeature nearestFeature = null;
            double nearestDistance = Double.POSITIVE_INFINITY;

            Coordinate[] coords = new Coordinate[latlons.length];
            for (int i = 0; i < coords.length; i++) {
                coords[i] = new Coordinate(latlons[i][1], latlons[i][0]);
            }
            Polygon polygon = gf.createPolygon(coords);

            // get all features that are within maxSearchDistance of the query
            SimpleFeatureCollection candidates = getIntersectingFeatures(polygon);
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
                results.add(new AbstractMap.SimpleImmutableEntry<>(nearestFeature.getID(), distance));
            } else {
                results.add(new AbstractMap.SimpleImmutableEntry<>("0", -1.0));
            }
        }

        candidateCounts.add(nearestCount == 0 ? 0 : (candidateCount / nearestCount));
        nearestCounts.add(nearestCount);
    }

    @TearDown
    public void teardown() {
        super.teardown();
    }

    private static final FilterFactory ff = CommonFactoryFinder.getFilterFactory();

    private SimpleFeatureCollection getIntersectingFeatures(Geometry geometry) {
        SimpleFeatureType schema = index.getSchema();
        Filter filter = ff.intersects(ff.property(schema.getGeometryDescriptor().getName()), ff.literal(geometry));
        return getFeatures(filter);
    }

    private SimpleFeatureCollection getFeatures(Filter filter) {
        return index.subCollection(filter);
    }

    public static void main(String[] args) {
//        for (Iterator<DataStoreFactorySpi> i = DataStoreFinder.getAvailableDataStores(); i.hasNext(); ) {
//            DataStoreFactorySpi factory = i.next();
//            System.out.println(factory.getDisplayName());
//        }

        try (GeotoolsBenchmark benchmark = new GeotoolsBenchmark()) {
            benchmark.setup();
            benchmark.pointIntersectsQuery();
            benchmark.teardown();
            benchmark.polygonIntersectsQuery();
            benchmark.teardown();
        } catch (Exception e) {
            logger.error("", e);
        }
    }

    @Override
    public void close() {
        index.clear();
    }
}

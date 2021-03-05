package uk.ac.shef.wit.geo.benchmark;

import org.geotools.data.DataUtilities;
import org.geotools.data.collection.SpatialIndexFeatureCollection;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.factory.CommonFactoryFinder;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.geometry.jts.JTS;
import org.geotools.referencing.GeodeticCalculator;
import org.geotools.referencing.crs.DefaultGeographicCRS;
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

import java.lang.invoke.MethodHandles;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

@State(Scope.Thread)
public class GeotoolsBenchmark
        extends AbstractBenchmark {

    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    private static final FilterFactory2 ff = CommonFactoryFinder.getFilterFactory2();
    private SpatialIndexFeatureCollection index;

    @Setup
    public void setup() {
        SimpleFeatureTypeBuilder polyTypeBuilder = new SimpleFeatureTypeBuilder();
        polyTypeBuilder.setName("Polygon");
        polyTypeBuilder.setNamespaceURI("Polygon");
        polyTypeBuilder.setCRS(DefaultGeographicCRS.WGS84);
        polyTypeBuilder.add("polyGeom", Polygon.class);
        polyTypeBuilder.setDefaultGeometry("polyGeom");
        polyTypeBuilder.add("id", Integer.class);
        SimpleFeatureType polygonFeature = polyTypeBuilder.buildFeatureType();

        ArrayList<SimpleFeature> features = new ArrayList<>();
        int id = 0;
        for (double[][] latlons : getIndexPolygons()) {
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
        }

        SimpleFeatureCollection featureCollection = DataUtilities.collection(features);
        index = new SpatialIndexFeatureCollection(featureCollection.getSchema());
        index.addAll(features);
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
        results.clear();
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
        results.clear();
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
    public void teardown() {
        super.teardown();
    }

    public static SimpleFeature createSimpleFeature(SimpleFeatureType schema, Geometry geometry) {
        if (geometry != null && geometry.isValid()) {
            SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(schema);
            featureBuilder.add(geometry);
            return featureBuilder.buildFeature(null);
        }
        return null;
    }

    private SimpleFeatureCollection getCandidateFeatures(Geometry geometry) {
        SimpleFeatureType schema = index.getSchema();
        Filter filter = ff.intersects(ff.property(schema.getGeometryDescriptor().getName()), ff.literal(geometry));
        return index.subCollection(filter);
    }
}

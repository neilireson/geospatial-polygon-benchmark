package uk.ac.shef.wit.geo.benchmark;

import org.locationtech.jts.geom.GeometryFactory;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.results.RunResult;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import static java.lang.Math.*;

@State(Scope.Benchmark)
public abstract class AbstractBenchmark {

    static final GeometryFactory gf = new GeometryFactory();
    static final String outputDirectoryName = "out";
    private static final Random random = new Random();

    //        @Param({"10000", "100000"})
    int numberOfIndexPolygons = 1234;

    int numberOfQueryPoints = 10000;

    //    @Param({"1000", "10000", "100000"})
    int queryRadiusMetres = 1000;

    // roughly the bounding box for England
    int minLat = 50;
    int maxLat = 56;
    int minLon = -2;
    int maxLon = 2;

    // random polygon parameters
    int minVertices = 5;
    int maxVertices = 12;
    // the beta distribution determines the likelihood of polygon radius
    // 1,1 is uniform; 1,2 linear decrease; 2,1 linear increase; 5,5 approx normal pdf;
    // 8,2
//    AbstractRealDistribution radiusBetaPDF = new BetaDistribution(1.0, 1.0);
    float irregularity = 1f;
    float spikiness = 1f;
    // increasing lambda increases the likelihood of polygons with a minRadius
    // -log(1 - (1 - exp(-lambda)) * U) / lambda;
    // basically the exponential distribution (-log(1-U)/2) bounded [0,1]
    double lambda = 10;
    float minRadius = 400;
    float maxRadius = 4000;

    final List<Long> candidateCounts = new ArrayList<>();
    final List<Long> nearestCounts = new ArrayList<>();
    final List<Map.Entry<Integer, Double>> results = new ArrayList<>();

    synchronized List<double[][]> getIndexPolygons() {
        return getPolygons("benchmark-index-polygons", numberOfIndexPolygons);
    }

    synchronized List<double[]> getQueryPoints() {
        return getPoints("benchmark-query-points", numberOfQueryPoints);
    }

    synchronized List<double[][]> getQueryPolygons() {
        return getPolygons("benchmark-query-polygons", numberOfQueryPoints);
    }

    private List<double[]> getPoints(String prefix, int numberOfPoints) {
        createDirectory(outputDirectoryName);
        String filename = prefix + "-" + numberOfPoints + ".csv";
        Path path = Paths.get(outputDirectoryName, filename);
        List<double[]> indexPoints = new ArrayList<>();
        if (Files.exists(path)) {
            try (BufferedReader reader = new BufferedReader(new FileReader(path.toFile()))) {
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
            for (int i = 0; i < numberOfPoints; i++) {
                indexPoints.add(createRandomLatLon());
            }
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(path.toFile()))) {
                for (double[] latlon : indexPoints) {
                    writer.write(latlon[0] + "," + latlon[1]);
                    writer.newLine();
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        return indexPoints;
    }

    private List<double[][]> getPolygons(String prefix, int numberOfPolygons) {
        createDirectory(outputDirectoryName);
        String filename = prefix + "-" + numberOfPolygons + ".csv";
        Path path = Paths.get(outputDirectoryName, filename);
        List<double[][]> indexPolygons = new ArrayList<>();
        if (Files.exists(path)) {
            try (BufferedReader reader = new BufferedReader(new FileReader(path.toFile()))) {
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
            for (int i = 0; i < numberOfPolygons; i++) {
                indexPolygons.add(createPolygon());
            }
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(path.toFile()))) {
                for (double[][] latlons : indexPolygons) {
                    writer.write(latlons[0][1] + "," + latlons[0][0]);
                    for (int i = 1; i < latlons.length; i++) {
                        writer.write(" " + latlons[i][1] + "," + latlons[i][0]);
                    }
                    writer.newLine();
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        return indexPolygons;
    }

    void createDirectory(String directoryPath) {
        Path outputDirectory = Paths.get(directoryPath);
        if (!Files.exists(outputDirectory)) {
            if (!outputDirectory.toFile().mkdirs()) {
                throw new RuntimeException("Failed to create output directory: " + outputDirectory.toAbsolutePath());
            }
        } else if (!Files.isDirectory(outputDirectory)) {
            throw new RuntimeException("Output directory is not a directory: " + outputDirectory.toAbsolutePath());
        }
    }


    private double[] createRandomLatLon() {
        final double latitude = ThreadLocalRandom.current().nextDouble(minLat, maxLat);
        final double longitude = ThreadLocalRandom.current().nextDouble(minLon, maxLon);
        return new double[]{latitude, longitude};
    }

    double[][] createPolygon() {
        final double latitude = ThreadLocalRandom.current().nextDouble(minLat, maxLat);
        final double longitude = ThreadLocalRandom.current().nextDouble(minLon, maxLon);
        return createPolygon(latitude, longitude);
    }

    double[][] createPolygon(double lat, double lon) {
        return createPolygon((float) lat, (float) lon,
                minVertices + random.nextInt(maxVertices - minVertices),
                irregularity, spikiness, minRadius, maxRadius, lambda);
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

    protected void teardown() throws IOException {
        System.out.format("%n%s: average number of candidates per successful query = %.0f, " +
                        "number of intersecting polygons found = %.0f/%d%n",
                getClass().getSimpleName(),
                candidateCounts.stream().mapToLong(x -> x).average().orElse(0),
                nearestCounts.stream().mapToLong(x -> x).average().orElse(0),
                numberOfQueryPoints);

        // write results
        String filename = "results" +
                "-" + this.getClass().getSimpleName() +
                "-" + numberOfIndexPolygons +
                "-" + queryRadiusMetres +
                ".csv";
        Path path = Paths.get(outputDirectoryName, filename);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(path.toFile()))) {
            for (Map.Entry<Integer, Double> result : results) {
                writer.write(String.valueOf(result.getKey()));
                writer.write('\t');
                writer.write(String.valueOf(result.getValue()));
                writer.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static void main(String[] args) throws RunnerException, IOException {
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
//                .include(GeotoolsBenchmark.class.getSimpleName())
                .include(LuceneBenchmark.class.getSimpleName())
//                .include(MongoDbBenchmark.class.getSimpleName())
                .build();

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
    }
}
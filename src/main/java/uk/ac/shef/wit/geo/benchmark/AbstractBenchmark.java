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
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import static java.lang.Math.*;

@State(Scope.Benchmark)
public abstract class AbstractBenchmark {

    static final GeometryFactory gf = new GeometryFactory();
    private final Random random = new Random();

    //        @Param({"10000", "100000"})
    int numberOfIndexPolygons = 100000;

    int numberOfQueryPoints = 1000;

    //    @Param({"1000", "10000", "100000"})
    int queryRadiusMetres = 1000;

    // roughly the UK
    int minLat = 52;
    int maxLat = 53;
    int minLon = -1;
    int maxLon = 1;

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
        createDirectory("out");
        String filename = prefix + "-" + numberOfPoints + ".csv";
        Path path = Paths.get("out", filename);
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
        createDirectory("out");
        String filename = prefix + "-" + numberOfPolygons + ".csv";
        Path path = Paths.get("out", filename);
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

    private void createDirectory(String directoryPath) {
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
//        return createPolygon(lat, lon, 4 + random.nextInt(8), 0.0005, 0.005);
        return createPolygon(lat, lon, 5 + random.nextInt(7), 100, 0.9, 0.1);
    }

    /**
     * @param lat           The centre latitude of the polygon
     * @param lon           The centre longitude of the polygon
     * @param numberOfNodes Number of vertices
     * @param minRadius     Minimum radius degrees
     * @param maxRadius     Maximum radius degrees
     * @return polygon
     */
    double[][] createPolygon(double lat, double lon, int numberOfNodes, double minRadius, double maxRadius) {
        double[] angles = new double[numberOfNodes];
        for (int i = 0; i < numberOfNodes; i++) {
            angles[i] = random.nextDouble() * Math.PI * 2;
        }
        Arrays.sort(angles);

        double[][] latlons = new double[numberOfNodes + 1][];
        for (int i = 0; i < numberOfNodes; i++) {
            double radius = minRadius + random.nextDouble() * (maxRadius - minRadius);
            double coordLat = Math.toRadians(lat) + (Math.sin(angles[i]) * radius);
            double coordLon = Math.toRadians(lon) + (Math.cos(angles[i]) * radius);
            latlons[i] = new double[]{Math.toDegrees(coordLat), Math.toDegrees(coordLon)};
        }
        // close the ring
        latlons[numberOfNodes] = latlons[0];
        return latlons;
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
     * @param numVerts     number of vertices
     * @param aveRadius    in px, the average radius of this polygon, this roughly controls how large the polygon is, really only useful for order of magnitude.
     * @param irregularity [0,1] indicating how much variance there is in the angular spacing of vertices. [0,1] will map to [0, 2pi/numberOfVerts]
     * @param spikeyness   [0,1] indicating how much variance there is in each vertex from the circle of radius aveRadius. [0,1] will map to [0, aveRadius]
     * @return a list of vertices forming a polygon.
     */
    public static double[][] createPolygon(double latitude, double longitude,
                                           int numVerts, double aveRadius,
                                           double irregularity, double spikeyness) {
        double lat = toRadians(latitude);
        double lon = toRadians(longitude);

        // Radius of Earth at given latitude
        double latRadius = WGS84EarthRadius(lat);

        Random random = new Random();
        irregularity = clip(irregularity, 0, 1) * 2 * PI / numVerts;
        spikeyness = clip(spikeyness, 0, 1) * aveRadius;

        // generate n angle steps
        double[] angleSteps = new double[numVerts];
        double lower = (2 * PI / numVerts) - irregularity;
        double upper = (2 * PI / numVerts) + irregularity;
        double sum = 0;
        for (int i = 0; i < numVerts; i++) {
            double tmp = lower + ((upper - lower) * random.nextDouble());
            angleSteps[i] = tmp;
            sum = sum + tmp;
        }

        // normalize the steps so that point 0 and point n+1 are the same
        double k = sum / (2 * PI);
        for (int i = 0; i < numVerts; i++) {
            angleSteps[i] = angleSteps[i] / k;
        }
        // now generate the points
        double[][] points = new double[numVerts + 1][];
        double angle = 0;
        for (int i = 0; i < numVerts; i++) {
            angle += angleSteps[i];
            double d =
                    clip(random.nextGaussian() * spikeyness + aveRadius, 0, 2 * aveRadius)
                            / latRadius;
            double vertLat = lat + (sin(angle) * d);
            double vertLon = lon + (cos(angle) * d);
            points[i] = new double[]{toDegrees(vertLat), toDegrees(vertLon)};
        }
        points[numVerts] = points[0];

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

    protected void teardown() {
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
        Path path = Paths.get("out", filename);
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
                .include(GeotoolsBenchmark.class.getSimpleName())
                .include(LuceneBenchmark.class.getSimpleName())
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
package uk.ac.shef.wit.geo.benchmark;

import org.locationtech.jts.geom.Envelope;
import org.springframework.context.support.ClassPathXmlApplicationContext;
import org.springframework.context.support.FileSystemXmlApplicationContext;

import java.net.URL;

public class TrialConfiguration {

    public enum SimplificationType {
        BoundingBox,
        ConvexHull,
        ConcaveHull,
        DouglasPeucker,
        TopologyPreserving,
        None;
    }

    private String shapeFile = "/data/osm/gb/shp/wales/gis_osm_buildings_a_free_1.shp";
    private String typeName = "gis_osm_buildings_a_free_1";

    // The maximum difference between the number of features in the shapefile and index.
    // Used for cases where the shapefile contains geometries deemed invalid by the indexes to prevent re-indexing.
    private long missingDataThreshold = 0;

    // Need to specify either shapefile or number of polygons to generate and index
    private Integer numberOfIndexPolygons = null;
    private int numberOfQueryPoints = 10000;
    private int numberOfQueryPolygons = 10000;

    // Simplification
    private boolean removeHoles = false;
    private SimplificationType simplificationType = SimplificationType.None;
    /**
     * The (non-negative) distance tolerance for the simplification.
     * All vertices in the simplified geometry will be within this distance of the original geometry.
     * Therefore if the units are latitude/longitude degrees at the Equator the distances are approximately:
     * 1°       = 111 km  (60 nautical miles)
     * 0.1°     = 11.1 km
     * 0.01°    = 1.11 km
     * 0.001°   = 111 m
     * 0.0001°  = 11.1 m
     * 0.00001° = 1.11 m
     * As latitude moves away from the Equator the longitude distance values will decrease:
     * degrees      equator     lat=23N/S     lat=45N/S    lat=67N/S    lat=90N/S
     * 1°           111 km      102 km        78 km        43 km        0 km
     */
    private Double SimplificationThreshold = 0.0;

    // the boundingBox used to calculate the location of the random points and polygons
    // if the index is created from a shapefile and the boundingBox is null:
    // the boundingBox is calculated from the approximate convexhull of the features
    // Default: roughly the bounding box for England
    private Envelope boundingBox;// = new Envelope(-2, 2, 50, 56);

    // random polygon parameters
    private int minVertices = 5;
    private int maxVertices = 20;
    // beta distribution used to determine the likelihood of polygon radius
    // 1,1 is uniform; 1,2 linear decrease; 2,1 linear increase; 5,5 approx normal pdf;
    // 8,2
//    AbstractRealDistribution radiusBetaPDF = new BetaDistribution(1.0, 1.0);
    private float irregularity = 1f;
    private float spikiness = 0.5f;
    // exponential distribution used to determine the likelihood of polygon radius
    // increasing lambda (>0) increases the likelihood of polygons with a smaller radius
    // -log(1 - (1 - exp(-lambda)) * U) / lambda;
    // this is a bounded ([0,1]) exponential distribution (-log(1-U)/2)
    private double lambda = 10;
    private float minRadius = 400;
    private float maxRadius = 4000;

    public String getShapeFile() {
        return shapeFile;
    }

    public void setShapeFile(String shapeFile) {
        this.shapeFile = shapeFile;
    }

    public String getTypeName() {
        return typeName;
    }

    public void setTypeName(String typeName) {
        this.typeName = typeName;
    }

    public Integer getNumberOfIndexPolygons() {
        return numberOfIndexPolygons;
    }

    public void setNumberOfIndexPolygons(Integer numberOfIndexPolygons) {
        this.numberOfIndexPolygons = numberOfIndexPolygons;
    }

    public long getMissingDataThreshold() {
        return missingDataThreshold;
    }

    public void setMissingDataThreshold(long missingDataThreshold) {
        this.missingDataThreshold = missingDataThreshold;
    }

    public int getNumberOfQueryPoints() {
        return numberOfQueryPoints;
    }

    public void setNumberOfQueryPoints(int numberOfQueryPoints) {
        this.numberOfQueryPoints = numberOfQueryPoints;
    }

    public int getNumberOfQueryPolygons() {
        return numberOfQueryPolygons;
    }

    public void setNumberOfQueryPolygons(int numberOfIndexPolygons) {
        this.numberOfIndexPolygons = numberOfIndexPolygons;
    }

    public Envelope getBoundingBox() {
        return boundingBox;
    }

    public void setBoundingBox(Envelope boundingBox) {
        this.boundingBox = boundingBox;
    }

    public int getMinVertices() {
        return minVertices;
    }

    public void setMinVertices(int minVertices) {
        this.minVertices = minVertices;
    }

    public int getMaxVertices() {
        return maxVertices;
    }

    public void setMaxVertices(int maxVertices) {
        this.maxVertices = maxVertices;
    }

    public float getIrregularity() {
        return irregularity;
    }

    public void setIrregularity(float irregularity) {
        this.irregularity = irregularity;
    }

    public float getSpikiness() {
        return spikiness;
    }

    public void setSpikiness(float spikiness) {
        this.spikiness = spikiness;
    }

    public double getLambda() {
        return lambda;
    }

    public void setLambda(double lambda) {
        this.lambda = lambda;
    }

    public float getMinRadius() {
        return minRadius;
    }

    public void setMinRadius(float minRadius) {
        this.minRadius = minRadius;
    }

    public float getMaxRadius() {
        return maxRadius;
    }

    public void setMaxRadius(float maxRadius) {
        this.maxRadius = maxRadius;
    }

    public boolean getRemoveHoles() {
        return removeHoles;
    }

    public void setRemoveHoles(boolean removeHoles) {
        this.removeHoles = removeHoles;
    }

    public SimplificationType getSimplificationType() {
        return simplificationType;
    }

    public void setSimplificationType(SimplificationType simplificationType) {
        this.simplificationType = simplificationType;
    }

    public Double getSimplificationThreshold() {
        return SimplificationThreshold;
    }

    public void setSimplificationThreshold(Double simplificationThreshold) {
        SimplificationThreshold = simplificationThreshold;
    }

    public static TrialConfiguration create() {
        URL url = TrialConfiguration.class.getResource("/default-trial-configuration.xml");
        return (url != null) ? create(url.toString()) : new TrialConfiguration();
    }

    public static TrialConfiguration create(String path) {
        if (path.startsWith("classpath:")) {
            return new ClassPathXmlApplicationContext(path).getBean(TrialConfiguration.class);
        } else {
            return new FileSystemXmlApplicationContext(path).getBean(TrialConfiguration.class);
        }
    }
}

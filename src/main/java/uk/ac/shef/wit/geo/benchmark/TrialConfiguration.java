package uk.ac.shef.wit.geo.benchmark;

import org.locationtech.jts.geom.Envelope;
import org.springframework.context.support.ClassPathXmlApplicationContext;
import org.springframework.context.support.FileSystemXmlApplicationContext;

import java.net.URL;

public class TrialConfiguration {

    private String shapeFile = "/data/osm/gb/shp/wales/gis_osm_buildings_a_free_1.shp";
    private String typeName = "gis_osm_buildings_a_free_1";

    private Integer numberOfIndexPoints = null;
    private long missingDataThreshold = 0;
    private int numberOfQueryPoints = 10000;

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

    public Integer getNumberOfIndexPoints() {
        return numberOfIndexPoints;
    }

    public void setNumberOfIndexPoints(Integer numberOfIndexPoints) {
        this.numberOfIndexPoints = numberOfIndexPoints;
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

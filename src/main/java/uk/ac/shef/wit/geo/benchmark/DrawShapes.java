package uk.ac.shef.wit.geo.benchmark;

import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.geometry.jts.WKTReader2;
import org.geotools.util.factory.FactoryFinder;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.hull.ConcaveHull;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.simplify.TopologyPreservingSimplifier;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Path2D;
import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

public class DrawShapes extends JPanel {

    private static final Logger logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
    private static final ThreadLocal<WKTReader2> wktReader = ThreadLocal.withInitial(WKTReader2::new);

    private final Map<Color, List<Path2D.Double>> polygonList = new HashMap<>();
    private double minX;
    private double minY;
    private double rangeX;
    private double rangeY;

    private final int width = 600;
    private final int height = 600;

    public static void draw(Geometry geometry) {
        new DrawShapes(geometry);
    }

    private DrawShapes(Geometry geometry) {
        init(geometry);
        addGeometry(geometry);
    }

    private void init(Geometry geometry) {
        Geometry envelope = geometry.getEnvelope();
        if (geometry instanceof Point || geometry instanceof LineString) {
            logger.error("Geometry in not a 2d object");
            return;
        }
        Coordinate[] boundary = envelope.getCoordinates();
        double minX = boundary[0].x;
        double minY = boundary[0].y;
        double maxX = boundary[2].x;
        double maxY = boundary[2].y;
        init(minX, maxX, minY, maxY);
    }

    private void init(double minX, double maxX, double minY, double maxY) {
        this.minX = minX;
        this.minY = minY;
        rangeX = maxX - minX;
        rangeY = maxY - minY;

        JFrame frame = new JFrame();
        JScrollPane scroll = new JScrollPane(this);
        frame.setBounds(10, 10, 600, 600); //(10, 10, 1100, 900);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setContentPane(scroll);
        frame.setVisible(true);
    }

    public void addPolygons(List<double[][]> polygonLatlons, Color color) {
        List<Path2D.Double> polygons = polygonList.computeIfAbsent(color, list -> new ArrayList<>());
        logger.info("initialising polygons: {}", polygonLatlons.size());

        for (double[][] latlons : polygonLatlons) {
            Path2D.Double polygon = new Path2D.Double();
            for (int i = 0; i < latlons.length; i++) {
                double x = width * (latlons[i][1] - minX) / rangeX;
                double y = height * (latlons[i][0] - minY) / rangeY;
                if (i == 0)
                    polygon.moveTo(x, y);
                else
                    polygon.lineTo(x, y);
            }
            polygon.closePath();
            polygons.add(polygon);
        }
    }

    public void addGeometry(Geometry geometry) {
        addGeometry(geometry, Color.BLACK);
    }

    public void addGeometry(Geometry geometry, Color color) {
        if (geometry.getNumGeometries() > 1) {
            for (int i = 0; i < geometry.getNumGeometries(); i++) {
                addGeometry(geometry.getGeometryN(i));
            }
            return;
        }

        List<Path2D.Double> paths = polygonList.computeIfAbsent(color, list -> new ArrayList<>());

        Path2D.Double path = new Path2D.Double();
        Coordinate[] coordinates = geometry.getCoordinates();
        for (int i = 0; i < coordinates.length; i++) {
            double x = width * (coordinates[i].x - minX) / rangeX;
            double y = height * (coordinates[i].y - minY) / rangeY;
            if (i == 0)
                path.moveTo(x, y);
            else
                path.lineTo(x, y);
        }
        path.closePath();
        paths.add(path);
    }

    @Override
    public Dimension getPreferredSize() {
        return new Dimension(width, height);
    }

    @Override
    public void paintComponent(Graphics g1) {
        super.paintComponent(g1);

        Graphics2D g = (Graphics2D) g1;
        this.setOpaque(true);
        this.setBackground(Color.WHITE);

        for (Map.Entry<Color, List<Path2D.Double>> entry : polygonList.entrySet()) {
            Color color = entry.getKey();
            if (Color.WHITE.equals(color)) {
                throw new IllegalArgumentException("Polygon color " + color + " is same as background");
            }
            g.setColor(entry.getKey());
            for (Path2D.Double polygon : entry.getValue()) {
                g.draw(polygon);
            }
        }
    }


    public static void main(String[] args) throws ParseException {
        String wkt = "MULTIPOLYGON (((-1.3123997 53.3445653, -1.3112742 53.3447766, -1.3102408 53.3450808, -1.3102832 53.345309, -1.3094834 53.3454189, -1.30891 53.345571, -1.3087755 53.345571, -1.3079403 53.3453555, -1.3067582 53.3448695, -1.3066852 53.3446706, -1.3066478 53.3446656, -1.3066041 53.3446598, -1.3066874 53.3449625, -1.3071192 53.3460189, -1.3050806 53.346133, -1.3051514 53.3468683, -1.3049815 53.3475021, -1.3062769 53.3476627, -1.3074236 53.3472105, -1.3072537 53.346002, -1.3072962 53.3458795, -1.3078129 53.3459429, -1.3080182 53.3471514, -1.3079261 53.3473289, -1.3072254 53.3475782, -1.3065761 53.3484642, -1.3065461 53.3485051, -1.3056648 53.3484987, -1.3034714 53.3483922, -1.3031857 53.348303, -1.3028176 53.3482385, -1.3023765 53.3481735, -1.3019692 53.3480892, -1.3044904 53.3482685, -1.3049303 53.3474616, -1.3050591 53.3468595, -1.3046406 53.3439516, -1.3028165 53.3440318, -1.3007899 53.3440809, -1.3006507 53.3439055, -1.3027094 53.3438257, -1.301775 53.3418015, -1.3044947 53.3419902, -1.3050457 53.3420285, -1.3052003 53.3421385, -1.3026527 53.3420508, -1.3023767 53.3424016, -1.3029076 53.3434665, -1.3041746 53.3438765, -1.3052788 53.343251, -1.3056186 53.3433313, -1.3048116 53.3443497, -1.3058026 53.3445526, -1.308995 53.3449794, -1.3097878 53.3437497, -1.3096281 53.3433802, -1.3096015 53.3433188, -1.3097192 53.3433277, -1.3098642 53.3433678, -1.3108708 53.3433989, -1.3111327 53.3431919, -1.3101912 53.3422072, -1.3092852 53.3418395, -1.3070697 53.3420339, -1.307105 53.3422495, -1.3075864 53.3427439, -1.3093347 53.3427228, -1.3094857 53.3430557, -1.3077067 53.3432383, -1.3061019 53.3430768, -1.3057175 53.3430381, -1.3055832 53.3424692, -1.3052021 53.3421349, -1.3050481 53.3420287, -1.3093772 53.3416705, -1.3097807 53.3417297, -1.3113521 53.3427439, -1.311407 53.3434852, -1.3114229 53.3436271, -1.3119042 53.3443878, -1.3123997 53.3445653)))";
        Geometry geometry = wktReader.get().read(wkt);
        DrawShapes drawShapes = new DrawShapes(geometry);

        TopologyPreservingSimplifier simplifier = new TopologyPreservingSimplifier(geometry.convexHull());
        simplifier.setDistanceTolerance(0.00);
        Geometry simplifiedGeometry = simplifier.getResultGeometry();
        drawShapes.addGeometry(simplifiedGeometry, Color.RED);

        Geometry concaveHull = ConcaveHull.compute(geometry);
        drawShapes.addGeometry(concaveHull, Color.BLUE);
    }
}




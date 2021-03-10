package uk.ac.shef.wit.geo.benchmark;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
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
    private final Map<Color, List<Path2D.Double>> polygonList = new HashMap<>();
    private double minX;
    private double minY;
    private double rangeX;
    private double rangeY;

    private final int width = 600;
    private final int height = 600;

    public DrawShapes(double minX, double maxX, double minY, double maxY) {
        init(minX, maxX, minY, maxY);
    }

    private DrawShapes(Geometry geometry) {
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
        addGeometry(geometry);
    }

    public static void draw(Geometry geometry) {
        new DrawShapes(geometry);
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
        if (geometry.getNumGeometries() > 1) {
            for (int i = 0; i < geometry.getNumGeometries(); i++) {
                addGeometry(geometry.getGeometryN(i));
            }
            return;
        }
        Color color = Color.BLACK;
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
}




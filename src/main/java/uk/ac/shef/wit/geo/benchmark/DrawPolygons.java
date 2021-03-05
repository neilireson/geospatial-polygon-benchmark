package uk.ac.shef.wit.geo.benchmark;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Path2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

public class DrawPolygons extends JPanel {

    private final Map<Color, List<Path2D.Double>> polygonList = new HashMap<>();
    private final double minX;
    private final double minY;
    private final double rangeX;
    private final double rangeY;

    private final int width = 600;
    private final int height = 600;

    public DrawPolygons(double minX, double maxX, double minY, double maxY) {
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
        System.out.println("initialising polygons: " + polygonLatlons.size());

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




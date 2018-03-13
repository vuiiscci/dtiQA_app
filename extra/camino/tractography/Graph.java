package tractography;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;


public class Graph {
		
	private final int MAXTHICKNESS = 10;
	private final int MAXRADIUS = 20;
	private final int FRAMEDIMENSION = 600;
	ArrayList<Edge> computedEdges;
	HashMap<Integer, Integer> verticesX;
	HashMap<Integer, Integer> verticesY;
	ArrayList<Integer> roiIntensities;
	HashMap<Integer, String> indexLabelMap;
	HashMap<Integer, Integer> nodeStrength;
	int threshold = 0;
	
	public Graph(ArrayList<Integer> aRoiIntensities,
			HashMap<Integer, Integer> aVerticesX,
			HashMap<Integer, Integer> aVerticesY,
			ArrayList<Edge> aComputedEdges,
			HashMap<Integer, String> aIndexLabelMap,
			HashMap<Integer, Integer> aNodeStrength,
			int aThreshold)
	{
		roiIntensities = aRoiIntensities;
		verticesX = aVerticesX;
		verticesY = aVerticesY;
		computedEdges = aComputedEdges;
		indexLabelMap = aIndexLabelMap;
		nodeStrength = aNodeStrength;
		threshold = aThreshold;
	}
	
	public void drawGraph()
	{
		BufferedImage image = new BufferedImage(FRAMEDIMENSION,FRAMEDIMENSION,BufferedImage.TYPE_INT_RGB);
		
		Graphics graphics = image.getGraphics();
		Graphics2D g2 = (Graphics2D) graphics;
		g2.setBackground(Color.WHITE);
		g2.setColor(Color.WHITE);
		g2.fillRect(0, 0, FRAMEDIMENSION, FRAMEDIMENSION);
		
    	double temp;
    	double min=0,max=1.0f;
    	if (computedEdges.size() > 0)
    	{
    		min = max = computedEdges.get(0).getWeight();
    	}
    	
    	for (int i=1; i < computedEdges.size(); i++)
		{
    		temp = computedEdges.get(i).getWeight();
    		if (min > temp)
    			min = temp;
    		if (max < temp)
    			max = temp;
		}
    	//System.err.println("Minimum : " + min + ",  "+ max);
    	
    	for (int i=0; i < computedEdges.size(); i++)
		{	
    		int x1= verticesX.get(computedEdges.get(i).getStartNode());
    		int y1= verticesY.get(computedEdges.get(i).getStartNode()); 
    		int x2= verticesX.get(computedEdges.get(i).getEndNode());
    		int y2= verticesY.get(computedEdges.get(i).getEndNode()); 
    		float thickness = (float)((computedEdges.get(i)
				.getWeight()- min) /(max - min) * MAXTHICKNESS);
    		if (computedEdges.get(i).getCount() > threshold)
    		{
    			g2.setStroke(new BasicStroke(thickness));	
    			g2.setPaint(Color.BLUE);
        		g2.drawLine(x1, y1, x2, y2);
    		}
		} 
    	
    	
    	if (nodeStrength.size() > 0)
    	{
    		min = max = nodeStrength.get(roiIntensities.get(0));
    	}
    	
    	for (int i=1; i < nodeStrength.size(); i++)
		{
    		temp = nodeStrength.get(roiIntensities.get(i));
    		if (min > temp)
    			min = temp;
    		if (max < temp)
    			max = temp;
		}
    	//System.err.println("nodeStrength size:" + nodeStrength.size());
    	
 		for (int i=0; i < verticesX.size(); i++)
     	{
 			if (verticesX.get(roiIntensities.get(i)) != null)
 			{
 				int x = verticesX.get(roiIntensities.get(i));
 				int y = verticesY.get(roiIntensities.get(i));
 				
                                // See message from Ian Malone, who found a problem with this line.
 				int thickness = (int)Math.round((nodeStrength.get(roiIntensities.get(i))- min)
 						/(max - min) * MAXRADIUS);
 				g2.drawOval(x - thickness/2, y - thickness/2, thickness, thickness);    			
 				
 	    		g2.setPaint(Color.RED);
 				g2.fillOval(x - thickness/2, y - thickness/2, thickness, thickness);
 				
 				if ((indexLabelMap.size() > 0) && 
 						(indexLabelMap.get(roiIntensities.get(i))) != null)
 				{	
 					Font font = new Font("Arial", Font.BOLD, 11);
 					g2.setFont(font);
 					g2.setPaint(Color.BLACK);
 					g2.drawString(indexLabelMap.get(roiIntensities.get(i)), x+thickness/2+2, y+thickness/2+2);
 				}
 			}
 			
     	}
 		try {
 			System.err.println("Written");
			ImageIO.write( image, "PNG", new File( "Connectivity.png" ));
		} 
 		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
 		
 		JLabel jlabel=new JLabel(new ImageIcon(image));
 		JFrame frame = new JFrame("Connectivity Mapping");
    	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	frame.getContentPane().add(jlabel);
    	frame.pack();
    	frame.setVisible(true);
	}
}
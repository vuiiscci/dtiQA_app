package tractography;

import imaging.ImageHeader;
import misc.*;
import numerics.*;

import java.util.ArrayList;
import java.io.*;


/**
 * An ROI for tractography. The ROI consists of a list of points in mm space. There is no restriction
 * on the number of points per voxel.
 * 
 *
 * @author Philip Cook
 * @version $Id$
 */
public final class PointListROI implements RegionOfInterest {

    private final Point3D[] points;
    
    private final int label;


    /**
     * Create an ROI from a list of points in physical space, associated with a specific label.
     *  
     */
    public PointListROI(Point3D[] pointList, int lab) {    
        
        points = new Point3D[pointList.length];
        
	System.arraycopy(pointList, 0, points, 0, points.length);

        label = lab;
    }


    /**
     * For testing. Assumes identity rotation to physical space. 
     * Use VoxelROI get ROIs for all labels in a volume from an image file.
     *
     */
    protected PointListROI(int[][][] seedVol, int index, double xVoxelDim, double yVoxelDim, double zVoxelDim) {
        	
        int xDataDim = seedVol.length;
        int yDataDim = seedVol[0].length;
        int zDataDim = seedVol[0][0].length;

        ArrayList<Point3D> pointList = new ArrayList<Point3D>(100);
        
        for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
                    
                    if (seedVol[i][j][k] == index) {

                        // Camino space is indexed from voxel corner, hence the offset
                        Point3D seed = new Point3D((i + 0.5) * xVoxelDim, (j + 0.5) * yVoxelDim, (k + 0.5) * zVoxelDim);

                        pointList.add(seed);
                    }
                    
                }
	    }
	}

        points = pointList.toArray(new Point3D[] {});

        label = index;

    }


    /**
     * Create an ROI from a list of points, in physical space, stored as text in the format
     * <code>x y z\nx y z\n...</code>. 
     *  
     */
    public static PointListROI readPoints(String file, ImageHeader diffusionSpace) {

	
	Point3D[] pointList = null;

	ArrayList<Point3D> tmpPoints = new ArrayList<Point3D>(100);
	
	try {

	    InputStreamReader isr = new InputStreamReader(new FileInputStream(file));
	    
	    BufferedReader reader = new BufferedReader(isr);
	    
	    String line = reader.readLine();

	    readLines: 
	    while (line != null) {

		String[] elements = line.split("\\s+");

		if (elements.length < 3) {
		    //blank line
		    break readLines;
		}

		Point3D point = new Point3D(Double.parseDouble(elements[0]), 
					    Double.parseDouble(elements[1]), 
					    Double.parseDouble(elements[2]));

		tmpPoints.add(point);

		line = reader.readLine();
	    }

	    pointList = new Point3D[tmpPoints.size()];
	
	    tmpPoints.toArray(pointList);

	    isr.close();
	}
	catch (IOException e) {
	    throw new LoggedException(e);
	}

        RealMatrix physToVox = diffusionSpace.getPhysicalToVoxelTransform();

        double[] voxelDims = diffusionSpace.getVoxelDims();
        
        for (int i = 0; i < pointList.length; i++) {

            // point in voxel space
            Point3D tmp = pointList[i].transform(physToVox);
            
            // point in Camino space
            pointList[i] = new Point3D( (tmp.x + 0.5) * voxelDims[0] , (tmp.y + 0.5) * voxelDims[1], (tmp.z + 0.5) * voxelDims[2]);
        }
        
	return new PointListROI(pointList, 1);
	
    }




    /**
     * @return the list of seed points
     */
    public Point3D[] getSeedPoints() {

        Point3D[] defCopy = new Point3D[points.length];
        
	System.arraycopy(points, 0, defCopy, 0, points.length);

        return defCopy;

    }


    /**
     *
     * @return the voxel indices corresponding to the points in the ROI, in order
     *
     */
    public Voxel[] getSeedVoxels(double xVoxelDim, double yVoxelDim, double zVoxelDim) {
        Voxel[] voxels = new Voxel[points.length];

        for (int i = 0; i < voxels.length; i++) {
            voxels[i] = new Voxel( (int)(points[i].x / xVoxelDim), (int)(points[i].y / yVoxelDim),  (int)(points[i].z / zVoxelDim));
        }

        return voxels;
    }


   
    public int getRegionLabel() {
        return label;
    }

}

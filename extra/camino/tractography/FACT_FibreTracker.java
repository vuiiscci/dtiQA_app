package tractography;

import numerics.*;

import java.util.Arrays;

/**
 *
 * Implements the FACT algorithm. No interpolation is applied.
 *
 * @version $Id: FACT_FibreTracker.java 1155 2012-12-11 21:43:30Z ucacpco $
 * @author  Philip Cook
 * 
 */
public class FACT_FibreTracker extends FibreTracker {

    private final double[] planeIntersections = new double[6];


    /** 
     * @param data the dataset within which the tracking will take place.
     */
    public FACT_FibreTracker(TractographyImage data) {

	super(data);
    }


    protected final Vector3D getFirstStep(Point3D seedPos, int pdIndex, boolean direction) {
	
	Point3D currentPos = seedPos;

	int i,j,k;

	i = (int)(seedPos.x / xVoxelDim);
	j = (int)(seedPos.y / yVoxelDim);
	k = (int)(seedPos.z / zVoxelDim);

	Vector3D[] pds = image.getPDs(i,j,k);

	Vector3D trackingDirection = pds[pdIndex];

        if (!direction) {
	    return getNextStep(currentPos, pds[pdIndex].negated());
	}
        
        return getNextStep(currentPos, pds[pdIndex]);

    }

    
    protected final Vector3D getNextStep(Point3D currentPos, Vector3D previousDirection) {


	int i,j,k;
        
	i = (int)(currentPos.x / xVoxelDim);
	j = (int)(currentPos.y / yVoxelDim);
	k = (int)(currentPos.z / zVoxelDim);
        
	Vector3D[] voxelPDs = image.getPDs(i,j,k);
	
	double maxProd = 0.0;
	int pdIndex = 0;
        
	if (voxelPDs.length > 1) {
            
	    for (int p = 0; p < voxelPDs.length; p++) {
		double dotProd = Math.abs(voxelPDs[p].dot(previousDirection));
		if (dotProd > maxProd) {
		    maxProd = dotProd;
		    pdIndex = p;
		}
	    }
	}
	
	Vector3D trackingDirection = voxelPDs[pdIndex];
        
	if ( trackingDirection.dot(previousDirection) < 0.0) {
	    // move along -e1
	    trackingDirection = trackingDirection.negated();
	}
        
        
        double frontPlane = yVoxelDim * (1 + j);
        double backPlane = yVoxelDim * j;
        
        double topPlane = zVoxelDim * (1 + k);
        double bottomPlane = zVoxelDim * k;
        
        double rightPlane = xVoxelDim * (1 + i);
        double leftPlane = xVoxelDim * i;
        
        if (trackingDirection.y != 0.0) {
            
            planeIntersections[0] = (frontPlane - currentPos.y) / trackingDirection.y;
            planeIntersections[1] = (backPlane - currentPos.y) / trackingDirection.y;

        }
        else {
            // meets at infinity...
            planeIntersections[0] = Double.MAX_VALUE;
            planeIntersections[1] = Double.MAX_VALUE;
        }		
        
        if (trackingDirection.z != 0.0) {
            planeIntersections[2] = (topPlane - currentPos.z) / trackingDirection.z;
            planeIntersections[3] = (bottomPlane - currentPos.z) / trackingDirection.z;
            
        }
        else {
            planeIntersections[2] = Double.MAX_VALUE;
            planeIntersections[3] = Double.MAX_VALUE;
        }		
        
        if (trackingDirection.x != 0.0 ) {
            planeIntersections[4] = (leftPlane - currentPos.x) / trackingDirection.x;
            planeIntersections[5] = (rightPlane - currentPos.x) / trackingDirection.x;
        }
        else {
            planeIntersections[4] = Double.MAX_VALUE;
            planeIntersections[5] = Double.MAX_VALUE;
        }
	
        // now sort
        Arrays.sort(planeIntersections);
        
        // scale by smallest positive element
        // allow value of 0.0; occurs only if point sits precisely on a voxel boundary
        int planeCounter = 0;
        while (planeIntersections[planeCounter] < 0.0) {
            planeCounter++;
            
            // this could happen if point was on boundary between voxels and planeIntersection had 
            // to be > 0. No longer a problem
            if (planeCounter == 6) {
                throw new IllegalStateException("No plane intersects tracking vector.\nPosition " + currentPos + "\nTracking vector " + 
                                                trackingDirection);
            }
        }
        
        double displacement = planeIntersections[planeCounter] + 0.001;

        // hack for approximate backward compatibility

        // Alternative solution: would be a global minimum step size, but need to do that carefully so as to not 
        // mess up tracking in high-resolution data

        Point3D newPos = currentPos.displace(trackingDirection.scaled(displacement));

        int a = (int)(newPos.x / xVoxelDim);
        int b = (int)(newPos.y / yVoxelDim);
        int c = (int)(newPos.z / zVoxelDim);
        
        if (a == i && b == j && c == k) {
            displacement += 0.001;
        }
        else if (inBounds(newPos) && visitedVoxel[a][b][c] > 0) {
            displacement = displacement < 0.1 ? 0.1 : displacement;
        }
        
        return trackingDirection.scaled(displacement);
        
    }




}

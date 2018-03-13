package tractography;

import numerics.*;

/**
 * Tracks using the fourth-order Runge-Kutta method. Tracts will use the the resultant
 * vector of the four sub-steps, rather than taking four independent steps. This keeps
 * the number of points in a streamline from getting too large.
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class RK4FibreTracker extends FibreTracker {

    private final double stepSize;
   
    
    protected ImageInterpolator interpolator;


    /** 
     *
     * @param interp the interpolator that provides the tracking direction at any point in the image.
     * @param stepLength the distance (in mm) to move before performing the next interpolation. 
     *
     */
    public RK4FibreTracker(ImageInterpolator interp, double stepLength) {
	super(interp.getImage());
	stepSize = stepLength;

        interpolator = interp;

	if (stepSize <= 0.0) {
	    throw new IllegalArgumentException("Can't track with step size " + stepSize);
	}
    }


    /**
     * @return the distance in mm between points of Tracts from this tracker.
     */ 
    public double stepSize() {
	return stepSize;
    }

    protected Vector3D[] getPDs(int i, int j, int k) {
	throw new UnsupportedOperationException("Interpolated trackers should take PDs from interpolator");
    }


    protected final Vector3D getFirstStep(Point3D seedPos, int pdIndex, boolean direction) {
        
        Vector3D initDirection = interpolator.getTrackingDirection(seedPos, pdIndex, direction);

        return getNextStep(seedPos, initDirection);
    }
    

    protected final Vector3D getNextStep(Point3D currentPos, Vector3D previousDirection) {
            
	Vector3D v1, v2, v3, v4;
	
        Point3D v2Pos, v3Pos, v4Pos;

        double[] stepDirection = new double[3];
        
        v1 = interpolator.getTrackingDirection(currentPos, previousDirection);
        
        // half steps along v1 and v2
        v2Pos = currentPos.displace(v1.scaled(stepSize * 0.5));

        if (!inBounds(v2Pos) || isBackground(v2Pos) ) {
            return v1.scaled(stepSize);
        }
        
        v2 = interpolator.getTrackingDirection(v2Pos, v1);
        
        v3Pos = currentPos.displace(v2.scaled(stepSize * 0.5));

        if (!inBounds(v3Pos) || isBackground(v3Pos) ) {
            return v1.scaled(stepSize);
        }
        
        v3 = interpolator.getTrackingDirection(v3Pos, v1);
        
        // take full step along v3 and v4
        v4Pos = currentPos.displace(v3.scaled(stepSize));
        
        if (!inBounds(v4Pos) || isBackground(v4Pos) ) {
            return v1.scaled(stepSize);
        }

        v4 = interpolator.getTrackingDirection(v4Pos, v1);
        
        
        stepDirection[0] = (1.0 / 6.0) * ( v1.x + 2.0 * v2.x + 2.0 * v3.x + v4.x );
        stepDirection[1] = (1.0 / 6.0) * ( v1.y + 2.0 * v2.y + 2.0 * v3.y + v4.y );
        stepDirection[2] = (1.0 / 6.0) * ( v1.z + 2.0 * v2.z + 2.0 * v3.z + v4.z );
        
        Vector3D trackingDirection = new Vector3D(stepDirection[0], stepDirection[1], stepDirection[2]).normalized();
        
        return trackingDirection.scaled(stepSize);
    }



    public final boolean isBackground(Point3D currentPos) {
        
	int i = (int)(currentPos.x / xVoxelDim);
	int j = (int)(currentPos.y / yVoxelDim);
	int k = (int)(currentPos.z / zVoxelDim);
        
        return isotropic[i][j][k];
    }

    
}

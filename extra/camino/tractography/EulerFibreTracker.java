package tractography;

import numerics.*;

/**
 * Uses Euler's method to track through the data, ie a fixed step size taking some interpolated
 * value at each point. The user assigns the interpolation method.
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class EulerFibreTracker extends FibreTracker {

    private final double stepSize;
   
    protected ImageInterpolator interpolator;


    /**
     *
     * @param interp the interpolator that provides the tracking direction at any point in the image.
     * @param stepLength the distance (in mm) to move before performing the next interpolation. 
     *
     */
    public EulerFibreTracker(ImageInterpolator interp, double stepLength) {

	super(interp.getImage());

	stepSize = stepLength;

        interpolator = interp;
        
	if (stepSize <= 0.0) {
	    throw new IllegalArgumentException("Can't track with step size " + stepSize);
	}
    }


    /**
     *
     * @return the distance in mm between points of Tracts from this tracker.
     *
     */ 
    public double stepSize() {
	return stepSize;
    }

    protected final Vector3D getFirstStep(Point3D seedPos, int pdIndex, boolean direction) {
        
        return interpolator.getTrackingDirection(seedPos, pdIndex, direction).scaled(stepSize);
    }


    protected final Vector3D getNextStep(Point3D currentPos, Vector3D previousDirection) {
        
        return interpolator.getTrackingDirection(currentPos, previousDirection).scaled(stepSize);
        
    }
    
}

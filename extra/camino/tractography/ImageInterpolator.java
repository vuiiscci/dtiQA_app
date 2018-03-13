package tractography;

import numerics.*;

/**
 * 
 * Interface for interpolators
 *
 * @author Philip Cook
 * @version $Id$
 */
public interface ImageInterpolator {


    /** 
     * Gets tracking direction at some point. The interpolator decides the direction 
     * to follow, and returns it.
     * 
     * @param point the point in mm to interpolate at. 
     * @param previousDirection the direction of the previous tracking step.
     * @return the tracking direction for this point.
     * 
     */
    public Vector3D getTrackingDirection(Point3D point, Vector3D previousDirection);


    /** 
     * Gets the initial tracking direction, given a pdIndex and a seed point. Implementing classes should
     * ensure that specifying the same point and pdIndex gets you the same vector, which differs only in sign
     * depending on the boolean parameter. 
     *
     * @param point the seed point in mm
     * @param pdIndex the numerical index of the PD to follow, from 0 to (number of PDs) - 1.
     * @param direction if true, the direction will be the PD, if false, it will be the negated PD.
     * @return the tracking direction for this point. 
     * 
     */
    public Vector3D getTrackingDirection(Point3D point, int pdIndex, boolean direction);



    /**
     * @return the image for this interpolator.
     *
     */
    public TractographyImage getImage();

}

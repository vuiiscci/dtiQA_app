package tractography;

import numerics.Point3D;

/**
 * Interface for a tractography region of interest. The ROI is in Camino space, ie a grid defined by the voxel
 * space, with origin at the corner of voxel (0,0,0) and distances measured in mm from that point, ie
 * (right, anterior, superior) is a positive vector.
 *
 * @author Philip Cook
 * @version $Id$
 */
public interface RegionOfInterest {

    /**
     * @return All the seed points for this region.
     */
    public Point3D[] getSeedPoints();

    
    /**
     *
     * @return the voxel indices corresponding to the points in the ROI, in order
     *
     */
    public Voxel[] getSeedVoxels(double xVoxelDim, double yVoxelDim, double zVoxelDim);


    /**
     * Get the integer label for this ROI
     *
     */
    public int getRegionLabel();


}

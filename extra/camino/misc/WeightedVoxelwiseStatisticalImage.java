package misc;

import numerics.*;

/**
 * Interface for images that compute statistics on the values in each voxel.
 * 
 * @author Philip Cook
 * @version $Id$
 *
 */
public interface WeightedVoxelwiseStatisticalImage extends VoxelwiseStatisticalImage {


 
    /**
     * Add a value to the voxel (i,j,k), with an associated weight.
     *
     */
    public void addValue(int i, int j, int k, double v, double weight);


    /**
     * Add an image of values. 
     * 
     * @param values must have the same dimensions as the image.
     * @param weights must have the same dimensions as the image.
     */
    public void addValues(double[][][] values, double[][][] weights);


}

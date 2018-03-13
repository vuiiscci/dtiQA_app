package misc;

import numerics.*;

/**
 * Interface for images that compute statistics on the values in each voxel.
 * 
 * @author Philip Cook
 * @version $Id$
 *
 */
public interface VoxelwiseStatisticalImage {


    /**
     * Add a value to the voxel (i,j,k).
     *
     */
    public void addValue(int i, int j, int k, double v);


    /**
     * Computes a 3D image, where each voxel intensity is some statistic of the values in the voxel.
     *
     * @param stat a statistic supported by the image.
     */
    public double[][][] getVoxelStatistic(String stat);


    /**
     * @return the data dimensions of the image.
     */
    public int[] getDataDims();

    /**
     * @return the voxel dimensions of the image.
     */
    public double[] getVoxelDims();


    

}

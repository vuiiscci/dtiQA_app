package tractography;

import numerics.*;


/**
 * Interface for tractography images
 *
 * @version $Id$
 * @author Philip Cook
 * 
 */
public interface TractographyImage {

    /**
     *
     * @return the principal direction or directions in the voxel. Probabilistic images will 
     * return a sample principal direction, successive calls to this method will return
     * the same direction until <code>resetRandomization</code> is called.
     *
     */
    public Vector3D[] getPDs(int i, int j, int k);


    /**
     * @return the princpal direction or directions in the voxel, given an estimated local fiber orientation.
     * If the PDs in the voxel are not dependent on the fiber orientation, then this method returns the
     * same thing as getPDs(int, int, int).
     *
     * @return a sample fibre orientation from this voxel. Note that this will <b>not</b> correct the orientation 
     * for agreement with the previous tracking direction so there is no guarantee that the returned vector 
     * will have a positive dot product with <code>fibreOrientation</code>.
     *
     */
    public Vector3D[] getPDs(int i, int j, int k, Vector3D fibreOrientation);


    /**
     * @return the number of PDs in this voxel.
     *
     */
    public int numberOfPDs(int i, int j, int k);
  

    /** Size of image in x dimension. */
    public int xDataDim();


    /** Size (mm) of x voxel length. */
    public double xVoxelDim();


    /** Size of image in y dimension. */
    public int yDataDim();


    /** Size (mm) of y voxel length. */
    public double yVoxelDim();


    /** Size of image in z dimension. */
    public int zDataDim();


    /** Size (mm) of z voxel length. */
    public double zVoxelDim();


    /** Array of voxel dimensions in mm. */
    public double[] getVoxelDims();


    /** Array of data dimensions. */
    public int[] getDataDims();


    /**
     * Gets a boolean mask image for tracking. A voxel in the mask is <code>true</code> 
     * if streamlines should terminate on entry to the voxel.
     */
    public boolean[][][] getIsotropicMask();


    /**
     * Computes mask for tracking. A voxel in the mask is <code>true</code> 
     * if streamlines should terminate on entry to the voxel. Background voxels are always 
     * <code>true</code>.
     *
     * @param anisMap an image containing the some quantitative measurement of the 
     * anisotropy of the diffusion in the voxel. 
     * @param threshold voxels with anisotropy below this value will be set 
     * to <code>true</code>.
     * 
     */
    public void computeIsotropicMask(double[][][] scalarMap, double threshold);

    /**
     * For probabilistic images, this method resets the cache so that the next call to 
     * <code>getPDs</code> returns a new sample. Does nothing for deterministic images.
     *
     */
    public void resetRandomization();

}

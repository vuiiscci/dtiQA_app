package tractography;

import numerics.*;

import java.util.Random;

/**
 * Superclass to probabilistic images. 
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public abstract class ProbabilisticTractographyImage extends PD_TractographyImage {
 
    private boolean[][][] randomized;

    protected final Random ran;

    /** line of dimension zDataDim containing false values */
    private final boolean[] line;


    /**
     *
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param vc the voxel classification.
     * @param r a source of random numbers. The java.util.Random class is not recommended, use
     * tools.MTRandom instead.
     *
     */
    public ProbabilisticTractographyImage(int[] dataDims, double[] voxelDims, Random r) {

	super(dataDims, voxelDims);

	ran = r;

        randomized = new boolean[xDataDim][yDataDim][zDataDim];

        line = new boolean[zDataDim];

    }


    /**
     * Returns the principal directions in this voxel. Multiple calls will return the same
     * vectors until resetRandomization() is called.
     *
     * @return the principal directions in this voxel.
     */
    public final Vector3D[] getPDs(int i, int j, int k) {
        
        if (!randomized[i][j][k]) {
            setVectors(i,j,k);
            randomized[i][j][k] = true;
        }
        
        return super.getPDs(i,j,k);
    }


    /**
     * Returns the principal directions in this voxel. Multiple calls will return the same
     * vectors until resetRandomization() is called.
     *
     * @return the principal directions in this voxel.
     */
    public final Vector3D[] getPDs(int i, int j, int k, Vector3D fibreOrientation) {  
        if (!randomized[i][j][k]) {
            setVectors(i,j,k, fibreOrientation);
            randomized[i][j][k] = true;
        }
        
        return super.getPDs(i,j,k);
   
    }


    /**
     *
     * Sets the vectors for position i, j, k. Called only if a voxel is not already randomized.
     *
     */
    protected abstract void setVectors(int i, int j, int k);


    /**
     *
     * Sets the vectors for position i, j, k, given a previous direction. Unless overriden, the previous direction is not used
     * and the method simply calls setVectors(int, int, int).
     *
     */
    protected void setVectors(int i, int j, int k, Vector3D previousDirection) {
        setVectors(i,j,k);
    }
   



    /**
     * Resets the cache so that the next call to <code>getPDs</code> returns a new sample.
     *
     */
    public final void resetRandomization() {
        
        // doing it this way will make this method have a higher percentage in the profiler
        // but the total execution time will be less than if we allocate a new volume
        for (int j = 0; j < yDataDim; j++) {
            for (int i = 0; i < xDataDim; i++) {
                System.arraycopy(line, 0, randomized[i][j], 0, zDataDim);
            }
        }
    }

}

package tractography;

import data.*;
import numerics.*;


/**
 * Superclass for tractography images. Images store the data for tracking, including PICo PDF parameters.
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class PD_TractographyImage implements TractographyImage {

    
    /** Size of image in x dimension. */
    protected final int xDataDim;

    /** Size of image in y dimension. */
    protected final int yDataDim;

    /** Size of image in z dimension. */
    protected final int zDataDim;

    /** Size (mm) of x voxel length. */
    protected final double xVoxelDim;

    /** Size (mm) of y voxel length. */
    protected final double yVoxelDim;

    /** Size (mm) of z voxel length. */
    protected final double zVoxelDim;
    
    /** Number of PDs in each voxel */
    protected final int[][][] numPDs;

    /** 
     * List of unit vectors containing the principal directions, if any, in each voxel.
     *
     * May be lazily initialized by subclasses.
     */
    protected final Vector3D[][][][] vectors; 

    protected boolean[][][] isotropicMask = null;

    /**
     * Basic constructor initialises dimensions only.
     *
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public PD_TractographyImage(int[] dataDims, double[] voxelDims) {
	
	xVoxelDim = voxelDims[0];
	yVoxelDim = voxelDims[1];
	zVoxelDim = voxelDims[2];

	xDataDim = dataDims[0];
	yDataDim = dataDims[1];
	zDataDim = dataDims[2];

	vectors = new Vector3D[xDataDim][yDataDim][zDataDim][];
	numPDs = new int[xDataDim][yDataDim][zDataDim];

        isotropicMask = new boolean[xDataDim][yDataDim][zDataDim];

    }


    /**
     * Matlab constructor.
     *
     * @param vectors vector data, should be one vector per PD. The number of PDs per voxel is inferred 
     * from this array.
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public PD_TractographyImage(Vector3D[][][][] vectors, double[] voxelDims) {
	
	xVoxelDim = voxelDims[0];
	yVoxelDim = voxelDims[1];
	zVoxelDim = voxelDims[2];

	xDataDim = vectors.length;
	yDataDim = vectors[0].length;
	zDataDim = vectors[0][0].length;
	
	this.vectors = vectors;
	numPDs = new int[xDataDim][yDataDim][zDataDim];

        for (int k = 0; k < zDataDim; k++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
		    numPDs[i][j][k] = vectors[i][j][k].length;
		}
	    }
	}	
	
	computeIsotropicMask();
    }



    /**
     *
     * @return the princpal direction or directions in the voxel. Probabilistic images will 
     * return a sample principal direction, successive calls to this method will return
     * the same direction until <code>resetRandomization</code> is called.
     *
     */
    public Vector3D[] getPDs(int i, int j, int k) {
	

	if (numPDs[i][j][k] == 1) {
	    return new Vector3D[] {vectors[i][j][k][0]};
	}
        else if (numPDs[i][j][k] == 0) {
            return new Vector3D[] {};
        }
	else { 

            Vector3D[] pds = new Vector3D[numPDs[i][j][k]];
            
            for (int x = 0; x < numPDs[i][j][k]; x++) {
                pds[x] = vectors[i][j][k][x];
            }

            return pds;
	}
 
    }


    /**
     * @return the princpal direction or directions in the voxel, given an estimated local fiber orientation.
     * If the PDs in the voxel are not dependent on the fiber orientation, then this method returns the
     * same thing as getPDs(int, int, int).
     *
     * @return a sample fibre orientation from this voxel. Note that this will <b>not</b> correct the orientation 
     * for agreement with the previous tracking direction so there is no guarantee that the returned vector 
     *will have a positive dot product with <code>fibreOrientation</code>.
     *
     */
    public Vector3D[] getPDs(int i, int j, int k, Vector3D fibreOrientation) {
	return getPDs(i,j,k);
    }



    /**
     * @return the number of PDs in this voxel.
     *
     */
    public final int numberOfPDs(int i, int j, int k) {
        return numPDs[i][j][k];
    }
   

    /** Size of image in x dimension. */
    public final int xDataDim() {
	return xDataDim;
    }

    /** Size (mm) of x voxel length. */
    public final double xVoxelDim() {
	return xVoxelDim;
    }

    /** Size of image in y dimension. */
    public final int yDataDim() {
	return yDataDim;
    }

    /** Size (mm) of y voxel length. */
    public final double yVoxelDim() {
	return yVoxelDim;
    }

    /** Size of image in z dimension. */
    public final int zDataDim() {
	return zDataDim;
    }

    /** Size (mm) of z voxel length. */
    public final double zVoxelDim() {
	return zVoxelDim;
    }

    /** Array of voxel dimensions in mm. */
    public final double[] getVoxelDims() {
	return new double[] {xVoxelDim, yVoxelDim, zVoxelDim};
    }

    /** Array of data dimensions. */
    public final int[] getDataDims() {
	return new int[] {xDataDim, yDataDim, zDataDim};
    }


    /**
     * Computes mask for tracking. A voxel in the mask is <code>true</code> 
     * if streamlines should terminate on entry to the voxel. This method creates a brain / background
     * mask. A voxel is isotropic if the number of PDs in the voxel is zero.
     * 
     *
     */
    public final void computeIsotropicMask() {

	isotropicMask = new boolean[xDataDim][yDataDim][zDataDim];

	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
			
                    if (numPDs[i][j][k] == 0) {
                        isotropicMask[i][j][k] = true;
                    }
                    
                }
            }
        }
		
    }


    public final void computeIsotropicMask(double[][][] anisMap, double threshold) {

        isotropicMask = new boolean[xDataDim][yDataDim][zDataDim];

	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
					
		    if (numPDs[i][j][k] == 0) {
			isotropicMask[i][j][k] = true;
		    }
		    else if (anisMap[i][j][k] < threshold) {
			isotropicMask[i][j][k] = true;
		    }
		    
		}
	    }
	}
	
    }


    /**
     * Gets a boolean mask image for tracking. A voxel in the mask is <code>true</code> 
     * if streamlines should terminate on entry to the voxel.
     */
    public final boolean[][][] getIsotropicMask() {

        boolean[][][] copy = new boolean[xDataDim][yDataDim][zDataDim];

        for (int j = 0; j < yDataDim; j++) {
            for (int i = 0; i < xDataDim; i++) {
                System.arraycopy(isotropicMask[i][j], 0, copy[i][j], 0, zDataDim);
            }
        }

	return copy;
    }


    /**
     * For probabilistic images, this method resets the cache so that the next call to 
     * <code>getPDs</code> returns a new sample. Does nothing for deterministic images.
     *
     */
    public void resetRandomization() {
        
    }


    /**
     * Gets an image from the data files. 
     * 
     *
     * @param inputFiles the data files, each containing one vector per voxel. Usually
     * some voxels will have fewer than the maximum number of vectors, in which case the last
     * vector or vectors will be zero. The input files should be ordered such that the non-zero 
     * vectors are always first. Data files may be images or raw files with the specified type
     * and dimensions.
     * 
     * @param dataType the data type of the data files and <code>anisMapFile</code>. Used for
     * raw data only. 
     * 
     * @param anisMap the anisotropy map, which is used to create the tract mask.
     * May be <code>null</code> if not required, in which case voxels will be marked isotropic
     * only if they contain zero vectors.
     *
     * @param anisThresh threshold for the anisotropy, used only if anisMap is not null.
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public static final PD_TractographyImage getTractographyImage(String[] inputFiles,
                                                                  String dataType, 
                                                                  double[][][] anisMap,
                                                                  double anisThresh,
                                                                  int[] dataDims,
                                                                  double[] voxelDims) {

        int numComponents = 3;

        // Read vectors into an array and then construct image

	int maxPDs = inputFiles.length;
	
	Vector3D[][][][] vectors = new Vector3D[dataDims[0]][dataDims[1]][dataDims[2]][];

	DataSource[] sources = new DataSource[maxPDs];

	for (int i = 0; i < maxPDs; i++) {
	    sources[i] = ExternalDataSource.getDataSource(inputFiles[i], numComponents, dataType);
	}

	for (int k = 0; k < dataDims[2]; k++) {
	    for (int j = 0; j < dataDims[1]; j++) {
		for (int i = 0; i < dataDims[0]; i++) {
		    
		    // determine how many PDs exist for this voxel
		    Vector3D[] tmpVec = new Vector3D[maxPDs];
		    
		    int numPDs = 0;

		    // catch PDs in wrong order with this flag
		    boolean gotLastPD = false;

		    for (int n = 0; n < maxPDs; n++) {
			double[] components = sources[n].nextVoxel();
			
			if (Math.abs(components[0]) + Math.abs(components[1]) + Math.abs(components[2]) > 0.0) {
			    
			    if (gotLastPD) {
				throw new IllegalArgumentException("Vector input in wrong order, non-zero vectors " + 
								   "must come first");
			    }

			    tmpVec[n] = new Vector3D(components);
			    numPDs++;
			}
			else {
			    gotLastPD = true;
			    // don't break here or input files will get out of sync
			}
		    }
		    
		    vectors[i][j][k] = new Vector3D[numPDs];

		    System.arraycopy(tmpVec, 0, vectors[i][j][k], 0, numPDs);

		}
	    }
	}
	

        PD_TractographyImage image = new PD_TractographyImage(vectors, voxelDims);

        if (anisMap != null) {
            image.computeIsotropicMask(anisMap, anisThresh);
        }

        return image;

    }

 
}

package tractography;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;


/**
 * Defines image for tractography from tensor data. 
 * 
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class DT_TractographyImage extends PD_TractographyImage implements TensorTractographyImage {

    
    private final DT[][][][] tensors;
    private final double[][][][] mix;
    

    /**
     * Constructs an image directly from an array of DTs and an array
     * of mixing parameters.  This allows direct access from matlab.
     *
     * @param dts The array of tensors.
     *
     * @param mps The array of mixing parameters
     *
     * @param nPDs The array of numbers of directions
     *
     * @param dataDims array of data dimensions {xDataDim,
     * yDataDim, zDataDim}.
     *
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim,
     * yVoxelDim, zVoxelDim}.
     */
    public DT_TractographyImage(DT[][][][] dts, double[][][][] mps, int[][][] nPDs, int[] dataDims, 
				double[] voxelDims) {
	
        super(dataDims, voxelDims);
        tensors = dts;
        mix = mps;


	for (int j = 0; j < yDataDim; j++) {
            for (int i = 0; i < xDataDim; i++) {
                System.arraycopy(nPDs[i][j], 0, numPDs[i][j], 0, zDataDim);
            }
        }
	

        computeIsotropicMask();
    }


    /**
     * Constructs an image from the data sources. 
     *
     * @param tensorSource sources for the diffusion tensor volumes.
     * @param maxPDs the maximum number of PDs in a voxel.
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     *
     */
    public DT_TractographyImage(DataSource tensorSource, int maxPDs, 
                                int[] dataDims, double[] voxelDims) {
	

	super(dataDims, voxelDims);

	tensors = new DT[xDataDim][yDataDim][zDataDim][];
        mix = new double[xDataDim][yDataDim][zDataDim][];
        
        for (int k = 0; k < zDataDim; k++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
                    
                    double[] voxel = tensorSource.nextVoxel();
                    
                    if (maxPDs == 1) {
                        numPDs[i][j][k] = (voxel[0] < 0.0) ? 0 : 1;
                        
                        tensors[i][j][k] = new DT[numPDs[i][j][k]];

                    
                        if (tensors[i][j][k].length == 1) {
                            tensors[i][j][k][0] = new DT(voxel[2], voxel[3], voxel[4],
                                                         voxel[5], voxel[6],voxel[7]);

                            mix[i][j][k] = new double[] {1.0};
                        }
                        else {
                             mix[i][j][k] = new double[] {};
                        }

                        

                    }
                    else {
                        numPDs[i][j][k] = (voxel[0] < 0.0) ? 0 : (int)voxel[2];

                        if (numPDs[i][j][k] < 0 || numPDs[i][j][k] > maxPDs) {
                            throw new LoggedException("Invalid number of components in input data. " + 
                                                      "Check -inputmodel and -numpds options");
                        }   

                        tensors[i][j][k] = new DT[numPDs[i][j][k]];
                        mix[i][j][k] = new double[numPDs[i][j][k]];
                        for (int p = 0; p < numPDs[i][j][k]; p++) {

                            // index of first tensor element
                            int start = 4 + p * 7;

                            mix[i][j][k][p] = voxel[start - 1];
                            tensors[i][j][k][p] = new DT(voxel[start], voxel[start+1], voxel[start+2],
                                                         voxel[start+3], voxel[start+4], voxel[start+5]);
                            
                        }
                        
                    }
                }
            }
        }
        
		
	computeIsotropicMask();

    }


    /**
     * Copy constructor. Does not copy cached or publically mutable data.
     *
     *
     */
    protected DT_TractographyImage(DT_TractographyImage im) {
	super(new int[] {im.xDataDim, im.yDataDim, im.zDataDim}, 
	      new double[] {im.xVoxelDim, im.yVoxelDim, im.zVoxelDim});

	tensors = new DT[xDataDim][yDataDim][zDataDim][];
        mix = new double[xDataDim][yDataDim][zDataDim][];

	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    
                    numPDs[i][j][k] = im.numPDs[i][j][k];
		    
		    tensors[i][j][k] = new DT[numPDs[i][j][k]];
		    mix[i][j][k] = new double[numPDs[i][j][k]];

		    for (int t = 0; t < numPDs[i][j][k]; t++) {
			tensors[i][j][k][t] = im.tensors[i][j][k][t];	
                        mix[i][j][k][t] = im.mix[i][j][k][t];
		    }
			    
		}
	    }
	    
	}
	computeIsotropicMask();

	
    }    
    
    
  
    /**
     *
     * Computes boolean mask, a voxel is true if the voxel has no PDs, or if the average 
     * FA over all compartments is below the threshold. Background voxels are always classified 
     * as isotropic.
     * 
     */
    public void computeIsotropicMask(double faThreshold) {

	isotropicMask = new boolean[xDataDim][yDataDim][zDataDim];

	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
			
		    if (numPDs[i][j][k] == 0) {
			isotropicMask[i][j][k] = true;
		    }
		    else if ( numPDs[i][j][k] == 1 && tensors[i][j][k][0].fa() < faThreshold) {
                        isotropicMask[i][j][k] = true;
                    }
                    else {
                        double sumFA = tensors[i][j][k][0].fa();

                        for (int d = 1; d < numPDs[i][j][k]; d++) {
                            sumFA += tensors[i][j][k][d].fa();
                        }

                        double avgFA = sumFA / numPDs[i][j][k];

                        if (avgFA < faThreshold) {
                            isotropicMask[i][j][k] = true;
                        }
                    }
		    
                } 
	    }
	}
	
    }


    /**
     * @return the principal directions in this voxel.
     */
    public Vector3D[] getPDs(int i, int j, int k) {

	if (vectors[i][j][k] == null) {
	    
	    vectors[i][j][k] = new Vector3D[tensors[i][j][k].length];

	    for (int t = 0; t < tensors[i][j][k].length; t++) {

		double[][] eig = tensors[i][j][k][t].sortedEigenSystem();

		vectors[i][j][k][t] = new Vector3D(eig[1][0], eig[2][0], eig[3][0]);
	    }
	    
	}

	return super.getPDs(i,j,k);
        
    }




    /**
     * @return the tensors in this voxel. 
     */
    public final DT[] getDTs(int i, int j, int k) {
     
        DT[] defCopy = new DT[numPDs[i][j][k]];
        
        System.arraycopy(tensors[i][j][k], 0, defCopy, 0, numPDs[i][j][k]);

        return defCopy;
    
    }

    
    /**
     * @return the mixing parameters in this voxel.
     */
    public final double[] getMix(int i, int j, int k) {

        double[] defCopy = new double[numPDs[i][j][k]];
        
        System.arraycopy(mix[i][j][k], 0, defCopy, 0, numPDs[i][j][k]);

        return defCopy;
    }

   
    /**
     * Gets an image from the data file. 
     * If <code>anisMapFile</code> is not <code>null</code>, it is read and used 
     * for isotropic masking.
     * 
     *
     * @param inputFile the data file.
     * @param dataType the data type of the data file and <code>anisMapFile</code>.
     * @param maxPDs the maximum number of PDs in a voxel.
     * @param anisMap the anisotropy map, which is used to create the tract mask.
     * May be <code>null</code> if not required, in which case fractional anisotropy will be used.
     * @param anisThresh threshold for the anisotropy in the computation of the tract mask.
     * 
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public static final DT_TractographyImage getTractographyImage(String inputFile, 
                                                                  String dataType, int maxPDs, 
                                                                  double[][][] anisMap, 
                                                                  double anisThresh, 
                                                                  int[] dataDims, 
                                                                  double[] voxelDims) {


        int numComponents = 8;

        if (maxPDs > 1) {
            numComponents = 3 + maxPDs * 7;
        }
    
	DataSource dataSource = ExternalDataSource.getDataSource(inputFile, numComponents, dataType);
     
        DT_TractographyImage image = new DT_TractographyImage(dataSource, maxPDs, dataDims, voxelDims);

        if (anisMap != null) {
            image.computeIsotropicMask(anisMap, anisThresh);
        }
        else if (anisThresh > 0.0) {
            image.computeIsotropicMask(anisThresh);
        }

        return image;	    

    }

    
}

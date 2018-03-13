package tractography;

import data.*;
import imaging.*;
import inverters.*;
import misc.DT;
import numerics.*;

import java.util.Random;

/**
 * Superclass for bootstrap images that fit models to DWI data on the fly. 
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class RepBS_DWI_TractographyImage extends BS_DWI_TractographyImage {
  
    private final float[][][][][] data;
    
    private final int numRepeats; // number of repeats, ie number of times each measurement is made

    /**
     *
     * @param data the raw bootstrap data, order [x][y][z][repeat][measurement].
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param nPDs number of PDs per voxel.
     * @param fitter the method to fit a model to the data and extract the principal diretions. 
     * @param r a source of random numbers. The java.util.Random class is not recommended, use
     * tools.MTRandom instead.
     *
     */
    public RepBS_DWI_TractographyImage(float[][][][][] data, double[] voxelDims, int[][][] nPDs, 
                                       ModelFitter fitter, Random r) {

	super(new int[] {data.length, data[0].length, data[0][0].length}, voxelDims, fitter, r);

	this.data = data;

	numRepeats = data[0][0][0].length;

        for (int j = 0; j < yDataDim; j++) {
            for (int i = 0; i < xDataDim; i++) {
                System.arraycopy(nPDs[i][j], 0, numPDs[i][j], 0, zDataDim);
            }
        }

        computeIsotropicMask();
    }


     
       
    /**
     * @return the next bootstrap sample of data.
     *
     */
    protected final double[] getData(int i, int j, int k) {

	double[] sample = new double[numMeasurements];

	for (int m = 0; m < numMeasurements; m++) {
            
            int r = numRepeats > 1 ? ran.nextInt(numRepeats) : 0;

	    sample[m] = data[i][j][k][r][m];
	}

	return sample;
    }


    
   /**
     * Gets an image from the data file. 
     * If <code>anisMapFile</code> is not <code>null</code>, it is read and used 
     * for isotropic masking.
     * 
     *
     * @param inputFiles the data files.
     * @param dataType the data type of the data file and <code>anisMapFile</code>.
     * @param imPars the imaging scheme for the data.
     * @param inversionIndex the numerical index of the tensor inversion.
     * @param vc the voxel classification for this image, either the number of PDs (0, 1, 2...) or 
     * for tensor data the SH order (-1, 0, 2, 4). If <code>null</code>,  each voxel is assumed to contain one PD. 
     * @param anisMap the anisotropy map, which is used to create the isotropic mask.
     * May be <code>null</code> if not required.
     * @param anisThresh threshold for the anisotropy in the computation of the tract mask.
     * 
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param ran a source of random numbers.
     */
    public static final RepBS_DWI_TractographyImage getTractographyImage(String[] inputFiles, 
									    String dataType, 
									    DW_Scheme imPars,
									    ModelIndex[] indices,
									    int[][][] vc,
									    double[][][] anisMap, 
									    double anisThresh, 
									    int[] dataDims, 
									    double[] voxelDims, 
									    Random ran) {
        
	int xDataDim = dataDims[0];
	int yDataDim = dataDims[1];
	int zDataDim = dataDims[2];
	
	int numRepeats = inputFiles.length;
	
	int numMeasurements = imPars.numMeasurements();

	float[][][][][] data = new float[xDataDim][yDataDim][zDataDim][numRepeats][numMeasurements];
	
	for (int r = 0; r < numRepeats; r++) {

	    DataSource component = ExternalDataSource.getDataSource(inputFiles[r], numMeasurements, dataType);

	    for (int k = 0; k < zDataDim; k++) {
		for (int j = 0; j < yDataDim; j++) {
		    for (int i = 0; i < xDataDim; i++) {

			double[] d = component.nextVoxel();

			for (int c = 0; c < numMeasurements; c++) {
			    data[i][j][k][r][c] = (float)d[c];
			}
		    }
		}
	    }
	}

	int[][][] voxClass = vc;

	if (voxClass == null) {
	    voxClass = new int[xDataDim][yDataDim][zDataDim];
            
            for (int k = 0; k < zDataDim; k++) {
		for (int j = 0; j < yDataDim; j++) {
		    for (int i = 0; i < xDataDim; i++) {
                        voxClass[i][j][k] = 2;
                    }
                    
                }
            }
        }
        else {
            for (int k = 0; k < zDataDim; k++) {
		for (int j = 0; j < yDataDim; j++) {
		    for (int i = 0; i < xDataDim; i++) {
                        if (vc[i][j][k] < 0) {
                            voxClass[i][j][k] = 0;
                        }
                        else if (vc[i][j][k] > 3) {
                            voxClass[i][j][k] = 2;
                        }
                        else {
                            voxClass[i][j][k] = 1;
                        }
                        
                    }
                }
            }
        }
        
        ModelFitter fitter = ModelFitter.getModelFitter(imPars, indices);

        RepBS_DWI_TractographyImage image = null;

        if (fitter instanceof TensorModelFitter) {
            image = new DT_RepBS_DWI_TractographyImage(data, voxelDims, voxClass, fitter, ran);
            
        }
        else {
            image = new RepBS_DWI_TractographyImage(data, voxelDims, voxClass, fitter, ran);
        }

	
        if (anisMap != null) {
            image.computeIsotropicMask(anisMap, anisThresh);
        }

        return image;	    

    }

}

package tractography;

import data.*;
import imaging.*;
import inverters.*;
import misc.*;
import numerics.*;


/**
 * Superclass for images that take a single raw DWI data set, and fit some model on the fly.
 * 
 * 
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class DWI_TractographyImage extends PD_TractographyImage {

   
    private final float[][][][] data;

    private final ModelFitter fitter;

    private final int numMeasurements;
    
  
    /**
     * 
     *
     * @param dat DWI data array.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param nPDs array containing the number of PDs per voxel. Should be 0 for background.
     * @param fit the model fitter that provides principal directions given some data.
     * 
     */
    public DWI_TractographyImage(float[][][][] dat, double[] voxelDims, int[][][] nPDs, ModelFitter fit) { 
	
	super(new int[] {dat.length, dat[0].length, dat[0][0].length}, voxelDims);


        data = dat;

        fitter = fit;

	numMeasurements = fitter.getScheme().numMeasurements();

        for (int j = 0; j < yDataDim; j++) {
            for (int i = 0; i < xDataDim; i++) {
                System.arraycopy(nPDs[i][j], 0, numPDs[i][j], 0, zDataDim);
            }
        }
        
	
	computeIsotropicMask();
    }



    public Vector3D[] getPDs(int i, int j, int k) {
        
	if (vectors[i][j][k] == null) {
            fitter.fitModel(getData(i,j,k), numPDs[i][j][k]);
            
            vectors[i][j][k] = fitter.getPDs();

        }
        
        return super.getPDs(i,j,k);
    }


 
    public final double[] getData(int i, int j, int k) {
        
        double[] dat = new double[numMeasurements];

        for (int n = 0; n < numMeasurements; n++) {
            dat[n] = data[i][j][k][n];
        }

        return dat;
    }


    public final double getMeasurement(int i, int j, int k, int meas) {
        return data[i][j][k][meas];
    }

    
    public final int numMeasurements() {
        return numMeasurements;
    }


    public final Vector3D[] getPDs(double[] data, int numPDsInVox) {
        fitter.fitModel(data, numPDsInVox);

        return fitter.getPDs();
    }




    /**
     * Gets an image from the data file. 
     * If <code>anisMapFile</code> is not <code>null</code>, it is read and used 
     * for isotropic masking.
     * 
     *
     * @param inputFile the data file.
     * @param dataType the data type of the data file and <code>anisMapFile</code>.
     * @param imPars the imaging scheme for the data.
     * @param vc a voxel classification map containing SH orders
     * @param anisMap the anisotropy map, which is used to create the isotropic mask.
     * May be <code>null</code> if not required.
     * @param anisThresh threshold for the anisotropy in the computation of the tract mask.
     * 
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public static final DWI_TractographyImage getTractographyImage(String inputFile, 
                                                             String dataType, 
                                                             DW_Scheme imPars,
                                                             ModelIndex[] indices,
                                                             int[][][] vc,
                                                             double[][][] anisMap, 
                                                             double anisThresh, 
                                                             int[] dataDims, 
                                                             double[] voxelDims) {
      
	int xDataDim = dataDims[0];
	int yDataDim = dataDims[1];
	int zDataDim = dataDims[2];
        
        // components per voxel in the data
        int numComponents = imPars.numMeasurements();
        
        DataSource dataSource = ExternalDataSource.getDataSource(inputFile, numComponents, dataType);
        
	float[][][][] data = new float[xDataDim][yDataDim][zDataDim][numComponents];
        
        int[][][] numPDs = new int[xDataDim][yDataDim][zDataDim];
        
	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    data[i][j][k] = new float[numComponents];
                    
		    double[] tmp = dataSource.nextVoxel();
                    
		    for (int c = 0; c < numComponents; c++) {
			data[i][j][k][c] = (float)tmp[c];
		    }

                    if (vc[i][j][k] < 0) {
                        numPDs[i][j][k] = 0;
                    }
                    else if (vc[i][j][k] > 3) {
                        numPDs[i][j][k] = 2;
                    }
                    else {
                        numPDs[i][j][k] = 1;
                    }
		}
	    }
	}
        
        ModelFitter fitter = ModelFitter.getModelFitter(imPars, indices);
        
        DWI_TractographyImage image = null;

        if (fitter instanceof TensorModelFitter) {
            
            image = new DT_DWI_TractographyImage(data, voxelDims, numPDs, (TensorModelFitter)fitter);
	
        }
        else {
           image = new DWI_TractographyImage(data, voxelDims, numPDs, fitter); 
        }
        if (anisMap != null) {
            image.computeIsotropicMask(anisMap, anisThresh);
        }
        
        return image;	    
        
    }
    
    
    
    
}

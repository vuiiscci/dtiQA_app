package tractography;

import data.*;
import imaging.*;
import inverters.*;
import misc.DT;
import numerics.*;

import java.util.Random;

/**
 * Each call to #getPDs(int, int, int) returns a new bootstrap estimate of the PDs. 
 * The bootstrap data is inverted using the specified one and two tensor index. 
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class DT_WildBS_DWI_TractographyImage extends BS_DWI_TractographyImage implements TensorTractographyImage {
 
    private final float[][][][] data;

    private final DT[][][][] tensors;
  
    private final WildBS_DT_DataSynth dataSynth;
    
    private final TensorModelFitter tensorFitter;



    /**
     *
     * @param data raw diffusion-weighted data.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param vc the voxel classification.
     * @param imPars the imaging scheme for the data.
     * @param r a source of random numbers. The java.util.Random class is not recommended, use
     * tools.MTRandom instead.
     *
     *
     */
    public DT_WildBS_DWI_TractographyImage(float[][][][] data, double[] voxelDims, int[][][] nPDs, DW_Scheme imPars,  Random r) {

	super(new int[] {data.length, data[0].length, data[0][0].length}, 
	      voxelDims, new TensorModelFitter(imPars, new ModelIndex[] {ModelIndex.LDT}), r);

	this.data = data;

	dataSynth = new WildBS_DT_DataSynth(null, imPars, -1, r);

        tensors = new DT[xDataDim][yDataDim][zDataDim][];

        tensorFitter = (TensorModelFitter)fitter;

        
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
    protected double[] getData(int i, int j, int k) {
	dataSynth.setSourceData(data[i][j][k]);
	
	return dataSynth.nextVoxel();
    }


    /**
     *
     * Sets the vectors for position i, j, k. Called only if a voxel is not already randomized.
     *
     */
    protected void setVectors(int i, int j, int k) {
        
        double[] data = getData(i,j,k);

        tensorFitter.fitModel(data,numPDs[i][j][k]);

        vectors[i][j][k] = tensorFitter.getPDs();

        tensors[i][j][k] = tensorFitter.getDTs();
                
    }

    
    /**
     * @return the tensors in this voxel. 
     */
    public final DT[] getDTs(int i, int j, int k) {
       
        // Forces randomization if it hasn't already happened
        getPDs(i,j,k);
        
        return new DT[] {tensors[i][j][k][0]};
    
    }

    
    /**
     * @return the mixing parameters in this voxel.
     */
    public final double[] getMix(int i, int j, int k) {
        return new double[] {1.0};
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
     * @param vc the voxel classification for this image. This class only supports one PD per voxel. 
     * This volume should be used to distinguish brain from background.
     * @param anisMap the anisotropy map, which is used to create the isotropic mask.
     * May be <code>null</code> if not required.
     * @param anisThresh threshold for the anisotropy in the computation of the tract mask.
     * 
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public static final DT_WildBS_DWI_TractographyImage getTractographyImage(String inputFile, 
                                                                             String dataType, 
                                                                             DW_Scheme imPars,
                                                                             int[][][] vc,
                                                                             double[][][] anisMap, 
                                                                             double anisThresh, 
                                                                             int[] dataDims, 
                                                                             double[] voxelDims,
                                                                             Random ran) {

      
	int xDataDim = dataDims[0];
	int yDataDim = dataDims[1];
	int zDataDim = dataDims[2];
        
        // components per voxel in the data
        int numComponents = imPars.numMeasurements();

        DataSource dataSource = ExternalDataSource.getDataSource(inputFile, numComponents, dataType);

	float[][][][] data = new float[dataDims[0]][dataDims[1]][dataDims[2]][numComponents];


	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    data[i][j][k] = new float[numComponents];

		    double[] tmp = dataSource.nextVoxel();

		    for (int c = 0; c < numComponents; c++) {
			data[i][j][k][c] = (float)tmp[c];
		    }
		}
	    }
	}


        int[][][] voxClass = new int[xDataDim][yDataDim][zDataDim];

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
        

	DT_WildBS_DWI_TractographyImage image = new DT_WildBS_DWI_TractographyImage(data, voxelDims, voxClass, imPars, ran);
	
   
        if (anisMap != null) {
            image.computeIsotropicMask(anisMap, anisThresh);
        }

        return image;	    

    }

   

}

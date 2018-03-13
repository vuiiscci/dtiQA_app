package tractography;

import data.*;
import imaging.*;
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
public abstract class BS_DWI_TractographyImage extends ProbabilisticTractographyImage {


    protected final int numMeasurements; // number of measurements in each voxel

    protected final ModelFitter fitter;

 

    /**
     *
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param vc the voxel classification.
     * @param modelFitter the model fitter to use to get PDs from bootstrapped data.
     * @param r a source of random numbers. The java.util.Random class is not recommended, use
     * tools.MTRandom instead.
     */
    public BS_DWI_TractographyImage(int[] dataDims, double[] voxelDims, ModelFitter modelFitter, Random r) {
        
	super(dataDims, voxelDims, r);
           
        fitter = modelFitter;

	numMeasurements = fitter.numMeasurements();
     
    }



    /**
     *
     * Sets the vectors for position i, j, k. Called only if a voxel is not already randomized.
     *
     */
    protected void setVectors(int i, int j, int k) {
        
        double[] data = getData(i,j,k);

        fitter.fitModel(data,numPDs[i][j][k]);

        vectors[i][j][k] = fitter.getPDs();

    }




    /**
     * Gets a bootstrap sample of raw DWI data in some voxel. 
     *
     */
    protected abstract double[] getData(int i, int j, int k);




 
}

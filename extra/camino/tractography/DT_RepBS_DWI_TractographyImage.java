package tractography;

import data.*;
import imaging.*;
import inverters.*;
import misc.DT;
import numerics.*;

import java.util.Random;

/**
 * Each call to #setVectors(int, int, int) returns a new bootstrap estimate of the PDs. 
 * The bootstrap data is inverted using the specified one and two tensor index. 
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class DT_RepBS_DWI_TractographyImage extends RepBS_DWI_TractographyImage implements TensorTractographyImage {
 
   
    private final DT[][][][] tensors;
  
    private final double[][][][] mix;

    private final TensorModelFitter tensorFitter;


    /**
     *
     * @param data the raw bootstrap data, order [x][y][z][repeat][measurement].
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param nPDs number of PDs per voxel.
     * @param fitter the method to fit a model to the data and extract the principal diretions.
     * @param r a source of random numbers. The java.util.Random class is not recommended, use
     * tools.MTRandom instead.
     *
     * @see apps.ModelFit
     */
    public DT_RepBS_DWI_TractographyImage(float[][][][][] data, double[] voxelDims, int[][][] nPDs, 
                                          ModelFitter tensorFitter, Random r) {

	super(data, voxelDims, nPDs, tensorFitter, r);

        tensors = new DT[xDataDim][yDataDim][zDataDim][];
        mix = new double[xDataDim][yDataDim][zDataDim][];

        // need a reference cast to TensorModelFitter
        this.tensorFitter = (TensorModelFitter)tensorFitter;
    
                
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
        
        mix[i][j][k] = tensorFitter.getMix();

    }


     
     /**
     * @return the tensors in this voxel. 
     */
    public final DT[] getDTs(int i, int j, int k) {

        // Randomizes if necessary
        getPDs(i,j,k);

        DT[] defCopy = new DT[numPDs[i][j][k]];
        
        System.arraycopy(tensors[i][j][k], 0, defCopy, 0, numPDs[i][j][k]);

        return defCopy;
    
    }

    
    /**
     * @return the mixing parameters in this voxel.
     */
    public final double[] getMix(int i, int j, int k) {

        // Randomizes if necessary
        getPDs(i,j,k);

        double[] defCopy = new double[numPDs[i][j][k]];
        
        System.arraycopy(mix[i][j][k], 0, defCopy, 0, numPDs[i][j][k]);

        return defCopy;
    }

   


}

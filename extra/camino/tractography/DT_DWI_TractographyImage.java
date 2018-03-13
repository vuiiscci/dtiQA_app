package tractography;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;

import java.util.Random;


/**
 * Superclass for images that take a single raw DWI data set, and fit some model in the fly.
 * 
 * 
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class DT_DWI_TractographyImage extends DWI_TractographyImage implements TensorTractographyImage {

   
    private final DT[][][][] tensors;
    
    private final double[][][][] mix;
    
    private final TensorModelFitter tensorFitter;


    /**
     * 
     *
     * @param 
     * @param 
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public DT_DWI_TractographyImage(float[][][][] data, double[] voxelDims, int[][][] numPDs, 
                                    TensorModelFitter fitter) { 
	
	super(data, voxelDims, numPDs, fitter);

        tensorFitter = fitter;

        tensors = new DT[xDataDim][yDataDim][zDataDim][];
        mix = new double[xDataDim][yDataDim][zDataDim][];

    }



    public Vector3D[] getPDs(int i, int j, int k) {
        
	if (vectors[i][j][k] == null) {
            tensorFitter.fitModel(getData(i,j,k), numPDs[i][j][k]);
            
            vectors[i][j][k] = tensorFitter.getPDs();
            
            tensors[i][j][k] = tensorFitter.getDTs();
            
            mix[i][j][k] = tensorFitter.getMix();
        }

        return super.getPDs(i,j,k);

    }

     
    /**
     * @return the tensors in this voxel. 
     */
    public final DT[] getDTs(int i, int j, int k) {
       
        if (vectors[i][j][k] == null) {
            getPDs(i,j,k);
        }

        DT[] defCopy = new DT[numPDs[i][j][k]];
        
        System.arraycopy(tensors[i][j][k], 0, defCopy, 0, numPDs[i][j][k]);

        return defCopy;
    
    }

    
    /**
     * @return the mixing parameters in this voxel.
     */
    public final double[] getMix(int i, int j, int k) {
        
        if (vectors[i][j][k] == null) {
            getPDs(i,j,k);
        }

        double[] defCopy = new double[numPDs[i][j][k]];
        
        System.arraycopy(mix[i][j][k], 0, defCopy, 0, numPDs[i][j][k]);

        return defCopy;
    }

   

   
}

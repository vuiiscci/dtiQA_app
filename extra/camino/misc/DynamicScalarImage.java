package misc;

import tractography.*;
import numerics.*;


/**
 *
 * Scalar image that allows computation of the mean, min, and max values in each voxel. 
 * This is a subset of the capabilities of the SparseVectorImage, but this the above quantities 
 * can be generated over very large data sets without keeping individual voxel values, so this
 * class is used where possible to avoid running out of memory.
 *
 * @version $Id$
 * @author Philip Cook
 */
public class DynamicScalarImage implements WeightedVoxelwiseStatisticalImage {


    // inefficiently storing everything here, could subclass and only 
    // compute the things we actually need

    // all these quantities are weighted if weights are passed. Else they all have unit weights.
    private double[][][] sum;
    private double[][][] min;
    private double[][][] max;

  
    // norm is also weighted, if weights are not used then norm[i][j][k] == numValues[i][j][k]
    private double[][][] norm;
    
    private boolean[][][] touchedVoxel;

    private final int[] dataDims;
    private final double[] voxelDims;
    
    /**
     * Construct an image, default interpolation is nearest neighbour.
     *
     */
    public DynamicScalarImage(int[] dataDims, double[] voxelDims) {
	
	this.dataDims = dataDims;

	this.voxelDims = voxelDims;
	
	sum    = new double[dataDims[0]][dataDims[1]][dataDims[2]];
	min    = new double[dataDims[0]][dataDims[1]][dataDims[2]];
	max    = new double[dataDims[0]][dataDims[1]][dataDims[2]];
	norm   = new double[dataDims[0]][dataDims[1]][dataDims[2]];

	touchedVoxel = new boolean[dataDims[0]][dataDims[1]][dataDims[2]];
    }

    
    /**
     * Adds a value to the voxel (i,j,k), with unit weight.
     */
    public void addValue(int i, int j, int k, double v) {
	addValue(i,j,k,v,1.0);
    }


    /**
     * Add an image of values, with custom weight. Images should be in the same dimension 
     * as this scalar image.
     *
     * @param weights must be non-negative.
     *
     * 
     *
     */
    public void addValues(double[][][] values, double[][][] weights) {
        for (int i = 0; i < dataDims[0]; i++) {
	    for (int j = 0; j < dataDims[1]; j++) {	
		for (int k = 0; k < dataDims[2]; k++) {
                    addValue(i,j,k, values[i][j][k], weights[i][j][k]);
                }
            }
        }
    }

    
    /**
     * Add an image of values, with unit weight.
     *
     *
     */
    public void addValues(double[][][] values) {
        for (int i = 0; i < dataDims[0]; i++) {
	    for (int j = 0; j < dataDims[1]; j++) {	
		for (int k = 0; k < dataDims[2]; k++) {
                    addValue(i,j,k, values[i][j][k], 1.0);
                }
            }
        }
    }
    

    /**
     * Add a value to the voxel (i,j,k), with custom weight.
     *
     * @param weight must be non-negative.
     *
     */

    public void addValue(int i, int j, int k, double v, double weight) {

	if (weight < 0.0) {
	    throw new LoggedException("Weights must be non-negative");
	}

	if (!touchedVoxel[i][j][k]) {
	    touchedVoxel[i][j][k] = true;
	    min[i][j][k] = v;
	    max[i][j][k] = v;
	}

        double wtd = v * weight;

	sum[i][j][k] += wtd;

	norm[i][j][k] += weight;
	
	if (v < min[i][j][k]) {
	    min[i][j][k] = v;
	}
	
	if (v > max[i][j][k]) {
	    max[i][j][k] = v;
	}
	
    }


 
    /**
     * Computes a 3D image, where each voxel intensity is some statistic of the data in each voxel.
     *
     * @param stat one of "mean", "max", "min", "sum". If no values have been added to a voxel, all of these
     * statistics are zero. 
     */
    public double[][][] getVoxelStatistic(String stat) {
	
	double[][][] result = new double[dataDims[0]][dataDims[1]][dataDims[2]];

	for (int i = 0; i < dataDims[0]; i++) {
	    for (int j = 0; j < dataDims[1]; j++) {	
		for (int k = 0; k < dataDims[2]; k++) {
		    
		    if (stat.equals("mean")) {
			result[i][j][k] = norm[i][j][k] > 0.0 ? sum[i][j][k] / norm[i][j][k] : 0.0; 
		    }
		    else if (stat.equals("max")) {
			result[i][j][k] = touchedVoxel[i][j][k] ? max[i][j][k] : 0.0;
		    }
		    else if (stat.equals("sum")) {
			result[i][j][k] = touchedVoxel[i][j][k] ? sum[i][j][k] : 0.0;
		    }
		    else if (stat.equals("min")) {
			result[i][j][k] = touchedVoxel[i][j][k] ? min[i][j][k] : 0.0;
		    }
		    else { 
			throw new LoggedException("Unsupported statistic : " + stat);
		    }

		}
	    }
	}

	return result;
	
    }


    public int[] getDataDims() {
	return new int[] {dataDims[0], dataDims[1], dataDims[2]};
    }

 
    public double[] getVoxelDims() {
	return new double[] {voxelDims[0], voxelDims[1], voxelDims[2]};
    }
 
   
}

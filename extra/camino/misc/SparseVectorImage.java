package misc;

import tractography.*;
import tools.*;
import numerics.*;


/**
 * Sparsely populated 3D image. The image contains a vector of scalar values at each point
 * in 3D space. The vector length is variable, and is zero if no values are added to the voxel.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class SparseVectorImage implements VoxelwiseStatisticalImage {

    
    // number of zero values for each voxel. Because most values in a sparse image
    // are zero, storing them this way saves space 
    private final int[][][] zeros;

    private final double[][][][] data;

    private final int[] dataDims;

    private final double[] voxelDims;
    
    private final int[][][] vectorLengths;
    
    protected final int initialArrayLength;

    private final static double GROWTH_FACTOR = 2.0;
  
    
    /**
     * Construct an image, default interpolation is nearest neighbour.
     *
     * @param dataDims 3D dimensions of the image.
     * @param voxelDims 3D dimensions of the image, in mm.
     */
    public SparseVectorImage(int[] dataDims, double[] voxelDims) {
	this.voxelDims = voxelDims;
	this.dataDims = dataDims;

        initialArrayLength = 50;

	data = new double[dataDims[0]][dataDims[1]][dataDims[2]][];
	vectorLengths = new int[dataDims[0]][dataDims[1]][dataDims[2]];
	zeros = new int[dataDims[0]][dataDims[1]][dataDims[2]];
    }


    /**
     * Construct an image, default interpolation is nearest neighbour. This constructor
     * takes the total number of images in the input. When this is known, storage can be
     * pre-allocated.
     * 
     *
     * @param dataDims 3D dimensions of the image.
     * @param voxelDims 3D dimensions of the image, in mm.
     * @param numImages the number of images.
     */
    public SparseVectorImage(int[] dataDims, double[] voxelDims, int numImages) {
	this.voxelDims = voxelDims;
	this.dataDims = dataDims;

        initialArrayLength = numImages + 1;

	data = new double[dataDims[0]][dataDims[1]][dataDims[2]][];
	vectorLengths = new int[dataDims[0]][dataDims[1]][dataDims[2]];
	zeros = new int[dataDims[0]][dataDims[1]][dataDims[2]];
    }

 

    /**
     * Add a value to the voxel (i,j,k).
     *
     *
     */
    public void addValue(int i, int j, int k, double v) {

        if (v == 0.0) {
            zeros[i][j][k]++;
        }
        else {

            if (data[i][j][k] == null) {
                data[i][j][k] = new double[initialArrayLength];
            }
            
            data[i][j][k][vectorLengths[i][j][k]] = v;
            
            vectorLengths[i][j][k]++;
            
            if (vectorLengths[i][j][k] == data[i][j][k].length) {
                growDataVector(i,j,k, GROWTH_FACTOR);
            }
        }
	
    }





    /**
     * Computes a 3D image, where each voxel intensity is some statistic of the associated vector.
     *
     * @param stat one of "mean", "max", "min", "median", "sum". "var", "std".
     */
    public double[][][] getVoxelStatistic(String stat) {
	
	double[][][] result = new double[dataDims[0]][dataDims[1]][dataDims[2]];

	for (int i = 0; i < dataDims[0]; i++) {
	    for (int j = 0; j < dataDims[1]; j++) {	
		for (int k = 0; k < dataDims[2]; k++) {
		    
		    if (vectorLengths[i][j][k] > 0) {
			if (stat.equals("mean")) {
			    result[i][j][k] = ArrayOps.mean(valuesAt(data, i,j,k));
			}
			else if (stat.equals("max")) {
			    result[i][j][k] = ArrayOps.max(valuesAt(data, i,j,k));
			}
			else if (stat.equals("min")) {
			    result[i][j][k] = ArrayOps.min(valuesAt(data, i,j,k));
			}
			else if (stat.equals("median")) {
			    result[i][j][k] = ArrayOps.median(valuesAt(data, i,j,k));
			}
                        else if (stat.equals("sum")) {
                            result[i][j][k] = ArrayOps.sum(valuesAt(data, i,j,k));
                        }
			else if (stat.equals("var") || stat.equals("std")) {
			    
			    double[] dataVec = valuesAt(data, i,j,k);
			    
			    result[i][j][k] = ArrayOps.var( dataVec, ArrayOps.mean(dataVec) );

                            if (stat.equals("std")) {
                                result[i][j][k] = Math.sqrt(result[i][j][k]);
                            }
			}
                        else { 
                            throw new LoggedException("Unsupported statistic : " + stat);
                        }
		    }
			
		}
	    }
	}

	return result;
	
    }


    /**
     * Gets the subset of array data[i][j][k] that contains data.
     *
     * @return the data vector for voxel i,j,k. Prepends zeros, if there are any.
     *
     */
    private double[] valuesAt(double[][][][] vol, int i, int j, int k) {

	if (vectorLengths[i][j][k] == 0) {
	    return new double[zeros[i][j][k]];
	}

	double[] values = new double[zeros[i][j][k] + vectorLengths[i][j][k]];

	System.arraycopy(vol[i][j][k], 0, values, zeros[i][j][k], vectorLengths[i][j][k]);
	
	return values;
    }

 
    private void growDataVector(int i, int j, int k, double factor) {
	
	int oldLength = data[i][j][k].length;
	
	int newLength = (int)(oldLength * factor);
	
	double[] replacement = new double[newLength];
	
	System.arraycopy(data[i][j][k], 0, replacement, 0, oldLength);
	
	data[i][j][k] = replacement;
	
    }



    public int[] getDataDims() {
	return new int[] {dataDims[0], dataDims[1], dataDims[2]};
    }
    
    
    public double[] getVoxelDims() {
	return new double[] {voxelDims[0], voxelDims[1], voxelDims[2]};
    }
 
   
   
}

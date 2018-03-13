package tractography;

import numerics.*;

/**
 * <dl>
 * <dt>Purpose: To provide linearly interpolated measurements using eight-neighbour interpolation at any 
 * point within the dataset.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> Provides an implementation of some methods common to multiple 
 * breeds of interpolator that interpolate over the voxel containing the point and the eight 
 * neighbouring voxels.
 *
 *
 * @version $Id: EightNeighbourInterpolator.java 554 2008-07-03 09:08:26Z ucacdxa $
 * @author  Philip Cook
 * 
 */
public abstract class EightNeighbourLinearInterpolator extends EightNeighbourInterpolator {


    
    /** Construct an interpolator.
     *
     */
    public EightNeighbourLinearInterpolator( int xDataSize, int yDataSize, int zDataSize, 
                                             double xVoxelSize, double yVoxelSize, double zVoxelSize ) {
        super(xDataSize, yDataSize, zDataSize, xVoxelSize, yVoxelSize, zVoxelSize); 
    }
    
    /** Construct an interpolator.
     *
     */
    public EightNeighbourLinearInterpolator( int[] dataDims, double[] voxelDims ) { 
        super(dataDims, voxelDims);
    }



   /** 
     * Gets the interpolated vector at some point. Each element is 
     * interpolated independently. Exactly what gets interpolated is determined by a call
     * to <code>getVoxelDataAsArray</code>.
     * 
     * 
     * @param point the point in mm to interpolate at.
     * @return the interpolated vector.
     *
     */
    public final double[] interpolateVectorElements(Point3D point, Vector3D previousDir) {

	double[] interpFraction = new double[8];
	int[] dims = new int[6];
        
	// get the interpolation parameters
	int inVoxel = setInterpolationVoxels(point, interpFraction, dims);
	
        double[][] vectors = new double[8][];

        for (int i = 0; i < 8; i++) {

            int x = dims[i / 4];
            int y = dims[2 + ((i / 2) % 2)];
            int z = dims[4 + (i % 2)];
            
            vectors[i] = getVoxelDataAsArray(x,y,z, previousDir);
        }

        for (int i = 0; i < 8; i++) {

            if (vectors[i] == null) { // happens if one of the voxels is background
                vectors[i] = vectors[inVoxel];
            }

        }
        

        int components = vectors[0].length;

        double[] interpolated = new double[components];
        

	for (int i = 0; i < components; i++) {
	
	    interpolated[i] = 
		
		vectors[0][i] * interpFraction[0] +
		vectors[1][i] * interpFraction[1] +
		vectors[2][i] * interpFraction[2] +
		vectors[3][i] * interpFraction[3] +
		vectors[4][i] * interpFraction[4] +
		vectors[5][i] * interpFraction[5] +
		vectors[6][i] * interpFraction[6] +
		vectors[7][i] * interpFraction[7]; 
	
	}

				
	return interpolated;
	
    }



    protected abstract double[] getVoxelDataAsArray(int x, int y, int z, Vector3D previousDir);
    
 
}

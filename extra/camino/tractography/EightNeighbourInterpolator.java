package tractography;

import numerics.Point3D;

/**
 * Superclass for multiple breeds of interpolator that interpolate 
 * over the voxel containing the point and the eight neighbouring voxels.
 * 
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class EightNeighbourInterpolator {


    protected final int xDataDim;
    protected final int yDataDim;
    protected final int zDataDim;

    protected final double xVoxelDim;
    protected final double yVoxelDim;
    protected final double zVoxelDim;


   /** Construct an interpolator.
     *
     */
    public EightNeighbourInterpolator( int xDataSize, int yDataSize, int zDataSize, 
					double xVoxelSize, double yVoxelSize, double zVoxelSize ) {
	xDataDim = xDataSize;
	yDataDim = yDataSize;
	zDataDim = zDataSize;

	xVoxelDim = xVoxelSize;
	yVoxelDim = yVoxelSize;
	zVoxelDim = zVoxelSize;
	
    }

   /** Construct an interpolator.
     *
     */
    public EightNeighbourInterpolator( int[] dataDims, double[] voxelDims ) { 
					
	xDataDim = dataDims[0];
	yDataDim = dataDims[1];
	zDataDim = dataDims[2];

	xVoxelDim = voxelDims[0];
	yVoxelDim = voxelDims[1];
	zVoxelDim = voxelDims[2];
	
    }




    /** Sets values relating to the position of a point, for use in interpolatation.
     * @param point the point in mm where interpolation is to occur.
     *
     * @param interpComponents the linear interpolation components, in order 
     * {000, 001, 010, 011, 100, 101, 110, 111} where the indexes are xyz, so 000 is 
     * to the rear lower left (x0, y0, z0), and 011 is the front upper left (x0, y1, z1).
     *
     * @param indices a container for the dataset indices of the interpolation voxels. 
     * In cases where the point is on the edge of the dataset, the interpolation will be 
     * between the edge voxel and itself. The parameter should be of the format int[6]:
     * <p><code>
     * indices[0] == x0
     * indices[1] == x1
     * indices[2] == y0
     * indices[3] == y1
     * indices[4] == z0
     * indices[5] == z1
     * </code>
     * <BR><BR>where<BR>
     * <code>[x0][y0][z0]</code> == the voxel to the front lower left of 
     * <code>point</code>.
     * <code>[x1][y0][z0]</code> == the voxel to the front lower right of 
     * <code>point</code>.
     * <code>[x0][y1][z0]</code> == the voxel to the front upper left of 
     * <code>point</code>.
     * <code>[x0][y0][z1]</code> == the voxel to the rear lower left of 
     * <code>point</code>.
     *
     * @return the voxel containing the point.
     */
    public final int setInterpolationVoxels(Point3D point, double[] interpComponents, int[] indices) { 
	
	double[] relXYZ = new double[3];
	
        // Compute the indices and offsets.
        interpIndicesAndOffsets(point, indices, relXYZ);

	
	interpComponents[0] = (1.0 - relXYZ[0]) * (1.0 - relXYZ[1]) * (1.0 - relXYZ[2]);
	interpComponents[1] = (1.0 - relXYZ[0]) * (1.0 - relXYZ[1]) * relXYZ[2]; 
	interpComponents[2] = (1.0 - relXYZ[0]) * relXYZ[1] * (1.0 - relXYZ[2]);
	interpComponents[3] = (1.0 - relXYZ[0]) * relXYZ[1] * relXYZ[2];
	interpComponents[4] = relXYZ[0] * (1.0 - relXYZ[1]) * (1.0 - relXYZ[2]);
	interpComponents[5] = relXYZ[0] *  (1.0 - relXYZ[1]) * relXYZ[2];
	interpComponents[6] = relXYZ[0] * relXYZ[1] * (1.0 - relXYZ[2]);
	interpComponents[7] = relXYZ[0] * relXYZ[1] * relXYZ[2];
       
	return 4 * (relXYZ[0] < 0.5 ? 0 : 1) + 2 * (relXYZ[1] < 0.5 ? 0 : 1) + (relXYZ[2] < 0.5 ? 0 : 1);
	
    }
    

    /**
     * As setInterpolationVoxels, but the components are now
     * for computing the derivative with respect to x.
     */ 
    public final int setInterpolationVoxelsDX(Point3D point, double[] interpComponents, int[] indices) { 
	
	double[] relXYZ = new double[3];
	
        // Compute the indices and offsets.
        interpIndicesAndOffsets(point, indices, relXYZ);

	
	interpComponents[0] = (-1.0/xVoxelDim) * (1.0 - relXYZ[1]) * (1.0 - relXYZ[2]);
	interpComponents[1] = (-1.0/xVoxelDim) * (1.0 - relXYZ[1]) * relXYZ[2]; 
	interpComponents[2] = (-1.0/xVoxelDim) * relXYZ[1] * (1.0 - relXYZ[2]);
	interpComponents[3] = (-1.0/xVoxelDim) * relXYZ[1] * relXYZ[2];
	interpComponents[4] = (1.0/xVoxelDim) * (1.0 - relXYZ[1]) * (1.0 - relXYZ[2]);
	interpComponents[5] = (1.0/xVoxelDim) *  (1.0 - relXYZ[1]) * relXYZ[2];
	interpComponents[6] = (1.0/xVoxelDim) * relXYZ[1] * (1.0 - relXYZ[2]);
	interpComponents[7] = (1.0/xVoxelDim) * relXYZ[1] * relXYZ[2];
       
	return 4 * (relXYZ[0] < 0.5 ? 0 : 1) + 2 * (relXYZ[1] < 0.5 ? 0 : 1) + (relXYZ[2] < 0.5 ? 0 : 1);
	
    }
    

    /**
     * As setInterpolationVoxels, but the components are now
     * for computing the derivative with respect to y.
     */ 
    public final int setInterpolationVoxelsDY(Point3D point, double[] interpComponents, int[] indices) { 
	
	double[] relXYZ = new double[3];
	
        // Compute the indices and offsets.
        interpIndicesAndOffsets(point, indices, relXYZ);

	
	interpComponents[0] = (1.0 - relXYZ[0]) * (-1.0/yVoxelDim) * (1.0 - relXYZ[2]);
	interpComponents[1] = (1.0 - relXYZ[0]) * (-1.0/yVoxelDim) * relXYZ[2]; 
	interpComponents[2] = (1.0 - relXYZ[0]) * (1.0/yVoxelDim) * (1.0 - relXYZ[2]);
	interpComponents[3] = (1.0 - relXYZ[0]) * (1.0/yVoxelDim) * relXYZ[2];
	interpComponents[4] = relXYZ[0] * (-1.0/yVoxelDim) * (1.0 - relXYZ[2]);
	interpComponents[5] = relXYZ[0] * (-1.0/yVoxelDim) * relXYZ[2];
	interpComponents[6] = relXYZ[0] * (1.0/yVoxelDim) * (1.0 - relXYZ[2]);
	interpComponents[7] = relXYZ[0] * (1.0/yVoxelDim) * relXYZ[2];
       
	return 4 * (relXYZ[0] < 0.5 ? 0 : 1) + 2 * (relXYZ[1] < 0.5 ? 0 : 1) + (relXYZ[2] < 0.5 ? 0 : 1);
	
    }
    

    /**
     * As setInterpolationVoxels, but the components are now
     * for computing the derivative with respect to z.
     */ 
    public final int setInterpolationVoxelsDZ(Point3D point, double[] interpComponents, int[] indices) { 
	
	double[] relXYZ = new double[3];
	
        // Compute the indices and offsets.
        interpIndicesAndOffsets(point, indices, relXYZ);

	
	interpComponents[0] = (1.0 - relXYZ[0]) * (1.0 - relXYZ[1]) * (-1.0/zVoxelDim);
	interpComponents[1] = (1.0 - relXYZ[0]) * (1.0 - relXYZ[1]) * (1.0/zVoxelDim); 
	interpComponents[2] = (1.0 - relXYZ[0]) * relXYZ[1] * (-1.0/zVoxelDim);
	interpComponents[3] = (1.0 - relXYZ[0]) * relXYZ[1] * (1.0/zVoxelDim);
	interpComponents[4] = relXYZ[0] * (1.0 - relXYZ[1]) * (-1.0/zVoxelDim);
	interpComponents[5] = relXYZ[0] *  (1.0 - relXYZ[1]) * (1.0/zVoxelDim);
	interpComponents[6] = relXYZ[0] * relXYZ[1] * (-1.0/zVoxelDim);
	interpComponents[7] = relXYZ[0] * relXYZ[1] * (1.0/zVoxelDim);
       
	return 4 * (relXYZ[0] < 0.5 ? 0 : 1) + 2 * (relXYZ[1] < 0.5 ? 0 : 1) + (relXYZ[2] < 0.5 ? 0 : 1);
	
    }
    

    /**
     * Finds voxel indices and offsets from each in preparation
     * for setting up the weightings for various interpolations.
     */
    private void interpIndicesAndOffsets(Point3D point, int[] indices, double[] relXYZ) {
	// Find lowest index of voxel surrounding point
	int[] datasetIndex = new int[3];

   	// interpolate from voxel centres
	datasetIndex[0] = (int)( (point.x / xVoxelDim) - 0.5);
	datasetIndex[1] = (int)( (point.y / yVoxelDim) - 0.5);
	datasetIndex[2] = (int)( (point.z / zVoxelDim) - 0.5);

	datasetIndex[0] = datasetIndex[0] < 0 ? 0 : datasetIndex[0]; 
	datasetIndex[1] = datasetIndex[1] < 0 ? 0 : datasetIndex[1];
	datasetIndex[2] = datasetIndex[2] < 0 ? 0 : datasetIndex[2];
	
	// Need to know position relative to this voxel
	// IN VOXELS
	relXYZ[0] = (point.x / xVoxelDim) - 0.5 - datasetIndex[0];
	relXYZ[1] = (point.y / yVoxelDim) - 0.5 - datasetIndex[1];
	relXYZ[2] = (point.z / zVoxelDim) - 0.5 - datasetIndex[2];
	
	// Must make sure we don't go out of bounds
	// If we are along boundaries, interpolate boundary voxel with itself
	
	indices[0] = datasetIndex[0];
	indices[2] = datasetIndex[1];
	indices[4] = datasetIndex[2];
	
	
	// Catch special cases of points along boundaries
	
	// relative position of point to leftmost-lowermost voxel is negative only on boundaries
	// (elsewhere it is always a positive distance from a voxel to the lower left)
	// Example: relXYZ[0] == -0.2. In this case, the point is 0.2 voxels to the left of voxel 
	// (0,y,z). Negate relXYZ (interpolating voxel with itself, so relXYZ is irrelevant in this 
	// dimension)
	if (relXYZ[0] < 0) {
	    indices[1] = 0;
	    relXYZ[0] = -1.0 * relXYZ[0];
	}
	else if (indices[0] == xDataDim - 1) {
	    indices[1] = indices[0];
	}
	else {
	    indices[1] = indices[0] + 1;
	}
	
	if (relXYZ[1] < 0) {
	    indices[3] = 0;
	    relXYZ[1] = -1.0 * relXYZ[1];
	}
	else if (indices[2] == yDataDim - 1) {
	    indices[3] = indices[2];
	}
	else {
	    indices[3] = indices[2] + 1;
	}
	
	if (relXYZ[2] < 0) {
	    indices[5] = 0;
	    relXYZ[2] = -1.0 * relXYZ[2];
	}
	else if (indices[4] == zDataDim - 1) {
	    indices[5] = indices[4];
	}
	else {
	    indices[5] = indices[4] + 1;
	}
    }

    
 

}

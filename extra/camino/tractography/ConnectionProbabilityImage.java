package tractography;

import numerics.*;

import java.util.*;

/**
 * Connection probability image.
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class ConnectionProbabilityImage {


    // voxel dims of seed space
    private final double xVoxelDim;
    private final double yVoxelDim;
    private final double zVoxelDim;


    // dimensions of seed space
    private final int xDataDim;
    private final int yDataDim;
    private final int zDataDim;

    private final int[][][] streamlineCounts;


    // total number of streamlines used to make this image
    private int totalStreamlines = 0;


    /**
     * Initializes the filter with the dimensions of the seed space.
     *
     */
    public ConnectionProbabilityImage(int[] dataDims, double[] voxelDims) {

        this(dataDims[0], dataDims[1], dataDims[2], voxelDims[0], voxelDims[1], voxelDims[2]);
    }

    /**
     * Initializes the filter with the dimensions of the seed space.
     *
     */
    public ConnectionProbabilityImage(int xDataDim, int yDataDim, int zDataDim,
				double xVoxelDim, double yVoxelDim, double zVoxelDim) {

	this.xDataDim = xDataDim;
	this.yDataDim = yDataDim;
	this.zDataDim = zDataDim;

	this.xVoxelDim = xVoxelDim;
	this.yVoxelDim = yVoxelDim;
	this.zVoxelDim = zVoxelDim;

	streamlineCounts = new int[xDataDim][yDataDim][zDataDim];

    }

 
    /**
     * Add some streamlines to the image.
     *
     */
    public final void processTracts(TractCollection tc) {

	
	for (int tCounter = 0; tCounter < tc.numberOfTracts(); tCounter++) {

	    processTract(tc.getTract(tCounter));
	}
	    
    }

    
    /**
     * Add a single streamline to the image.
     *
     */
    public final void processTract(Tract t) {
	totalStreamlines++;
	    
	Voxel[] voxels = t.toVoxelList(xVoxelDim, yVoxelDim, zVoxelDim).getVoxels();

	int numVoxels = voxels.length;
	    
	// makes sure that degenerate voxels are processed together
	java.util.Arrays.sort(voxels);

	Voxel lastVoxel = voxels[0];

	int v = 1;

	streamlineCounts[voxels[0].x][voxels[0].y][voxels[0].z] += 1;

	while (v < numVoxels) {
		    
	    if (voxels[v].x != lastVoxel.x || voxels[v].y != lastVoxel.y || 
		voxels[v].z != lastVoxel.z) {
			
		streamlineCounts[voxels[v].x][voxels[v].y][voxels[v].z] += 1;
		lastVoxel = voxels[v];
	    }
	    else {
		// fibre has looped, is moving along a voxel boundary, or we have reached the 
		// end
		// would like to break here, but we can't be sure that the fibre is doing 
		// something stupid. It might be repeatedly crossing a boundary tangential 
		// to the fibre path, for example.
		// -----|
		//  b-> | ->
		// c----|d
		//  a-> | ->
		// -----|

		// if eigenvector a is (0.9999, 0.0001, 0.0) and b is (0.9999, -0.0001, 0.0)
		// then a path that enters a or b close to the boundary (c) might cross back
		// and forth between the voxels before reaching d
			
		// however we definitely don't want to count this fibre again
	    }


	    v++;
		    
                    

	} // end while

		
    }
		
    	


    /**
     * @return raw streamline counts.
     *
     */
    public double[][][] getStreamlineCounts() {

	double[][][] defCopy = new double[xDataDim][yDataDim][zDataDim];
	
	for (int i = 0; i < xDataDim; i++) {
	    for (int j = 0; j < yDataDim; j++) {
                for (int k = 0; k < zDataDim; k++) {
                    defCopy[i][j][k] = streamlineCounts[i][j][k];
                }
            }
        }

	return defCopy;
    }


    /**
     * @return streamline counts normalized by the total number of streamlines processed.
     *
     */
    public double[][][] getConnectionProbabilities() {
	double[][][] cp = new double[xDataDim][yDataDim][zDataDim];

	double norm = (double)(totalStreamlines);


        if (norm > 0) {
            for (int i = 0; i < xDataDim; i++) {
                for (int j = 0; j < yDataDim; j++) {
                    for (int k = 0; k < zDataDim; k++) {
                        cp[i][j][k] = ((double)streamlineCounts[i][j][k]) / norm;
                    }
                }
            }
        }

	return cp;

    }


    /**
     * @return the total number of streamlines processed by <code>processTracts</code>.
     *
     */
    public int totalStreamlines() {
	return totalStreamlines;
    }



}

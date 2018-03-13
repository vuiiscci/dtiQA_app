package tractography;

import numerics.*;

import java.util.Random;

/**
 * Provides interpolated tracking directions at any point within the dataset. 
 * Interpolates directions using a similar method to Behrens et al
 * (Magnetic Resonance in Medicine, 50:1077-1088, 2003). 
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class NeighbourChoiceInterpolator extends EightNeighbourInterpolator implements ImageInterpolator {
    
    private final double[] interpFraction = new double[8];
    private final int[] dims = new int[6];

    protected final TractographyImage image;

    protected final Random ran;

    private boolean[][][] background;



    public NeighbourChoiceInterpolator(TractographyImage image, Random r) {
	
	super( image.xDataDim(), image.yDataDim(),image.zDataDim(), 
	       image.xVoxelDim(), image.yVoxelDim(),image.zVoxelDim()
	       );

	this.image = image;
	ran = r;

        background = new boolean[xDataDim][yDataDim][zDataDim];
        
       	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
                    background[i][j][k] = (image.numberOfPDs(i,j,k) == 0);
                }
            }
        }

    }
    
  
    /** 
     * Get the tracking direction at some point. The direction will come from one of the 8 voxels 
     * surrounding the point. The probability of any voxel being chosen is the same as the trilinear 
     * interpolation fraction for that voxel. Note that successive calls with the same point may 
     * return different directions.
     *
     * @param point the point in mm to interpolate at. 
     * @return the list of tracking direction chosen from one of the neighbours.
     *
     */
    public Vector3D getTrackingDirection(Point3D point, Vector3D previousDirection) {
   
        int[] choice = getNeighborVoxelIndex(point);
	
        Vector3D[] voxelPDs = image.getPDs(choice[0], choice[1], choice[2], previousDirection);
        
        int vectorChoice = chooseVector(voxelPDs, previousDirection);

        if (voxelPDs[vectorChoice].dot(previousDirection) > 0.0) {
            return voxelPDs[vectorChoice];
        }
        else {
            return voxelPDs[vectorChoice].negated();
        }
        
      	
    }


    /**
     * Gets one of the eight voxels neighboring this point, with probabilities weighted by the position between the voxels.
     *
     */
    protected int[] getNeighborVoxelIndex(Point3D point) {
        
        setInterpolationVoxels(point, interpFraction, dims);
        
        double sumInterp = 0.0;
        
        // probability of choosing background voxel is zero
        for (int i = 0; i < 8; i++) {
            int x = dims[i / 4];
            int y = dims[2 + ((i / 2) % 2)];
            int z = dims[4 + (i % 2)];
            
            if (background[x][y][z]) {
                interpFraction[i] = 0.0;
            }
            
            sumInterp += interpFraction[i];
        }

        double random = ran.nextDouble() * sumInterp;
        
	double cumulSum = 0.0; // add up fractions as we go
        
	// loop over neighbours and see if random number is < cumulSum
	for (int i = 0; i < 8; i++) {

	    cumulSum = cumulSum + interpFraction[i];
            
	    // include a delta in case of precision errors 	    
	    if (random - cumulSum < 1E-9) {
		int x = dims[i / 4];
		int y = dims[2 + ((i / 2) % 2)];
		int z = dims[4 + (i % 2)];
                
                return new int[] {x, y, z};
                
            }
            
        }

        // should never get here
	throw new misc.LoggedException
	    ("Rejected all 8 neighbourhood voxels. cumulsum == " 
	     + cumulSum + ", interpFractions == {" + interpFraction[0] + ", " + interpFraction[1] + 
	     ", " + interpFraction[2] + ", " + interpFraction[3] + ", " + interpFraction[4] + ", " +  
	     interpFraction[5] + ", " + interpFraction[6] + ", " + interpFraction[7] + "}"
	     );

    }
    

    
    /** 
     * Get the initial tracking direction, given a pdIndex and a seed point.
     * 
     * @param pdIndex follow this pd for the first tracking step. 
     * @param direction if true, the direction will be the pd, if false, it will be the 
     * negated pd.
     * @return the tracking direction for this point. 
     * 
     */
    public Vector3D getTrackingDirection(Point3D point, int pdIndex, boolean direction) {
	
	int x = (int)(point.x / xVoxelDim);
	int y = (int)(point.y / yVoxelDim);
	int z = (int)(point.z / zVoxelDim);


        if (direction) {
            return getTrackingDirection(point,image.getPDs(x,y,z)[pdIndex]);
        }
        else {
            return getTrackingDirection(point,image.getPDs(x,y,z)[pdIndex].negated());
        }

    }



    /**
     * Chooses which vector to use for interpolation in this voxel, given the previous
     * direction. Will negate the chosen vector if necessary, so the dot product of the
     * previous direction and the returned vector is always non-negative.
     *
     */
    protected final int chooseVector(Vector3D[] vecs, Vector3D previousDir) {
	
	if (vecs.length == 0) {
	    return -1;
	}
        if (vecs.length == 1) {
            return 0;
	}
	
	int choice = 0;
	double maxDot = 0.0;
	
	for (int d = 0; d < vecs.length; d++) {
	    double dot = Math.abs(vecs[d].dot(previousDir));
	    if (dot > maxDot) {
		maxDot = dot;
		choice = d;
	    }
	}

        return choice;

    }


    public TractographyImage getImage() {
        return image;
    }

     
}

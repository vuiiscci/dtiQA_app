package tractography;

import numerics.*;

/**
 * 
 * Superclass for NN interpolators. 
 *
 * @author Philip Cook
 * @version $Id$
 */
public class NearestNeighbourInterpolator implements ImageInterpolator {

    private final TractographyImage image;

    protected final double xVoxelDim;
    protected final double yVoxelDim;
    protected final double zVoxelDim;


    public NearestNeighbourInterpolator(TractographyImage im) {
        image = im;
        
        xVoxelDim = image.xVoxelDim();
        yVoxelDim = image.yVoxelDim();
        zVoxelDim = image.zVoxelDim();
    }

 
    public final Vector3D getTrackingDirection(Point3D point, Vector3D previousDir) {
        
	int i = (int)(point.x / xVoxelDim);
	int j = (int)(point.y / yVoxelDim);
	int k = (int)(point.z / zVoxelDim);

        Vector3D[] voxelPDs = image.getPDs(i,j,k, previousDir);
	
	double maxProd = 0.0;

	int choice = chooseVector(voxelPDs, previousDir);

        if (voxelPDs[choice].dot(previousDir) > 0.0) {
            return voxelPDs[choice];
        }
        else {
            return voxelPDs[choice].negated();
        }

    }


    public final Vector3D getTrackingDirection(Point3D point, int pdIndex, boolean direction) {
              
	int i = (int)(point.x / xVoxelDim);
	int j = (int)(point.y / yVoxelDim);
	int k = (int)(point.z / zVoxelDim);

        Vector3D[] voxelPDs = image.getPDs(i,j,k);

	Vector3D[] pds = image.getPDs(i,j,k);
        
        if (direction) {
	    return pds[pdIndex];
	}
        else {
            return pds[pdIndex].negated();
        }

    }



    /**
     * Chooses which vector to use for interpolation in this voxel, given the previous
     * direction.
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

package tractography;

import numerics.*;

/**
 * Provides interpolated measurements of the principal direction
 * at any point within the dataset.
 *
 *
 * @version $Id: VectorLinearInterpolator.java 1147 2012-11-23 23:58:59Z ucacpco $
 * @author Philip Cook
 * 
 */
public final class DWI_LinearInterpolator extends EightNeighbourLinearInterpolator implements ImageInterpolator {


    private DWI_TractographyImage image;

    
    /** Construct an interpolator.
     * @param data the dataset to use for interpolation.
     */
    public DWI_LinearInterpolator(DWI_TractographyImage data) {
	
	super( data.xDataDim(), data.yDataDim(), data.zDataDim(), 
	       data.xVoxelDim(), data.yVoxelDim(), data.zVoxelDim()
	       );

	image = data;

    }
    


    /**
     * Returns null if voxel is background
     */
    protected double[] getVoxelDataAsArray(int i, int j, int k, Vector3D previousDir) {
        
        return image.numberOfPDs(i,j,k) > 0 ? image.getData(i,j,k) : null;
    } 


    /**
     * Chooses which vector to return given a list of choices and the previous direction
     *
     */
    private Vector3D chooseVector(Vector3D[] vecs, Vector3D previousDir) {
	
	if (vecs.length == 0) {
	    return null;
	}
        if (vecs.length == 1) {
	    if (vecs[0].dot(previousDir) > 0.0) {
		return vecs[0];
	    }
	    else {
		return vecs[0].negated();
	    }
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

	if (vecs[choice].dot(previousDir) > 0.0) {
	    return vecs[choice];
	}
	else {
	    return vecs[choice].negated();
	}
    }


  
    public Vector3D getTrackingDirection(Point3D point, Vector3D previousDir) {


	int i = (int)( (point.x / xVoxelDim) );
	int j = (int)( (point.y / yVoxelDim) );
	int k = (int)( (point.z / zVoxelDim) );

        double[] interpolated = interpolateVectorElements(point, previousDir);
	
	return chooseVector( image.getPDs(interpolated, image.numberOfPDs(i,j,k)), previousDir );

    }



    /** 
     * Get the initial tracking direction, given a pdIndex and a seed point.
     * 
     * @param direction if true, the direction will be the PD, if false, it will be the negated PD.
     * @return the tracking direction for this point. 
     * 
     */
    public Vector3D getTrackingDirection(Point3D point, int pdIndex, boolean direction) {

	int i = (int)( (point.x / xVoxelDim) );
	int j = (int)( (point.y / yVoxelDim) );
	int k = (int)( (point.z / zVoxelDim) );

	Vector3D[] pds = image.getPDs(i,j,k);

	if (direction) {
	    return getTrackingDirection(point, pds[pdIndex]);
	}
	else {
	    return getTrackingDirection(point, pds[pdIndex].negated());
	}

    }

    
    public TractographyImage getImage() {
        return image;
    }


}

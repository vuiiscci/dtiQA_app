package tractography;

import misc.DT;
import numerics.*;

/**
 * Provides interpolated measurements of the diffusion tensor 
 * at any point within the dataset.
 *
 *
 * @version $Id$
 * @author Philip Cook
 * 
 */
public final class DT_LinearInterpolator extends EightNeighbourLinearInterpolator 
    implements ImageInterpolator, TensorInterpolator {


    private TensorTractographyImage image;

    
    /** Construct an interpolator.
     * @param data the dataset to use for interpolation.
     */
    public DT_LinearInterpolator(TensorTractographyImage data) {
	
	super( data.xDataDim(), data.yDataDim(), data.zDataDim(), 
	       data.xVoxelDim(), data.yVoxelDim(), data.zVoxelDim()
	       );

	image = data;

    }
    
  
    /** 
     * Gets the interpolated tensor at some point. Each tensor component is 
     * interpolated independently. If there are multiple tensors in a voxel, 
     * the one aligned closest to the current tracking direction is used. If any 
     * neighbouring voxels are classified as background, their tensor is replaced by the
     * tensor in the voxel containing the point.
     * 
     * 
     * @param point the point in mm to interpolate at.
     * @return the interpolated tensor.
     *
     */
    public DT getDT(Point3D point, Vector3D previousDir) {

        double[] interpolated = interpolateVectorElements(point, previousDir);
				
	return new DT(interpolated);
	
    }


    protected double[] getVoxelDataAsArray(int i, int j, int k, Vector3D previousDir) {
        
        DT d = chooseTensor(i,j,k, previousDir);
        
        return d != null ? d.getComponents() : null;
    } 
    

    private DT chooseTensor(int i, int j, int k, Vector3D previousDir) {
	
	DT[] dts = image.getDTs(i,j,k);


        if (dts.length == 0) {
	    return null;
	}
        if (dts.length == 1) {
	    return dts[0];
	}
	
	Vector3D[] pds = image.getPDs(i,j,k);
	
	int choice = 0;
	double maxDot = 0.0;
	
	for (int d = 0; d < pds.length; d++) {
	    double dot = Math.abs(pds[d].dot(previousDir));
	    if (dot > maxDot) {
		maxDot = dot;
		choice = d;
	    }
	}
	
	return dts[choice];
    }

  
    public Vector3D getTrackingDirection(Point3D point, Vector3D previousDirection) {
	DT dt = getDT(point, previousDirection);
	double[][] seig = dt.sortedEigenSystem();
	
	Vector3D e1 = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);

	if (e1.dot(previousDirection) > 0.0) {
	    return e1;
	}
	else {
	    return e1.negated();
	}
    }


    /** 
     * Get the initial tracking direction, given a pdIndex and a seed point. 
     * 
     * @param direction toggles whether we go "forward" or "backward" from the seed point. Used to 
     * initiate tracking in opposite directions from the same seed.
     *
     * @return the tracking direction for this point. 
     * 
     */
    public Vector3D getTrackingDirection(Point3D point, int pdIndex, boolean direction) {

	int x = (int)( (point.x / xVoxelDim) );
	int y = (int)( (point.y / yVoxelDim) );
	int z = (int)( (point.z / zVoxelDim) );

	Vector3D[] pds = image.getPDs(x,y,z);

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

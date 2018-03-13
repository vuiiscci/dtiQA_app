package tractography;

import misc.*;
import numerics.*;


/**
 * 
 * Superclass for NN interpolators. 
 *
 * @author Philip Cook
 * @version $Id$
 */
public class DT_NN_Interpolator extends NearestNeighbourInterpolator implements TensorInterpolator {

    
    private TensorTractographyImage tensorImage;
    

    public DT_NN_Interpolator(TensorTractographyImage im) {
        super(im);

        tensorImage = im;
        
    }


    public DT getDT(Point3D point, Vector3D previousDirection) {

	int i = (int)( point.x / xVoxelDim );
	int j = (int)( point.y / yVoxelDim );
	int k = (int)( point.z / zVoxelDim );

        DT[] dts = tensorImage.getDTs(i,j,k);

        if (dts.length == 0) {
	    return null;
	}
        if (dts.length == 1) {
	    return dts[0];
	}
	
	int choice = chooseVector(tensorImage.getPDs(i,j,k), previousDirection);
	
	return dts[choice];
    }



}

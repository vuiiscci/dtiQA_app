package tractography;

import misc.*;
import numerics.*;

import java.util.Random;

/**
 * 
 * Superclass for NN interpolators. 
 *
 * @author Philip Cook
 * @version $Id$
 */
public class DT_NC_Interpolator extends NeighbourChoiceInterpolator implements TensorInterpolator {

    private final TensorTractographyImage tensorImage;


    public DT_NC_Interpolator(TensorTractographyImage im, Random r) {
        super(im, r);

        tensorImage = im;
    }


    public DT getDT(Point3D point, Vector3D previousDirection) {

        int[] voxel = getNeighborVoxelIndex(point);

        DT[] dts = tensorImage.getDTs(voxel[0], voxel[1], voxel[2]);

        if (dts.length == 0) {
	    return null;
	}
        if (dts.length == 1) {
	    return dts[0];
	}

	int i = (int)( (point.x / xVoxelDim) );
	int j = (int)( (point.y / yVoxelDim) );
	int k = (int)( (point.z / zVoxelDim) );
	
	int choice = chooseVector(tensorImage.getPDs(voxel[0], voxel[1], voxel[2]), previousDirection);
	
	return dts[choice];
    }



}

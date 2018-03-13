package tractography;

import misc.DT;
import numerics.*;


/**
 * 
 * Does tensor deflection. This is more regularization than interpolation, though
 * you may additionally interpolate the vector field. This class implements
 * ImageInterpolator because that allows it to work with the trackers.
 * 
 * e1, the previous tracking direction, and the input tracking direction multiplied
 * by the tensor:
 * <br><br>
 * v_{out}^* = f * e_1 + (1 - f)[(1 - g) * v_{in} + g * normalize(D * v_{in})]
 * <br><br>
 * v_{out} = v_{out}^* / | v_{out}^* .
 * <br><br>
 * This class does not interpolate the tensor field in any way.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class TendInterpolator implements ImageInterpolator {

    private final TensorTractographyImage image;

    private final TensorInterpolator interpolator;

    // relative weighting of previous direction and tend term
    // between 0.0 (ignore tend) and 1.0 (ignore previous direction)
    private double g = 0.0;

    private double xVoxelDim;
    private double yVoxelDim;
    private double zVoxelDim;

    private int xDataDim;
    private int yDataDim;
    private int zDataDim;

    // f is the extent to which we trust the 
    // local tensor e_1. Between 0 (not at all) and
    // 1.0 (totally). Even if this is set to zero, the 
    // local DT has some influence through the tend term.
    private double[][][] f;


    /**
     * Default constructor, doesn't use tend at all.
     *
     *
     */
    public TendInterpolator(TensorTractographyImage image) {
	this(image, new DT_NN_Interpolator(image), 1.0, 0.0);
    }


    /**
     * Weight the tend term by a constant.
     * 
     * @param tendF between 0.0 (ignore e_1) and 1.0 (ignore tend).
     * @param tendG between 0.0 (ignore deflection, tend term is previous direction) and 1.0 (trust deflection, 
     * tend term is D * v_{in}).
     *
     */
    public TendInterpolator(TensorTractographyImage image, double tendF, double tendG) {
        this(image, new DT_NN_Interpolator(image), tendF, tendG);
    }



    /**
     * Weight the tend term by a constant, and additionally interpolate the underlying tensor field.
     * 
     * @param tendF between 0.0 (ignore e_1) and 1.0 (ignore tend).
     * @param tendG between 0.0 (ignore deflection, tend term is previous direction) and 1.0 (trust deflection, 
     * tend term is D * v_{in}).
     *
     */
    public TendInterpolator(TensorTractographyImage image, TensorInterpolator interp, double tendF, double tendG) {

	this.image = image;

        interpolator = interp;

	xDataDim = image.xDataDim();
	yDataDim = image.yDataDim();
	zDataDim = image.zDataDim();

	xVoxelDim = image.xVoxelDim();
	yVoxelDim = image.yVoxelDim();
	zVoxelDim = image.zVoxelDim();
	
	f = new double[xDataDim][yDataDim][zDataDim];

	g = tendG;

	if (!(tendF >= 0.0 && tendF <= 1.0)) {
	    throw new misc.LoggedException("Tend parameters must be between 0 and 1. f is " + tendF);
	} 

	if (!(g >= 0.0 && g <= 1.0)) {
	    throw new misc.LoggedException("Tend parameters must be between 0 and 1. g is " + g);
	} 


	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    f[i][j][k] = tendF;
		}
	    }
	}

    }


    /**
     * Weight the tend term independently in each voxel, and additionally interpolate the underlying tensor field.
     *
     * @param tendF an image of tend parameters, of the same 
     * dimension as the tensor image and normalized between 0 and 1.
     *
     * @param tendG between 0.0 (ignore deflection, tend term is previous direction) and 1.0 (trust deflection, 
     * tend term is D * v_{in}).
     *
     */
    public TendInterpolator(TensorTractographyImage image, TensorInterpolator interp, double[][][] tendF, double tendG) {

	this.image = image;
    
        interpolator = interp;

	xDataDim = image.xDataDim();
	yDataDim = image.yDataDim();
	zDataDim = image.zDataDim();

	xVoxelDim = image.xVoxelDim();
	yVoxelDim = image.yVoxelDim();
	zVoxelDim = image.zVoxelDim();
	
	f = tendF;
	g = tendG;

	if (!(g >= 0.0 && g <= 1.0)) {
	    throw new misc.LoggedException("Tend parameters must be between 0 and 1. g is " + g);
	} 

	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    if (!(f[i][j][k] >= 0.0 && f[i][j][k] <= 1.0)) {
			throw new misc.LoggedException("Tend parameters must be between 0 " + 
						       " and 1. Value of f at voxel " + i + " " + j +
						       " " + k + " is " + f[i][j][k]);
		    }
		}
	    }
	}


    }

 
    /**
     * @return the tend tracking direction for this point. The returned value is the same
     * for all points within the same voxel, given the same previous direction. Ultimately
     * this will be replaced with an option to interpolate the tensors somehow.
     *
     */
    public Vector3D getTrackingDirection(Point3D point, Vector3D previousDirection) {

	int i = (int)( (point.x / xVoxelDim) );
	int j = (int)( (point.y / yVoxelDim) );
	int k = (int)( (point.z / zVoxelDim) );

        // The local fiber orientation from the nearest tensor, or interpolated tensor
	Vector3D e1 = interpolator.getTrackingDirection(point, previousDirection);

	DT dt = interpolator.getDT(point, previousDirection);

	double[] prev = new double[] {previousDirection.x, previousDirection.y,	previousDirection.z};

	double[] tend = dt.multiply(prev);

	Vector3D tendVec = new Vector3D(tend).normalized();

	Vector3D gTerm = previousDirection.scaled(1.0 - g).plus(tendVec.scaled(g));

	Vector3D vout = e1.scaled(f[i][j][k]).plus(gTerm.scaled(1.0 - f[i][j][k]));

	return vout.normalized();
    }



    /** 
     * Get the initial tracking direction, given a pdIndex and a seed point.
     * 
     * @param direction if true, the direction will be the PD, if false, it will be the negated PD.
     * @return the tracking direction for this point. 
     * 
     */
    public Vector3D getTrackingDirection(Point3D point, int pdIndex, boolean direction) {

	int x = (int)( (point.x / xVoxelDim) );
	int y = (int)( (point.y / yVoxelDim) );
	int z = (int)( (point.z / zVoxelDim) );

	Vector3D[] pds = image.getPDs(x,y,z);

	if (direction) {
	    return pds[pdIndex];
	}
	else {
	    return pds[pdIndex].negated();
	}

    }


    public TractographyImage getImage() {
        return image;
    }


}

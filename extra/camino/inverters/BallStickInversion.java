package inverters;

import imaging.*;
import misc.*;
import numerics.*;

/**
 * Fits the ball and stick partial volume model to diffusion-weighted data. 
 *
 * @see data.BallStick
 *
 * @author  Philip Cook
 * @version $Id$
 * 
 */
public class BallStickInversion extends DiffusionInversion {
	
    private BallStickFitter fitter;

    private LinearDT_Inversion linearInv;

    public BallStickInversion(DW_Scheme ip) {

	this.ip = ip;

	try {
	    fitter = new BallStickFitter(ip);
	}
	catch(Exception e) {
	    throw new LoggedException(e);
	}

	linearInv = new LinearDT_Inversion(ip);
    }


    /**
     * Does the inversion.
     *
     * @param data The MRI data. 
     *
     * @return {errorCode, ln(S_0), d, f, x, y, z}. 
     */
    public double[] invert(double[] data) {

	double exitCode = 0.0;

	// quick check for bad data
	for (int i = 0; i < data.length; i++) {
	    if (!(data[i] > 0.0)) {
		exitCode = 6.0;
	    }
	}

	try {
	    fitter.newDepVals(data);
	    fitter.minimise();
	}
	catch(Exception e) {
	    System.err.println(e);
	    
	    // If the fitter fails, the exit code is 2.

	    // return initial guess of parameters from DT
	    exitCode = 2.0;

	    double[] res = new double[7];
	    
	    double[] comps = linearInv.invert(data);
	    
	    DT dt = new DT(comps[2], comps[3], comps[4], comps[5], comps[6], comps[7]);

	    double[][] seig = dt.sortedEigenSystem();

	    double fa = dt.fa();
	    
	    res[0] = exitCode;
	    res[1] = comps[1];
	    res[2] = dt.trace() / 3.0;
	    res[3] = dt.fa();
	    res[4] = seig[1][0];
	    res[5] = seig[2][0];
	    res[6] = seig[3][0];

	    return res;

	}

	// [0.0, ln(S_0), d, f, theta, phi]
	double[] params = fitter.getParameters();
	
	// theta, phi give rotation of X-axis to the stick axis
	// Slightly odd formulation copied from Behrens et al 2003
	
	Vector3D v = Rotations.rotateVector(Rotations.X_AXIS, Rotations.Y_AXIS, params[4]);

	v = Rotations.rotateVector(v, Rotations.Z_AXIS, params[5]);

	// Construct the array to return.  The first element is
	// the exit code. 
	double[] res = new double[7];
	res[0] = exitCode;
	res[1] = params[1];
	res[2] = params[2];
	res[3] = params[3];
	res[4] = v.x;
	res[5] = v.y;
	res[6] = v.z;


	return res;

    }


    /**
     * Specifies the number of elements of the output array from the inversion.
     * 
     * @return The number of elements of the output array.
     */
    public int itemsPerVoxel() {
	return 7;
    }

}

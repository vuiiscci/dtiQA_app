package inverters;

import optimizers.*;
import misc.*;
import java.util.logging.Logger;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of diffusion tensor to DW-MR data.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of <code>MarquardtChiSqFitter</code> to provide
 * a Levenburg-Marquardt algorithm for fitting a diffusion tensor to DW-MR
 * measurements. The fitter fits the model to the normalized data directly
 * without taking logs, so that the noise statistics are less corrupted. The
 * diffusion tensor is constrained to be positive definite by optimizing the
 * parameters of its Cholesky decomposition.
 * 
 * </dl>
 * 
 * @version $Id$
 * @author Danny Alexander
 *  
 */
public class DiffTensorFitter extends TensorModelFitter {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.inversions.DiffTensorFitter");


    /**
     * Default constructor.
     */
    public DiffTensorFitter() {
    }

    /**
     * The number of unweighted acquisitions that are made (nob0s) is required to
     * estimate the noise levels of each data item.
     * 
     * @param indepVals
     *            The matrix of non-zero gradient directions without the zero measurements.
     * 
     * @param depVals
     *            The normalized measurements.
     * 
     * @param bValues
     *            The array of b-values.
     * 
     * @param nob0s
     *            The number of b=0 acquisitions.
     */
    public DiffTensorFitter(double[][] indepVals, double[] depVals, double[] bValues,
            int nob0s) throws MarquardtMinimiserException {
        noParams = 6;
        initialize(indepVals, depVals, bValues, nob0s);
    }

    /**
     * Sets the parameters to the starting point for the
     * optimization. We choose an isotropic diffusion tensor with
     * typical trace for brain data.
     */
    protected void initAs() {

        //Initial values of the parameters and the standard
        //deviations of each data point can be set here.

        //Initialise the off diagonal tensor elements to zero.
        a[2] = a[3] = a[5] = 0.0;

        //Initialise on diagonals to non-zero.
        double traceD = 21.0E-10;
        a[1] = Math.sqrt(traceD / 3.0);
        a[4] = Math.sqrt(traceD / 3.0);
        a[6] = Math.sqrt(traceD / 3.0);

    }

    /**
     * Sets the parameters to an alternative starting point.
     *
     * @param startDT An array with the format output by all DT
     * inversions: [exitcode lnA(0) Dxx Dxy Dxz Dyy Dyz Dzz].  The
     * optimization begins from the diffusion tensor it specifies.
     */
    public void setStartPoint(double[] dtParams) {

        DT startDT = new DT(dtParams[2], dtParams[3], dtParams[4], dtParams[5], dtParams[6], dtParams[7]);
        double[] u = getCholParams(startDT);

        // Look out for failure of the Cholesky decomposition, which can
        // happen if the linear fit has a non-positive eigenvalue.
        boolean cholFailure = false;

        for(int i=0; i<6; i++) {
            a[i+1] = u[i];
            if(u[i] == 0.0)
                cholFailure = true;
        }
        
        if(cholFailure) {

	    // If L1 and L2 are positive, maybe we can recover a tensor that is somewhat sensible

            double[][] seig = startDT.sortedEigenSystem();

            DT modDT = null;

            String ldtEig = "[ " + seig[0][0] + " , " + seig[0][1] + " , " + seig[0][2] + " ]";
            
            if (seig[0][0] < 0.0) {
		// All eigenvalues negative means something seriously wrong, so start from isotropic tensor
		// MD might be way off; probably doesn't matter since this voxel is most likely worthless
		// Can't guess MD more intelligently without knowing units of b-values
		double md = Math.abs(seig[0][0]) / 2.0;
		modDT = new DT(md, 0.0, 0.0, md, 0.0, md);

            }
	    else if (seig[0][1] < 0.0) {

		// probably means L1 >> L2 > L3

		double rd = seig[0][0] / 4.0;

		seig[0][1] = rd / 2.0;
		seig[0][2] = rd / 2.0;

		modDT = DT.dtFromEig(seig);
		
	    }
            else {
                // Could be L1 \approx L2 >> L3 or L1 >> L2 > L3
		seig[0][2] = Math.min(seig[0][1] * 0.9, Math.abs(seig[0][2])); 

		modDT = DT.dtFromEig(seig);
                
            }

            logger.fine("Chol failure.  Negative evals in linear fitted DT " + ldtEig + " . Using modified DT for starting point. (traceD, fa) = (" + modDT.trace() + " , " + modDT.fa() + ")" );
            
            u = getCholParams(modDT);

            for(int i=0; i<6; i++) {
                a[i+1] = u[i];
            }
       } 
    }



    /**
     * Compute the value of the model at the specified sample point with
     * parameter estimates in atry.
     * 
     * @param atry
     *            Parameter values to try.
     * 
     * @param i
     *            The index of the sample point.
     * 
     * @return The value of the model at the sample point.
     */
    protected double yfit(double[] atry, int i) {

        // Compute contribution from tensor.
        DT d = getDT_Chol(atry, 1);
        double[] g = getG(x, i);
        double yVal = Math.exp(-bValues[i] * d.contractBy(g));

        return yVal;
    }

    /**
     * Overrides the default to compute the derivatives analytically.
     */
    protected double[] dydas(double[] atry, int i) {

        double[] derivs = new double[ma + 1];

        double[] g = getG(x, i);
        double gDotRow1 = g[0] * atry[1] + g[1] * atry[2] + g[2] * atry[3];
        double gDotRow2 = g[1] * atry[4] + g[2] * atry[5];
        double gDotRow3 = g[2] * atry[6];

        // Get the tensor from the parameter values.
        DT d = getDT_Chol(atry, 1);

        // Estimate the measurement from the tensor
        double estim = Math.exp(-bValues[i] * d.contractBy(g));

        insertCholDerivs(derivs, bValues[i], 1, g, gDotRow1, gDotRow2, gDotRow3, 1.0, estim);

        //Check the derivatives against numerical derivatives.
        //   	double[] numDs2 = dydasNumerical(atry, i);
        //   	System.out.println("i = " + i);
        //   	for(int j=0; j<derivs.length; j++) {
        //   	   System.out.print(derivs[j] + " " + numDs2[j] + " : ");
        //  	}
        //   	System.out.println();

        return derivs;
    }

    /**
     * Returns the diffusion tensor represented by the current
     * parameter values.
     * 
     * @return The diffusion tensor specified by the parameter settings.
     */
    public DT getDiffTensor() {
        return getDT_Chol(a, 1);
    }


    

}

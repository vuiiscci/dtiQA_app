package inverters;

import optimizers.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of diffusion tensor to DW-MR data.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtFitter to provide
 * a Levenburg-Marquardt algorithm for fitting a diffusion tensor to DW-MR
 * measurements. The fitter fits the model to the normalized data directly
 * without taking logs, so that the noise statistics are less corrupted. The
 * diffusion tensor is unconstrained so may not end up positive definite.
 * 
 * </dl>
 * 
 * @version $Id: DiffTensorUnConFitter.java,v 1.2 2005/08/18 11:07:43 ucacmgh
 *          Exp $
 * @author Danny Alexander
 *  
 */
public class DiffTensorUnConFitter extends DiffTensorFitter {

    /**
     * Default constructor.
     */
    public DiffTensorUnConFitter() {
    }

    /**
     * The constructor requires a list of independent values (indepVals, the
     * gradient directions) and associated dependent values (depVals, the data). The
     * number of unweighted acquisitions that are made (nob0s) is required to
     * estimate the noise levels of each data item. 
     * 
     * @param indepVals
     *            The matrix of gradient directions without the zero gradients.
     * 
     * @param depVals
     *            The normalized measurements.
     * 
     * @param bValues
     *            The array of diffusion times.
     * 
     * @param nob0s
     *            The number of q=0 acquisitions.
     */
    public DiffTensorUnConFitter(double[][] indepVals, double[] depVals,
            double[] bValues, int nob0s) throws MarquardtMinimiserException {
        noParams = 6;
        initialize(indepVals, depVals, bValues, nob0s);
    }

    /**
     * Sets the parameters to the starting point for the optimization. We choose
     * an isotropic diffusion tensor with typical trace for brain data.
     */
    protected void initAs() {

        //Initial values of the parameters and the standard
        //deviations of each data point can be set here.

        //Initialise the off diagonal tensor elements to zero.
        a[2] = a[3] = a[5] = 0.0;

        //Initialise on diagonals to non-zero.
        double traceD = 21.0E-10;
        a[1] = traceD / 3.0;
        a[4] = traceD / 3.0;
        a[6] = traceD / 3.0;

    }


    public void setStartPoint(double[] dtParams) {

        for(int i=0; i<6; i++) {
            a[i+1] = dtParams[i+2];
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
        DT d = new DT(atry[1], atry[2], atry[3], atry[4], atry[5], atry[6]);
        double[] g = getG(x, i);
        double yVal = Math.exp(-bValues[i] * d.contractBy(g));

        return yVal;
    }

    /**
     * Overrides the default to compute the derivatives analytically.
     */
    protected double[] dydas(double[] atry, int i) {

        double[] derivs = new double[ma + 1];

        DT d = new DT(atry[1], atry[2], atry[3], atry[4], atry[5], atry[6]);
        double[] g = getG(x, i);
        double estim = Math.exp(-bValues[i] * d.contractBy(g));

        derivs[1] = -bValues[i] * g[0] * g[0] * estim;
        derivs[2] = -2.0 * bValues[i] * g[0] * g[1] * estim;
        derivs[3] = -2.0 * bValues[i] * g[0] * g[2] * estim;
        derivs[4] = -bValues[i] * g[1] * g[1] * estim;
        derivs[5] = -2.0 * bValues[i] * g[1] * g[2] * estim;
        derivs[6] = -bValues[i] * g[2] * g[2] * estim;

        //Check the derivatives against numerical derivatives.
        //    	double[] numDs2 = dydasNumerical(atry, i);
        //    	System.out.println("i = " + i);
        //    	for(int j=0; j<derivs.length; j++) {
        //    	   System.out.print(derivs[j] + " " + numDs2[j] + " : ");
        //   	}
        //    	System.out.println();

        return derivs;
    }

    /**
     * Returns the diffusion tensor represented by the current parameter values.
     * 
     * @return The diffusion tensor specified by the parameter settings.
     */
    public DT getDiffTensor() {
        return new DT(a[1], a[2], a[3], a[4], a[5], a[6]);
    }

}

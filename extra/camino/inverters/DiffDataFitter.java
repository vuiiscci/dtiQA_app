package inverters;

import optimizers.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Marquardt fitter especially for fitting models to diffusion data.
 * 
 * <dt>Description:
 * 
 * <dd>Contains methods that initialize the sigmas in the Levenberg-Marquardt
 * fitting routine and method that allow new measurements to be specified to run
 * repeated fits without creating new objects.
 * 
 * </dl>
 * 
 * @version $Id$
 * @author Danny Alexander
 *  
 */
abstract public class DiffDataFitter extends MarquardtChiSqFitter {

    /**
     * The non-zero b-values.
     */
    protected double[] bValues;

    /**
     * The number of b=0 images.
     */
    protected double M;

    /**
     * The number of parameters in the model to fit.
     */
    protected int noParams;

    /**
     * Initializes the instance variables and data structures.
     * 
     * @param indepVals
     *            The matrix of gradient directions g without the zero directions.
     * 
     * @param depVals
     *            The normalized measurements.
     * 
     * @param bVals
     *            The array of non-zero b-values.
     * 
     * @param nob0s
     *            The number of b=0 measurements.
     */
    protected void initialize(double[][] indepVals, double[] depVals,
            double[] bVals, int nob0s) throws MarquardtMinimiserException {

	// index diffusion times from 1 like y and g
        bValues = new double[bVals.length + 1];

	for (int i = 0; i < bVals.length; i++) {
	    bValues[i+1] = bVals[i];
	}

        M = (double) nob0s;

        // Pass sampled points to initialiser, as well as the
        // number of parameters of the model.
        initData(indepVals, depVals, noParams);

        // Initialize the parameters to the starting point and
        // set the standard deviations of the measurements.
        initASigs();
    }

    /**
     * Initialises a new fitting procedure with the same independent variables
     * (b-vectors), but a new set of measurements (dependent variables).
     * 
     * @param depVals
     *            The new set of measurements.
     */
    public void newDepVals(double[] depVals) throws MarquardtMinimiserException {
        if (depVals.length != y.length - 1) {
            System.err.println(depVals.length + " " + (y.length - 1));
            throw new MarquardtMinimiserException(
                    "New data contains the wrong number of values.");
        }
        for (int i = 0; i < ndata; i++) {
            y[i + 1] = depVals[i];
        }

        // Reinitialize the starting parameters and expected variances
        // of the measurements.
        initASigs();

    }

    /**
     * Initialises the parameter values and their standard deviations.
     */
    protected void initASigs() {

        // Initialize the parameters
        initAs();

        // Initialize the standard deviations.
        initSigs();

    }

    /**
     * Initializes the parameters to a default starting position before
     * optimization.
     */
    abstract protected void initAs();

    /**
     * Initializes the standard deviations of the samples.
     */
    protected void initSigs() {

        //Initialise standard deviations. Assuming constant
        //sd, then var(s/s0) = (sig(s)/s0)^2 (1 + (s/s0)^2/M)
        //Where M is the number of s0 acquisitions, sig(s) is
        //the standard deviation of each MR measurement. Here
        //we can ignore the constant sig(s)/s0 scaling.
        for (int i = 1; i <= ndata; i++) {
            if (M > 0) {
                sig[i] = Math.sqrt(1.0 + y[i] * y[i] / M);
            }
        }
    }

    /**
     * Returns the indexed gradient direction.
     * 
     * @param x
     *            The independent variable array.
     * 
     * @param i
     *            The index of the gradient direction required.
     * 
     * @return The gradient direction of the i-th sample.
     */
    public static double[] getG(double[][] x, int i) {
        double[] g = new double[3];

        g[0] = x[i][0];
        g[1] = x[i][1];
        g[2] = x[i][2];

        return g;
    }

}

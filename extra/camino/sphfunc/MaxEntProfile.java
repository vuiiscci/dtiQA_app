package sphfunc;

import mesd.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Function of the sphere represented using the maximum entropy
 * multiplicative basis.
 * 
 * <dt>Description:
 * 
 * <dd>The function is:
 * 
 * f(x) = exp(\lambda_0 + \sum_{i=1}^N \lambda_i R(q_i; x)),
 * 
 * where R is the response function, which must be specified by the user in
 * {@link mesd.SphDeconvKernels SphDeconvKernels} and q_i, i = 1, ..., N, is a set of unit vectors on the
 * sphere, which must be specified before using this class by calling
 * <code>setReconDirs</code>.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 */
public class MaxEntProfile extends SphericalFunction {

    /**
     * Array containing the lambda_i.
     */
    private double[] lambda;

    /**
     * Array of parameters of the response function, R.
     */
    private double[] kernelParams;

    /**
     * Array containing the q_i.
     */
    private static double[][] reconDirs;

    /**
     * number of parameters in kernal summation
     */
    private static int noParams;

    /**
     * Sets the set of reconstruction directions.
     * 
     * @param newKs
     *            {{x1, y1, z1}, {x2, y2, z2}, ..., {xN, yN, zN}}.
     */
    public static void setReconDirs(double[][] newKs) {
        reconDirs = newKs;
        noParams = reconDirs.length + 1;
    }

    /**
     * Returns the array of reconstruction directions (directions for the
     * lambdas), as distinct from gradient directions.
     * 
     * @return the reconDirs array
     */
    public static double[][] getLambdaDirs() {
        return reconDirs;
    }

    /**
     * Gets the number of lambdas to fit
     * 
     * @return number of parameters
     */
    public static int getNumParams() {
        return noParams;
    }

    /**
     * Creates a max. ent. profile given the list of lambdas and the parameters
     * of the response function.
     * 
     * @param coeffs
     *            {exitcode, ln A^\star(0), lambda_0, ..., lambda_N}
     * 
     * @param params
     *            The specification of the response function.
     */
    public MaxEntProfile(double[] coeffs, double[] params) {

        kernelParams = new double[params.length];
        for (int i = 0; i < params.length; i++) {
            kernelParams[i] = params[i];
        }

        lambda = new double[coeffs.length - 2];

        if(lambda.length != noParams) {
            throw new LoggedException("Number of lambdas provided is " + lambda.length + ".  Expected " + noParams + ".");
        }

        for (int i = 0; i < lambda.length; i++) {
            lambda[i] = coeffs[i + 2];
        }
    }


    public double getRadius(double x, double y, double z) {
        double expo = lambda[0];
        for (int i = 1; i < lambda.length; i++) {

            double[] qvec = reconDirs[i - 1];
            double[] xvec = new double[3];
            xvec[0] = x;
            xvec[1] = y;
            xvec[2] = z;

            expo += lambda[i] * SphDeconvKernels.kernel(qvec, xvec, kernelParams);
        }

        return Math.exp(expo);
    }

    /**
     * Returns the number of lambdas defining the profile.
     * 
     * @return N + 1.
     */
    public static int numLambdas() {
        return reconDirs.length + 1;
    }

}


package inverters;

import optimizers.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of a two tensor model to DW-MR data.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtChiSqFitter to
 * provide a Levenburg-Marquardt algorithm for fitting a two diffusion tensor
 * model to DW-MR measurements on a sphere in Fourier space. This class uses the
 * Cholesky decomposition to ensure that the fitted diffusion tensors are
 * positive definite
 * 
 * </dl>
 * 
 * $Id$
 * 
 * @author Danny Alexander
 *  
 */
public class TwoTensorCholFitter extends TwoTensorFitter {

    // Each diffusion tensor is positive definite so that
    // D = U^T U, where U is upper triangular. We optimize the elements
    // of U in each component tensor.
    // We use p = sin^2(\theta) to constrain the mixing parameter to
    // the interval [0, 1].
    // The optimized parameters are enumerated as follows:
    // a1 = u11, a2 = u12, a3 = u13, a4 = u14. a5 = u15, a6 = u16
    // a7 = u21, a8 = u22, a9 = u23, a10 = u24. a11 = u25, a12 = u26
    // a13 = \theta.
    // U_i = {{ui1, ui2, ui3},{0, ui4, ui5},{0, 0, ui6}}.

    /**
     * Default constructor does nothing.
     */
    public TwoTensorCholFitter() {
    }

    /**
     * The constructor requires a list of independent values (indepVals) and
     * associated dependent values (depVals) - these are the data. Together with
     * the number of unweighted acquisitions that are made (nob0s), which is
     * required to estimate the noise levels of each data item. The diffusion
     * time used in the imaging sequence is also required.
     * 
     * @param indepVals
     *            The matrix of gradient directions g without the zero entries.
     * 
     * @param depVals
     *            The normalized measurements.
     * 
     * @param bValues
     *            The array of non-zero b-values.
     * 
     * @param nob0s
     *            The number of b=0 measurements.
     */
    public TwoTensorCholFitter(double[][] indepVals, double[] depVals,
            double[] bValues, int nob0s) throws MarquardtMinimiserException {

        //This model has 13 parameters.
        noParams = 13;
        dt1StartIndex = 1;
        dt2StartIndex = 7;
        mixParIndex = 13;
        initialize(indepVals, depVals, bValues, nob0s);
    }


    /**
     * Overridden to use the Cholesky decomposition.
     */
    protected void insertTensorParamDerivs(double[] derivs, double[] atry, double[] g,
					   double d1Contrib, double d2Contrib, double p, double b) {

        double gDotRow11 = g[0] * atry[dt1StartIndex] + g[1] * atry[dt1StartIndex + 1]
                + g[2] * atry[dt1StartIndex + 2];
        double gDotRow12 = g[1] * atry[dt1StartIndex + 3] + g[2]
                * atry[dt1StartIndex + 4];
        double gDotRow13 = g[2] * atry[dt1StartIndex + 5];

        double gDotRow21 = g[0] * atry[dt2StartIndex] + g[1] * atry[dt2StartIndex + 1]
                + g[2] * atry[dt2StartIndex + 2];
        double gDotRow22 = g[1] * atry[dt2StartIndex + 3] + g[2]
                * atry[dt2StartIndex + 4];
        double gDotRow23 = g[2] * atry[dt2StartIndex + 5];

        insertCholDerivs(derivs, b, dt1StartIndex, g, gDotRow11, gDotRow12, gDotRow13,
                p, d1Contrib);
        insertCholDerivs(derivs, b, dt2StartIndex, g, gDotRow21, gDotRow22, gDotRow23,
                1 - p, d2Contrib);

    }

    /**
     * Overridden to use p = sin^2 a.
     */
    protected void insertOtherDerivs(double[] derivs, double[] atry, double d1Contrib,
            double d2Contrib, double p) {

        //dp/d\theta is 2cos(\theta)sin(\theta).
        derivs[mixParIndex] = (2.0 * Math.cos(atry[mixParIndex]) * Math
                .sin(atry[mixParIndex]))
                * (d1Contrib - d2Contrib);
    }

    /**
     * Overridden to use Cholesky decomposition.
     */
    protected DT getDT1(double[] a) {
        return getDT_Chol(a, dt1StartIndex);
    }

    /**
     * Overridden to use Cholesky decomposition.
     */
    protected DT getDT2(double[] a) {
        return getDT_Chol(a, dt2StartIndex);
    }

    /**
     * Overridden to use p = sin^2(a)
     */
    protected double getP(double[] a) {
        double rootP = Math.sin(a[mixParIndex]);
        return rootP * rootP;
    }

    /**
     * Overridden to use Cholesky decomposition
     */
    protected double[] dt1ToParams(DT dt1) {
        return getCholParams(dt1);
    }

    /**
     * Overridden to use Cholesky decomposition
     */
    protected double[] dt2ToParams(DT dt2) {
        return getCholParams(dt2);
    }

    /**
     * Overridden to use p = sin^2(a)
     */
    protected double mixParToParam(double mix) {
        return Math.asin(Math.sqrt(mix));
    }

}

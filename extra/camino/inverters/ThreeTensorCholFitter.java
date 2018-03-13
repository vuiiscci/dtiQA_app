package inverters;

import optimizers.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of a three-tensor model to DW-MR data.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtChiSqFitter to
 * provide a Levenburg-Marquardt algorithm for fitting a three diffusion tensor
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
public class ThreeTensorCholFitter extends ThreeTensorFitter {

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
    public ThreeTensorCholFitter() {
    }


    /**
     * The constructor requires a list of independent values (indepVals) and
     * associated dependent values (depVals) - these are the data. Together with
     * the number of unweighted acquisitions that are made (nob0s), which is
     * required to estimate the noise levels of each data item. The b-values used
     * in the imaging sequence are also required.
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
    public ThreeTensorCholFitter(double[][] indepVals, double[] depVals,
            double[] diffusionTimes, int nob0s) throws MarquardtMinimiserException {

        //This model has 20 parameters.
        noParams = 20;
        dt1StartIndex = 1;
        dt2StartIndex = 7;
        dt3StartIndex = 13;
        mixPar1Index = 19;
        mixPar2Index = 20;
        initialize(indepVals, depVals, diffusionTimes, nob0s);
    }

    /**
     * Overridden to use p = sin^2(\theta) for both constraints.
     */
    protected void insertMixPars(double[] startPoint, double mix1, double mix2) {
        startPoint[mixPar1Index] = Math.asin(Math.sqrt(mix1));
        startPoint[mixPar2Index] = Math.asin(Math.sqrt(mix2 / (1.0 - mix1)));
    }

    /**
     * Overridden to use the Cholesky decomposition.
     */
    protected void insertTensorParamDerivs(double[] derivs, double[] atry, double[] g,
            double d1Contrib, double d2Contrib, double d3Contrib, double p1, double p2, double b) {

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

        double gDotRow31 = g[0] * atry[dt3StartIndex] + g[1] * atry[dt3StartIndex + 1]
                + g[2] * atry[dt3StartIndex + 2];
        double gDotRow32 = g[1] * atry[dt3StartIndex + 3] + g[2]
                * atry[dt3StartIndex + 4];
        double gDotRow33 = g[2] * atry[dt3StartIndex + 5];

        insertCholDerivs(derivs, b, dt1StartIndex, g, gDotRow11, gDotRow12, gDotRow13,
                p1, d1Contrib);
        insertCholDerivs(derivs, b, dt2StartIndex, g, gDotRow21, gDotRow22, gDotRow23,
                p2, d2Contrib);
        insertCholDerivs(derivs, b, dt3StartIndex, g, gDotRow31, gDotRow32, gDotRow33,
                (1 - p1 - p2), d3Contrib);

    }

    /**
     * Overridden to use p = sin^2 a.
     */
    protected void insertOtherDerivs(double[] derivs, double[] atry, double d1Contrib,
            double d2Contrib, double d3Contrib, double p1, double p2) {

        //dp/d\theta is 2cos(\theta)sin(\theta).
        double dp1da1 = (2.0 * Math.cos(atry[mixPar1Index]) * Math
                .sin(atry[mixPar1Index]));
        double mix23 = p2 / (1 - p1);
        derivs[mixPar1Index] = dp1da1
                * (d1Contrib - mix23 * d2Contrib - (1 - mix23) * d3Contrib);

        // Add the derivative with respect to the second mixing parameter.
        double dp2da2 = (1 - p1)
                * (2.0 * Math.cos(atry[mixPar2Index]) * Math.sin(atry[mixPar2Index]));
        derivs[mixPar2Index] = dp2da2 * (d2Contrib - d3Contrib);
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
     * Overridden to use Cholesky decomposition.
     */
    protected DT getDT3(double[] a) {
        return getDT_Chol(a, dt3StartIndex);
    }

    /**
     * Overridden to use p1 = sin^2(a)
     */
    protected double getP1(double[] a) {
        double rootP = Math.sin(a[mixPar1Index]);
        return rootP * rootP;
    }

    /**
     * Overridden to use p2 = (1-p1) sin^2(a)
     */
    protected double getP2(double[] a) {
        double rootP = Math.sin(a[mixPar2Index]);
        double p1 = getP1(a);
        return (1.0 - p1) * rootP * rootP;
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
     * Overridden to use Cholesky decomposition
     */
    protected double[] dt3ToParams(DT dt3) {
        return getCholParams(dt3);
    }

}

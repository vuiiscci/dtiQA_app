package inverters;

import optimizers.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of a two tensor model with one prolate
 * axisymmetric tensor to DW-MR data.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtChiSqFitter to
 * provide a Levenburg-Marquardt algorithm for fitting a two tensor model with
 * one prolate axisymmetric tensor to DW-MR measurements.
 * 
 * </dl>
 * 
 * $Id$
 * 
 * @author Danny Alexander
 *  
 */
public class TwoTensorOneAxiSymFitter extends TwoTensorCholFitter {

    /**
     * Default constructor does nothing.
     */
    public TwoTensorOneAxiSymFitter() {
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
    public TwoTensorOneAxiSymFitter(double[][] indepVals, double[] depVals,
            double[] bValues, int nob0s) throws MarquardtMinimiserException {
        noParams = 11;
        dt1StartIndex = 1;
        dt2StartIndex = 5;
        mixParIndex = 11;
        initialize(indepVals, depVals, bValues, nob0s);
    }

    /**
     * Overridden as there is no need to compute the axisymmetric diffusion
     * tensor.
     */
    protected double[] initParams1() {
        double traceD1 = 2.1E-9;
        double[] params = new double[4];
        for (int i = 0; i < 3; i++) {
            params[i] = 0.0;
        }
        params[3] = Math.sqrt(traceD1 / 3.0);

        return params;
    }

    /**
     * Overridden as there is no need to compute the axisymmetric diffusion
     * tensor.
     */
    protected double[] fAndBetaToParams1(double[] fBeta) {
        double[] params = new double[4];

        params[0] = fBeta[0];
        params[1] = fBeta[1];
        params[2] = fBeta[2];
        params[3] = Math.sqrt(fBeta[3]);

        return params;
    }

    /**
     * Overridden to replace the first DT by the axisymmetric representation.
     */
    protected void insertTensorParamDerivs(double[] derivs, double[] atry, double[] g,
            double d1Contrib, double d2Contrib, double p, double b) {

        double s = atry[dt1StartIndex + 3];

        double gDotF1 = g[0] * atry[dt1StartIndex] + g[1] * atry[dt1StartIndex + 1]
                + g[2] * atry[dt1StartIndex + 2];

        double gDotRow21 = g[0] * atry[dt2StartIndex] + g[1] * atry[dt2StartIndex + 1]
                + g[2] * atry[dt2StartIndex + 2];
        double gDotRow22 = g[1] * atry[dt2StartIndex + 3] + g[2]
                * atry[dt2StartIndex + 4];
        double gDotRow23 = g[2] * atry[dt2StartIndex + 5];

        insertAxiSymDerivs(derivs, b, dt1StartIndex, g, s, gDotF1, p, d1Contrib);
        insertCholDerivs(derivs, b, dt2StartIndex, g, gDotRow21, gDotRow22, gDotRow23,
                1 - p, d2Contrib);

    }

    /**
     * Overridden to use the axisymmetric representation.
     */
    protected DT getDT1(double[] a) {
        return getDT_AxiSym(a, dt1StartIndex);
    }

    /**
     * Overridden to use the axisymmetric representation.
     */
    protected double[] dt1ToParams(DT dt1) {
        return getAxiSymParams(dt1);
    }

}

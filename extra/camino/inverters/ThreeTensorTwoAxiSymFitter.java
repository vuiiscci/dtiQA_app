package inverters;

import optimizers.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of a three tensor model with two prolate
 * axisymmetric tensors to DW-MR data.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtChiSqFitter to
 * provide a Levenburg-Marquardt algorithm for fitting a three tensor model with
 * two prolate axisymmetric tensors to DW-MR measurements.
 * 
 * </dl>
 * 
 * $Id$
 * 
 * @author Danny Alexander
 *  
 */
public class ThreeTensorTwoAxiSymFitter extends ThreeTensorOneAxiSymFitter {

    /**
     * Default constructor does nothing.
     */
    public ThreeTensorTwoAxiSymFitter() {
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
    public ThreeTensorTwoAxiSymFitter(double[][] indepVals, double[] depVals,
            double[] bValues, int nob0s) throws MarquardtMinimiserException {
        noParams = 16;
        dt1StartIndex = 1;
        dt2StartIndex = 5;
        dt3StartIndex = 9;
        mixPar1Index = 15;
        mixPar2Index = 16;
        initialize(indepVals, depVals, bValues, nob0s);
    }

    /**
     * Overridden as there is no need to compute the axisymmetric diffusion
     * tensor.
     */
    protected double[] initParams2() {
        double traceD2 = 2.4E-9;
        double[] params = new double[4];
        for (int i = 0; i < 3; i++) {
            params[i] = 0.0;
        }
        params[3] = Math.sqrt(traceD2 / 3.0);

        return params;
    }

    /**
     * Overridden as there is no need to compute the axisymmetric diffusion
     * tensor.
     */
    protected double[] fAndBetaToParams2(double[] fBeta) {
        double[] params = new double[4];

        params[0] = fBeta[0];
        params[1] = fBeta[1];
        params[2] = fBeta[2];
        params[3] = Math.sqrt(fBeta[3]);

        return params;
    }

    /**
     * Overridden to replace the second DT by the axisymmetric representation.
     */
    protected void insertTensorParamDerivs(double[] derivs, double[] atry, double[] g,
            double d1Contrib, double d2Contrib, double d3Contrib, double p1, double p2, double b) {

        double s1 = atry[dt1StartIndex + 3];
        double gDotF1 = g[0] * atry[dt1StartIndex] + g[1] * atry[dt1StartIndex + 1]
                + g[2] * atry[dt1StartIndex + 2];

        double s2 = atry[dt2StartIndex + 3];
        double gDotF2 = g[0] * atry[dt2StartIndex] + g[1] * atry[dt2StartIndex + 1]
                + g[2] * atry[dt2StartIndex + 2];

        double gDotRow31 = g[0] * atry[dt3StartIndex] + g[1] * atry[dt3StartIndex + 1]
                + g[2] * atry[dt3StartIndex + 2];
        double gDotRow32 = g[1] * atry[dt3StartIndex + 3] + g[2]
                * atry[dt3StartIndex + 4];
        double gDotRow33 = g[2] * atry[dt3StartIndex + 5];

        insertAxiSymDerivs(derivs, b, dt1StartIndex, g, s1, gDotF1, p1,
                d1Contrib);
        insertAxiSymDerivs(derivs, b, dt2StartIndex, g, s2, gDotF2, p2,
                d2Contrib);
        insertCholDerivs(derivs, b, dt3StartIndex, g, gDotRow31, gDotRow32, gDotRow33,
                1 - p1 - p2, d3Contrib);

    }

    /**
     * Overridden to use the axisymmetric representation.
     */
    protected DT getDT2(double[] a) {
        return getDT_AxiSym(a, dt2StartIndex);
    }

    /**
     * Overridden to use the axisymmetric representation.
     */
    protected double[] dt2ToParams(DT dt2) {
        return getAxiSymParams(dt2);
    }

}

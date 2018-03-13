package inverters;

import optimizers.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of a three-tensor model to DW-MR data
 * with fixed mixing parameters.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtChiSqFitter to
 * provide a Levenburg-Marquardt algorithm for fitting a three diffusion tensor
 * model to DW-MR measurements with fixed mixing parameters.
 * 
 * </dl>
 * 
 * $Id$
 * 
 * @author Danny Alexander
 *  
 */
public class ThreeTensorCholFixMP_Fitter extends ThreeTensorCholFitter {

    /**
     * Default constructor does nothing.
     */
    public ThreeTensorCholFixMP_Fitter() {
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
    public ThreeTensorCholFixMP_Fitter(double[][] indepVals, double[] depVals,
            double[] bValues, int nob0s) throws MarquardtMinimiserException {
        noParams = 18;
        dt1StartIndex = 1;
        dt2StartIndex = 7;
        dt3StartIndex = 13;
        mixPar1Index = 19;
        mixPar2Index = 20;
        initialize(indepVals, depVals, bValues, nob0s);
    }

    /**
     * Overridden to do nothing as mixing parameter is fixed.
     */
    protected void insertMixPars(double[] startPoint, double mix1, double mix2) {
    }

    /**
     * Overridden to do nothing as mixing parameter is fixed.
     */
    protected void insertOtherDerivs(double[] derivs, double[] atry, double d1Contrib,
            double d2Contrib, double d3Contrib, double p1, double p2) {

    }

    /**
     * Overridden to return a constant value.
     */
    protected double getP1(double[] a) {

        //The mixing parameter is fixed at 1/3.
        return 1.0 / 3.0;
    }

    /**
     * Overridden to return a constant value.
     */
    protected double getP2(double[] a) {

        //The mixing parameter is fixed at 1/3.
        return 1.0 / 3.0;
    }

}

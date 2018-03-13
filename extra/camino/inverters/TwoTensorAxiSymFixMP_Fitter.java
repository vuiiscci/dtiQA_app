package inverters;

import optimizers.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of a two prolate axisymmetric tensor
 * model to DW-MR data with fixed mixing parameter.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtChiSqFitter to
 * provide a Levenburg-Marquardt algorithm for fitting a two prolate
 * axisymmetric diffusion tensor model to DW-MR measurements with fixed mixing
 * parameter.
 * 
 * </dl>
 * 
 * $Id$
 * 
 * @author Danny Alexander
 *  
 */
public class TwoTensorAxiSymFixMP_Fitter extends TwoTensorAxiSymFitter {

    /**
     * Default constructor does nothing.
     */
    public TwoTensorAxiSymFixMP_Fitter() {
    }

    /**
     * The constructor requires a list of independent values (indepVals) and
     * associated dependent values (depVals) - these are the data. Together with
     * the number of unweighted acquisitions that are made (nob0s), which is
     * required to estimate the noise levels of each data item. 
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
     *            The number of bq=0 measurements.
     */
    public TwoTensorAxiSymFixMP_Fitter(double[][] indepVals, double[] depVals,
            double[] bValues, int nob0s) throws MarquardtMinimiserException {
        noParams = 8;
        dt1StartIndex = 1;
        dt2StartIndex = 5;
        mixParIndex = 9;
        initialize(indepVals, depVals, bValues, nob0s);
    }

    /**
     * Overridden to do nothing as mixing parameter is fixed.
     */
    protected void insertMixPar(double[] startPoint, double mix) {
    }

    /**
     * Overridden to do nothing as mixing parameter is fixed.
     */
    protected void insertOtherDerivs(double[] derivs, double[] atry, double d1Contrib,
            double d2Contrib, double p) {

    }

    /**
     * Computes the mixing parameter from the array of parameters.
     */
    protected double getP(double[] a) {

        //The mixing parameter is fixed at 0.5.
        return 0.5;
    }

}

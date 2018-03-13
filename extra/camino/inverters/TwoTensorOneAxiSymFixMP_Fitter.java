package inverters;

import optimizers.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of a two tensor model with one prolate
 * axisymmetric tensor to DW-MR data with fixed mixing parameter.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtChiSqFitter to
 * provide a Levenburg-Marquardt algorithm for fitting a two tensor model with
 * one prolate axisymmetric tensor to DW-MR measurements with fixed mixing
 * parameter.
 * 
 * </dl>
 * 
 * $Id: TwoTensorOneAxiSymFixMP_Fitter.java,v 1.3 2005/08/18 11:07:43 ucacmgh
 * Exp $
 * 
 * @author Danny Alexander
 *  
 */
public class TwoTensorOneAxiSymFixMP_Fitter extends TwoTensorOneAxiSymFitter {

    /**
     * Default constructor does nothing.
     */
    public TwoTensorOneAxiSymFixMP_Fitter() {
    }

    /**
     * The constructor requires a list of independent values (indepVals) and
     * associated dependent values (depVals) - these are the data. Together with
     * the number of unweighted acquisitions that are made (nob0s), which is
     * required to estimate the noise levels of each data item. 
     * 
     * @param indepVals
     *            The matrix of gradient directions g without b=0 entries.
     * 
     * @param depVals
     *            The normalized measurements.
     * 
     * @param bValues
     *            The array of non-zero b-values.
     * 
     * @param nob0s
     *            The number of q=0 measurements.
     */
    public TwoTensorOneAxiSymFixMP_Fitter(double[][] indepVals, double[] depVals,
            double[] bValues, int nob0s) throws MarquardtMinimiserException {
        noParams = 10;
        dt1StartIndex = 1;
        dt2StartIndex = 5;
        mixParIndex = 11;
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

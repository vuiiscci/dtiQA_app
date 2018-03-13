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
 * <dd>This class implements the abstract methods of MarquardtchiSqFitter to
 * provide a Levenburg-Marquardt algorithm for fitting a two diffusion tensor
 * model to DW-MR measurements on a sphere in Fourier space.
 * 
 * </dl>
 * 
 * @see inverters.ModelIndex 
 * 
 *
 * @author Danny Alexander
 * @version $Id$
 *  
 * 
 */
public class TwoTensorFitter extends TensorModelFitter {

    // With the two tensor model,
    // f(x) = (1-a)e^(x^T D1 x) + a e^(x^T D2 x)
    // and D1 and D2 are each specified by 6 parameters:
    // dixx, diyy, dizz, dixy, dixz, diyz.
    // The mixing parameter a is a 13th.

    // This class constrains the
    // values of the tensor diagonal elements to be positive by setting
    // d1xx = a2^2, d2yy = a9^2, etc. Note that this does NOT ensure
    // that the fitted tensor is positive definite.
    // The mixing parameter is constrained to lie in the range 0 to 1
    // by expressing it as:
    // a = 1 / (a1^2 + 1).

    /**
     * First index in the parameter array of the parameters defining the first
     * component tensor.
     */
    protected int dt1StartIndex;

    /**
     * First index in the parameter array of the parameters defining the second
     * component tensor.
     */
    protected int dt2StartIndex;

    /**
     * Index in the parameter array of mixing parameter.
     */
    protected int mixParIndex;

    /**
     * Default constructor does nothing.
     */
    public TwoTensorFitter() {
    }

    /**
     * The constructor requires a list of independent values (indepVals) and
     * associated dependent values (depVals) - these are the data. Together with
     * the number of unweighted acquisitions that are made (nob0s), which is
     * required to estimate the noise levels of each data item. The b values
     *  used in the imaging sequence is also required.
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
    public TwoTensorFitter(double[][] indepVals, double[] depVals, double[] bValues,
            int nob0s) throws MarquardtMinimiserException {

        //This model has 13 parameters.
        noParams = 13;
        dt1StartIndex = 2;
        dt2StartIndex = 8;
        mixParIndex = 1;

        initialize(indepVals, depVals, bValues, nob0s);
    }

    /**
     * Set the default starting configuration for the minimization. Both
     * component tensors are isotropic, but with different size.
     */
    protected void initAs() {

        // Initialise the proportions to equal
        double startMix = 0.5;

        // Get the parameter values.
        double[] startPoint = getParamArray(initParams1(), initParams2(), startMix);

        // Copy the starting point into the parameters array.
        for (int i = 0; i < startPoint.length; i++) {
            a[i] = startPoint[i];
        }
    }

    /**
     * Returns the default set of starting parameters for the first component
     * diffusion tensor.
     * 
     * @return An array of default starting parameter values.
     */
    protected double[] initParams1() {
        double traceD1 = 2.1E-9;
        DT dt1 = new DT(traceD1 / 3.0, 0.0, 0.0, traceD1 / 3.0, 0.0, traceD1 / 3.0);
        return dt1ToParams(dt1);
    }

    /**
     * Returns the default set of starting parameters for the second component
     * diffusion tensor.
     * 
     * @return An array of default starting parameter values.
     */
    protected double[] initParams2() {
        double traceD2 = 2.4E-9;
        DT dt2 = new DT(traceD2 / 3.0, 0.0, 0.0, traceD2 / 3.0, 0.0, traceD2 / 3.0);
        return dt2ToParams(dt2);
    }

    /**
     * Choose the starting point from a single diffusion tensor.
     */
    public void setStartFromSingleDT(DT singleDT) {

        double[] startPoint = getStartParamsFromSingleDT(singleDT);
        for (int i = 1; i <= noParams; i++) {
            a[i] = startPoint[i];
        }

    }

    /**
     * Computes the array of starting parameters from a single diffusion tensor.
     * 
     * @param singleDT
     *            The single DT fitted to the data.
     * 
     * @return An array of starting values for the model parameters.
     */
    protected double[] getStartParamsFromSingleDT(DT singleDT) {

        //Call a separate function that computes F and beta
        //for each component.
        double[] fsAndBetas = fsAndBetasFromSingleDT(singleDT);

        //The elements in fsAndBetas are f1, f2, beta1, beta2.

        //Convert f and beta to DiffTensors
        double[] fBeta = new double[4];

        fBeta[0] = fsAndBetas[0];
        fBeta[1] = fsAndBetas[1];
        fBeta[2] = fsAndBetas[2];
        fBeta[3] = fsAndBetas[6];
        double[] p1 = fAndBetaToParams1(fBeta);

        fBeta[0] = fsAndBetas[3];
        fBeta[1] = fsAndBetas[4];
        fBeta[2] = fsAndBetas[5];
        fBeta[3] = fsAndBetas[7];
        double[] p2 = fAndBetaToParams2(fBeta);

        //Set the initial mixing parameter to 0.5.
        double startMix = 0.5;

        return getParamArray(p1, p2, startMix);
    }

    /**
     * Computes the parameters for the first component DT from its axisymmetric
     * representation.
     * 
     * @param {F1,
     *            F2, F3, beta}
     * 
     * @return The array of parameters for the first component DT.
     */
    protected double[] fAndBetaToParams1(double[] fBeta) {
        DT d = fAndBetaToDT(fBeta);
        return dt1ToParams(d);
    }

    /**
     * Computes the parameters for the second component DT from its axisymmetric
     * representation.
     * 
     * @param {F1,
     *            F2, F3, beta}
     * 
     * @return The array of parameters for the second component DT.
     */
    protected double[] fAndBetaToParams2(double[] fBeta) {
        DT d = fAndBetaToDT(fBeta);
        return dt2ToParams(d);
    }

    /**
     * Computes the optimization parameters array corresponding to the parameter
     * arrays for the two diffusion tensors and the mixing parameter.
     * 
     * @param p1
     *            The parameters for the first diffusion tensor.
     * 
     * @param p2
     *            The parameters for second diffusion tensor.
     * 
     * @param mix
     *            The mixing parameter.
     * 
     * @return The array of parameter values.
     */
    protected double[] getParamArray(double[] p1, double[] p2, double mix) {

        //Construct the starting parameter array.
        double[] startPoint = new double[noParams + 1];

        for (int i = 0; i < p1.length; i++) {
            startPoint[dt1StartIndex + i] = p1[i];
        }
        for (int i = 0; i < p2.length; i++) {
            startPoint[dt2StartIndex + i] = p2[i];
        }

        insertMixPar(startPoint, mix);

        return startPoint;
    }

    /**
     * Inserts the mixing parameter in the starting point array.
     * 
     * @params startPoint The array of starting parameter.
     * 
     * @param mix
     *            The mixing parameter value.
     */
    protected void insertMixPar(double[] startPoint, double mix) {
        startPoint[mixParIndex] = mixParToParam(mix);
    }

    /**
     * Converts a DT object into a set of parameter values for the first
     * component diffusion tensor.
     * 
     * @param dt1
     *            The diffusion tensor.
     * 
     * @return The array of parameter values.
     */
    protected double[] dt1ToParams(DT dt1) {
        return getPosDiagParams(dt1);
    }

    /**
     * Converts a DT object into a set of parameter values for the second
     * component diffusion tensor.
     * 
     * @param dt2
     *            The diffusion tensor.
     * 
     * @return The array of parameter values.
     */
    protected double[] dt2ToParams(DT dt2) {
        return getPosDiagParams(dt2);
    }

    /**
     * Converts a mixing parameter value to the value in the representation used
     * in the optimization.
     * 
     * @param mix
     *            The mixing parameter value.
     * 
     * @return The optimized parameter value.
     */
    protected double mixParToParam(double mix) {
        return Math.sqrt((1.0 - mix) / mix);
    }

    /**
     * Returns the first component diffusion tensor.
     * 
     * @return The first component diffusion tensor.
     */
    public DT getDT1() {
        return getDT1(a);
    }

    /**
     * Returns the second component diffusion tensor.
     * 
     * @return The second component diffusion tensor.
     */
    public DT getDT2() {
        return getDT2(a);
    }

    /**
     * Returns the mixing parameter.
     * 
     * @return The mixing parameter.
     */
    public double getMix() {
        return getP(a);
    }

    protected double yfit(double[] atry, int i) {

        double[] g = getG(x, i);

        //Compute contribution from first tensor.
        DT d1 = getDT1(atry);
        double d1Contrib = Math.exp(-bValues[i] * d1.contractBy(g));

        //Compute contribution from second tensor.
        DT d2 = getDT2(atry);
        double d2Contrib = Math.exp(-bValues[i] * d2.contractBy(g));

        double p = getP(atry);

        double yVal = p * d1Contrib + (1 - p) * d2Contrib;

        return yVal;
    }

    /**
     * Overrides the default to compute the derivatives analytically.
     */
    protected double[] dydas(double[] atry, int i) {
        double[] derivs = new double[ma + 1];

        double[] g = getG(x, i);

        //Compute contribution from first tensor.
        DT d1 = getDT1(atry);
        double d1Contrib = Math.exp(-bValues[i] * d1.contractBy(g));

        //Compute contribution from second tensor.
        DT d2 = getDT2(atry);
        double d2Contrib = Math.exp(-bValues[i] * d2.contractBy(g));

        double p = getP(atry);

        // Insert the derivatives of the parameters of the diffusion
        // tensors.
        insertTensorParamDerivs(derivs, atry, g, d1Contrib, d2Contrib, p, bValues[i]);

        // Insert the derivative with respect to other parameters,
        // such as the mixing parameter, directly into the derivs array.
        insertOtherDerivs(derivs, atry, d1Contrib, d2Contrib, p);

        //Check the derivatives against numerical derivatives.
        boolean checkDerivatives = false;
        if (checkDerivatives) {
            double[] numDs2 = dydasNumerical(atry, i);
            System.out.println("i = " + i);
            for (int j = 0; j < derivs.length; j++) {
                System.out.print(atry[j] + " " + derivs[j] + " " + numDs2[j] + "  :  ");
            }
            System.out.println();
        }

        return derivs;
    }

    /**
     * Computes the derivates with respect to the parameters of the two
     * diffusion tensors.
     * 
     * @param derivs
     *            The derivatives array to put the results in.
     * 
     * @param atry
     *            The parameter values.
     * 
     * @param g
     *            The gradient direction.
     * 
     * @param d1Contrib
     *            The contribution from the first component.
     * 
     * @param d2Contrib
     *            The contribution from the second component.
     * 
     * @param p
     *            The mixing parameter.
     * 
     * @param b
     *            The b-value
     */
    protected void insertTensorParamDerivs(double[] derivs, double[] atry, double[] g,
            double d1Contrib, double d2Contrib, double p, double b) {

        // To compute the derivatives, we need q2 and the parameters
        // representing the diagonal elements of the DTs.
        double[] g2 = new double[6];
        g2[0] = g[0] * g[0];
        g2[1] = 2.0 * g[0] * g[1];
        g2[2] = 2.0 * g[0] * g[2];
        g2[3] = g[1] * g[1];
        g2[4] = 2.0 * g[1] * g[2];
        g2[5] = g[2] * g[2];

        insertPosDiagDerivs(derivs, b, dt1StartIndex, g2, atry[dt1StartIndex],
                atry[dt1StartIndex + 3], atry[dt1StartIndex + 5], p, d1Contrib);

        insertPosDiagDerivs(derivs, b, dt2StartIndex, g2, atry[dt2StartIndex],
                atry[dt2StartIndex + 3], atry[dt2StartIndex + 5], 1 - p, d2Contrib);

    }

    /**
     * Adds derivatives for parameters other than those defining the diffusion
     * tensors to the derivs array. The only other parameters are usually the
     * mixing parameter if it is variable.
     * 
     * @param derivs
     *            The derivatives array to put the results in.
     * 
     * @param atry
     *            The parameter values.
     * 
     * @param d1Contrib
     *            The contribution from the first component.
     * 
     * @param d2Contrib
     *            The contribution from the second component.
     * 
     * @param p
     *            The mixing parameter.
     */
    protected void insertOtherDerivs(double[] derivs, double[] atry, double d1Contrib,
            double d2Contrib, double p) {

        // Add the derivative with respect to the mixing parameter.
        derivs[mixParIndex] = (-2.0 * atry[mixParIndex] * (p * p))
                * (d1Contrib - d2Contrib);
    }

    /**
     * Constructs the DT for the first component from the array of parameters.
     * Assumes positive diagonals representation.
     * 
     * @param a
     *            The full list of model parameters.
     * 
     * @return The first diffusion tensor.
     */
    protected DT getDT1(double[] a) {
        return getDT_PosDiag(a, dt1StartIndex);
    }

    /**
     * Constructs the DT for the second component from the array of parameters.
     * Assumes positive diagonals representation.
     * 
     * @param a
     *            The full list of model parameters.
     * 
     * @return The second diffusion tensor.
     */
    protected DT getDT2(double[] a) {
        return getDT_PosDiag(a, dt2StartIndex);
    }

    /**
     * Computes the mixing parameter from the array of parameters. Assumes p =
     * 1/(a^2+1).
     * 
     * @param a
     *            The full list of model parameters.
     * 
     * @return The mixing parameter.
     */
    protected double getP(double[] a) {
        return invX2P1(a[mixParIndex]);
    }

    /**
     * Returns a TwoTensorFitter with type specified by the index and
     * initialized with the passed parameters. 
     *
     * @param indepVals
     *            The matrix of gradient directions g without the zero directions.
     * 
     * @param depVals
     *            The normalized measurements.
     * 
     * @param bValues
     *            The b values for each measurement.
     * 
     * @param nob0s
     *            The number of q=0 measurements.
     * 
     * @param index
     *            The index of the required model fitter.
     * 
     * @return The model fitter.
     *
     * @see inverters.ModelIndex
     */
    public static TwoTensorFitter getIndexedTwoTensorFitter(double[][] indepVals,
            double[] depVals, double[] bValues, int nob0s, ModelIndex index)
            throws MarquardtMinimiserException {

        if (index == ModelIndex.CYLCYL) {
            return new TwoTensorAxiSymFitter(indepVals, depVals, bValues, nob0s);
        }
	if (index == ModelIndex.CYLCYL_EQ) {
            return new TwoTensorAxiSymFixMP_Fitter(indepVals, depVals, bValues,
                    nob0s);
        }
	if (index == ModelIndex.POSPOS) {
            return new TwoTensorCholFitter(indepVals, depVals, bValues, nob0s);
        }
	if (index == ModelIndex.POSPOS_EQ) {
            return new TwoTensorCholFixMP_Fitter(indepVals, depVals, bValues, nob0s);
        }
	if (index == ModelIndex.POSCYL) {
            return new TwoTensorOneAxiSymFitter(indepVals, depVals, bValues, nob0s);
        }
	if (index == ModelIndex.POSCYL_EQ) {
            return new TwoTensorOneAxiSymFixMP_Fitter(indepVals, depVals, bValues,
                    nob0s);
        }
        else {
            return new TwoTensorFitter(indepVals, depVals, bValues, nob0s);
        }
    }

}

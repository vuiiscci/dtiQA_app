package inverters;

import optimizers.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of a three tensor model to DW-MR data.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtchiSqFitter to
 * provide a Levenburg-Marquardt algorithm for fitting a three diffusion tensor
 * model to DW-MR measurements on a sphere in Fourier space.
 * 
 * </dl>
 * 
 * @see inverters.ModelIndex 
 * 
 * @author Danny Alexander
 * @version $Id$  
 *
 */
public class ThreeTensorFitter extends TensorModelFitter {

    // With the three tensor model,
    // f(x) = a1 e^(x^T D1 x) + a2 e^(x^T D2 x) + a3 e^(x^T D3 x)
    // where D1, D2 and D3 are each specified by 6 parameters:
    // dixx, diyy, dizz, dixy, dixz, diyz, and
    // the mixing parameters a1, a2 and a3 sum to one.

    // This class constrains the
    // values of the tensor diagonal elements to be positive by setting
    // d1xx = a2^2, d2yy = a9^2, etc. Note that this does NOT ensure
    // that the fitted tensor is positive definite.
    // The mixing parameters are constrained to sum to one using
    // p_1 = 1/(a_1^2 + 1), p2 = (1-p1)/(a_2^2 + 1) and p3 = 1-p1-p2.

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
     * First index in the parameter array of the parameters defining the third
     * component tensor.
     */
    protected int dt3StartIndex;

    /**
     * Index in the parameter array of mixing parameter a1.
     */
    protected int mixPar1Index;

    /**
     * Index in the parameter array of mixing parameter a2.
     */
    protected int mixPar2Index;

    /**
     * Default constructor does nothing.
     */
    public ThreeTensorFitter() {
    }

    /**
     * The constructor requires a list of independent values (indepVals) and
     * associated dependent values (depVals) - these are the data. Together with
     * the number of unweighted acquisitions that are made (nob0s), which is
     * required to estimate the noise levels of each data item.
     * 
     * @param indepVals
     *             The matrix of gradient directions g without b=0 entries.
     * 
     * @param depVals
     *            The normalized measurements.
     * 
     * @param bValues
     *            The array of diffusion times.
     * 
     * @param nob0s
     *            The number of b=0 measurements.
     */
    public ThreeTensorFitter(double[][] indepVals, double[] depVals,
            double[] bValues, int nob0s) throws MarquardtMinimiserException {

        //This model has 20 parameters.
        noParams = 20;
        dt1StartIndex = 1;
        dt2StartIndex = 7;
        dt3StartIndex = 13;
        mixPar1Index = 19;
        mixPar2Index = 20;

        initialize(indepVals, depVals, bValues, nob0s);
    }

    /**
     * Set the default starting configuration for the minimization. Both
     * component tensors are isotropic, but with different size.
     */
    protected void initAs() {

        // Initialise the proportions approximately equal
        double startMix1 = 0.333;
        double startMix2 = 0.333;

        // Get the parameter values.
        double[] startPoint = getParamArray(initParams1(), initParams2(), initParams3(),
                startMix1, startMix2);

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
     * Returns the default set of starting parameters for the second component
     * diffusion tensor.
     * 
     * @return An array of default starting parameter values.
     */
    protected double[] initParams3() {
        double traceD3 = 1.8E-9;
        DT dt3 = new DT(traceD3 / 3.0, 0.0, 0.0, traceD3 / 3.0, 0.0, traceD3 / 3.0);
        return dt3ToParams(dt3);
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
        double[] fsAndBetas = fsAndBetasThreeFromSingleDT(singleDT);

        //The elements in fsAndBetas are f1, f2, beta1, beta2.

        //Convert f and beta to DiffTensors
        double[] fBeta = new double[4];

        fBeta[0] = fsAndBetas[0];
        fBeta[1] = fsAndBetas[1];
        fBeta[2] = fsAndBetas[2];
        fBeta[3] = fsAndBetas[9];
        double[] p1 = fAndBetaToParams1(fBeta);

        fBeta[0] = fsAndBetas[3];
        fBeta[1] = fsAndBetas[4];
        fBeta[2] = fsAndBetas[5];
        fBeta[3] = fsAndBetas[10];
        double[] p2 = fAndBetaToParams2(fBeta);

        fBeta[0] = fsAndBetas[6];
        fBeta[1] = fsAndBetas[7];
        fBeta[2] = fsAndBetas[8];
        fBeta[3] = fsAndBetas[11];
        double[] p3 = fAndBetaToParams3(fBeta);

        //Set the initial mixing parameters equal.
        double startMix1 = 0.333;
        double startMix2 = 0.333;

        return getParamArray(p1, p2, p3, startMix1, startMix2);
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
     * Computes the parameters for the third component DT from its axisymmetric
     * representation.
     * 
     * @param {F1,
     *            F2, F3, beta}
     * 
     * @return The array of parameters for the third component DT.
     */
    protected double[] fAndBetaToParams3(double[] fBeta) {
        DT d = fAndBetaToDT(fBeta);
        return dt3ToParams(d);
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
     * @param p3
     *            The parameters for third diffusion tensor.
     * 
     * @param mix1
     *            The mixing parameter for the first component.
     * 
     * @param mix2
     *            The mixing parameter for the second component.
     * 
     * @return The array of parameter values.
     */
    protected double[] getParamArray(double[] p1, double[] p2, double[] p3, double mix1,
            double mix2) {

        //Construct the starting parameter array.
        double[] startPoint = new double[noParams + 1];

        for (int i = 0; i < p1.length; i++) {
            startPoint[dt1StartIndex + i] = p1[i];
        }
        for (int i = 0; i < p2.length; i++) {
            startPoint[dt2StartIndex + i] = p2[i];
        }
        for (int i = 0; i < p3.length; i++) {
            startPoint[dt3StartIndex + i] = p3[i];
        }

        insertMixPars(startPoint, mix1, mix2);

        return startPoint;
    }

    /**
     * Inserts the mixing parameters in the starting point array.
     * 
     * @params startPoint The array of starting parameter.
     * 
     * @param mix1
     *            The mixing parameter for the first component.
     * 
     * @param mix2
     *            The mixing parameter for the second component.
     */
    protected void insertMixPars(double[] startPoint, double mix1, double mix2) {
        startPoint[mixPar1Index] = Math.sqrt((1.0 - mix1) / mix1);
        startPoint[mixPar2Index] = Math.sqrt((1.0 - mix1 - mix2) / mix2);
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
     * Converts a DT object into a set of parameter values for the third
     * component diffusion tensor.
     * 
     * @param dt3
     *            The diffusion tensor.
     * 
     * @return The array of parameter values.
     */
    protected double[] dt3ToParams(DT dt3) {
        return getPosDiagParams(dt3);
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
     * Returns the third component diffusion tensor.
     * 
     * @return The third component diffusion tensor.
     */
    public DT getDT3() {
        return getDT3(a);
    }

    /**
     * Returns the mixing parameter for the first component.
     * 
     * @return The mixing parameter.
     */
    public double getMix1() {
        return getP1(a);
    }

    /**
     * Returns the mixing parameter for the second component.
     * 
     * @return The mixing parameter.
     */
    public double getMix2() {
        return getP2(a);
    }

    protected double yfit(double[] atry, int i) {

        double[] g = getG(x, i);

        //Compute contribution from first tensor.
        DT d1 = getDT1(atry);
        double d1Contrib = Math.exp(-bValues[i] * d1.contractBy(g));

        //Compute contribution from second tensor.
        DT d2 = getDT2(atry);
        double d2Contrib = Math.exp(-bValues[i] * d2.contractBy(g));

        //Compute contribution from third tensor.
        DT d3 = getDT3(atry);
        double d3Contrib = Math.exp(-bValues[i] * d3.contractBy(g));

        double p1 = getP1(atry);
        double p2 = getP2(atry);

        double yVal = p1 * d1Contrib + p2 * d2Contrib + (1 - p1 - p2) * d3Contrib;

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

        //Compute contribution from third tensor.
        DT d3 = getDT3(atry);
        double d3Contrib = Math.exp(-bValues[i] * d3.contractBy(g));

        double p1 = getP1(atry);
        double p2 = getP2(atry);

        // Insert the derivatives of the parameters of the diffusion
        // tensors.
        insertTensorParamDerivs(derivs, atry, g, d1Contrib, d2Contrib, d3Contrib, p1, p2, bValues[i]);

        // Insert the derivative with respect to other parameters,
        // such as the mixing parameter, directly into the derivs array.
        insertOtherDerivs(derivs, atry, d1Contrib, d2Contrib, d3Contrib, p1, p2);

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
     * Computes the derivates with respect to the parameters of the three
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
     * @param d3Contrib
     *            The contribution from the third component.
     * 
     * @param p1
     *            The mixing parameter for the first component.
     * 
     * @param p2
     *            The mixing parameter for the second component.
     * 
     * @param b
     *            The b-value
     */
    protected void insertTensorParamDerivs(double[] derivs, double[] atry, double[] g,
					   double d1Contrib, double d2Contrib, double d3Contrib, 
					   double p1, double p2, double b) {

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
                atry[dt1StartIndex + 3], atry[dt1StartIndex + 5], p1, d1Contrib);

        insertPosDiagDerivs(derivs, b, dt2StartIndex, g2, atry[dt2StartIndex],
                atry[dt2StartIndex + 3], atry[dt2StartIndex + 5], p2, d2Contrib);

        insertPosDiagDerivs(derivs, b, dt3StartIndex, g2, atry[dt3StartIndex],
                atry[dt3StartIndex + 3], atry[dt3StartIndex + 5], 1 - p1 - p2, d3Contrib);

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
            double d2Contrib, double d3Contrib, double p1, double p2) {

        // Add the derivative with respect to the first mixing parameter.
        double dp1da1 = (-2.0 * atry[mixPar1Index] * (p1 * p1));
        double mix23 = p2 / (1 - p1);
        derivs[mixPar1Index] = dp1da1
                * (d1Contrib - mix23 * d2Contrib - (1 - mix23) * d3Contrib);

        // Add the derivative with respect to the second mixing parameter.
        double dp2da2 = (1 - p1) * (-2.0 * atry[mixPar2Index] * (mix23 * mix23));
        derivs[mixPar2Index] = dp2da2 * (d2Contrib - d3Contrib);
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
     * Constructs the DT for the third component from the array of parameters.
     * Assumes positive diagonals representation.
     * 
     * @param a
     *            The full list of model parameters.
     * 
     * @return The third diffusion tensor.
     */
    protected DT getDT3(double[] a) {
        return getDT_PosDiag(a, dt3StartIndex);
    }

    /**
     * Computes the first mixing parameter from the array of parameters. Assumes
     * p = 1/(a^2+1).
     * 
     * @param a
     *            The full list of model parameters.
     * 
     * @return The mixing parameter.
     */
    protected double getP1(double[] a) {
        return invX2P1(a[mixPar1Index]);
    }

    /**
     * Computes the secocnd mixing parameter from the array of parameters.
     * Assumes p = 1/(a^2+1).
     * 
     * @param a
     *            The full list of model parameters.
     * 
     * @return The mixing parameter.
     */
    protected double getP2(double[] a) {
        return (1.0 - getP1(a)) * invX2P1(a[mixPar2Index]);
    }

    /**
     * Returns a ThreeTensorFitter with type specified by the index and
     * initialized with the passed parameters. The indexes are:
     * 
     * 0 - ThreeTensorFitter. 1 - ThreeTensorAxiSymFitter. 2 -
     * ThreeTensorAxiSymFixMP_Fitter. 3 - ThreeTensorCholFitter. 4 -
     * ThreeTensorCholFixMP_Fitter. 5 - ThreeTensorOneAxiSymFitter. 6 -
     * ThreeTensorOneAxiSymFixMP_Fitter. 7 - ThreeTensorTwoAxiSymFitter. 8 -
     * ThreeTensorTwoAxiSymFixMP_Fitter. 9 - ThreeTensorFitter.
     * 
     * @param indepVals
     *            The matrix of gradient directions g without the zero entries.
     * 
     * @param depVals
     *            The normalized measurements.
     * 
     * @param bValues
     *            The b-values.
     * 
     * @param nob0s
     *            The number of b=0 measurements.
     * 
     * @param index
     *            The index of the required model fitter.
     * 
     * @return The model fitter.
     */
    public static ThreeTensorFitter getIndexedThreeTensorFitter(double[][] indepVals,
            double[] depVals, double[] bValues, int nob0s, ModelIndex index)
            throws MarquardtMinimiserException {

        if (index == ModelIndex.CYLCYLCYL) {
            return new ThreeTensorAxiSymFitter(indepVals, depVals, bValues, nob0s);
        }
        if (index == ModelIndex.CYLCYLCYL_EQ) {
            return new ThreeTensorAxiSymFixMP_Fitter(indepVals, depVals, bValues, nob0s);
        }
        if (index == ModelIndex.POSPOSPOS) {
            return new ThreeTensorCholFitter(indepVals, depVals, bValues, nob0s);
        }
        if (index == ModelIndex.POSPOSPOS_EQ) {
            return new ThreeTensorCholFixMP_Fitter(indepVals, depVals, bValues, nob0s);
        }
        if (index == ModelIndex.POSPOSCYL) {
            return new ThreeTensorOneAxiSymFitter(indepVals, depVals, bValues, nob0s);
        }
        if (index == ModelIndex.POSPOSCYL_EQ) {
            return new ThreeTensorOneAxiSymFixMP_Fitter(indepVals, depVals, bValues, nob0s);
        }
        if (index == ModelIndex.POSCYLCYL) {
            return new ThreeTensorTwoAxiSymFitter(indepVals, depVals, bValues, nob0s);
        }
        if (index == ModelIndex.POSCYLCYL_EQ) {
            return new ThreeTensorTwoAxiSymFixMP_Fitter(indepVals, depVals, bValues, nob0s);
        }
        else {
            return new ThreeTensorFitter(indepVals, depVals, bValues, nob0s);
        }
    }

}

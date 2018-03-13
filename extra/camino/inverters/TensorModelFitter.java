package inverters;

import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General abstract class for a simple-diffusion based model fitter.
 * 
 * <dt>Description:
 * 
 * <dd>Contains methods specifically for manipulation of data comprising
 * diffusion tensors.
 * 
 * </dl>
 * 
 * @version $Id$
 * @author Danny Alexander
 *  
 */
abstract public class TensorModelFitter extends DiffDataFitter {

    /**
     * Constructs a DT assuming the axisymmetric decomposition.
     * 
     * @param a
     *            The parameter values.
     * 
     * @param start
     *            The first index of the parameters of the diffusion tensor.
     * 
     * @return The diffusion tensor.
     */
    protected DT getDT_AxiSym(double[] a, int start) {

        double[] fAndBeta = new double[4];
        for (int i = 0; i < 3; i++) {
            fAndBeta[i] = a[i + start];
        }
        fAndBeta[3] = a[start + 3] * a[start + 3];

        return fAndBetaToDT(fAndBeta);
    }

    /**
     * Finds the axisymmetric parameters from a DT. This method assumes that the
     * DT is prolate.
     * 
     * @param d
     *            The diffusion tensor.
     * 
     * @return {F1, F2, F3, beta}
     */
    protected static double[] getAxiSymParams(DT d) {

        double[][] eigenSys = d.sortedEigenSystem();

        // Set beta to the smallest DT eigenvalue. If the DT
        // has a negative eigenvalue, we use the absolute value.
        double beta = Math.abs(eigenSys[0][2]);

        // Alpha is difference between the largest and smallest
        // eigenvalues.
        double alpha = eigenSys[0][0] - eigenSys[0][2];

        //Construct the parameter array.
        double[] params = new double[4];

        // n is the first eigenvector.

        double modF = Math.sqrt(alpha);
        params[0] = modF * eigenSys[1][0];
        params[1] = modF * eigenSys[2][0];
        params[2] = modF * eigenSys[3][0];
        params[3] = Math.sqrt(beta);

        return params;
    }

    /**
     * Constructs the DT from the parameters of its Cholesky decomposition.
     * 
     * @param a
     *            The parameter values.
     * 
     * @param start
     *            The first index of the parameters of the diffusion tensor.
     * 
     * @return The diffusion tensor.
     */
    protected static DT getDT_Chol(double[] a, int start) {
        return new DT(a[start] * a[start], a[start] * a[start + 1], a[start]
                * a[start + 2],
                a[start + 1] * a[start + 1] + a[start + 3] * a[start + 3], a[start + 1]
                        * a[start + 2] + a[start + 3] * a[start + 4], a[start + 2]
                        * a[start + 2] + a[start + 4] * a[start + 4] + a[start + 5]
                        * a[start + 5]);
    }

    /**
     * Returns the elements of the Cholesky decomposition U of the diffusion
     * tensor.
     * 
     * @param d
     *            The diffusion tensor.
     * 
     * @return {U11, U12, U13, U22, U23, U33}.
     */
    protected static double[] getCholParams(DT d) {
        double[] dData = d.getComponents();
        double d11 = dData[0];
        double d12 = dData[1];
        double d13 = dData[2];
        double d22 = dData[3];
        double d23 = dData[4];
        double d33 = dData[5];

        double[] u = new double[6];

        u[0] = (d11 > 0.0) ? Math.sqrt(d11) : 0.0;
        u[1] = (u[0] > 0) ? d12 / u[0] : 0.0;
        u[2] = (u[0] > 0) ? d13 / u[0] : 0.0;
        double eq2 = d22 - u[1] * u[1];
        u[3] = (eq2 > 0.0) ? Math.sqrt(eq2) : 0.0;
        u[4] = (u[3] > 0) ? (d23 - u[1] * u[2]) / u[3] : 0.0;
        double eq3 = d33 - u[2] * u[2] - u[4] * u[4];
        u[5] = (eq3 > 0.0) ? Math.sqrt(eq3) : 0.0;

        return u;
    }

    /**
     * Constructs the diffusion tensor from the values in array a starting from
     * position start assuming that they describe a diffusion tensor with
     * positive diagonal elements.
     * 
     * @param a
     *            Full array of model parameters.
     * 
     * @param start
     *            The first index of the parameters that specify the diffusion
     *            tensor.
     * 
     * @return The diffusion tensor.
     */
    protected static DT getDT_PosDiag(double[] a, int start) {
        return new DT(a[start] * a[start], a[start + 1], a[start + 2], a[start + 3]
                * a[start + 3], a[start + 4], a[start + 5] * a[start + 5]);
    }

    /**
     * Constructs an array of parameter values from a DT assuming the positive
     * diagonal elements representation.
     * 
     * @param dt
     *            The diffusion tensor.
     * 
     * @return The array of parameter values.
     */
    protected static double[] getPosDiagParams(DT dt) {

        double[] dData = dt.getComponents();

        dData[0] = Math.sqrt(dData[0]);
        dData[3] = Math.sqrt(dData[3]);
        dData[5] = Math.sqrt(dData[5]);

        return dData;
    }

    /**
     * Inserts the derivatives with respect to the parameters of the diffusion
     * tensor assuming the parameters are the Cholesky decomposition.
     * 
     * @param derivs
     *            The array containing the derivatives.
     * 
     * @param b
     *            The b-value.
     * 
     * @param start
     *            The first index in derivs array at which to put the
     *            derivatives with respect to the parameters of the diffusion
     *            tensor.
     * 
     * @param g
     *            The gradient direction.
     * 
     * @param gDotRow1
     *            The dot product of the gradient direction and the first row of the
     *            Cholesky decomposition of the diffusion tensor.
     * 
     * @param gDotRow2
     *            The dot product of the gradient direction and the second row of the
     *            Cholesky decomposition of the diffusion tensor.
     * 
     * @param gDotRow3
     *            The dot product of the gradient direction and the third row of the
     *            Cholesky decomposition of the diffusion tensor.
     * 
     * @param mix
     *            The weighting of this diffusion tensor.
     * 
     * @param estim
     *            The estimate of the measurement from this diffusion tensor at
     *            measurement [g, b].
     */
    protected static void insertCholDerivs(double[] derivs, double b, int start,
            double[] g, double gDotRow1, double gDotRow2, double gDotRow3, double mix,
            double estim) {

        derivs[start] = -2.0 * g[0] * gDotRow1 * mix * b * estim;
        derivs[start + 1] = -2.0 * g[1] * gDotRow1 * mix * b * estim;
        derivs[start + 2] = -2.0 * g[2] * gDotRow1 * mix * b * estim;
        derivs[start + 3] = -2.0 * g[1] * gDotRow2 * mix * b * estim;
        derivs[start + 4] = -2.0 * g[2] * gDotRow2 * mix * b * estim;
        derivs[start + 5] = -2.0 * g[2] * gDotRow3 * mix * b * estim;

    }

    /**
     * Inserts the derivatives with respect to the parameters of the diffusion
     * tensor assuming the parameters are the Cholesky decomposition.
     * 
     * @param derivs
     *            The array containing the derivatives.
     * 
     * @param b
     *            The b-value.
     * 
     * @param start
     *            The first index in derivs array at which to put the
     *            derivatives with respect to the parameters of the diffusion
     *            tensor.
     * 
     * @param g2
     *            The elements of g.g^T.
     * 
     * @param rootDxx
     *            The square root of the xx element of D.
     * 
     * @param rootDyy
     *            The square root of the yy element of D.
     * 
     * @param rootDzz
     *            The square root of the zz element of D.
     * 
     * @param mix
     *            The weighting of this diffusion tensor.
     * 
     * @param estim
     *            The estimate of the measurement from this diffusion tensor at
     *            measurement [g, b].
     */
    protected static void insertPosDiagDerivs(double[] derivs, double b, int start,
            double[] g2, double rootDxx, double rootDyy, double rootDzz, double mix,
            double estim) {
        derivs[start] = -2.0 * rootDxx * mix * b * g2[0] * estim;
        derivs[start + 1] = -mix * b * g2[1] * estim;
        derivs[start + 2] = -mix * b * g2[2] * estim;
        derivs[start + 3] = -2.0 * rootDyy * mix * b * g2[3] * estim;
        derivs[start + 4] = -mix * b * g2[4] * estim;
        derivs[start + 5] = -2.0 * rootDzz * mix * b * g2[5] * estim;
    }

    /**
     * Inserts the derivatives assuming the parameters are F and beta of the
     * axisymmetric diffusion tensor.
     * 
     * @param derivs
     *            The array containing the derivatives.
     * 
     * @param b
     *            The b-value.
     * 
     * @param start
     *            The first index in derivs array at which to put the
     *            derivatives with respect to the parameters of the diffusion
     *            tensor.
     * 
     * @param g
     *            The gradient direction.
     * 
     * @param s
     *            The square root of beta.
     * 
     * @param gDotF
     *            The dot product of g and F
     * 
     * @param mix
     *            The weighting of this diffusion tensor.
     * 
     * @param estim
     *            The estimate of the measurement from this diffusion tensor at
     *            measurement [g, b].
     */
    protected void insertAxiSymDerivs(double[] derivs, double b, int start, double[] g,
            double s, double gDotF, double mix, double estim) {

        derivs[start] = -2.0 * b * mix * g[0] * gDotF * estim;
        derivs[start + 1] = -2.0 * b * mix * g[1] * gDotF * estim;
        derivs[start + 2] = -2.0 * b * mix * g[2] * gDotF * estim;

        derivs[start + 3] = -2.0 * b * mix * s * estim;

    }

    /**
     * Computes the Fs and betas defining two axisymmetric tensors derived from
     * a single DT to use as a starting point for fitting the two-tensor model.
     * 
     * @param singleDT
     *            The single DT to use.
     * 
     * @return An array containing the Fs and betas for the two components.
     */
    public static double[] fsAndBetasFromSingleDT(DT singleDT) {

        double[][] eigenSys = singleDT.sortedEigenSystem();

        // Set both betas to the smallest DT eigenvalue. If the DT
        // has a negative eigenvalue, we use the absolute value.
        double beta1 = Math.abs(eigenSys[0][2]);
        double beta2 = Math.abs(eigenSys[0][2]);

        //Set alpha so that the averages of the eigenvalues of
        //the pair of diffusion tensors along corresponding axes
        //are the eigenvalues of the single diffusion tensor.
        double alpha1 = 2.0 * (eigenSys[0][0] - eigenSys[0][2]);
        double alpha2 = 2.0 * (eigenSys[0][1] - eigenSys[0][2]);

        //Construct the starting parameter array.
        double[] fsAndBetas = new double[8];

        //n1 and n2 are the first and second eigenvectors of the
        //single diffusion tensor, respectively.

        double modF1 = Math.sqrt(alpha1);
        fsAndBetas[0] = modF1 * eigenSys[1][0];
        fsAndBetas[1] = modF1 * eigenSys[2][0];
        fsAndBetas[2] = modF1 * eigenSys[3][0];
        fsAndBetas[6] = beta1;

        double modF2 = Math.sqrt(alpha2);
        fsAndBetas[3] = modF2 * eigenSys[1][1];
        fsAndBetas[4] = modF2 * eigenSys[2][1];
        fsAndBetas[5] = modF2 * eigenSys[3][1];
        fsAndBetas[7] = beta2;

        return fsAndBetas;
    }

    /**
     * Computes the Fs and betas defining three axisymmetric tensors derived
     * from a single DT to use as a starting point for fitting the three tensor
     * model. The three tensors have fixed anisotropy and are oriented along
     * each of the eigenvectors of the single DT.
     * 
     * @param singleDT
     *            The single DT to use.
     * 
     * @return An array containing the Fs and betas for the three components.
     */
    public static double[] fsAndBetasThreeFromSingleDT(DT singleDT) {

        double[][] eigenSys = singleDT.sortedEigenSystem();

        // Each starting DT is axisymmetric with the same
        // anisotropy so that \alpha_i = n \beta_i. One
        // points along each eigenvector of the single DT and
        // the average eigenvalue along each is the corresponding
        // eigenvalue of the single DT.
        double n = 7.0;
        double l1 = Math.abs(eigenSys[0][0]);
        double l2 = Math.abs(eigenSys[0][1]);
        double l3 = Math.abs(eigenSys[0][2]);
        double scal = 3.0 / (n * (n + 3.0));

        double beta1 = scal * (l1 * (n + 2.0) - l2 - l3);
        double beta2 = scal * (l2 * (n + 2.0) - l1 - l3);
        double beta3 = scal * (l3 * (n + 2.0) - l1 - l2);

        //Set alpha so that the averages of the eigenvalues of
        //the diffusion tensors along corresponding axes
        //are the eigenvalues of the single diffusion tensor.
        double alpha1 = n * beta1;
        double alpha2 = n * beta2;
        double alpha3 = n * beta3;

        //Construct the starting parameter array.
        double[] fsAndBetas = new double[12];

        //n1 and n2 are the first and second eigenvectors of the
        //single diffusion tensor, respectively.

        double modF1 = Math.sqrt(alpha1);
        fsAndBetas[0] = modF1 * eigenSys[1][0];
        fsAndBetas[1] = modF1 * eigenSys[2][0];
        fsAndBetas[2] = modF1 * eigenSys[3][0];
        fsAndBetas[9] = beta1;

        double modF2 = Math.sqrt(alpha2);
        fsAndBetas[3] = modF2 * eigenSys[1][1];
        fsAndBetas[4] = modF2 * eigenSys[2][1];
        fsAndBetas[5] = modF2 * eigenSys[3][1];
        fsAndBetas[10] = beta2;

        double modF3 = Math.sqrt(alpha3);
        fsAndBetas[6] = modF3 * eigenSys[1][2];
        fsAndBetas[7] = modF3 * eigenSys[2][2];
        fsAndBetas[8] = modF3 * eigenSys[3][2];
        fsAndBetas[11] = beta3;

        return fsAndBetas;
    }

    /**
     * Computes the DT from F and Beta.
     * 
     * @param fAndBeta
     *            The array containing F and Beta.
     * 
     * @return The diffusion tensor.
     */
    public static DT fAndBetaToDT(double[] fAndBeta) {
        return new DT(fAndBeta[0] * fAndBeta[0] + fAndBeta[3], fAndBeta[0] * fAndBeta[1],
                fAndBeta[0] * fAndBeta[2], fAndBeta[1] * fAndBeta[1] + fAndBeta[3],
                fAndBeta[1] * fAndBeta[2], fAndBeta[2] * fAndBeta[2] + fAndBeta[3]);
    }

    /**
     * One of the standard functions for constraining optimization parameters
     * (eg a mixing parameter) to [0, 1]. It is (x^2 + 1)^(-1), where x has any
     * real value.
     * 
     * @param The
     *            optimization parameter.
     * 
     * @return The number in [0, 1]
     */
    public static double invX2P1(double x) {
        return 1.0 / (x * x + 1.0);
    }

}

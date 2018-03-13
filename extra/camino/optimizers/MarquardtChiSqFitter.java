package optimizers;

import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Model fitting algorithm using a Levenburg Marquardt algorithm.
 * 
 * <dt>Description:
 * 
 * <dd>Extends the basic Levenburg Marquardt algorithm in <code>MarquardtMinimiser</code> to
 * minimise a chi^2 goodness of fit measure between a model and some data. yfit
 * and dydas need to be implemented in subclasses, which encapsulate the model
 * for the data and its derivatives. The second derivatives are approximated by
 * multiples of first derivatives avoiding the need for their computation.
 * 
 * </dl>
 * 
 * @version $Id$
 * @author Danny Alexander
 *  
 */
abstract public class MarquardtChiSqFitter extends MarquardtMinimiser {

    //Number of sampled data points.
    protected int ndata;

    //Independent values of sampled points.
    protected double[][] x;

    //Dependent values at sampled points.
    protected double[] y;

    //Standard variations for each sampled point.
    protected double[] sig;

    // Specifies whether the optimization should use the full Hessian of
    // the objective function or the Levenburg-Marquardt approximation.
    protected boolean useHessian = false;

    /**
     * Initialises data and working arrays.
     * 
     * @param indepVals
     *            The values of the independent variables (sample points).
     * 
     * @param depVals
     *            The values of the dependent variables (data).
     * 
     * @param noParams
     *            The number of parameters to fit.
     */
    protected void initData(double[][] indepVals, double[] depVals, int noParams)
            throws MarquardtMinimiserException {

        if (depVals.length != indepVals.length) {
            throw new MarquardtMinimiserException(
                    "Numbers of dependent and independent values differ in MarquardtFitter.initData");
        }

        x = new double[indepVals.length + 1][indepVals[0].length];
        y = new double[depVals.length + 1];
        ndata = y.length - 1;
        sig = new double[ndata + 1];
        for (int i = 0; i < ndata; i++) {
            for (int j = 0; j < x[0].length; j++) {
                x[i + 1][j] = indepVals[i][j];
            }
            y[i + 1] = depVals[i];
            sig[i + 1] = 1.0;
        }

        init(noParams);
    }

    /**
     * Returns the value of the model function at x[i] with parameters atry.
     * 
     * @param atry
     *            The parameters of the model to use.
     * 
     * @param i
     *            The index of the data point to use.
     * 
     * @return The error between the fit and the data point.
     */
    abstract protected double yfit(double[] atry, int i);

    /**
     * Returns an array of values of the derivatives of the model at x[i] with
     * respect to each of the parameters in a, using values in atry. This
     * default computes the derivatives numerically, but this method can be
     * overridden to compute the derivatives analytically.
     * 
     * @param atry
     *            The parameters of the model to use.
     * 
     * @param i
     *            The index of the data point to use.
     * 
     * @return The derivates of the model at atry for data point i.
     */
    protected double[] dydas(double[] atry, int i) {
        return dydasNumerical(atry, i);
    }

    /**
     * Compute the derivatives numerically.
     * 
     * @param atry
     *            The parameters of the model to use.
     * 
     * @param i
     *            The index of the data point to use.
     * 
     * @return The derivates of the model at atry for data point i.
     */
    protected double[] dydasNumerical(final double[] atry, final int i) {
        double[] derivs = new double[ma + 1];

        final double[] a = new double[atry.length];
        for (int k = 1; k < atry.length; k++) {
            a[k] = atry[k];
        }

        for (int j = 1; j <= ma; j++) {

            final int index = j;

            NumDeriv df = new NumDeriv() {
                protected float func(float arg) {

                    a[index] = arg;
                    float t = (float) yfit(a, i);
                    return t;
                }
            };

            float[] err = new float[1];
            derivs[j] = df.dfridr((float) atry[j],
                    (atry[j] != 0.0) ? 0.1f * (float) atry[j] : 0.1f, err);

            //Replace the original value in a[index] ready to compute
            //the next derivative.
            a[index] = atry[index];
        }

        return derivs;
    }

    /**
     * Returns the matrix of second derivatives of the model at x[i] with
     * respect to each pair of parameters in a, using values in atry. This
     * default computes the derivatives numerically, but this method can be
     * overridden to compute the derivatives analytically.
     * 
     * @param atry
     *            The parameters of the model to use.
     * 
     * @param i
     *            The index of the data point to use.
     * 
     * @return The second derivates of the model at atry for data point i.
     */
    protected double[][] d2yda2s(double[] atry, int i) {
        return d2yda2sNumerical(atry, i);
    }

    /**
     * Compute the second derivatives numerically.
     * 
     * @param atry
     *            The parameters of the model to use.
     * 
     * @param i
     *            The index of the data point to use.
     * 
     * @return The second derivates of the model at atry for data point i.
     */
    protected double[][] d2yda2sNumerical(final double[] atry, final int i) {
        double[][] derivs2 = new double[ma + 1][ma + 1];

        final double[] a = new double[atry.length];
        for (int k = 1; k < atry.length; k++) {
            a[k] = atry[k];
        }

        for (int j = 1; j <= ma; j++) {
            for (int j2 = j; j2 <= ma; j2++) {

                final int index1 = j;
                final int index2 = j2;

                NumDeriv df = new NumDeriv() {
                    protected float func(float arg) {

                        a[index1] = arg;
                        double[] ts = dydas(a, i);

                        return (float) ts[index2];
                    }
                };

                float[] err = new float[1];
                derivs2[j][j2] = derivs2[j2][j] = df.dfridr((float) atry[j2],
                        (atry[j2] != 0.0) ? 0.1f * (float) atry[j2] : 0.1f, err);

                //Replace the original value in a[index] ready to compute
                //the next derivative.
                a[index1] = atry[index1];
            }
        }

        return derivs2;
    }

    /**
     * Implements the chi-squared function and its derivatives. The second
     * derivative matrix is replaced by values proportional to the products of
     * first derivatives: <code>d2fda2[i][j] = dfda[i]*dfda[j]*sig2i</code>, which avoids
     * computation and can improve stability.
     * 
     * @param params
     *            The current model parameter settings.
     * 
     * @param dfda
     *            The first derivatives of the model with respect to each
     *            parameter.
     * 
     * @param d2fda2
     *            The second derivatives of the model with respect to each pair
     *            of parameters.
     */
    protected double fObj(double[] params, double[] dfda, double[][] d2fda2) {

        //Initialise values and arrays.
        double chisq = 0.0;
        for (int i = 1; i <= ma; i++) {
            for (int j = 1; j <= ma; j++) {
                d2fda2[i][j] = 0.0;
            }
            dfda[i] = 0.0;
        }

        for (int i = 1; i <= ndata; i++) {
            double ymod = yfit(params, i);

            double[] dyda = dydas(params, i);

            double[][] d2yda2 = null;
            if (useHessian) {
                d2yda2 = d2yda2s(params, i);
            }

            double sig2i = 1.0 / (sig[i] * sig[i]);
            double dy = y[i] - ymod;
            for (int l = 1; l <= ma; l++) {
                double wt = 2.0 * dyda[l] * sig2i;
                for (int m = 1; m <= l; m++) {
                    d2fda2[m][l] += wt * dyda[m];
                    if (useHessian) {
                        d2fda2[m][l] += 2.0 * d2yda2[m][l] * dy * sig2i;
                    }
                }
                dfda[l] -= dy * wt;
            }
            chisq += dy * dy * sig2i;
        }

        //Fill in other off-diagonal of second derivatives assuming
        //symmetry.
        for (int i = 2; i <= ma; i++) {
            for (int j = 1; j < i; j++) {
                d2fda2[i][j] = d2fda2[j][i];
            }
        }

        return chisq;
    }

    /**
     * Allows the standard deviations of the data points to be specified.
     * 
     * @param index
     *            The index of the standard deviation to set.
     * 
     * @param std
     *            The new value of the standard deviation.
     */
    public void setSig(int index, double std) {
        if (index > 0 && index < sig.length) {
            sig[index] = std;
        }
    }

    /**
     * Returns the value of chi-squared.
     * 
     * @return chi-squared.
     */
    public double getChiSq() {
        return getFObjVal();
    }

    /**
     * Returns the average residual squared error given the current values in a.
     * 
     * @return Average residual squared error.
     */
    public double getResidual() {
        double residual = 0.0;
        for (int i = 1; i <= ndata; i++) {
            double ymod = yfit(a, i);
            residual += (ymod - y[i]) * (ymod - y[i]);
        }

        residual /= (double) ndata;

        return residual;
    }

    /**
     * Returns the residual error (absolute difference from fitted value) for
     * each data point.
     * 
     * @return Array of residual errors.
     */
    public double[] getResiduals() {
        double[] residuals = new double[ndata];
        for (int i = 1; i <= ndata; i++) {
            double ymod = yfit(a, i);
            residuals[i - 1] = Math.abs(ymod - y[i]);
        }

        return residuals;
    }

    /**
     * Returns the mean squared relative error given the current values in a.
     * 
     * @return The average relative residual squared error.
     */
    public double getRelativeResidual() {
        double residual = 0.0;
        for (int i = 1; i <= ndata; i++) {
            double ymod = yfit(a, i);
            residual += (ymod - y[i]) * (ymod - y[i]) / (y[i] * y[i]);
        }

        residual /= (double) ndata;

        return residual;
    }

}

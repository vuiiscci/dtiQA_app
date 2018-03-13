
package mesd;

import java.util.Vector;
import java.util.logging.Logger;

import optimizers.*;
import sphfunc.*;
import misc.*;
import tools.CL_Initializer;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Non-linear least squares fitter of the maximum entropy representation of
 * the spherical deconvolution.
 * 
 * <dt>Description:
 * 
 * <dd>This class implements the abstract methods of MarquardtChiSqFitter to
 * provide a Levenburg--Marquardt algorithm for fitting the max. ent. form: p(x) =
 * \exp(\lambda_0 \exp(\sum_j \lambda_j R(q_j; x) to DW-MR data which is related
 * via \int_{|x|=r} p(x) R(q; x) dx = A(q)
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id: MESD_Fitter.java,v 1.8 2005/08/16 11:50:57
 *         ucacmgh Exp $
 *  
 */
public class MESD_Fitter extends MarquardtChiSqFitter {

    Logger logger = Logger.getLogger("camino." + this.getClass().getName());

    //Variables containing values of parameters for which the function
    //and derivatives were last calculated in case they don't need
    //to be recomputed.
    protected double[] lastATry;

    protected int lastI;

    protected double[] derivs;

    protected double yVal;

    //Arrays storing various values that remain unchanged during the
    //minimisation, but are continually used to compute the integral over
    //the sphere.
    //Point sets used to do integration during minimisation.
    protected double[][] xvecs;

    protected Vector<SphericalPointSet> pointSets;
    //J1.4 protected Vector pointSets;

    protected int psIndex = 0;

    //Test point set.
    protected double[][] testPS;

    protected double[] pxsTS;

    protected double ptAreaTS;

    //The test point set is the Archimedean point set of size
    //350^2. This should not be changed unless a new regression
    //line is computed to estimate integration errors with the new
    //point set.
    final protected int testPSResolution = 350;

    /**
     * Parameters of the deconvolution kernel. The first entry specifies the
     * function itself. The subsequent entries are parameter settings of the
     * kernel.
     */
    protected double[] kernelParams;

    //Other parameters.
    protected double ptArea;


    //We precompute the response function for each measurement and each x,
    //which are used in each computation of the numerical integral. First
    //index is k (measurement), second is x.
    protected double[][] ckdx;

    protected Vector<DPSet> ckdxs;
    //J1.4 protected Vector ckdxs;

    protected double[][] ckdxTS;

    //Sampled unit vectors - equivalent to indepVals but expressed in
    //x,y,z coordinates.
    protected double[][] kvecs;


    //We also need the response function for each representation
    //parameter vector and each x, which again are used in each
    //computation of the numerical integral. First index is p
    //(parameter vector), second is x.
    protected double[][] cpdx;

    protected Vector<DPSet> cpdxs;
    //J1.4 protected Vector cpdxs;

    protected double[][] cpdxTS;

    //direction set for parameter summation (by default, same as sample
    // directions)
    protected double[][] pvecs;

    //For estimation of standard deviations, we need to know the number of
    //b0 images acquired.
    protected double n;

    //We store the value of p at each x
    //in the current point set to avoid repeating
    //computation unnecessarily for a fixed set of parameter values.
    //The values of p(x) are stored in this array.
    protected double[] pxs;

    /**
     * Constructor.
     * 
     * @param indepVals
     *            Normalized gradient directions with zeros removed.
     * 
     * @param nob0s
     *            The number of b=0 measurements in the acquisition.
     * 
     * @param params
     *            Specification of the deconvolution kernel.
     */
    public MESD_Fitter(double[][] indepVals, int nob0s, double[] params, int numLambdas)
            throws MarquardtMinimiserException {
        n = (double) nob0s;

        kernelParams = new double[params.length];
        for (int i = 0; i < params.length; i++) {
            kernelParams[i] = params[i];
        }

        // Add a point at zero at the start of the list of independent
        // values.
        double[][] indepValsWithZero = new double[indepVals.length + 1][indepVals[0].length];
        indepValsWithZero[0][0] = 0.0;
        indepValsWithZero[0][1] = 0.0;
        indepValsWithZero[0][2] = 0.0;
        int numSkipped=0;
        for (int i = 0; i < indepVals.length; i++) {
            indepValsWithZero[i + 1][0] = indepVals[i][0];
            indepValsWithZero[i + 1][1] = indepVals[i][1];
            indepValsWithZero[i + 1][2] = indepVals[i][2];
        }
        
        initialise(indepValsWithZero, numLambdas);
                
    }

    /**
     * Specifies the index of the point set to be used to compute the numerical
     * integrals.
     * 
     * @param i
     *            The point set index.
     */
    public void setPointSet(int i) {
        psIndex = i;
        xvecs = ((SphericalPointSet) pointSets.elementAt(i)).data;
        ckdx = ((DPSet) ckdxs.elementAt(i)).getData();
        cpdx = ((DPSet) cpdxs.elementAt(i)).getData();

        //Initialise the array that stores the values of p(x) for each
        //iteration.
        pxs = new double[xvecs.length];

        //The contribution of each point is assumed equal.
        // This is the weighting factor for each point.
        ptArea = computePT_Area(xvecs.length);

    }

    /**
     * Reinitialises the parameters of the function to be minimised to the
     * starting position.
     */
    public void reInit() {
        initAs();
    }

    /**
     * Initialises a new fitting procedure with the same independent variables
     * (sampled directions), but a new set of measurements (dependent
     * variables).
     * 
     * @param depVals
     *            The set of normalized diffusion weighted measurements.
     */
    public void newDepVals(double[] depVals) throws MarquardtMinimiserException {

        if (depVals.length != y.length - 2) {
            logger.severe("New data contains the wrong number of values."
                    + depVals.length + " " + (y.length - 2));
            System.err.println(depVals.length + " " + (y.length - 2));
            throw new MarquardtMinimiserException(
                    "New data contains the wrong number of values.");
        }

        // The first dependent value is always 1.0 corresponding to
        // the measurement at b=0.
        y[1] = 1.0;
        for (int i = 0; i < ndata - 1; i++) {
            y[i + 2] = depVals[i];
        }

        initASigs();

        setPointSet(0);

    }

    /**
     * Called by constructor to initialize the data structures for the
     * optimization.
     * 
     * @param indepVals
     *            Normalized gradient directions with zeros removed.
     */
    protected void initialise(double[][] indepVals, int noParams)
            throws MarquardtMinimiserException {

        // Initialize the independent values. Pass a dummy set of
        // measurements to the initializer here.
        //initData(indepVals, new double[indepVals.length], indepVals.length);
        initData(indepVals, new double[indepVals.length], noParams);

        // Initialise the array of k-vectors.
        initKVecs(indepVals);

        // Initialise parameter k-vectors
        initParamKVecs(MaxEntProfile.getLambdaDirs());

        //Initialise the test point set.
        testPS = SphericalPoints.getZPhiXs(testPSResolution);
        pxsTS = new double[testPS.length];

        //Precalculate the values required for the summed approx.
        //to the integral over the sphere.
        initPointSets();

        //Initial values of the parameters and the standard
        //deviations of each data point can be set here.
        initASigs();

        //Start with the smallest point set.
        setPointSet(0);

        //Initialise remaining working variables.
        lastATry = new double[ma + 1];
        derivs = new double[ma + 1];

    }

    /**
     * Initialises the array of parameter values to be optimised and their
     * standard deviations.
     */
    protected void initASigs() {

        // Initialize the lambda_j
        initAs();

        // Initialize the standard deviations.
        initSigs();

    }

    /**
     * Initializes the lambda_j before optimization. All lambda_j for j>=1 are
     * initially zero. Lambda_0 is initially log(4 pi) to ensure that the
     * initial profile integrates to one.
     */
    protected void initAs() {

        a[1] = -Math.log(4.0 * Math.PI);

        //Initialise all the lambdas to zero.
        for (int i = 2; i <= ma; i++) {
            a[i] = 0.0;
        }

    }

    /**
     * Initializes the standard deviations of the normalized measurements.
     */
    protected void initSigs() {

        // Initialise standard deviations. Assuming constant
        // sd, then var(s/s0) = (sig(s)/s0)^2 (1 + (s/s0)^2/n)
        // Where n is the number of s0 acquisitions, sig(s) is
        // the standard deviation of each MR measurement. Here
        // we can ignore the constant sig(s)/s0 scaling.
        for (int i = 1; i <= ndata; i++) {
            sig[i] = Math.sqrt(1.0 + y[i] * y[i] / n);
        }

        // The standard deviation of the normalized measurement at
        // zero depends on the number of repeats. This makes the
        // normalization much stricter and the reconstructed functions
        // less spiky. Not used by default.
        //sig[1] = Math.sqrt(1.0/n);
    }

    /**
     * Computes the value of the model at the specified sample point with
     * parameter estimates in atry.
     * 
     * @param atry
     *            The set of objective function parameters at which to evaluate
     *            the function.
     * 
     * @param i
     *            The index of the independent variable at which to evaluate the
     *            objective function.
     * 
     * @return The value of the objective function.
     */
    protected double yfit(double[] atry, int i) {

        if (!sameAsLast(atry, i)) {

            //Compute the value and derivates.
            evaluateFunction(atry, i);
        }

        return yVal;
    }

    /**
     * Computes the derivatives of the model at the specified sample point with
     * parameter estimates in atry.
     * 
     * @param atry
     *            The set of objective function parameters at which to evaluate
     *            the function.
     * 
     * @param i
     *            The index of the independent variable at which to evaluate the
     *            objective function.
     * 
     * @return [dy/da_1, dy/da_2...dy/da_N] all evaluated at independent
     *         variable i.
     */
    protected double[] dydas(double[] atry, int i) {

        if (!sameAsLast(atry, i)) {
            //Compute the value and derivates.
            evaluateFunction(atry, i);
        }

        return derivs;
    }

    /**
     * Overrides the original version in MarquardtChiSqFitter so that we can
     * only compute the required values of p(x) once rather than repeating
     * computation.
     */
    protected void mrqcof(double[] atry, double[][] alpha, double[] beta)
            throws MarquardtMinimiserException {

        //The function and gradients can be computed more quickly if
        //we compute the value of p(x) at each x location outside the
        //loop below - they are the same each time round.
        computePXs(atry);

        //Now run mrqcof as usual.
        super.mrqcof(atry, alpha, beta);
    }

    /**
     * Computes the values of the function at all sampled locations as well as
     * all the gradients in one go.
     * 
     * @param atry
     *            The set of objective function parameters at which to evaluate
     *            the function.
     */
    protected void computePXs(double[] atry) {

        //Compute the value of p(x) at every x vector.
        for (int x = 0; x < xvecs.length; x++) {

            double expo = 0.0;

            //Add the other lambda contributions.
            for (int j = 1; j <= ma; j++) {
                expo += cpdx[j - 1][x] * atry[j];
            }

            //Note that we absorb the normalisation by the
            //area of the sphere in here.
            pxs[x] = ptArea * Math.exp(expo);
        }
    }

    /**
     * Checks to see if the parameters passed are the same as they were the last
     * time the function was evaluated.
     * 
     * @param atry
     *            The set of objective function parameters at which to evaluate
     *            the function.
     * 
     * @param i
     *            The index of the independent variable at which to evaluate the
     *            objective function.
     * 
     * @return test result.
     */
    protected boolean sameAsLast(double[] atry, int i) {
        boolean same = true;

        same &= (i == lastI);

        for (int j = 1; j < atry.length; j++) {
            same &= (atry[j] == lastATry[j]);
        }

        return same;
    }

    /**
     * Evaluates the function and its derivatives at the specified point.
     * 
     * @param atry
     *            The set of objective function parameters at which to evaluate
     *            the function.
     * 
     * @param i
     *            The index of the independent variable at which to evaluate the
     *            objective function.
     */
    protected void evaluateFunction(double[] atry, int i) {

        for (int j = 0; j < atry.length; j++) {
            lastATry[j] = atry[j];
        }
        lastI = i;

        for (int j = 0; j < derivs.length; j++) {
            derivs[j] = 0.0;
        }
        yVal = 0.0;

        for (int x = 0; x < xvecs.length; x++) {

            //p(x) for each x unit vector is stored in pxs.
            //R(k_i; x) for each k_i and x
            //is stored in ckdx.
            double contrib = pxs[x] * ckdx[i - 1][x];
            yVal += contrib;
            for (int j = 1; j <= ma; j++) {

                //Note that the indexes to the ks and ps start
                //counting at zero, whereas the corresponding lambdas
                //start from a[1].
                derivs[j] += cpdx[j - 1][x] * contrib;
            }
        }
    }

    /**
     * This method initialises the array of k-vectors - sampled locations or
     * gradient directions.
     * 
     * @param indepVals
     *            Normalized gradient directions with zeros removed.
     */
    protected void initKVecs(double[][] indepVals) {
        kvecs = new double[ndata][3];
        for (int i = 0; i < ndata; i++) {
            kvecs[i][0] = indepVals[i][0];
            kvecs[i][1] = indepVals[i][1];
            kvecs[i][2] = indepVals[i][2];

        }
    }

    protected void initParamKVecs(double[][] pVecVals) {
        if (pVecVals.length != ma - 1) {
            logger.severe("length of lambda summation array passed to initialiser"
                    + "is not of proscribed length! ma = " + ma + " but length is "
                    + pVecVals.length);
            throw new RuntimeException(
                    "length of lambda summation array passed to initialiser"
                            + "is not of proscribed length! ma = " + ma
                            + " but length is " + pVecVals.length);

        }
        pvecs = new double[ma][3];
        // Note that we need the zero at the start like kvecs.
        for (int i = 0; i < pvecs.length - 1; i++) {
            pvecs[i + 1][0] = pVecVals[i][0];
            pvecs[i + 1][1] = pVecVals[i][1];
            pvecs[i + 1][2] = pVecVals[i][2];
        }
    }

    /**
     * Member class used to store the sets of dot products between the k-vectors
     * and the x-vectors.
     */
    protected class DPSet {
        private double[][] data;

        public DPSet(double[][] d) {
            data = d;
        }

        public double[][] getData() {
            return data;
        }
    }

    /**
     * Initialises the database of point sets and the test point set.
     */
    protected void initPointSets() {
        pointSets = new Vector<SphericalPointSet>();
        ckdxs = new Vector<DPSet>();
        cpdxs = new Vector<DPSet>();
        // J1.4 pointSets = new Vector();
        // J1.4 ckdxs = new Vector();
        // J1.4 cpdxs = new Vector();
        //System.err.print("initialising point sets... ");
        for (int p = 0; p < ISCodes.getNoPointSetsForMaxEnt(); p++)
            try {
                //Add the i-th point set in the database to a vector.
                SphericalPointSet sps = ISCodes.getPointSetForMaxEnt(p);
                pointSets.addElement(sps);

                //Create the array of R(q, x)s
                double[][] dpsk = new double[ndata][sps.data.length];

                // Similarly for the R(p, x)s
                double[][] dpsp = new double[ma][sps.data.length];
                
                for (int x = 0; x < sps.data.length; x++) {
                    for(int i=0; i<ndata; i++) {

                        double[] qvec = kvecs[i];
                        double[] xvec = sps.data[x];

                        dpsk[i][x] = SphDeconvKernels.kernel(qvec, xvec, kernelParams);
                    }
                    for (int i = 0; i < ma; i++) {

                        double[] qvec = pvecs[i];
                        double[] xvec = sps.data[x];

                        dpsp[i][x] = SphDeconvKernels.kernel(qvec, xvec, kernelParams);
                    }
                }
                ckdxs.addElement(new DPSet(dpsk));
                cpdxs.addElement(new DPSet(dpsp));
            }
            catch (Exception e) {
                throw new RuntimeException(e);
            }
        //System.err.println("done.");

        //Set up the array of dot products for the test point set.
        ckdxTS = new double[ndata][testPS.length];
        cpdxTS = new double[ma][testPS.length];
        //System.err.print("initialising kernel... ");
        for (int x = 0; x < testPS.length; x++) {
            for (int i = 0; i < ndata; i++) {
                double[] qvec = kvecs[i];
                double[] xvec = testPS[x];

                ckdxTS[i][x] = SphDeconvKernels.kernel(qvec, xvec, kernelParams);
            }
            for(int i = 0; i<ma; i++){
                double[] qvec = pvecs[i];
                double[] xvec = testPS[x];

                cpdxTS[i][x] = SphDeconvKernels.kernel(qvec, xvec, kernelParams);
            }
        }
        //System.err.println("done.");
        // This is the weighting factor for numerical integrals
        // computed over the sphere using the test point set.
        ptAreaTS = computePT_Area(testPS.length);

    }

    /**
     * Computes the weighting factor for the numerical integrals over the sphere
     * given the number of points in the point set used.
     * 
     * @param numSampledPoints
     *            The number of points in the point set used to evaluate the
     *            numerical integrals.
     * 
     * @return 4 \pi / numSampledPoints.
     */
    protected double computePT_Area(int numSampledPoints) {
        return 4.0 * Math.PI / ((double) numSampledPoints);
    }

    /**
     * This method looks at the max. ent. function with the current set of
     * lambdas and computes the errors of each numerical integral computed with
     * the current integration point set and the test point set. The maximum
     * relative error is the first element of the array returned. The function
     * also computes the estimated value of the error based on the peakedness of
     * the function and comparison to a linear model of error tables of
     * integration of \exp(H.x). The estimated max. rel. error is the second
     * element of the array returned.
     * 
     * @return {Max. relative error compared to test point set, Max estimated
     *         error of test point set integral}
     */
    public double[] integrationErrorStats() {

        //Compute the value of p(x) at every x vector in the test set.
        for (int x = 0; x < testPS.length; x++) {
            double expo = 0.0;
            for (int j = 1; j <= ma; j++) {
                expo += cpdxTS[j - 1][x] * a[j];
            }
            pxsTS[x] = ptAreaTS * Math.exp(expo);
        }

        //Compute each integral with the new set of x's
        double[][] tsIntegrals = new double[ndata][ndata + 1];
        double[][] peakLengthScale = new double[ndata][ndata + 1];
        for (int i = 1; i <= ndata; i++) {

            //Compute the integral and the max of the function.
            double[] maxVals = new double[ndata + 1];
            for (int x = 0; x < testPS.length; x++) {
                double contrib = pxsTS[x] * ckdxTS[i - 1][x];

                //This is the integral of the function itself.
                tsIntegrals[i - 1][0] += contrib;
                maxVals[0] = (maxVals[0] > contrib) ? maxVals[0] : contrib;

                for (int j = 1; j <= ma; j++) {

                    //These are the integrals of the derivatives.
                    double dcontrib = cpdxTS[j - 1][x] * contrib;
                    tsIntegrals[i - 1][j] += dcontrib;
                    maxVals[j] = (maxVals[j] > dcontrib) ? maxVals[j] : dcontrib;
                }
            }

            //maxVal includes the scaling by ptAreaTS, which we need to remove.
            for (int j = 0; j <= ndata; j++) {
                maxVals[j] /= ptAreaTS;
            }

            //The lengthscale of the peak of each function is obtained
            //by looking at its integral divided by the peak value.
            //This is used to predict the expected error.
            //Note this is differs from the superclass as we now have
            //the PAS defined on the unit sphere.
            for (int j = 0; j <= ndata; j++) {
                peakLengthScale[i - 1][j] = Math.sqrt(tsIntegrals[i - 1][j]
                        / (maxVals[j]));
            }
        }

        //Compute each integral with the current set of x's.
        computePXs(a);
        double[][] curIntegrals = new double[ndata][ndata + 1];
        for (int i = 1; i <= ndata; i++) {
        //for (int i = 1; i <= ma; i++) {
            double[] numints = computeFIntegralDerivs(i);
            for (int j = 0; j < numints.length; j++) {
                curIntegrals[i - 1][j] = numints[j];
            }
        }

        //Compute the maximum relative error of the old integrals, assuming
        //the new ones are accurate.
        double maxRelErr = 0.0;

        //Also compute the maximum expected error.
        double expMRE = 0.0;

        for (int i = 0; i < curIntegrals.length; i++) {
            for (int j = 0; j < curIntegrals[0].length; j++) {
                double relErr = Math.abs(curIntegrals[i][j] - tsIntegrals[i][j])
                        / tsIntegrals[i][j];

                //Estimate the expected errors on the integrals given
                //the lengthscale of the function peak and the number
                //of points.
                double expErr = estimateErrors350(peakLengthScale[i][j]);

                maxRelErr = (relErr > maxRelErr) ? relErr : maxRelErr;
                expMRE = (expErr > expMRE) ? expErr : expMRE;
            }
        }

        double[] results = new double[2];
        results[0] = maxRelErr;
        results[1] = expMRE;
        return results;
    }

    /**
     * Evaluates the numerical integrals of \exp(\lambda_0 + \sum_i (\lambda_i
     * cos(x \cdot k_i)) cos(x \cdot k) over the sphere |x|=r for k = k_i. The
     * derivatives wrt each lambda are also computed and the whole lot returned
     * in an array. computePXs must be called before calling this function.
     * 
     * @param i
     *            The index of the independent value for which to compute the
     *            integrals.
     * 
     * @return {\int dp/d\lambda_0 dx, \int dp/d\lambda_1 dx, ..., \int
     *         dp/d\lambda_N dx}
     */
    protected double[] computeFIntegralDerivs(int i) {

        double[] numints = new double[ndata + 1];

        for (int x = 0; x < xvecs.length; x++) {

            //p(x) for each x unit vector is stored in pxs.
            //cos(k_i . x) for each k_i and x
            //is stored in ckdx.
            double contrib = pxs[x] * ckdx[i - 1][x];

            //This is the value of the function itself.
            numints[0] += contrib;

            //These are the derivatives.
            for (int j = 1; j <= ma; j++) {
                numints[j] += contrib * cpdx[j - 1][x];
            }

        }

        return numints;
    }

    /**
     * This function applies the parameters of a linear regression to esimate
     * the maximum error in the estimate of the integral of a function with a
     * specified peak length scale using a z-theta 3:1 point set with 350^2
     * elements.
     * 
     * @param lengthscale
     *            \sqrt(\int f(x) dx /max(f(x)))
     * 
     * @return Estimate of numerical integration error.
     */
    protected double estimateErrors350(double lengthscale) {

        //Parameters obtained by linear regression of the lengthscale
        //of the peak of \exp(H.x) with the maximum relative error in
        //numerical integration.
        double cLS = -3.06354;
        double cBase = -12.3513;

        //Note that |H| = lengthscale ^ -2.
        double predErr = Math.exp(cLS * Math.log(lengthscale) + cBase);

        return predErr;
    }

}

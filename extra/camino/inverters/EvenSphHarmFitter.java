    package inverters;

    import numerics.*;
    import imaging.*;
    import misc.*;
    import tools.ArrayOps;
    import Jama.Matrix;

    /**
     * <dl>
     * 
     * <dt>Purpose:
     * 
     * <dd>Fits an even spherical harmonic series to diffusion MRI data.
     * 
     * <dt>Description:
     * 
     * <dd>Computes and stores the complex matrix for linear mapping of a set of
     * measurements in a single voxel to the coefficients of a spherical harmonic
     * series. This class also contains methods for finding a suitable order of
     * truncation of the series using the ANOVA test for deletion of variables.
     * 
     * </dl>
     * 
     * @author Danny Alexander
     * @version $Id$  
     */
    public class EvenSphHarmFitter {

        /**
         * This is the inversion matrix.
         */
        protected RealMatrix invX;

        /**
         * This is the forward mapping matrix.
         */
        protected RealMatrix X;

        /**
         * This is the order that the fitted series goes up to.
         */
        protected int maxOrder;

        /**
         * This is the number of free parameters in the series.
         */
        protected int numParams;

        /**
         * These are the parameters of the imaging sequence.
         */
        protected DW_Scheme ip;


        protected double meanNonZeroB;

        /**
         * Zero, two and four-level thresholds in F-test for series truncation.
         */
        private double FDT1 = 1.0E-20;

        private double FDT2 = 1.0E-7;

        private double FDT3 = 1.0E-7;

        /**
         * This is the threshold on the singular values. If a singular value has
         * value less that the largest singular value divided by this value, it is
         * inverted to zero.
         */
        protected static final double SVTHRESH = 1.0E12;


        /**
         * Default constructor required for inheritance.
         */
        public EvenSphHarmFitter() {
        }


        /**
         * The constructor computes the inversion matrix from the imaging parameters
         * for the sequence used in the data acquisition.
         * 
         * @param imParams
         *            The imaging parameters of the acquisition sequence.
         * 
         * @param order
         *            The maximum order of the spherical harmonic series to fit.
         */
        public EvenSphHarmFitter(DW_Scheme imParams, int order) {

            maxOrder = order;
            ip = imParams;

            meanNonZeroB = ArrayOps.mean(ip.getNonZeroB_Values());

            // Compute the number of free parameters in the spherical
            // harmonic coefficients for a real even series up to the
            // specified order (see MRM02).
            numParams = (maxOrder + 1) * (maxOrder / 2 + 1);

            // Set up the matrix X.
            X = new RealMatrix(ip.numMeasurements(), numParams + 1);

            for (int i = 0; i < X.rows(); i++) {
                Vector3D g = new Vector3D(ip.getG_Dir(i));
                
     	        double[] tp = Vector3D.thetaPhi(g);

     	        double theta = tp[0];
     	        double phi = tp[1];

    //             double[] sphPolQ = SphericalPoints.getSphPolars(new double[] {g.x, g.y, g.z});

    //             double theta = sphPolQ[1];
    //             double phi = sphPolQ[2];

                double b = ip.getB_Value(i);

                X.setEntry(i, 0, 1.0);

                int nextEntry = 1;

                for (int l = 0; l <= maxOrder; l += 2) {
                    try {
                        Complex c = SphericalHarmonics.Y(l, 0, theta, phi);
                        X.setEntry(i, nextEntry, -b * c.real());
                        nextEntry += 1;

                        for (int m = 1; m <= l; m++) {
                            c = SphericalHarmonics.Y(l, m, theta, phi);
                            X.setEntry(i, nextEntry, -2.0 * b * c.real());
                            X.setEntry(i, nextEntry + 1, 2.0 * b * c.imag());
                            nextEntry += 2;
                        }
                    }
                    catch (Exception e) {
                        // This should never happen.
                        throw new RuntimeException(e);
                    }
                }    
            }

            // Now get the singular value decomposition.
            RealMatrix[] svd = null;
            try {
                svd = X.svd();
            }
            catch (Exception e) {
                throw new RuntimeException(e);
            }

            // Find the maximum singular value.
            double maxSV = svd[1].entry(0, 0);
            for (int i = 0; (i < svd[1].rows() && i < svd[1].columns()); i++) {
                if (svd[1].entry(i, i) > maxSV) {
                    maxSV = svd[1].entry(i, i);
                }
            }

            // Invert the singular values.
            for (int i = 0; (i < svd[1].rows() && i < svd[1].columns()); i++) {
                double curSV = svd[1].entry(i, i);
                svd[1].setEntry(i, i, (curSV > maxSV / SVTHRESH) ? (1.0 / curSV) : 0.0);
            }

            // Get the pseudo inverse.
            invX = svd[2].product(svd[1].transpose()).product(svd[0].transpose());

        }


        /**
         * Fits the spherical harmonic series using the inverse matrix.
         * 
         * @param data
         *            The MRI data.
         * 
         * @return The fitted series. The array contains [exitCode, log A^\star(0),
         *         c00, c20, Re(c21), Im(c21), Re(c22), Im(c22), c40, Re(c41),
         *         Im(c41), ...]
         */
        public double[] fit(double[] data) {

            // The failure free exit code is 0.
            double exitCode = 0.0;

            // Make the data matrix.
            RealMatrix A = new RealMatrix(data.length, 1);
            for (int i = 0; i < data.length; i++) {
                if (data[i] > 0) {
                    A.setEntry(i, 0, Math.log(data[i]));
                }
                else {
                    A.setEntry(i, 0, 0.0);

                    // The exit code is 6 if the data is bad and has to
                    // be changed to perform the inversion.
                    exitCode = 6.0;
                }
            }

            // Do the inversion.
            RealMatrix S = invX.product(A);

            // Create the return array.
            double[] res = new double[numParams + 2];
            res[0] = exitCode;
            for (int i = 0; i < numParams + 1; i++) {
                res[i + 1] = S.entry(i, 0);
            }

            return res;
        }

       /**
         * Fits the spherical harmonic series after applying gradient corrections to 
         * the inverse matrix. This is for use with Human Connectome Project data.
         * 
         * @param data
         *            The MRI data.
         * @param  gradAdj specifies the gradient adjustment per voxel. 
         * 
         * @return The fitted series. The array contains [exitCode, log A^\star(0),
         *         c00, c20, Re(c21), Im(c21), Re(c22), Im(c22), c40, Re(c41),
         *         Im(c41), ...]
         */
        public double[] fit(double[] data, double[] gradAdj) {
            // Gradient adjustment matrix.
            Matrix gmat = new Matrix(gradAdj,3);
            Matrix eye = Matrix.identity(3,3);
            Matrix iPlusG = eye.plus(gmat);
            // Apply gradient adjustment to forward matrix looping through rows.
            for ( int i = 0; i < X.rows(); i++ ) {

                double[] bvec = ip.getG_Dir(i);
                Matrix bmat = new Matrix(bvec,bvec.length);
                // compute the new bvec
                Matrix v = iPlusG.times(bmat);
                double[][] a = v.getArray();
                // magnitude of new bvec
                double n = Math.sqrt((a[0][0]*a[0][0])+(a[1][0]*a[1][0])+(a[2][0]*a[2][0]));
                // normalise and set new bvecs
                double[] new_gDir = new double[]{a[0][0]/n, a[1][0]/n, a[2][0]/n};
                // Modify bValue
                double bval = ip.getB_Value(i);
                double b = (n*n)*bval;;

                // same as constructor from here
                Vector3D g = new Vector3D(new_gDir);
                
                double[] tp = Vector3D.thetaPhi(g);

                double theta = tp[0];
                double phi = tp[1];
            
                // Set X values same way as in constructor.
                X.setEntry(i, 0, 1.0);

                int nextEntry = 1;

                for (int l = 0; l <= maxOrder; l += 2) {
                    try {
                        Complex c = SphericalHarmonics.Y(l, 0, theta, phi);
                        X.setEntry(i, nextEntry, -b * c.real());
                        nextEntry += 1;

                        for (int m = 1; m <= l; m++) {
                            c = SphericalHarmonics.Y(l, m, theta, phi);
                            X.setEntry(i, nextEntry, -2.0 * b * c.real());
                            X.setEntry(i, nextEntry + 1, 2.0 * b * c.imag());
                            nextEntry += 2;
                        }
                    }
                    catch (Exception e) {
                        // This should never happen.
                        throw new RuntimeException(e);
                    }
                }
            }

            // Now get the singular value decomposition.
            RealMatrix[] svd = null;
            try {
                svd = X.svd();
            }
            catch (Exception e) {
                throw new RuntimeException(e);
            }

            // Find the maximum singular value.
            double maxSV = svd[1].entry(0, 0);
            for (int i = 0; (i < svd[1].rows() && i < svd[1].columns()); i++) {
                if (svd[1].entry(i, i) > maxSV) {
                    maxSV = svd[1].entry(i, i);
                }
            }

            // Invert the singular values.
            for (int i = 0; (i < svd[1].rows() && i < svd[1].columns()); i++) {
                double curSV = svd[1].entry(i, i);
                svd[1].setEntry(i, i, (curSV > maxSV / SVTHRESH) ? (1.0 / curSV) : 0.0);
            }

            // Get the pseudo inverse.
            invX = svd[2].product(svd[1].transpose()).product(svd[0].transpose());

            // Now call the regular fit method. 
            return fit(data);

        }

        /**
         * Specifies the number of elements in the output array.
         * 
         * @return The number of elements in the output array.
         */
        public int itemsPerVoxel() {
            return numParams + 2;
        }

        /**
         * Sets the f-test thresholds.
         * 
         * @param f1
         *            The 0 versus higher order threshold.
         * 
         * @param f2
         *            The 2 versus higher order threshold.
         * 
         * @param f3
         *            The 4+ versus higher order threshold.
         */
        public void setF_TestThresholds(double f1, double f2, double f3) {
            FDT1 = f1;
            FDT2 = f2;
            FDT3 = f3;
        }

        /**
         * Returns the maximum spherical harmonic order that the fitter uses.
         */
        public int getMaxOrder() {
            return maxOrder;
        }

        /**
         * Determines the level of truncation of the series using a series of
         * f-tests between models from truncations at different orders.
         * 
         * @param params
         *            The fitted spherical harmonic parameters. The array is the
         *            output of EvenSphHarmFitter.fit.
         * 
         * @param data
         *            The measurements.
         * 
         * @return The order of truncation of the series.
         */
        public int truncate(double[] params, double[] data) {

            // Compute the array of probabilities of equivalence between
            // each pair of models.
            double[] p = getF_TestProbabilities(params, data);

            // Choose the truncation level
            return selectModel(p, FDT1, FDT2, FDT3);

            // 	// Perform F-test.
            // 	int dataPoints = data.length;
            // 	int trunc = runRepeatedF_Tests(ssqdtrs, ssqars, dataPoints);
            // 	//int trunc2 = runRepeatedF_Tests(ssqdtrs2, ssqars, dataPoints);

            // 	return trunc;
        }

        /**
         * Computes the set of probabilities of each pair of truncatated spherical
         * harmonic models being equivalent. These are the probabilities thresholded
         * in the f-test for model selection.
         * 
         * @param params
         *            The fitted spherical harmonic parameters. The array is the
         *            output of EvenSphHarmFitter.fit.
         * 
         * @param data
         *            The measurements.
         * 
         * @return The array of probabilities. The structure of the array is p02,
         *         p04, ..., p0n, p24, ..., p2n, ..., p(n-2)n, where n is the
         *         largest order of the series.
         */
        public double[] getF_TestProbabilities(double[] params, double[] data) {

            // Compute a matrix of fitted values for every non-zero sample
            // point at each level of truncation of the series.
            RealMatrix fittedValues = getFittedValues(params);

            // The first two elements of params are the exit code
            // and the estimated log(A^\star(0)).
            double ls0 = params[1];

            // Need the log normalized measurements to compute the
            // least squares differences.
            double[] logNormData = new double[data.length];
            for (int i = 0; i < logNormData.length; i++) {
                double logData = (data[i] > 0) ? Math.log(data[i]) : 0.0;
                logNormData[i] = logData - ls0;
            }

            // Compute the statistics required for the f-test.
            double[] ssqars = getSumSqsAboutRegression(logNormData, fittedValues);
            double[] ssqdtrs = getSumSqsDueToRegression(params, logNormData.length);

            //	for(int i=0; i<ssqdtrs.length; i++) {
            //	    System.out.println("SSQDTR " + i + ": " + ssqdtrs[i]);
            //	    System.out.println("SSQDTR2 " + i + ": " + ssqdtrs2[i]);
            //	    System.out.println("SSQAR " + i + ": " + ssqars[i]);
            //	}

            double[] p = modelProbabilities(ssqdtrs, ssqars, logNormData.length);
            return p;
        }

        /**
         * Constructs matrix of measurements approximations from the model obtained
         * from truncating the even spherical harmonic series at each even level.
         * 
         * @param params
         *            The spherical harmonic coefficients. The structure of the
         *            array is as output by EvenSphHarmFitter.fit.
         * 
         * @return The matrix of approximations from each model.
         */
        private RealMatrix getFittedValues(double[] params) {

            // The first two elements of params are the exit code
            // and the estimated log(A^\star(0)).
            double ls0 = params[1];

            RealMatrix fittedValues = new RealMatrix(X.rows(), maxOrder / 2 + 1);

            for (int i = 0; i < X.rows(); i++) {
                double bdx = 0.0;
                int nextInd = 0;
                for (int j = 0; j < maxOrder + 1; j += 2) {

                    //The number of parameters at each level is
                    //1, 5, 9, 13, ...
                    for (int k = 0; k < 2 * j + 1; k++) {
                        bdx = bdx + params[nextInd + 2] * X.entries[i][nextInd + 1];
                        nextInd += 1;
                    }

                    //Note log S(x) = log S_0 - b D(x). bdx contains the
                    //diffusion coefficient scaled by -b (in eshMat).
                    fittedValues.entries[i][j / 2] = bdx;
                }
            }

            return fittedValues;
        }

        /**
         * Computes the sums of squares due to and about each regression. These are
         * statistics required for the f-tests in the ANOVA model selection.
         * 
         * This is redundant now and is replaced by getSumSqsAboutRegression and
         * getSumSqsDueToRegression. The latter computes the variance of the model
         * analytically rather than numerically, which gives a slight improvement in
         * overall performance.
         * 
         * @param logNormData
         *            Array of log normalized measurements.
         * 
         * @param fittedValues
         *            Matrix of approximations to each measurement from each model.
         * 
         * @param ssqdtrs
         *            Array in which to place the sums of squares due to regression.
         * 
         * @param ssqars
         *            Array in which to place the sums of squares about the
         *            regression.
         */
        private void getRegressionStats(double[] logNormData, RealMatrix fittedValues,
                double[] ssqdtrs, double[] ssqars) {

            int totalPoints = logNormData.length;
            int noModels = fittedValues.columns();

            for (int j = 0; j < noModels; j++) {
                ssqdtrs[j] = 0.0;
                ssqars[j] = 0.0;
            }

            for (int i = 0; i < totalPoints; i++) {
                for (int j = 0; j < noModels; j++) {

                    // Compute the variance of the model (sum of squares
                    // "due to" the regression.
                    double edt = fittedValues.entries[i][0] - fittedValues.entries[i][j];

                    // Compute the sums of squared errors for each model
                    // ("about" the regression).
                    double ea = logNormData[i] - fittedValues.entries[i][j];

                    ssqdtrs[j] += edt * edt;
                    ssqars[j] += ea * ea;

                }
            }
        }

        /**
         * Computes the sums of squares about each regression. These are statistics
         * required for the f-tests in the ANOVA model selection.
         * 
         * @param logNormData
         *            Array of log normalized measurements.
         * 
         * @param fittedValues
         *            Matrix of approximations to each measurement from each model.
         */
        private double[] getSumSqsAboutRegression(double[] logNormData,
                RealMatrix fittedValues) {

            int totalPoints = logNormData.length;
            int noModels = fittedValues.columns();

            double[] ssqars = new double[noModels];

            for (int i = 0; i < totalPoints; i++) {
                for (int j = 0; j < noModels; j++) {

                    // Compute the sums of squared errors for each model
                    // ("about" the regression).
                    double ea = logNormData[i] - fittedValues.entries[i][j];

                    ssqars[j] += ea * ea;

                }
            }

            return ssqars;
        }

        /**
         * This method computes the sum of squares due to regression for the
         * spherical harmonic models analytically.
         * 
         * @param params
         *            The list of spherical harmonic coefficients computed by
         *            EvenSphHarmFitter.fit.
         * 
         * @param numPoints
         *            The number of measurements.
         * 
         * @return A list containing the sums of squares due to regression for each
         *         model.
         */
        private double[] getSumSqsDueToRegression(double[] params, int numPoints) {

            double[] results = new double[maxOrder / 2 + 1];

            // First index is 2, as the first two elements of params
            // are the exit code and log zero measurement.
            int index = 2;

            // Compute the mean of the models (all the same).
            double M0 = params[index] / (2.0 * Math.sqrt(Math.PI));

            // Initialize the mean square of the model.
            double V = 0.0;

            // We will need to know the number of zero measurements.
            int numZeros = ip.numZeroMeasurements();

            for (int l = 0; l <= maxOrder; l += 2) {

                // Just one of the first parameter (c_l0)
                V += params[index] * params[index];
                index += 1;

                // Note we need to add both the real and imaginary
                // components of the spherical harmonic terms.
                for (int m = 1; m <= 2 * l; m++) {

                    // Two of all the others.
                    V += 2.0 * params[index] * params[index];
                    index += 1;
                }

                results[l / 2] = ((double) (numPoints - numZeros))
                        * (V / (4.0 * Math.PI) - M0 * M0);
            }

            for (int i = 0; i < results.length; i++) {

                results[i] = results[i] * meanNonZeroB * meanNonZeroB;
            }

            return results;
        }

        /**
         * Computes the f-statistic probabilities from the arrays of regression
         * statistics output by the methods above.
         * 
         * @param ssqdtrs
         *            A list of sum of squares due to regression (variances of the
         *            models).
         * 
         * @param ssqars
         *            A list of sum of squares about the regression (sum of squared
         *            errors).
         * 
         * @param dataPoints
         *            The number of measurements modelled.
         * 
         * @return The array of probabilities. The structure of the array is p02,
         *         p04, ..., p0n, p24, ..., p2n, ..., p(n-2)n, where n is the
         *         largest order of the series.
         */
        private double[] modelProbabilities(double[] ssqdtrs, double[] ssqars, int dataPoints) {

            int numModels = ssqars.length;
            int numModelComparisons = numModels * (numModels - 1) / 2;
            double[] probs = new double[numModelComparisons];

            // Start with zero-th order model;
            int params1 = 0;

            // Keeps track of the next index of the probs array.
            int nextIndex = 0;

            // Compare best with each subsequent even order model.
            for (int i = 0; i < numModels - 1; i++)
                try {

                    // Update the parameter count for model 1.
                    params1 += 4 * i + 1;

                    // For the second model, start with the same as the first.
                    int params2 = params1;

                    for (int j = i + 1; j < numModels; j++) {

                        params2 += 4 * j + 1;

                        // Compute F-statistic.
                        double dof1 = (double) (params2 - params1);
                        double dof2 = (double) (dataPoints - params2 - 1);
                        double f = ((ssqdtrs[j] - ssqdtrs[i]) / dof1) / (ssqars[j] / dof2);

                        // Compute probability. Wary of outlying values of f,
                        // which can result from the higher order model giving a
                        // higher error, or from very little improvement.
                        double betaArg = dof2 / (f * (dof1 + dof2));
                        double p = (betaArg > 0.0 && betaArg < 1.0) ? IncompleteBeta.betai(
                                dof2 / 2.0, dof1 / 2.0, betaArg) : 1.0;

                        probs[nextIndex] = p;
                        nextIndex += 1;
                    }

                }
                catch (Exception e) {
                    throw new RuntimeException(e);
                }

            return probs;
        }

        /**
         * Performs a series of F-tests given the probabilities of equivalence of
         * each pair of spherical harmonic models.
         * 
         * @param p
         *            The array of probabilities output by getF_TestProbabilities.
         * 
         * @param T1
         *            The f-test threshold for current model 0.
         * 
         * @param T2
         *            The f-test threshold for current model 2.
         * 
         * @param T3
         *            The f-test threshold for current model 4 or greater.
         * 
         * @return The selected truncation order.
         */
        public int selectModel(double[] p, double T1, double T2, double T3) {
            return selectModel(p, maxOrder, T1, T2, T3);
        }

        /**
         * Performs a series of F-tests given the probabilities of equivalence of
         * each pair of spherical harmonic models.
         * 
         * @param p
         *            The array of probabilities output by getF_TestProbabilities.
         * 
         * @param order
         *            The maximum order in the series used to construct the array p.
         * 
         * @param T1
         *            The f-test threshold for current model 0.
         * 
         * @param T2
         *            The f-test threshold for current model 2.
         * 
         * @param T3
         *            The f-test threshold for current model 4 or greater.
         * 
         * @return The selected truncation order.
         */
        public static int selectModel(double[] p, int order, double T1, double T2, double T3) {

            int numModels = order / 2 + 1;

            // Start with zero-th order model;
            int currentBest = 0;
            int bestBaseIndex = 0;

            // Compare best with each subsequent model.
            for (int i = 1; i < numModels; i++) {

                //Compare to decision threshold and update best model if
                //appropriate.
                int pIndex = bestBaseIndex + i - currentBest / 2 - 1;
                if (currentBest == 0 && p[pIndex] < T1) {
                    currentBest = 2 * i;
                    bestBaseIndex = numModels * i - i * (i + 1) / 2;
                }
                else if (currentBest == 2 && p[pIndex] < T2) {
                    currentBest = 2 * i;
                    bestBaseIndex = numModels * i - i * (i + 1) / 2;
                }
                else if (currentBest > 2 && p[pIndex] < T3) {
                    currentBest = 2 * i;
                    bestBaseIndex = numModels * i - i * (i + 1) / 2;
                }

            }

            return currentBest;
        }

    }

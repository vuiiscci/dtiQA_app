package inverters;

import imaging.*;
import misc.*;
import tools.*;
import java.io.*;

import optimizers.MarquardtMinimiserException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fit the diffusion tensor with outlier rejection.
 * 
 * <dt>Description:
 * 
 * <dd>Fits the diffusion tensor using the RESTORE method in Chang, Jones and
 * Pierpaoli MRM 53 2005.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class RestoreDT_Inversion extends DT_Inversion {

    /**
     * The non-linear diffusion tensor fitter.
     */
    protected DiffTensorFitter fitter;

    /**
     * Linear diffusion tensor fitter to provide starting points for non-linear
     * optimization.
     */
    protected LinearDT_Inversion ldti;

    /**
     * The standard deviation of the noise in the data must be provided for the
     * restore algorithm.
     */
    protected double sigma;

    /**
     * Need the full set of qs to be global.
     */
    protected double[][] gNoZeros;

    /**
     * The b-values.
     */
    protected double[] bValues;

    /**
     * The number of zero measurements.
     */
    protected int M;

    /**
     * The index of the fitter used throughout.
     */
    protected int fitterIndex;

    /**
     * The factor for the MAD estimator of C in the Geman-McClure M-estimator.
     * Value taken from Chang et al.
     */
    private final double MADFACTOR = 1.4826;

    /**
     * Number of standard deviations from the mean for outlier threshold.
     */
    private final double MAXRESTHRESH = 3.0;


    /**
     *
     * Maximum iterations of the iterative reweighting procedure
     */ 
    private final int MAXITERATIONS = 100;


    /**
     * Convergence threshold on the change in the max residual in the RESTORE
     * iterative reweighting scheme.
     */
    private final double RESCHANGETHRESH = 0.001;

    /**
     * Stores the histogram of outliers.
     */
    private int[] outlierHist;

    /**
     * Output stream for the outlier map if specified.
     */
    private DataOutputStream outlierMap;

  
  /**
     * Default constructor.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param noiseLevel
     *            The standard deviation of the noise in the data.
     */
    public RestoreDT_Inversion(DW_Scheme imParams, double noiseLevel) {
        ip = imParams;
        sigma = noiseLevel;
        init(ip, 4);
    }

    /**
     * This constructor takes the imaging sequence parameters and an index
     * specifying the type of non-linear fitting to use.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param noiseLevel
     *            The standard deviation of the noise in the data.
     * 
     * @param index
     *            The fitter index.
     */
    public RestoreDT_Inversion(DW_Scheme imParams, double noiseLevel, int index) {
        ip = imParams;
        sigma = noiseLevel;
        init(ip, index);
    }


    /**
     * Sets up the fitter to use for inversion.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param index
     *            The fitter index.
     */
    protected void init(DW_Scheme imParams, int index) {

        ip = imParams;
        fitterIndex = index;

        // Set up the linear inversion.
        ldti = new LinearDT_Inversion(ip);

        // Set up the non-linear fitter.
        M = ip.numZeroMeasurements();
        
	bValues = ip.getNonZeroB_Values();

        gNoZeros = ip.getNonZeroG_Dirs();

        // Initialize it with dummy data.
        double[] dummyData = new double[ip.numMeasurements() - M];

        fitter = getFitter(gNoZeros, dummyData, bValues, M, fitterIndex);

        // Initialize the histogram
        outlierHist = new int[gNoZeros.length];

        // Initialize the outlier map output stream.
        if (CL_Initializer.outlierMapFile != null){
            try {
                outlierMap = new DataOutputStream(new FileOutputStream(
                        CL_Initializer.outlierMapFile));
            }
            catch (FileNotFoundException fnfe) {
                logger.warning(""+fnfe);
            }
        }
    }

    /**
     * Fits the diffusion tensor.
     * 
     * @param data
     *            The MRI data.
     * 
     * @return {exitcode, log A*(0), Dxx, Dxy, Dxz, Dyy, Dyz, Dzz}
     */
    public double[] invert(double[] data) {

        double exitCode = 0.0;

        // First do the linear inversion.
        double[] linDT = ldti.invert(data);

        if (linDT[0] != 0.0) {
            // problems in the data, probably can't do restore on this voxel

            // Make sure the outlier map is complete
            if (outlierMap != null) {
                outputOutlierMap(new int[0]);
            }

            return linDT;
        }

        // Set the starting point for the non-linear fitting.
        double[] nlStart = new double[6];
        for (int i = 0; i < nlStart.length; i++) {
            nlStart[i] = linDT[i + 2];
        }

        // Get the normalized data.
        double[] normData = ip.normalizeData(data);

        // Compute noise std for normalized data.
        double aStarZero = ip.geoMeanZeroMeas(data);
        double normSigma = sigma / aStarZero;

        // Initialize to the linear fitted DT.
        DT d = new DT(linDT[2], linDT[3], linDT[4], linDT[5], linDT[6], linDT[7]);

        try {

            // Initialize the data and starting point.
            fitter.newDepVals(normData);
            fitter.setInitParams(nlStart);

            // Initialize the weights in the non-linear fitter.
            for (int i = 1; i <= normData.length; i++) {
                fitter.setSig(i, normSigma);
            }

            fitter.minimise();
            d = fitter.getDiffTensor();

           // Check the residuals for potential outliers.
            double[] residuals = fitter.getResiduals();
            double maxResidual = ArrayOps.max(residuals);

            if (maxResidual > MAXRESTHRESH * normSigma) {

                // Invoke the iterative reweighting procedure to
                // highlight outliers.
                // Update weights until convergence, which we
                // detect as insignificant change in the
                // largest residual.
                double lastMaxRes = 0.0;

                int iterationCounter = 0;

                while (Math.abs(lastMaxRes - maxResidual) > RESCHANGETHRESH * maxResidual && 
                       iterationCounter < MAXITERATIONS) {

                    lastMaxRes = maxResidual;

                    // Compute the Geman-McClure weighting factor C
                    double medRes = ArrayOps.median(residuals);
                    double[] resSep = new double[residuals.length];
                    for (int i = 0; i < resSep.length; i++) {
                        resSep[i] = Math.abs(residuals[i] - medRes);
                    }
                    double C = MADFACTOR * ArrayOps.median(resSep);

                    fitter.newDepVals(normData);
                    fitter.setInitParams(nlStart);

                    // Set the robust weights in the SSD
                    for (int i = 1; i <= normData.length; i++) {
                        fitter.setSig(i, Math.sqrt(residuals[i - 1] * residuals[i - 1]
                                + C * C));
                    }

                    // Rerun the fitting process
                    fitter.minimise();

                    // Recompute the residuals
                    residuals = fitter.getResiduals();
                    maxResidual = ArrayOps.max(residuals);

                    iterationCounter += 1;
                }

                boolean converged = Math.abs(lastMaxRes - maxResidual) < RESCHANGETHRESH * maxResidual;

                if (!converged) {
                    logger.info("Maximum iterations reached without convergence. " + 
                                "Final tensor:\n\t" + d);
                    exitCode = 2.0;
                }

                // Get the list of outliers.
                int[] outlierIndices = getOutliers(residuals, normSigma);

                // Check there is enough data left to fit the tensor
                if(outlierIndices.length > residuals.length - 6) {
                    logger.info("Too many outliers.  Outputting linearly fitted DT.");
                    d = new DT(linDT[2], linDT[3], linDT[4], linDT[5], linDT[6], linDT[7]);
                    exitCode = 2.0;
                }
                else {

                    // Finally exclude the outliers and refit the
                    // diffusion tensor.
                    double[][] gNoZerosRobust = new double[gNoZeros.length
                                                           - outlierIndices.length][3];
                    double[] normDataRobust = new double[normData.length
                                                         - outlierIndices.length];

		    double[] bRobust = new double[normData.length - outlierIndices.length];

                    int nextGNZ = 0;
                    int nextOutlier = 0;
                    for (int i = 0; i < residuals.length; i++) {
                        if (nextOutlier < outlierIndices.length
                            && i == outlierIndices[nextOutlier]) {
                            nextOutlier++;
                        }
                        else {
                            gNoZerosRobust[nextGNZ][0] = gNoZeros[i][0];
                            gNoZerosRobust[nextGNZ][1] = gNoZeros[i][1];
                            gNoZerosRobust[nextGNZ][2] = gNoZeros[i][2];
                            normDataRobust[nextGNZ] = normData[i];
			    bRobust[nextGNZ] = bValues[i];
                            nextGNZ += 1;
                        }
                    }
                    DiffTensorFitter tempFitter = getFitter(gNoZerosRobust, normDataRobust,
                                                            bRobust, M, fitterIndex);
                    tempFitter.setInitParams(nlStart);
                    tempFitter.minimise();
                    d = tempFitter.getDiffTensor();

                }

                // Output the outlier map.
                if (outlierMap != null) {
                    outputOutlierMap(outlierIndices);
                }
                
                // set exit code to usual RESTORE value if fitting has not failed somehow
                if (exitCode == 0.0) {
                    // The error code is 1000 + the number of removed outliers.
                    exitCode = 1000 + outlierIndices.length;
                }

                // Update the outlier histogram
                for (int i = 0; i < outlierIndices.length; i++) {
                    outlierHist[outlierIndices[i]] += 1;
                }

            }
            else {
                if (outlierMap != null) {
                    outputOutlierMap(new int[0]);
                }
            }
        }
        catch (MarquardtMinimiserException e) {
            logger.info(e.toString() + "Fitting failed.  Outputting best DT found, which may be affected by outliers.");

            // If the fitter fails, the exit code is 2.
            exitCode = 2.0;

            // Make sure the outlier map is complete
            if (outlierMap != null) {
                outputOutlierMap(new int[0]);
            }
        }

        double[] dArr = d.getComponents();

        // Construct the array to return. The first element is
        // the exit code. The second element is ln A^\star(0).
        double[] res = new double[8];
        res[0] = exitCode;
        double geoMean = ip.geoMeanZeroMeas(data);
        res[1] = (geoMean>0)?Math.log(geoMean):0.0;
        for (int i = 0; i < 6; i++) {
            res[i + 2] = dArr[i];
        }

        return res;
    }

    /**
     * Overrides the default to output the outlier histogram.
     */
    public void close() {
        System.err.println("Outlier histogram:");
        for (int i = 0; i < outlierHist.length; i++) {
            System.err.println(i + " : " + outlierHist[i]);
        }
        if (outlierMap != null) try {
            outlierMap.close();
        }
        catch (Exception e) {
            logger.warning("" + e);
        }
    }

    /**
     * In background voxels this inverter needs to output the outlier map.
     */
    public void background() {
        if (outlierMap != null) {
            outputOutlierMap(new int[0]);
        }
    }

    /**
     * Finds the outliers by thresholding on residuals.
     * 
     * @param residuals
     *            Array of residuals of each data point.
     * 
     * @return A list of indices of outliers.
     */
    protected int[] getOutliers(double[] residuals, double normSigma) {
        int numOutliers = 0;
        for (int i = 0; i < residuals.length; i++) {
            if (residuals[i] > MAXRESTHRESH * normSigma) numOutliers += 1;
        }

        int[] outliers = new int[numOutliers];
        int next = 0;
        for (int i = 0; i < residuals.length; i++) {
            if (residuals[i] > MAXRESTHRESH * normSigma) {
                outliers[next] = i;
                next += 1;
            }
        }

        return outliers;
    }

    /**
     * Outputs the outlier map for one voxel.
     * 
     * @param outlierIndices
     *            The list of indices of the outliers.
     */
    protected void outputOutlierMap(int[] outlierIndices) {

        int nextOutlier = 0;
        try {
            for (int i = 0; i < gNoZeros.length; i++) {
                if (nextOutlier < outlierIndices.length
                        && i == outlierIndices[nextOutlier]) {
                    outlierMap.writeByte(1);
                    nextOutlier++;
                }
                else {
                    outlierMap.writeByte(0);
                }
            }
        }
        catch (FileNotFoundException e) {
            throw new LoggedException(e);
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }
    }

    /**
     * Returns a type of DiffTensorFitter specified by the index and initialized
     * with the other parameters.
     * 
     * @param q
     *            Independent variables.
     * 
     * @param data
     *            Normalized measurements.
     * 
     * @param bValues
     *            The array of b-values.
     * 
     * @param M
     *            The number of zero measurements in the acquisition.
     * 
     * @param index
     *            The index specifying the type of fitter.
     * 
     * @return The fitter.
     */
    protected DiffTensorFitter getFitter(double[][] g, double[] data, double[] bValues, int M,
            int index) {

        DiffTensorFitter dtf = null;
        try {
            if (index == 4) {
                dtf = new DiffTensorUnConFitter(g, data, bValues, M);
            }
            else {
                dtf = new DiffTensorFitter(g, data, bValues, M);
            }
        }
        catch (MarquardtMinimiserException e) {
            LoggedException.logExceptionSevere(e, Thread.currentThread().getName());
            System.exit(1);
        }

        return dtf;
    }
}

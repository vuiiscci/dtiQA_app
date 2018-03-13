package mesd;

import inverters.*;
import imaging.*;
import sphfunc.*;
import tools.CL_Initializer;

import java.util.Date;
import java.util.logging.Logger;

import optimizers.MarquardtMinimiserException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Provides methods for fitting max. ent. PAS to data.
 * 
 * <dt>Description:
 * 
 * <dd>This class provides an interface to MaxEntPDispFitter to allow the PAS
 * to be calculated from samples of the FT of a function restricted to a sphere.
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id: MESD_Inversion.java,v 1.6 2005/08/18 11:08:50
 *         ucacmgh Exp $
 *  
 */
public class MESD_Inversion extends DiffusionInversion {

    Logger logger = Logger.getLogger("mesd.MESD_Inversion");
    
    /**
     * Specifies whether to do the numerical integration tests and
     * increase point set sizes automatically or whether to just stick
     * with the result from the base point set.
     */
    public static boolean doNumIntTests = true;

    /**
     * Specifies the index of the base point set, ie the one used for
     * the numerical integration during the first minimization.
     */
    public static int basePointSet = 0;

    /**
     * The threshold on the error of the integral above which the minimisation
     * needs to be performed again.
     */
    public static double intErrThresh = 0.01;

    /**
     * A second threshold is used to semi-acceptance of a result once
     * integration has been performed using the largest point set.
     */
    public static double intErrThreshUpper = 0.05;

    /**
     * The threshold on the predicted error of the numerical integrals obtained
     * using the test point set.
     */
    public static double predErrThresh = 0.05;

    /**
     * The threshold on the RMS relative fitting error.
     */
    public static double relResidThresh = 0.5;

    /**
     * Convergence threshold in Lev. Mar. algorithm.
     */
    public static double convt = 1.0E-5;

    /**
     * Does the MESD fitting.
     */
    public MESD_Fitter mesdf;

    /**
     * This constructor sets up the default fitter.
     * 
     * @param imPars
     *            The details of the imaging sequence.
     * 
     * @param params
     *            The specification of the deconvolution kernel.
     */
    public MESD_Inversion(DW_Scheme imPars, double[] params, int numParams) {

        ip = imPars;
        int M = ip.numZeroMeasurements();
        
        
        // Get the list of gradient directions. Replaces getMeanNonZeroQ, 
	// which is the same thing unless there are multiple shells
        double[][] gNoZeros = ip.getNonZeroG_Dirs();

        try {
            mesdf = new MESD_Fitter(gNoZeros, M, params, numParams);
        }
        catch (MarquardtMinimiserException e) {
            throw new RuntimeException(e);
        }
                
    }

    /**
     * Computes the MESD.
     */
    public double[] invert(double[] data) {

        // The MESD fitter needs normalized data.
        double[] normData = ip.normalizeData(data);

        try {

            // Set up the data and parameters for the minimization.
            mesdf.newDepVals(normData);
            mesdf.setPointSet(basePointSet);
            mesdf.setConvergence(convt);

            long start = new Date().getTime();
            long orig = start;

            // Do the minimization.
            mesdf.minimise();

            // If nothing went wrong...
            int exitcode = 0;

            // Perform error checks to select the appropriate point
            // set for numerical integration.
            if(doNumIntTests) {

                //System.err.println("\n******** initial minimisation took "+(new
                // Date().getTime()-start)+" milliseconds ******* ");

                //Test the integration to see if it has become
                //unstable. If so then use a larger point set
                //and reminimise.
                int nextIcos = basePointSet + 1;
                start = new Date().getTime();
                double[] intStats = mesdf.integrationErrorStats();
                //System.err.println("*********** inital error stats error stats
                // calculation took "+(new Date().getTime()-start));
                double mRelErr = intStats[0];
                double mExpErr = intStats[1];
                while ((nextIcos < ISCodes.getNoPointSetsForMaxEnt())
                       && (mRelErr > intErrThresh)) {
                    
                    mesdf.setPointSet(nextIcos);

                    //First try reminimising from where
                    //the last minimisation converged.
                    start = new Date().getTime();
                    mesdf.minimise();
                    //System.err.println("*********** second minimisation took
                    // "+(new Date().getTime()-start));
                    start = new Date().getTime();
                    intStats = mesdf.integrationErrorStats();
                    //System.err.println("*********** second integration error
                    // stats took "+(new Date().getTime()-start));
                    mRelErr = intStats[0];
                    mExpErr = intStats[1];
                    if (mRelErr < intErrThresh) {
                        nextIcos++;
                        break;
                    }
                    
                    //If that failed, go back to the starting point and
                    //minimise with this larger point set.
                    //System.err.println("reinitialising.");
                    mesdf.reInit();
                    start = new Date().getTime();
                    mesdf.minimise();
                    //System.err.println("*********** third minimisation took
                    // "+(new Date().getTime()-start));
                    
                    //Do the same tests again.
                    start = new Date().getTime();
                    intStats = mesdf.integrationErrorStats();
                    //System.err.println("*********** third error stats calcualtion
                    // took "+(new Date().getTime()-start));
                    mRelErr = intStats[0];
                    mExpErr = intStats[1];
                    nextIcos++;
                }

                // Construct the exit code based on various
                // diagnostics.

                // Do a test on the maximum relative error
                // between the test point set and the point set
                // used to find the \lambda_j.
                int mreClass = getMREClass(mRelErr);

                // Do a test on the predicted error of the
                // numerical integrals evaluated using the test
                // point set.
                int predErrClass = getPredErrClass(mExpErr);

                // Do a test on the RMS relative fitting error.
                int relResidClass = getPredErrClass(Math.sqrt(mesdf.getRelativeResidual()));

                exitcode = mreClass + 10 * predErrClass + 100 * relResidClass;
            
            }

            // Construct the output array.
            double[] lambdas = mesdf.getParameters();

            double[] output = new double[lambdas.length + 1];
            output[0] = exitcode;
            output[1] = Math.log(ip.geoMeanZeroMeas(data));
            for (int i = 1; i < lambdas.length; i++) {
                output[i + 1] = lambdas[i];
            }

            // spew out some evaluations from the lambdas

            return output;

        }
        catch (Exception e) {
            logger.warning("Fitting failed in MESD_Inversion.");
            logger.warning(e.toString());

            double[] output = new double[MaxEntProfile.getNumParams() + 2];

            // Indicate failure of fitting with exitcode of 5.
            output[0] = 5;

            return output;
        }

    }

    /**
     * Performs a threshold test on the maximum relative error of the numerical
     * integrals obtained using two point sets. A double threshold is applied.
     * The result is 2 if the MRE is greater than both thresholds, 1 if between
     * the two and zero if less than both.
     * 
     * @param mre
     *            The maximum relative error.
     * 
     * @return Error code 0, 1 or 2.
     */
    public int getMREClass(double mre) {
        if (mre > intErrThreshUpper) {
            return 2;
        }
        else if (mre > intErrThresh) {
            return 1;
        }
        return 0;
    }

    /**
     * Performs a threshold test on the predicted error of the numerical
     * integrals obtained using the test point set. The result is 1 if the error
     * is greater than the threshold, otherwise zero.
     * 
     * @param mre
     *            The predicted error.
     * 
     * @return Error code 0 or 1.
     */
    public int getPredErrClass(double mre) {
        if (mre > predErrThresh) {
            return 1;
        }
        return 0;
    }

    /**
     * Performs a threshold test on the RMS fitting error. The result 1 if the
     * error is greater than the threshold, otherwise zero.
     * 
     * @param relResid
     *            The RMS fitting error.
     * 
     * @return Error code 0 or 1.
     */
    public int getResidErrClass(double relResid) {
        if (relResid > relResidThresh) {
            return 1;
        }
        return 0;
    }

    public int itemsPerVoxel() {
        return MaxEntProfile.getNumParams() + 2;
    }

}

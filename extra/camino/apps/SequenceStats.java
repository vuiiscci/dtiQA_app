package apps;

import java.io.*;
import java.util.Random;
import java.util.logging.Logger;

import tools.CL_Initializer;
import misc.LoggedException;
import numerics.*;
import data.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Computes statistics of a number of similar results from different random
 * rotations.
 * 
 * <dt>Description:
 * 
 * <dd>For single tensor inversions, this outputs: {success rate (mean, max,
 * min), deflection angle (mean, max, min), dyad max eigenvalue (mean, max,
 * min), trace mean (mean, max, min), trace std (mean, max, min), fa mean (mean,
 * max, min), fa std (mean, max, min)}.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class SequenceStats {

    private static Logger logger = Logger.getLogger("camino.apps.SequenceStats");

    /**
     * The number of statistics in single DT trials.
     */
    public static final int NUMSINGLEDTSTATS = 12;

    /**
     * The number of statistics in two tensor trials.
     */
    public static final int NUMTWOTENSORSTATS = 24;

    /**
     * The number of statistics in three tensor trials.
     */
    public static final int NUMTHREETENSORSTATS = 36;

    /**
     * The number of statistics in MFR trials.
     */
    public static final int NUMMFRSTATS = 54;

    public static void main(String[] args) {

        // The number of rotations to process.
        int numRotations = 0;

        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-rotations")) {
                numRotations = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
        }

        // Set up the required objects.

        // Set up the input and output streams.
        DataInputStream in = null;
        try {
            in = new DataInputStream(System.in);
        }
        catch (Exception e1) {

            LoggedException.logExceptionWarning(e1, Thread.currentThread().getName());
        }

        OutputManager om = new OutputManager();

        // Allocate the array to store the input.
	// default to MFR
        int numStats = NUMMFRSTATS;

	// test for tensor inversions
        if (CL_Initializer.inversionIndices[0].numDTs == 1) {
            numStats = NUMSINGLEDTSTATS;
        }
        if (CL_Initializer.inversionIndices[0].numDTs == 2) {
            numStats = NUMTWOTENSORSTATS;
        }
        if (CL_Initializer.inversionIndices[0].numDTs == 3) {
            numStats = NUMTHREETENSORSTATS;
        }

        double[][] rotStats = new double[numRotations][numStats];

        // Read in the trial statistics.
        for (int i = 0; i < numRotations; i++) {
            for (int j = 0; j < numStats; j++)
                try {
                    rotStats[i][j] = in.readDouble();
                }
                catch (Exception e) {
                    StackTraceElement[] stackTrace = e.getStackTrace();
                    String stString = new String();

                    // log the exceptions message
                    logger.severe(e.toString());

                    // log the stack trace
                    for (int jj = 0; jj < stackTrace.length; jj++) {
                        stString += stackTrace[jj] + "\n";
                    }
                    logger.severe(stString);

                    throw new RuntimeException(e);
                }
        }

        // Compute and output the statistics.
        double[] results = null;
        if (CL_Initializer.inversionIndices[0].numDTs == 1) {
            results = singleDT_Stats(rotStats);
        }
        else if (CL_Initializer.inversionIndices[0].numDTs == 2) {
            results = twoTensorStats(rotStats, CL_Initializer.dt2rotangle);
        }
        else if (CL_Initializer.inversionIndices[0].numDTs == 3) {
            results = threeTensorStats(rotStats);
        }
        else {
            results = mfrStats(rotStats, CL_Initializer.testFunction, CL_Initializer.dt2rotangle);
        }

        om.output(results);
        om.close();
    }

    /**
     * Computes the statistics for single DT data and outputs them.
     * 
     * @param rotStats
     *            The input statistics.
     * 
     * @return {success rate (mean, max, min), deflection angle (mean, max,
     *         min), dyad max eigenvalue (mean, max, min), trace mean (mean,
     *         max, min), trace std (mean, max, min), fa mean (mean, max, min),
     *         fa std (mean, max, min)}
     */
    public static double[] singleDT_Stats(double[][] rotStats) {

        // Compute the statistics.
        double[] failureStats = InversionStats.meanSTDMaxAndMin(rotStats, 1);
        double[] meanDirStats = dirStatsWithRotations(rotStats, 2, 1, 0.0);
        double[] dyadStats = InversionStats.meanSTDMaxAndMin(rotStats, 5);
        double[] traceMeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 8);
        double[] traceStdStats = InversionStats.meanSTDMaxAndMin(rotStats, 9);
        double[] faMeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 10);
        double[] faStdStats = InversionStats.meanSTDMaxAndMin(rotStats, 11);

        double[] results = new double[21];
        results[0] = failureStats[0];
        results[1] = failureStats[2];
        results[2] = failureStats[3];
        results[3] = meanDirStats[0];
        results[4] = meanDirStats[2];
        results[5] = meanDirStats[3];
        results[6] = dyadStats[0];
        results[7] = dyadStats[2];
        results[8] = dyadStats[3];
        results[9] = traceMeanStats[0];
        results[10] = traceMeanStats[2];
        results[11] = traceMeanStats[3];
        results[12] = traceStdStats[0];
        results[13] = traceStdStats[2];
        results[14] = traceStdStats[3];
        results[15] = faMeanStats[0];
        results[16] = faMeanStats[2];
        results[17] = faMeanStats[3];
        results[18] = faStdStats[0];
        results[19] = faStdStats[2];
        results[20] = faStdStats[3];

        return results;
    }

    /**
     * Computes the statistics for two tensors and outputs them.
     * 
     * @param rotStats
     *            The input statistics.
     * 
     * @param dt2rotangle
     *            The rotation angle of the second component
     * 
     * @return {success rate (mean, max, min), dyad1 max eigenvalue (mean, max,
     *         min), trace1 mean (mean, max, min), trace1 std (mean, max, min),
     *         fa1 mean (mean, max, min), fa1 std (mean, max, min), dyad2 max
     *         eigenvalue (mean, max, min), trace2 mean (mean, max, min), trace2
     *         std (mean, max, min), fa2 mean (mean, max, min), fa2 std (mean,
     *         max, min), prop mean (mean, max, min), prop std (mean, max, min)}
     */
    public static double[] twoTensorStats(double[][] rotStats, double dt2rotangle) {
        // Compute the statistics.
        double[] failureStats = InversionStats.meanSTDMaxAndMin(rotStats, 1);

        double[] dyad1Stats = InversionStats.meanSTDMaxAndMin(rotStats, 5);
        double[] meanDir1Stats = dirStatsWithRotations(rotStats, 2, 3, dt2rotangle);
        double[] trace1MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 8);
        double[] trace1StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 9);
        double[] fa1MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 10);
        double[] fa1StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 11);

        double[] dyad2Stats = InversionStats.meanSTDMaxAndMin(rotStats, 15);
        double[] meanDir2Stats = dirStatsWithRotations(rotStats, 12, 3, dt2rotangle);
        double[] trace2MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 18);
        double[] trace2StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 19);
        double[] fa2MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 20);
        double[] fa2StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 21);

        double[] propsMeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 22);
        double[] propsStdStats = InversionStats.meanSTDMaxAndMin(rotStats, 23);

        double[] results = new double[45];
        results[0] = failureStats[0];
        results[1] = failureStats[2];
        results[2] = failureStats[3];
        results[3] = meanDir1Stats[0];
        results[4] = meanDir1Stats[2];
        results[5] = meanDir1Stats[3];
        results[6] = dyad1Stats[0];
        results[7] = dyad1Stats[2];
        results[8] = dyad1Stats[3];
        results[9] = trace1MeanStats[0];
        results[10] = trace1MeanStats[2];
        results[11] = trace1MeanStats[3];
        results[12] = trace1StdStats[0];
        results[13] = trace1StdStats[2];
        results[14] = trace1StdStats[3];
        results[15] = fa1MeanStats[0];
        results[16] = fa1MeanStats[2];
        results[17] = fa1MeanStats[3];
        results[18] = fa1StdStats[0];
        results[19] = fa1StdStats[2];
        results[20] = fa1StdStats[3];
        results[21] = meanDir2Stats[0];
        results[22] = meanDir2Stats[2];
        results[23] = meanDir2Stats[3];
        results[24] = dyad2Stats[0];
        results[25] = dyad2Stats[2];
        results[26] = dyad2Stats[3];
        results[27] = trace2MeanStats[0];
        results[28] = trace2MeanStats[2];
        results[29] = trace2MeanStats[3];
        results[30] = trace2StdStats[0];
        results[31] = trace2StdStats[2];
        results[32] = trace2StdStats[3];
        results[33] = fa2MeanStats[0];
        results[34] = fa2MeanStats[2];
        results[35] = fa2MeanStats[3];
        results[36] = fa2StdStats[0];
        results[37] = fa2StdStats[2];
        results[38] = fa2StdStats[3];
        results[39] = propsMeanStats[0];
        results[40] = propsMeanStats[2];
        results[41] = propsMeanStats[3];
        results[42] = propsStdStats[0];
        results[43] = propsStdStats[2];
        results[44] = propsStdStats[3];

        return results;
    }

    /**
     * Computes the statistics for three tensors and outputs them.
     * 
     * @param rotStats
     *            The input statistics.
     * 
     * @return {success rate (mean, max, min), dyad1 max eigenvalue (mean, max,
     *         min), trace1 mean (mean, max, min), trace1 std (mean, max, min),
     *         fa1 mean (mean, max, min), fa1 std (mean, max, min), dyad2 max
     *         eigenvalue (mean, max, min), trace2 mean (mean, max, min), trace2
     *         std (mean, max, min), fa2 mean (mean, max, min), fa2 std (mean,
     *         max, min), prop mean (mean, max, min), dyad3 max eigenvalue
     *         (mean, max, min), trace3 mean (mean, max, min), trace3 std (mean,
     *         max, min), fa3 mean (mean, max, min), fa3 std (mean, max, min),
     *         prop1 mean (mean, max, min), prop1 std (mean, max, min), prop2
     *         mean (mean, max, min), prop2 std (mean, max, min)}
     */
    public static double[] threeTensorStats(double[][] rotStats) {
        // Compute the statistics.
        double[] failureStats = InversionStats.meanSTDMaxAndMin(rotStats, 1);

        double[] dyad1Stats = InversionStats.meanSTDMaxAndMin(rotStats, 5);
        double[] meanDir1Stats = dirStatsWithRotations(rotStats, 2, 3, 0.0);
        double[] trace1MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 8);
        double[] trace1StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 9);
        double[] fa1MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 10);
        double[] fa1StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 11);

        double[] dyad2Stats = InversionStats.meanSTDMaxAndMin(rotStats, 15);
        double[] meanDir2Stats = dirStatsWithRotations(rotStats, 12, 3, 0.0);
        double[] trace2MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 18);
        double[] trace2StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 19);
        double[] fa2MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 20);
        double[] fa2StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 21);

        double[] dyad3Stats = InversionStats.meanSTDMaxAndMin(rotStats, 25);
        double[] meanDir3Stats = dirStatsWithRotations(rotStats, 22, 3, 0.0);
        double[] trace3MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 28);
        double[] trace3StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 29);
        double[] fa3MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 30);
        double[] fa3StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 31);

        double[] prop1MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 32);
        double[] prop1StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 33);

        double[] prop2MeanStats = InversionStats.meanSTDMaxAndMin(rotStats, 34);
        double[] prop2StdStats = InversionStats.meanSTDMaxAndMin(rotStats, 35);

        double[] results = new double[69];
        results[0] = failureStats[0];
        results[1] = failureStats[2];
        results[2] = failureStats[3];
        results[3] = meanDir1Stats[0];
        results[4] = meanDir1Stats[2];
        results[5] = meanDir1Stats[3];
        results[6] = dyad1Stats[0];
        results[7] = dyad1Stats[2];
        results[8] = dyad1Stats[3];
        results[9] = trace1MeanStats[0];
        results[10] = trace1MeanStats[2];
        results[11] = trace1MeanStats[3];
        results[12] = trace1StdStats[0];
        results[13] = trace1StdStats[2];
        results[14] = trace1StdStats[3];
        results[15] = fa1MeanStats[0];
        results[16] = fa1MeanStats[2];
        results[17] = fa1MeanStats[3];
        results[18] = fa1StdStats[0];
        results[19] = fa1StdStats[2];
        results[20] = fa1StdStats[3];
        results[21] = meanDir2Stats[0];
        results[22] = meanDir2Stats[2];
        results[23] = meanDir2Stats[3];
        results[24] = dyad2Stats[0];
        results[25] = dyad2Stats[2];
        results[26] = dyad2Stats[3];
        results[27] = trace2MeanStats[0];
        results[28] = trace2MeanStats[2];
        results[29] = trace2MeanStats[3];
        results[30] = trace2StdStats[0];
        results[31] = trace2StdStats[2];
        results[32] = trace2StdStats[3];
        results[33] = fa2MeanStats[0];
        results[34] = fa2MeanStats[2];
        results[35] = fa2MeanStats[3];
        results[36] = fa2StdStats[0];
        results[37] = fa2StdStats[2];
        results[38] = fa2StdStats[3];
        results[39] = meanDir3Stats[0];
        results[40] = meanDir3Stats[2];
        results[41] = meanDir3Stats[3];
        results[42] = dyad3Stats[0];
        results[43] = dyad3Stats[2];
        results[44] = dyad3Stats[3];
        results[45] = trace3MeanStats[0];
        results[46] = trace3MeanStats[2];
        results[47] = trace3MeanStats[3];
        results[48] = trace3StdStats[0];
        results[49] = trace3StdStats[2];
        results[50] = trace3StdStats[3];
        results[51] = fa3MeanStats[0];
        results[52] = fa3MeanStats[2];
        results[53] = fa3MeanStats[3];
        results[54] = fa3StdStats[0];
        results[55] = fa3StdStats[2];
        results[56] = fa3StdStats[3];
        results[57] = prop1MeanStats[0];
        results[58] = prop1MeanStats[2];
        results[59] = prop1MeanStats[3];
        results[60] = prop1StdStats[0];
        results[61] = prop1StdStats[2];
        results[62] = prop1StdStats[3];
        results[63] = prop2MeanStats[0];
        results[64] = prop2MeanStats[2];
        results[65] = prop2MeanStats[3];
        results[66] = prop2StdStats[0];
        results[67] = prop2StdStats[2];
        results[68] = prop2StdStats[3];

        return results;
    }

    /**
     * Computes the statistics for data from a multiple-fibre reconstruction and
     * outputs them.
     * 
     * @param rotStats
     *            The input statistics.
     * 
     * @param testFunc
     *            The test function used in the simulation
     * 
     * @param dt2rotangle
     *            The rotation of the second component in test function 3
     * 
     * @return {success rate (mean, max, min), dyad max eigenvalue (mean, max,
     *         min), trace mean (mean, max, min), trace std (mean, max, min), fa
     *         mean (mean, max, min), fa std (mean, max, min)}
     */
    public static double[] mfrStats(double[][] rotStats, int testFunc, double dt2rotangle) {

        // Compute the statistics.
        double[] failureStats = InversionStats.meanSTDMaxAndMin(rotStats, 1);
        double[] meanmeanF_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 2);
        double[] stdmeanF_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 3);
        double[] meanstdF_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 4);
        double[] stdstdF_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 5);

        double[] mean1_DirStats = dirStatsWithRotations(rotStats, 6, testFunc,
                dt2rotangle);
        double[] dyad11_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 9);
        double[] meanF1_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 12);
        double[] stdF1_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 13);
        double[] meanH11_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 14);
        double[] meanH12_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 16);

        double[] mean2_DirStats = dirStatsWithRotations(rotStats, 18, testFunc,
                dt2rotangle);
        double[] dyad21_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 21);
        double[] meanF2_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 24);
        double[] stdF2_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 25);
        double[] meanH21_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 26);
        double[] meanH22_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 28);

        double[] mean3_DirStats = dirStatsWithRotations(rotStats, 30, testFunc,
                dt2rotangle);
        double[] dyad31_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 33);
        double[] meanF3_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 36);
        double[] stdF3_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 37);
        double[] meanH31_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 38);
        double[] meanH32_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 40);

        double[] mean4_DirStats = dirStatsWithRotations(rotStats, 42, testFunc,
                dt2rotangle);
        double[] dyad41_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 45);
        double[] meanF4_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 48);
        double[] stdF4_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 49);
        double[] meanH41_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 50);
        double[] meanH42_Stats = InversionStats.meanSTDMaxAndMin(rotStats, 52);

        double[] results = new double[87];
        results[0] = failureStats[0];
        results[1] = failureStats[2];
        results[2] = failureStats[3];
        results[3] = meanmeanF_Stats[0];
        results[4] = meanmeanF_Stats[2];
        results[5] = meanmeanF_Stats[3];
        results[6] = stdmeanF_Stats[0];
        results[7] = stdmeanF_Stats[2];
        results[8] = stdmeanF_Stats[3];
        results[9] = meanstdF_Stats[0];
        results[10] = meanstdF_Stats[2];
        results[11] = meanstdF_Stats[3];
        results[12] = stdstdF_Stats[0];
        results[13] = stdstdF_Stats[2];
        results[14] = stdstdF_Stats[3];

        results[15] = mean1_DirStats[0];
        results[16] = mean1_DirStats[2];
        results[17] = mean1_DirStats[3];
        results[18] = dyad11_Stats[0];
        results[19] = dyad11_Stats[2];
        results[20] = dyad11_Stats[3];
        results[21] = meanF1_Stats[0];
        results[22] = meanF1_Stats[2];
        results[23] = meanF1_Stats[3];
        results[24] = stdF1_Stats[0];
        results[25] = stdF1_Stats[2];
        results[26] = stdF1_Stats[3];
        results[27] = meanH11_Stats[0];
        results[28] = meanH11_Stats[2];
        results[29] = meanH11_Stats[3];
        results[30] = meanH12_Stats[0];
        results[31] = meanH12_Stats[2];
        results[32] = meanH12_Stats[3];

        results[33] = mean2_DirStats[0];
        results[34] = mean2_DirStats[2];
        results[35] = mean2_DirStats[3];
        results[36] = dyad21_Stats[0];
        results[37] = dyad21_Stats[2];
        results[38] = dyad21_Stats[3];
        results[39] = meanF2_Stats[0];
        results[40] = meanF2_Stats[2];
        results[41] = meanF2_Stats[3];
        results[42] = stdF2_Stats[0];
        results[43] = stdF2_Stats[2];
        results[44] = stdF2_Stats[3];
        results[45] = meanH21_Stats[0];
        results[46] = meanH21_Stats[2];
        results[47] = meanH21_Stats[3];
        results[48] = meanH22_Stats[0];
        results[49] = meanH22_Stats[2];
        results[50] = meanH22_Stats[3];

        results[51] = mean3_DirStats[0];
        results[52] = mean3_DirStats[2];
        results[53] = mean3_DirStats[3];
        results[54] = dyad31_Stats[0];
        results[55] = dyad31_Stats[2];
        results[56] = dyad31_Stats[3];
        results[57] = meanF3_Stats[0];
        results[58] = meanF3_Stats[2];
        results[59] = meanF3_Stats[3];
        results[60] = stdF3_Stats[0];
        results[61] = stdF3_Stats[2];
        results[62] = stdF3_Stats[3];
        results[63] = meanH31_Stats[0];
        results[64] = meanH31_Stats[2];
        results[65] = meanH31_Stats[3];
        results[66] = meanH32_Stats[0];
        results[67] = meanH32_Stats[2];
        results[68] = meanH32_Stats[3];

        results[69] = mean4_DirStats[0];
        results[70] = mean4_DirStats[2];
        results[71] = mean4_DirStats[3];
        results[72] = dyad41_Stats[0];
        results[73] = dyad41_Stats[2];
        results[74] = dyad41_Stats[3];
        results[75] = meanF4_Stats[0];
        results[76] = meanF4_Stats[2];
        results[77] = meanF4_Stats[3];
        results[78] = stdF4_Stats[0];
        results[79] = stdF4_Stats[2];
        results[80] = stdF4_Stats[3];
        results[81] = meanH41_Stats[0];
        results[82] = meanH41_Stats[2];
        results[83] = meanH41_Stats[3];
        results[84] = meanH42_Stats[0];
        results[85] = meanH42_Stats[2];
        results[86] = meanH42_Stats[3];

        return results;
    }

    /**
     * Computes the mean, standard deviation, maximum and minimum of the angle
     * between directions and the closest matching principal direction of the
     * test function.
     * 
     * @param props
     *            The array.
     * 
     * @param col
     *            The column.
     * 
     * @param functionIndex
     *            The test function used in the simulation
     * 
     * @param dt2rotangle
     *            The rotation of the second component in test function 3
     * 
     * @return The mean, standard deviation, max and min: {mean, std, max, min}.
     */
    public static double[] dirStatsWithRotations(double[][] props, int col,
            int functionIndex, double dt2rotangle) {

        double sum = 0.0;
        double sumSq = 0.0;
        double max = -12.0;
        double min = 12.0;

        for (int i = 0; i < props.length; i++) {

            // Construct the test function.
            StandardTestFunctions.setTransformation(Rotations.randomRotMat(new Random(
                    i + 1)));
            StandardTestFunctions.setDT2RotationAngle(dt2rotangle);
            ModelPDF p = StandardTestFunctions.getFunction(functionIndex);
            double[][] pds = p.getPDs();

            // Find the maximum absolute dot product.
            double maxDP = -1.0;
            for (int j = 0; j < pds.length; j++) {
                double adp = Math.abs(props[i][col] * pds[j][0] + props[i][col + 1]
                        * pds[j][1] + props[i][col + 2] * pds[j][2]);
                maxDP = (adp > maxDP) ? adp : maxDP;
            }

            // Take the inverse cosine to turn it into an angle.
            double sepAng = Math.acos(maxDP);

            sum += sepAng;
            sumSq += sepAng * sepAng;
            max = (sepAng > max) ? sepAng : max;
            min = (sepAng < min) ? sepAng : min;
        }

        double[] results = new double[4];
        results[0] = sum / (double) props.length;
        results[1] = Math.sqrt(sumSq / (double) props.length - results[0] * results[0]);
        results[2] = max;
        results[3] = min;

        return results;
    }

}


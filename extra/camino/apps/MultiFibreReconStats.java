package apps;

import java.util.*;

import tools.*;
import numerics.*;
import inverters.ModelIndex;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Computes statistics of properties of the multiple-fibre
 * reconstruction outputs.
 * 
 * <dt>Description:
 * 
 * <dd>Reads in the output of SphFuncPD_Stats and computes statistics
 * of the distribution of principal directions and shape properties in
 * the results.  Used mostly for simulations.
 * 
 * The output is: - number of trials over which statistics computed, -
 * fraction of successful trials - mean mean(SF) - std mean(SF) - mean
 * std(SF) - std std(SF) - for each direction - mean direction (x, y,
 * z), - mean dyadic eigenvalues (l1, l2, l3), - mean SF value, - std
 * SF value, - mean SF largest Hessian eigenvalue - std SF largest
 * Hessian eigenvalue - mean SF smallest Hessian eigenvalue - std SF
 * smallest Hessian eigenvalue
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class MultiFibreReconStats extends InversionStats {

    public static void main(String[] args) {
        MultiFibreReconStats is = new MultiFibreReconStats(args);
    }


    /**
     * Runs the main program given the command line arguments.
     * 
     * @param args
     *            The command line arguments.
     */
    protected MultiFibreReconStats(String[] args) {
        super(args);
    }


    /**
     * Tests whether the results of an inversion should be counted as
     * successful. This version checks that the principal direction
     * extraction procedure was successful and that the number of
     * directions in the output is at least the number of expected
     * directions.
     * 
     * @param fittedData
     *            The output data for the inversion.
     * 
     * @param numExpectedDirections
     *            The expected number of directions.
     * 
     * @return A boolean indicating whether the inversion was successful.
     */
    protected boolean successfulInversion(double[] fittedData, int numExpectedDirections) {
        return ((int) fittedData[3] == 1 && (fittedData[2] >= numExpectedDirections));
    }


    /**
     * Returns the number of parameters per voxel given the inversion
     * index.
     * 
     * @param index
     *            The inversion index.
     * 
     * @return The number of parameters per voxel the inversion provides.
     */
    protected int getNumParams(ModelIndex index) {
        return SphFuncPD_Stats.GLOBALSTATS + SphFuncPD_Stats.STATSPERPD * CL_Initializer.numPDsIO;
    }


    /**
     * Returns the number of properties per direction given the
     * inversion index.
     * 
     * @param index
     *            The inversion index.
     * 
     * @return The number of properties the inversion provides.
     */
    protected int getNumProperties(ModelIndex index) {
        return getNumParams(index);
    }


    /**
     * Computes the properties of an inversion from the data
     * values. Usually for multiple-fibre reconstructions, no extra
     * calculation is required and the input array is the properties.
     * 
     * @param data
     *            The data values.
     * 
     * @param inversion
     *            The index of the inversion that generated the data.
     * 
     * @return The array of properties.
     */
    protected double[] getProperties(double[] data, ModelIndex inversion) {
        return data;
    }


    /**
     * Computes some statistics from a list of inversion properties.
     * 
     * @param props
     *            The array of properties.
     * 
     * @param inversion
     *            The inversion index.
     * 
     * @return The list of statistics: {mean(mean(SF)), std(mean(SF)),
     *         mean(std(SF)), std(std(SF)), (i=1:3: meaniX, meaniY, meaniZ,
     *         mean(SFi), std(SFi), mean(Hi1), std(Hi1), mean(Hi2), std(Hi2))}
     */
    protected double[] computeStats(double[][] props, ModelIndex inversion) {

        // Sort the corresponding components.
        int firstDirIndex = SphFuncPD_Stats.GLOBALSTATS;
        int propsPerComponent = SphFuncPD_Stats.STATSPERPD;
        orderComponentsByDirection(props, CL_Initializer.numPDsIO, propsPerComponent,
                firstDirIndex);

        // Now reorder the peaks by their mean strength.
        orderPeaksByMeanStrength(props, CL_Initializer.numPDsIO, propsPerComponent,
                firstDirIndex);

        // Compute the stats of the other properties.
        double[] meanSF_Stats = meanSTDMaxAndMin(props, 4);
        double[] stdSF_Stats = meanSTDMaxAndMin(props, 5);

        // Gather all the results together.
        // Four for the mean and std (mean and std)
        // Each direction has DIRSTATS and 3 associated values each
        // with a mean and std for a total of DIRSTATS+6.
        double[] stats = new double[4 + CL_Initializer.numPDsIO * (DIRSTATS + 6)];

        stats[0] = meanSF_Stats[0];
        stats[1] = meanSF_Stats[1];
        stats[2] = stdSF_Stats[0];
        stats[3] = stdSF_Stats[1];

        for (int j = 0; j < CL_Initializer.numPDsIO; j++) {
            int dirInd = j * propsPerComponent + firstDirIndex;
            double[] meanDirStats = oneDirectionStats(props, dirInd);
            for (int i = 0; i < meanDirStats.length; i++) {
                stats[4 + j * (DIRSTATS + 6) + i] = meanDirStats[i];
            }
            double[] peakStats = meanSTDMaxAndMin(props, firstDirIndex + j * propsPerComponent + 3);

            double[][] hessEvals = getHessianEigenVals(props, firstDirIndex + j * propsPerComponent + 4);
            double[] h1_Stats = meanSTDMaxAndMin(hessEvals, 0);
            double[] h2_Stats = meanSTDMaxAndMin(hessEvals, 1);

            stats[4 + j * (DIRSTATS + 6) + DIRSTATS] = peakStats[0];
            stats[4 + j * (DIRSTATS + 6) + DIRSTATS + 1] = peakStats[1];
            stats[4 + j * (DIRSTATS + 6) + DIRSTATS + 2] = h1_Stats[0];
            stats[4 + j * (DIRSTATS + 6) + DIRSTATS + 3] = h1_Stats[1];
            stats[4 + j * (DIRSTATS + 6) + DIRSTATS + 4] = h2_Stats[0];
            stats[4 + j * (DIRSTATS + 6) + DIRSTATS + 5] = h2_Stats[1];
        }

        return stats;
    }


    /**
     * Reorders the array of inversion properties when it contains
     * multiple directions by the mean strength of the peak directions
     * across all trials.
     * 
     * @param props
     *            The array of properties.
     * 
     * @param numComponents
     *            The number of directions in the props array
     * 
     * @param propsPerComponent
     *            The number of values in the props array for each direction
     * 
     * @param firstDirIndex
     *            The index of props that is the start of the first direction
     */
    public void orderPeaksByMeanStrength(double[][] props, int numComponents,
            int propsPerComponent, int firstDirIndex) {

        // If no trials were successful, nothing to do.
        if (props.length == 0) return;

        // Get the array of mean peak strengths
        double[] meanPeakStrengths = new double[numComponents];
        double[] sortedMeanPeakStrengths = new double[numComponents];
        for (int j = 0; j < numComponents; j++) {
            int peakStrengthInd = j * propsPerComponent + firstDirIndex + 3;
            for (int i = 0; i < props.length; i++) {
                meanPeakStrengths[j] += props[i][peakStrengthInd];
            }
            meanPeakStrengths[j] /= (double) props.length;
            sortedMeanPeakStrengths[j] = meanPeakStrengths[j];
        }

        Arrays.sort(sortedMeanPeakStrengths);

        // Copy the props array
        double[][] propsCopy = new double[props.length][props[0].length];
        for (int i = 0; i < props.length; i++) {
            for (int j = 0; j < props[0].length; j++) {
                propsCopy[i][j] = props[i][j];
            }
        }

        // Refill in the right order
        for (int i = 0; i < numComponents; i++) {

            // Find the component with the i-th highest mean peak strength
            int thisComp = 0;
            for (int j = 0; j < numComponents; j++) {
                if (meanPeakStrengths[j] == sortedMeanPeakStrengths[numComponents - i - 1]) {
                    thisComp = j;
                }
            }

            // Replace the contents of the props array
            if (thisComp != i) {
                for (int j = 0; j < props.length; j++) {
                    for (int k = 0; k < propsPerComponent; k++) {
                        props[j][i * propsPerComponent + firstDirIndex + k] = propsCopy[j][thisComp
                                * propsPerComponent + firstDirIndex + k];
                    }
                }
            }
        }
    }


    /**
     * Returns the eigenvalues of the Hessian matrix the four elements
     * of which are in each row of the the props array starting from
     * element firstInd.
     *
     * @param props The array of properties.
     *
     * @param firstInd The first index of the Hessian matrix.
     *
     * @return A list of Hessian eigenvalues, larger then smaller,
     * from the Hessian in each row of props.
     */
    public double[][] getHessianEigenVals(double[][] props, int firstInd) {

        double[][] evals = new double[props.length][2];
        for(int i=0; i<props.length; i++) {

            // Reconstruct the Hessian matrix.
            RealMatrix H = new RealMatrix(2, 2);
            H.entries[0][0] = props[i][firstInd];
            H.entries[0][1] = props[i][firstInd + 1];
            H.entries[1][0] = props[i][firstInd + 2];
            H.entries[1][1] = props[i][firstInd + 3];
            
            RealMatrix[] eig = H.jacobi();

            evals[i][0] = eig[0].entries[0][0];
            evals[i][1] = eig[0].entries[1][1];
        }

        return evals;
    }

}

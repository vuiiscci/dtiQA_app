package apps;

import java.io.*;
import java.util.logging.Logger;

import inverters.ModelIndex;
import numerics.*;
import data.*;
import misc.*;
import tools.*;

/**
 * <dl>
 *
 * <dt>Purpose:
 *
 * <dd>Computes statistics of properties of the inversion outputs.
 *
 * <dt>Description:
 *
 * <dd>Reads in the output of ModelFit and computes statistics of the
 * distribution of principal directions and shape properties in the
 * results.  Used mostly for simulations.
 *
 * For a single tensor inversion, the output is: {number of trials
 * over which statistics computed, fraction of successful trials, mean
 * direction (x, y, z), mean dyadic eigenvalues (l1, l2, l3), mean
 * trace, std trace, mean FA, std FA}.
 *
 * For two tensors: {number of trials over which statistics computed,
 * fraction of successful trials, mean direction 1 (x, y, z), mean
 * dyadic eigenvalues 1 (l1, l2, l3), mean trace 1, std trace 1, mean
 * FA 1, std FA 1, mean direction 2 (x, y, z), mean dyadic eigenvalues
 * 2 (l1, l2, l3), mean trace 2, std trace 2, mean FA 2, std FA 2,
 * mean prop, std prop}.
 *
 * For three tensors: {number of trials over which statistics
 * computed, fraction of successful trials, mean direction 1 (x, y,
 * z), mean dyadic eigenvalues 1 (l1, l2, l3), mean trace 1, std trace
 * 1, mean FA 1, std FA 1, mean direction 2 (x, y, z), mean dyadic
 * eigenvalues 2 (l1, l2, l3), mean trace 2, std trace 2, mean FA 2,
 * std FA 2, mean direction 3 (x, y, z), mean dyadic eigenvalues 3
 * (l1, l2, l3), mean trace 3, std trace 3, mean FA 3, std FA 3, mean
 * prop1, std prop1, mean prop2, std prop2}.
 *
 * </dl>
 *
 * @author Danny Alexander
 * @version $Id$
 *
 */
public class InversionStats extends Executable {
    
    public InversionStats(String[] args) {
        super(args);
    }
    
    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.InversionStats");
    
    /**
     * The number of properties computed for one DT.
     */
    protected final int NUMDTPROPERTIES = 5;
    
    /**
     * The number of properties computed for two tensors.
     */
    protected final int NUMTWOTENSORPROPERTIES = 12;
    
    /**
     * The number of properties computed for three tensors.
     */
    protected final int NUMTHREETENSORPROPERTIES = 18;
    
    /**
     * The number of statistics computed from each list of directions.
     */
    protected static final int DIRSTATS = 6;
    
    /**
     * The number of statistics computed from a set of single DTs.
     */
    protected final int NUMDTSTATS = DIRSTATS + 4;
    
    /**
     * The number of statistics computed from a set of DT pairs.
     */
    protected final int NUMTWOTENSORSTATS = 2 * NUMDTSTATS + 2;
    
    /**
     * The number of statistics computed from a set of DT triples.
     */
    protected final int NUMTHREETENSORSTATS = 3 * NUMDTSTATS + 4;
    
    
    /**
     * The threshold on the dot product between directions used by
     * ConsistencyFraction.
     */
    protected static double dpThresh = 0.95;
    
    // The expected number of direction in the reconstructions.
    // This is ignored in the basic one and two-tensor inversions.
    private static int numExpectedDirections;// = -1;
	private DataInputStream in;
    
    public void initDefaultVals() {
        //final int NUMDTPROPERTIES = 5;
        //final int NUMTWOTENSORPROPERTIES = 12;
        //final int NUMTHREETENSORPROPERTIES = 18;
        //final int DIRSTATS = 6;
        //final int NUMDTSTATS = DIRSTATS + 4;
        //final int NUMTWOTENSORSTATS = 2 * NUMDTSTATS + 2;
        //final int NUMTHREETENSORSTATS = 3 * NUMDTSTATS + 4;
        dpThresh = 0.95;
        int numExpectedDirections = -1;
        in = null;
    }
    
/*    public static void main(String[] args) {
        InversionStats is = new InversionStats(args);
    }
 */
    /**
     * Runs the main program given the command line arguments.
     *
     * @param args
     *            The command line arguments.
     */
//    protected InversionStats(String[] args) {
//    public static void main(String[] args) {
	public void initOptions(String[] args) {        
        CL_Initializer.inputDataType = "double";
        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-expect")) {
                numExpectedDirections = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
            if (args[i].equals("-threshold")) {
                dpThresh = Double.parseDouble(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
        }
        CL_Initializer.checkParsing(args);
    }
    
    public void initVariables() {
        
        // Set up the required objects.
        
        // Set up the input and output streams.
        int numParams = getNumParams(CL_Initializer.inversionIndices[0]);
        CL_Initializer.data = ExternalDataSource.getDataSource(CL_Initializer.inputFile, numParams, CL_Initializer.inputDataType);
    }
    
	public void execute(OutputManager om) {	
        
        //    OutputManager om = new OutputManager();
        
        
        // Array for containing the salient properties of each
        // inversion.

	// we use the first inversion index. This might be the only one specified on the 
	// command line, or the user may have specified a multi-tensor inversion and
	// a fallback DT inversion

        double[][] allProperties = 
	    getPropertiesList(CL_Initializer.numVoxels, CL_Initializer.inversionIndices[0],
			      numExpectedDirections);

        // Compute the property stats.
        double[] stats = computeStats(allProperties, CL_Initializer.inversionIndices[0]);
        
        // Output the statistics.
        double[] o = new double[stats.length + 2];
        o[0] = (double) allProperties.length;
        o[1] = (double) allProperties.length / (double) CL_Initializer.numVoxels;
        for (int i = 0; i < stats.length; i++) {
            o[i + 2] = stats[i];
        }
        om.output(o);
        
        
        // Finish up.
        om.close();
    }
    
    
    /**
     * Computes the list of inversion properties from data on the standard
     * input.
     *
     * @param numVoxels
     *            The number of trials in the input data.
     *
     * @param inversion
     *            The index of the inversion that gave rise to the data.
     *
     * @param numExpectedDirections
     *            The expected number of directions.
     *
     * @return The list of properties computed from each trial.
     */
    protected double[][] getPropertiesList(int numVoxels,
    ModelIndex inversion, int numExpectedDirections) {
        
        int numParams = getNumParams(CL_Initializer.inversionIndices[0]);
        int numProperties = getNumProperties(inversion);
        
        // Read in all the data.
        double[][] data = new double[numVoxels][numParams];
        int validSamples = 0;
        for (int i = 0; i < numVoxels; i++) {
            double[] fittedData = CL_Initializer.data.nextVoxel();
            for (int j = 0; j < numParams; j++) {
                data[i][j] = fittedData[j];
            }
            if (successfulInversion(fittedData, numExpectedDirections)) {
                validSamples += 1;
            }
        }
        
        double[][] allProperties = new double[validSamples][numProperties];
        
        // Loop over the data computing the properties of each
        // inversion. Note that we discard any for which the exit
        // code is non-zero.
        int sample = 0;
        for (int i = 0; i < numVoxels; i++) {
            if (successfulInversion(data[i], numExpectedDirections)) {
                
                double[] props = getProperties(data[i], inversion);
                for (int j = 0; j < props.length; j++) {
                    allProperties[sample][j] = props[j];
                }
                sample += 1;
            }
        }
        
        return allProperties;
    }
    
    
    /**
     * Tests whether the results of an inversion should be counted as
     * successful. The default implementation checks the exit code.
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
        return (fittedData[0] == 0.0);
    }
    
    
    /**
     * Returns the number of parameters per voxel given the inversion index.
     *
     * @param index
     *            The inversion index.
     *
     * @return The number of parameters per voxel the inversion provides.
     */
    protected int getNumParams(ModelIndex index) {
        if (index.numDTs == 1) {
            // One tensor
            return 8;
        }
        else if (index.numDTs == 2) {
            // Two tensors
            return 17;
        }
        else if (index.numDTs == 3) {
            // Three tensors
            return 24;
        }
        else {
            throw new LoggedException("Inversion " + index + " is not supported");
        }
    }
    
    
    /**
     * Returns the number of properties given the inversion index.
     *
     * @param index
     *            The inversion index.
     *
     * @return The number of properties the inversion provides.
     */
    protected int getNumProperties(ModelIndex index) {
        if (index.numDTs == 1) {
            return NUMDTPROPERTIES;
        }
        if (index.numDTs == 2) {
            return NUMTWOTENSORPROPERTIES;
        }
        if (index.numDTs == 3) {
            return NUMTHREETENSORPROPERTIES;
        }
        else {
            throw new LoggedException("Inversion " + index + " is not supported");
        }
    }
    
    
    /**
     * Computes the properties of an inversion from the data values.
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
        if (inversion.numDTs == 1) {
            return getDT_Properties(data);
        }
        if (inversion.numDTs == 2) {
            return getTwoTensorProperties(data);
        }
        else {
            return getThreeTensorProperties(data);
        }
    }
    
    
    /**
     * Computes properties of the DT, in particular, the principal direction,
     * trace and fractional anisotropy, from the array of values output for a
     * single DT by <code>ModelFit</code>.
     *
     * @param data
     *            The data values.
     *
     * @return The properties: {x, y, z, trace, fa}.
     */
    protected double[] getDT_Properties(double[] data) {
        double[] properties = new double[NUMDTPROPERTIES];
        
        DT d = new DT(data[2], data[3], data[4], data[5], data[6], data[7]);
        
        double[] pd = d.getPD();
        for (int i = 0; i < pd.length; i++) {
            properties[i] = pd[i];
        }
        
        properties[3] = d.trace();
        properties[4] = d.fa();
        
        return properties;
    }
    
    
    /**
     * Computes properties of a pair of DTs, in particular, the principal
     * directions, traces and fractional anisotropies, from the array of values
     * output for a pair of DTs by <code>ModelFit</code>.
     *
     * @param data
     *            The data values.
     *
     * @return The properties: {x1, y1, z1, prop1, trace1, fa1, x2, y2, z2,
     *         prop2, trace2, fa2}.
     */
    protected double[] getTwoTensorProperties(double[] data) {
        double[] properties = new double[NUMTWOTENSORPROPERTIES];
        
        DT d1 = new DT(data[4], data[5], data[6], data[7], data[8], data[9]);
        DT d2 = new DT(data[11], data[12], data[13], data[14], data[15], data[16]);
        
        double[] pd1 = d1.getPD();
        double[] pd2 = d2.getPD();
        for (int i = 0; i < pd1.length; i++) {
            properties[i] = pd1[i];
            properties[i + 6] = pd2[i];
        }
        
        properties[3] = data[3];
        properties[4] = d1.trace();
        properties[5] = d1.fa();
        
        properties[9] = data[10];
        properties[10] = d2.trace();
        properties[11] = d2.fa();
        
        return properties;
    }
    
    
    /**
     * Computes properties of three mixed DTs, in particular, the principal
     * directions, traces and fractional anisotropies, from the array of values
     * output for a three tensor model by <code>ModelFit</code>.
     *
     * @param data
     *            The data values.
     *
     * @return The properties: {x1, y1, z1, prop1, trace1, fa1, x2, y2, z2,
     *         prop2, trace2, fa2, prop3, x3, y3, z3, trace3, fa3}.
     */
    protected double[] getThreeTensorProperties(double[] data) {
        double[] properties = new double[NUMTHREETENSORPROPERTIES];
        
        DT d1 = new DT(data[4], data[5], data[6], data[7], data[8], data[9]);
        DT d2 = new DT(data[11], data[12], data[13], data[14], data[15], data[16]);
        DT d3 = new DT(data[18], data[19], data[20], data[21], data[22], data[23]);
        
        double[] pd1 = d1.getPD();
        double[] pd2 = d2.getPD();
        double[] pd3 = d3.getPD();
        for (int i = 0; i < pd1.length; i++) {
            properties[i] = pd1[i];
            properties[i + 6] = pd2[i];
            properties[i + 12] = pd3[i];
        }
        
        properties[3] = data[3];
        properties[4] = d1.trace();
        properties[5] = d1.fa();
        
        properties[9] = data[10];
        properties[10] = d2.trace();
        properties[11] = d2.fa();
        
        properties[15] = data[17];
        properties[16] = d3.trace();
        properties[17] = d3.fa();
        
        return properties;
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
     * @return The list of statistics.
     */
    protected double[] computeStats(double[][] props, ModelIndex inversion) {
        
        if (inversion.numDTs == 1) {
            return computeStatsDT(props);
        }
        else if (inversion.numDTs == 2) {
            return computeStatsTwoTensor(props);
        }
        else {
            return computeStatsThreeTensor(props);
        }
    }
    
    
    /**
     * Computes statistics of a list of single tensor properties.
     *
     * @param props
     *            The array of properties.
     *
     * @return The list of statistics: {meanX, meanY, meanZ, dyad1, dyad2,
     *         dyad3, meanTrace, STD Trace, meanFA, stdFA}
     */
    protected double[] computeStatsDT(double[][] props) {
        
        double[] meanDir = oneDirectionStats(props, 0);
        double[] traceStats = meanSTDMaxAndMin(props, 3);
        double[] faStats = meanSTDMaxAndMin(props, 4);
        
        double[] stats = new double[NUMDTSTATS];
        for (int i = 0; i < meanDir.length; i++) {
            stats[i] = meanDir[i];
        }
        stats[DIRSTATS] = traceStats[0];
        stats[DIRSTATS + 1] = traceStats[1];
        stats[DIRSTATS + 2] = faStats[0];
        stats[DIRSTATS + 3] = faStats[1];
        
        return stats;
    }
    
    
    /**
     * Computes statistics of a list of two tensor properties.
     *
     * @param props
     *            The array of properties.
     *
     * @return The list of statistics: {mean1X, mean1Y, mean1Z, dyad11, dyad12,
     *         dyad13, meanTrace1, STD Trace1, meanFA1, stdFA1, mean2X, mean2Y,
     *         mean2Z, dyad21, dyad22, dyad23, meanTrace2, STD Trace2, meanFA2,
     *         stdFA2, meanProp, stdProp}
     */
    protected double[] computeStatsTwoTensor(double[][] props) {
        
        int firstDirIndex = 0;
        int propsPerComponent = 6;
        int numComponents = 2;
        
        // Sort the corresponding components.
        orderComponentsByDirection(props, numComponents, propsPerComponent, firstDirIndex);
        
        int dir1Index = 0;
        int dir2Index = propsPerComponent;
        
        double[] meanDir1 = oneDirectionStats(props, dir1Index);
        double[] meanDir2 = oneDirectionStats(props, dir2Index);
        
        // Compute the stats of the other properties.
        double[] propStats = meanSTDMaxAndMin(props, dir1Index + 3);
        
        double[] trace1Stats = meanSTDMaxAndMin(props, dir1Index + 4);
        double[] fa1Stats = meanSTDMaxAndMin(props, dir1Index + 5);
        
        double[] trace2Stats = meanSTDMaxAndMin(props, dir2Index + 4);
        double[] fa2Stats = meanSTDMaxAndMin(props, dir2Index + 5);
        
        // Gather all the results together.
        double[] stats = new double[NUMTWOTENSORSTATS];
        
        for (int i = 0; i < meanDir1.length; i++) {
            stats[i] = meanDir1[i];
        }
        stats[DIRSTATS] = trace1Stats[0];
        stats[DIRSTATS + 1] = trace1Stats[1];
        stats[DIRSTATS + 2] = fa1Stats[0];
        stats[DIRSTATS + 3] = fa1Stats[1];
        
        for (int i = 0; i < meanDir2.length; i++) {
            stats[i + NUMDTSTATS] = meanDir2[i];
        }
        stats[DIRSTATS + NUMDTSTATS] = trace2Stats[0];
        stats[1 + DIRSTATS + NUMDTSTATS] = trace2Stats[1];
        stats[2 + DIRSTATS + NUMDTSTATS] = fa2Stats[0];
        stats[3 + DIRSTATS + NUMDTSTATS] = fa2Stats[1];
        
        stats[2 * NUMDTSTATS] = propStats[0];
        stats[1 + 2 * NUMDTSTATS] = propStats[1];
        
        return stats;
    }
    
    
    /**
     * Computes statistics of a list of three tensor properties.
     *
     * @param props
     *            The array of properties.
     *
     * @return The list of statistics: {mean1X, mean1Y, mean1Z, dyad11, dyad12,
     *         dyad13, meanTrace1, STD Trace1, meanFA1, stdFA1, mean2X, mean2Y,
     *         mean2Z, dyad21, dyad22, dyad23, meanTrace2, STD Trace2, meanFA2,
     *         stdFA2, mean3X, mean3Y, mean3Z, dyad31, dyad32, dyad33,
     *         meanTrace3, STD Trace3, meanFA3, stdFA3, meanProp1, stdProp1,
     *         meanProp2, stdProp2}
     */
    protected double[] computeStatsThreeTensor(double[][] props) {
        
        int firstDirIndex = 0;
        int propsPerComponent = 6;
        int numComponents = 3;
        
        // Sort the corresponding components.
        orderComponentsByDirection(props, numComponents, propsPerComponent, firstDirIndex);
        
        int dir1Index = 0;
        int dir2Index = propsPerComponent;
        int dir3Index = propsPerComponent * 2;
        
        double[] meanDir1 = oneDirectionStats(props, dir1Index);
        double[] meanDir2 = oneDirectionStats(props, dir2Index);
        double[] meanDir3 = oneDirectionStats(props, dir3Index);
        
        // Compute the stats of the other properties.
        double[] prop1Stats = meanSTDMaxAndMin(props, dir1Index + 3);
        double[] prop2Stats = meanSTDMaxAndMin(props, dir2Index + 3);
        
        double[] trace1Stats = meanSTDMaxAndMin(props, dir1Index + 4);
        double[] fa1Stats = meanSTDMaxAndMin(props, dir1Index + 5);
        
        double[] trace2Stats = meanSTDMaxAndMin(props, dir2Index + 4);
        double[] fa2Stats = meanSTDMaxAndMin(props, dir2Index + 5);
        
        double[] trace3Stats = meanSTDMaxAndMin(props, dir3Index + 4);
        double[] fa3Stats = meanSTDMaxAndMin(props, dir3Index + 5);
        
        // Gather all the results together.
        double[] stats = new double[NUMTHREETENSORSTATS];
        
        for (int i = 0; i < meanDir1.length; i++) {
            stats[i] = meanDir1[i];
        }
        stats[DIRSTATS] = trace1Stats[0];
        stats[DIRSTATS + 1] = trace1Stats[1];
        stats[DIRSTATS + 2] = fa1Stats[0];
        stats[DIRSTATS + 3] = fa1Stats[1];
        
        for (int i = 0; i < meanDir2.length; i++) {
            stats[i + NUMDTSTATS] = meanDir2[i];
        }
        stats[DIRSTATS + NUMDTSTATS] = trace2Stats[0];
        stats[1 + DIRSTATS + NUMDTSTATS] = trace2Stats[1];
        stats[2 + DIRSTATS + NUMDTSTATS] = fa2Stats[0];
        stats[3 + DIRSTATS + NUMDTSTATS] = fa2Stats[1];
        
        for (int i = 0; i < meanDir3.length; i++) {
            stats[i + NUMDTSTATS * 2] = meanDir3[i];
        }
        stats[DIRSTATS + NUMDTSTATS * 2] = trace3Stats[0];
        stats[1 + DIRSTATS + NUMDTSTATS * 2] = trace3Stats[1];
        stats[2 + DIRSTATS + NUMDTSTATS * 2] = fa3Stats[0];
        stats[3 + DIRSTATS + NUMDTSTATS * 2] = fa3Stats[1];
        
        stats[3 * NUMDTSTATS] = prop1Stats[0];
        stats[1 + 3 * NUMDTSTATS] = prop1Stats[1];
        
        stats[2 + 3 * NUMDTSTATS] = prop2Stats[0];
        stats[3 + 3 * NUMDTSTATS] = prop2Stats[1];
        
        return stats;
    }
    
    
    /**
     * Reorders the array of inversion properties when it contains multiple
     * directions to maximize the similarity of the directions in each
     * component.
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
    public void orderComponentsByDirection(double[][] props, int numComponents,
    int propsPerComponent, int firstDirIndex) {
        
        for (int j = 0; j < numComponents - 1; j++) {
            
            int dirInd = j * propsPerComponent + firstDirIndex;
            double[] meanDir = oneDirectionStats(props, dirInd);
            
            boolean done = false;
            while (!done) {
                
                done = true;
                
                // Check each alternative direction and swap with the
                // existing one if it improves the match between the
                // current direction j and the mean direction for
                // component j.
                for (int i = 0; i < props.length; i++) {
                    
                    double dir1j = Math.abs(meanDir[0] * props[i][dirInd] + meanDir[1]
                    * props[i][dirInd + 1] + meanDir[2] * props[i][dirInd + 2]);
                    
                    for (int n = j + 1; n < numComponents; n++) {
                        int dirInd2 = n * propsPerComponent + firstDirIndex;
                        double dir1n = Math.abs(meanDir[0] * props[i][dirInd2]
                        + meanDir[1] * props[i][dirInd2 + 1] + meanDir[2]
                        * props[i][dirInd2 + 2]);
                        
                        if (dir1j < dir1n) {
                            
                            // Do the swap.
                            for (int k = 0; k < propsPerComponent; k++) {
                                double temp = props[i][dirInd + k];
                                props[i][dirInd + k] = props[i][dirInd2 + k];
                                props[i][dirInd2 + k] = temp;
                            }
                            
                            dir1j = dir1n;
                            
                            // We need another iteration if any swaps are made.
                            done = false;
                        }
                    }
                }
                
                if (!done) {
                    meanDir = oneDirectionStats(props, dirInd);
                }
            }
        }
    }
    
    
    /**
     * Computes the principal eigenvector (ex, ey, ez) and eigenvalues (l1, l2,
     * l3) of the mean dyadic of a list of vectors.
     *
     * @param props
     *            An array of properties from a number of trials.
     *
     * @param col
     *            The first column of the vector to compute the statistics of.
     *
     * @return The direction and length of the mean vector: {ex, ey, ez, l1, l2,
     *         l3}.
     */
    public static double[] oneDirectionStats(double[][] props, int col) {
        
        double[] dir = new double[3];
        
        RealMatrix meanDyadic = new RealMatrix(3, 3);
        
        int samples = 0;
        for (int i = 0; i < props.length; i++) {
            
            dir[0] = props[i][col + 0];
            dir[1] = props[i][col + 1];
            dir[2] = props[i][col + 2];
            
            // Check that a direction actually exists.
            if (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] > 0.0) {
                samples += 1;
                meanDyadic = meanDyadic.add(dyadic(dir));
            }
        }
        
        if (samples == 0) {
            // Just return zeros
            return new double[DIRSTATS];
        }
        
        if (samples == 1) {
            // The jacobi calculation occasionally goes wrong, but we know
            // that the mean direction is the single sample and the
            // eigenvalues must be 1 0 0.
            double[] st = new double[DIRSTATS];
            st[0] = dir[0];
            st[1] = dir[1];
            st[2] = dir[2];
            st[3] = 1.0;
            return st;
        }
        
        meanDyadic.scale(1.0 / samples);
        //System.err.println(meanDyadic);
        RealMatrix[] mdEig = meanDyadic.jacobi();
        
        // Find the indexes of the largest and smallest eigenvalues.
        int maxInd = 0;
        int minInd = 0;
        for (int j = 1; j < 3; j++) {
            if (mdEig[0].entry(j, j) >= mdEig[0].entry(maxInd, maxInd)) {
                maxInd = j;
            }
            if (mdEig[0].entry(j, j) < mdEig[0].entry(minInd, minInd)) {
                minInd = j;
            }
        }
        int midInd = 3 - minInd - maxInd;
        
        // Gather the results in a single array.
        double[] results = new double[DIRSTATS];
        for (int i = 0; i < 3; i++) {
            results[i] = mdEig[1].entry(i, maxInd);
        }
        results[3] = mdEig[0].entry(maxInd, maxInd);
        results[4] = mdEig[0].entry(midInd, midInd);
        results[5] = mdEig[0].entry(minInd, minInd);
        
        return results;
    }
    
    
    /**
     * Computes the dyadic of a direction.
     *
     * @param dir
     *            The direction: {x, y, z}.
     *
     * @return The dyadic matrix.
     */
    public static RealMatrix dyadic(double[] dir) {
        RealMatrix dy = new RealMatrix(dir.length, dir.length);
        
        for (int i = 0; i < dir.length; i++) {
            for (int j = 0; j < dir.length; j++) {
                dy.setEntry(i, j, dir[i] * dir[j]);
            }
        }
        
        return dy;
    }
    
    
    /**
     * Computes the mean, standard deviation, maximum and minimum of a column of
     * a 2D array.
     *
     * @param props
     *            The array.
     *
     * @param col
     *            The column.
     *
     * @return The mean, standard deviation, max and min: {mean, std, max, min}.
     */
    public static double[] meanSTDMaxAndMin(double[][] props, int col) {
        
        double[] stats = new double[4];
        
        // If no trials were successful, stats are meaningless so
        // just return zeros.
        if (props.length == 0) return stats;
        
        double sum = 0.0;
        double sumSq = 0.0;
        double max = props[0][col];
        double min = props[0][col];
        for (int i = 0; i < props.length; i++) {
            sum += props[i][col];
            sumSq += props[i][col] * props[i][col];
            max = (props[i][col] > max) ? props[i][col] : max;
            min = (props[i][col] < min) ? props[i][col] : min;
        }
        
        stats[0] = sum / (double) props.length;
        stats[1] = Math.sqrt(sumSq / (double) props.length - stats[0] * stats[0]);
        stats[2] = max;
        stats[3] = min;
        
        return stats;
    }
    
}

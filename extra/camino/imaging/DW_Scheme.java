package imaging;

import misc.*;
import numerics.*;

import tools.ArrayOps;

import java.awt.GradientPaint;
import java.io.*;
import java.text.*;
import java.util.regex.*;
import java.util.logging.Logger;
import java.util.*;



/**
 * Top level class for diffusion-weighted imaging schemes.
 *
 * @author Philip Cook, Daniel Alexander, Matt Hall
 * @version $Id$
 *  
 */
public abstract class DW_Scheme {

    /**
     * Gyromagnetic ratio of protons in water, in units of s^{-1} T^{-1}. 
     *
     */
    public static final double GAMMA = 2.6751525E8;

    /**
     * Logging object
     */
    private static final Logger logger = Logger.getLogger("camino.imaging.DW_Scheme");

    /** may be passed to gradOrder method to get X, Y, Z ordering of the gradient directions. */
    public static final int [] gradXYZ = {0,1,2};

    /** may be passed to gradOrder method to get X, Z, Y ordering of the gradient directions. */
    public static final int [] gradXZY = {0,2,1};

    /** may be passed to gradOrder method to get Y, X, Z ordering of the gradient directions. */
    public static final int [] gradYXZ = {1,0,2};

    /** may be passed to gradOrder method to get Y, Z, X ordering of the gradient directions. */
    public static final int [] gradYZX = {1,2,0};

    /** may be passed to gradOrder method to get Z, X, Y ordering of the gradient directions. */
    public static final int [] gradZXY = {2,0,1};

    /** may be passed to gradOrder method to get Z, Y, X ordering of the gradient directions. */
    public static final int [] gradZYX = {2,1,0};


    /**
     * This static method simply returns a default scheme for
     * various applications that require a scheme object without
     * actually using it.  Scheme objects returned from this
     * function should never be used other than as dummies.
     */
    public static DW_Scheme nullScheme() {
        return B_VectorScheme.elecPointSetScheme(1, 6, 1E9);
    }


    /**
     * Total number of measurements in this scheme.
     */
    protected final int numMeas;

    /**
     * Number of zero measurements. Zero measurements are defined as those that have zero gradient directions. They 
     * may have other non-zero parameters (such as echo time).
     */
    protected final int numZeroMeas;

    /** List of zero measurements.*/
    private final boolean[] zero;

    /** Gradient directions for this scheme. */
    protected double[][] gDir;


    /**
     * Default constructor is private to force initialization of gradient directions.
     */
    protected DW_Scheme() {

        numMeas = 0;
        numZeroMeas = 0;
        zero = null;
        gDir = null;
    }


    /**
     * Constructor for concrete subclasses. 
     *
     * @param gradDir gradient directions, one per measurement. Should be zero
     * for b=0 measurements. The remaining directions will be normalized to unit length if 
     * they are not already of unit length.
     *
     */
    protected DW_Scheme(double[][] gradDir) {

        gDir = normalizeGradDirs(gradDir);

        numMeas = gDir.length;

        zero = new boolean[numMeas];

        int zeroCounter = 0;

        for (int i = 0; i < numMeas; i++) {
            System.arraycopy(gradDir[i], 0, gDir[i], 0, 3);

            if (gDir[i][0] == 0.0 && gDir[i][1] == 0.0 && gDir[i][2] == 0.0) {
                zero[i] = true;
                zeroCounter += 1;
            }
        }

        numZeroMeas = zeroCounter;

    }



    /**
     * Gets the number of measurements in each voxel.
     *
     * @return the number of measurements in each voxel.
     */
    public int numMeasurements() {
        return numMeas;
    }


    /**
     * Gets the number of zero (b=0) measurements in each voxel.
     *
     * @return the number of b=0 measurements in each voxel.
     */
    public int numZeroMeasurements() {
        return numZeroMeas;
    }


    /**
     * Gets the number of  diffusion weighted (b > 0) measurements in each voxel.
     *
     * @return the number of diffusion weighted measurements in each voxel.
     */
    public int numNonZeroMeasurements() {
        return numMeasurements() - numZeroMeasurements();
    }


    /**
     * Gets a copy of the i-th gradient direction.
     *
     * @return a copy of the i-th gradient direction.
     */
    public final double[] getG_Dir(int i) {
        double[] copy = new double[3];

        copy[0] = gDir[i][0];
        copy[1] = gDir[i][1];
        copy[2] = gDir[i][2];

        return copy;
    }



    /**
     * @return gradient directions corresponding to non-zero measurements.
     *
     */
    public final double[][] getNonZeroG_Dirs() {
        double[][] nonZeroG_Dirs = new double[numMeas-numZeroMeas][3];

        int ind=0;

        for(int i = 0; i < numMeas; i++) {
            if(!zero[i]){
                nonZeroG_Dirs[ind][0] = gDir[i][0];
                nonZeroG_Dirs[ind][1] = gDir[i][1];
                nonZeroG_Dirs[ind][2] = gDir[i][2];
                ind++;
            }
        }

        if (ind != numMeas - numZeroMeas) {
            throw new LoggedException("Expected " + (numMeas - numZeroMeas) + " non-zero measurements but got " + 
                    ind);
        }

        return nonZeroG_Dirs;
    }


    /**
     * @return b-values corresponding to non-zero measurements.
     *
     */
    public final double[] getNonZeroB_Values() {
        double[] nonZeroB = new double[numMeas-numZeroMeas];

        int ind=0;

        for(int i = 0; i < numMeas; i++) {
            if(!zero[i]){
                nonZeroB[ind++] = getB_Value(i);
            }
        }

        if (ind != numMeas - numZeroMeas) {
            throw new LoggedException("Expected " + (numMeas - numZeroMeas) + " non-zero measurements but got " + 
                    ind);
        }

        return nonZeroB;

    }


    /**
     * Normalizes a voxel set of measurements acquired with this
     * sequence by the geometric mean of the b=0 measurements.  If the
     * mean is not positive, no normalization is performed.
     * 
     * @param data unnormalized set of measurements.
     * 
     * @return normalized measurements. Only data with b > 0 are returned.
     */
    public final double[] normalizeData(double[] data) {

        if (data.length != numMeas) {
            throw new LoggedException("Data to be normalized does not match scheme");
        }

        if (numZeroMeas == 0) {
            double[] notNormed = new double[numMeas];

            System.arraycopy(data, 0, notNormed, 0, numMeas);

            logger.warning("No b=0 data, cannot normalize");

            return notNormed;
        }

        return normalizeData(data, geoMeanZeroMeas(data));
    }


    /**
     * Normalizes a voxel set of measurements acquired with this sequence by the geometric mean 
     * of a specified set of measurements.
     * 
     * @param data unnormalized set of measurements.
     *
     * @param indices the integer indices of the measurements to use for the normalization.
     * 
     * @return normalized measurements. Data that are used to compute the normalization constant 
     * are not returned.
     */
    public final double[] normalizeData(double[] data, int[] indices) {

        if (data.length != numMeas) {
            throw new LoggedException("Data to be normalized does not match scheme");
        }

        double[] normed = new double[numMeas - indices.length];

        double[] dataForNormC = new double[indices.length];

        boolean[] skipMeas = new boolean[numMeas];

        for (int i = 0; i < indices.length; i++) {
            dataForNormC[i] = data[indices[i]];
            skipMeas[indices[i]] = true;
        }

        double normC = ArrayOps.geoMean(dataForNormC);

        int counter = 0;

        for (int i = 0; i < numMeas; i++) {

            if (!skipMeas[i]) {
                normed[counter++] = data[i] / normC;
            }
        }

        return normed;
    }



    /**
     * Normalizes a voxel set of measurements acquired with this
     * sequence. The normalization constant must be specified. This method is intended for use 
     * when there are no zero measurements and you want to normalize by an extrapolated zero signal.
     * 
     * @param data unnormalized set of measurements.
     * @param normC the normalization constant, which must be positive.
     * 
     * @return normalized measurements. Only data with b > 0 are returned.
     */
    public final double[] normalizeData(double[] data, double normC) {


        if (data.length != numMeas) {
            throw new LoggedException("Data to be normalized does not match scheme");
        }

        // if normC is zero, negative, or NaN
        if (!(normC > 0.0)) {
            logger.warning("Can't use normalization constant " + normC);
	    // do nothing if normC is bad, return original data with zero measurements removed
	    // but without normalizing anything
	    normC = 1.0;
	}

        double[] normed = new double[numMeas - numZeroMeas];

        int next = 0;

        for (int i = 0; i < numMeas; i++) {
            if (!zero[i]) {
                normed[next++] = data[i] / normC;
            }
        }

        return normed;
    }



    /**
     * Gets any imaging parameter that is stored in an array, for the non-zero measurements.
     *
     * @param value the sequence parameter in an array, where there is one parameter per measurement.
     *
     * @return the elements of the array <code>value</code> corresponding to non-zero measurements.
     */
    protected final double[] getNonZeroParam(double[] value) {
        double[] retArr = new double[numMeas-numZeroMeas];

        int ind=0;

        for(int i = 0; i < numMeas; i++) {
            if(!zero[i]){
                retArr[ind++]= value[i];
            }
        }

        if (ind != numMeas - numZeroMeas) {
            throw new LoggedException("Expected " + (numMeas - numZeroMeas) + " non-zero measurements but got " + 
                    ind);
        }

        return retArr;
    }


    /**
     * Computes the geometric mean of the b=0 measurements.
     * 
     * @param data unnormalized set of measurements.
     * 
     * @return the geometric mean, or 1.0 if the mean is undefined (no b=0 measurements).
     */
    public final double geoMeanZeroMeas(double[] data) {

        if (data.length != numMeas) {
            throw new LoggedException("Data to be normalized does not match scheme: data length is " + data.length);
        }
        if (numZeroMeas == 0) {
            return 1.0;
        }

        double gm = 1.0;

        for (int i = 0; i < numMeas; i++) {
            if (zero[i]) {
                gm = gm * data[i];
            }
        }

        return Math.pow(gm, 1.0 / numZeroMeas);
    }


    /**
     * Determines whether a measurement has zero diffusion weighting.
     *
     * @return true if measurement i is has zero diffusion weighting, false otherwise.
     *
     */
    public boolean zero(int i) {
        return zero[i];
    }


    /**
     * The B-matrix B is defined as follows. Given log DW data in a column vector Y, Y = BX, where
     * B is the B-matrix and X = [ln (A*(0)), dxx, dxy, dxz, dyy, dyz, dzz]^T. Each row i of B contains
     * {1, -b g_1^2, -2 b g_1 g_2, -2 b g_1 g_3, -b g_2^2, -2 b g_2 g_3, -b g_3^2}, where g is the gradient
     * direction for measurement i and b is the b-value for measurement i.
     *
     * @return the B-matrix. 
     */
    public final RealMatrix getB_Matrix() {
        RealMatrix B = new RealMatrix(numMeas, 7);

        for (int i = 0; i < numMeas; i++) {

            double bValue = getB_Value(i);

            B.setEntry(i, 0, 1.0);
            B.setEntry(i, 1, -bValue * gDir[i][0] * gDir[i][0]);
            B.setEntry(i, 2, -2.0 * bValue * gDir[i][0] * gDir[i][1]);
            B.setEntry(i, 3, -2.0 * bValue * gDir[i][0] * gDir[i][2]);
            B.setEntry(i, 4, -bValue * gDir[i][1] * gDir[i][1]);
            B.setEntry(i, 5, -2.0 * bValue * gDir[i][1] * gDir[i][2]);
            B.setEntry(i, 6, -bValue * gDir[i][2] * gDir[i][2]);

        }

        return B;
    }


    /**
     * Gets the b-value for measurement i. 
     *
     * @return the b-value for measurement i.
     */
    public abstract double getB_Value(int i);


    /**
     * Negates the X component of the gradient directions.
     *
     * @return a new scheme that is identical to this one except that the X component of the gradients are negated.
     */
    public abstract DW_Scheme flipX();


    /**
     * Negates the Y component of the gradient directions.
     *
     * @return a new scheme that is identical to this one except that the Y component of the gradients are negated.
     */
    public abstract DW_Scheme flipY();

    /**
     * Negates the Z component of the gradient directions.
     *
     * @return a new scheme that is identical to this one except that the Z component of the gradients are negated.
     */
    public abstract DW_Scheme flipZ();

    /**
     * Gets a scheme composed of a subset of the measurements in this scheme.
     *
     * @param indices integer indices of the measurements to include in the subset.
     * 
     * @return a new scheme consisting of all measurements in the indices list, in the
     * same order as they appear in that list.
     */
    public abstract DW_Scheme getSubsetScheme(int[] indices);


    /**
     *
     * Reorders the gradient directions.
     *
     * @param order the new order of the gradients, where x=0, y=1, z=2. The original order is
     * [0 1 2], the new order may be any combination of x,y,z. Convenience variables grad___ are 
     * provided by the <code>DW_Scheme</code> class.
     *
     *
     * @see imaging.DW_Scheme#gradXYZ gradXYZ and related variables.
     */
    public abstract DW_Scheme gradOrder(int [] order); 




    /**
     * Reads a scheme from a scheme file.
     *
     * @return a scheme object appropriate for the given scheme file.
     */
    public static final DW_Scheme readScheme(String filename) {

        try {

            Vector<String> lines = new Vector<String>();

            Scanner fileScanner = new Scanner(new File(filename));

            // expect UK formatting ie . for decimal
            fileScanner.useLocale(Locale.UK);

            // Read in the file line by line.
            fileScanner.useDelimiter("\r\n|\n");

            // Store all the lines in a vector.
            while(fileScanner.hasNext()) {
                // ignore empty lines and comment lines

                String next = fileScanner.next().trim();

                if (next.length() > 0) {
                    if (!next.startsWith("#")) {
                        lines.add(next);
                    } 
                }
            }
            fileScanner.close();

            // now figure out which subclass scheme to construct, based on first line
            // trim leading or trailing white space
            String version = lines.elementAt(0).trim();

            // legacy compatibility: SchemeV0 starts with diffusion time
            try {

                Double.parseDouble(version);

                // this may also throw a NumberFormatException if there are bad values in the scheme file
                return B_VectorScheme.readQ_VectorScheme(lines);
            }
            catch (NumberFormatException e) {
                // not a V0 scheme file
            }

            // if we get here, first line is a version ID string
            lines.removeElementAt(0);

            if (Pattern.matches("VERSION:\\s*1", version) || 
                    Pattern.matches("VERSION:\\s*" + RectGradSteTanScheme.VERSION, version)) {
                return RectGradSteTanScheme.readScheme(lines);
            }
            if (Pattern.matches("VERSION:\\s*2", version) || 
                    Pattern.matches("VERSION:\\s*" + B_VectorScheme.VERSION, version)) {
                return B_VectorScheme.readScheme(lines);
            }
            if (Pattern.matches("VERSION:\\s*3", version) || 
                    Pattern.matches("VERSION:\\s*" + RectGradTRSE_Scheme.VERSION, version)) {
                return RectGradTRSE_Scheme.readScheme(lines);
            }
            if (Pattern.matches("VERSION:\\s*" + GradientWaveform_Scheme.VERSION, version)) {
                return GradientWaveform_Scheme.readScheme(lines);
            }
            if (Pattern.matches("VERSION:\\s*" + RectGradDPFG_Scheme.VERSION, version)) {
                return RectGradDPFG_Scheme.readScheme(lines);
            }
            if (Pattern.matches("VERSION:\\s*" + RectQuadraticGradSteTanScheme.VERSION, version)) {
                return RectQuadraticGradSteTanScheme.readScheme(lines);
            }

            throw new LoggedException("Unrecognized scheme version string: " + version);

        } 
        catch(IOException e) {
            throw new LoggedException(e);
        }

    }



    /**
     * Returns gradient directions as vectors normalized to unit length if non-zero.
     * Will not modify directions that are within 1E-5 of unit length. 
     *
     *
     */
    public static final double[][] normalizeGradDirs(double[][] g) {

        double[][] gDir = new double[g.length][3];

        // set to true if we need to normalize a gradient direction
        boolean neededNorm = false;

        for (int i = 0; i < gDir.length; i++) {
            gDir[i][0] = g[i][0];
            gDir[i][1] = g[i][1];
            gDir[i][2] = g[i][2];

            if (gDir[i][0] == 0.0 && gDir[i][1] == 0.0 && gDir[i][2] == 0.0) {
                continue;
            }
            else {

                double modG_Dir = Math.sqrt(gDir[i][0] * gDir[i][0] + gDir[i][1] * gDir[i][1] + 
                        gDir[i][2] * gDir[i][2]);

                if (!(Math.abs(1.0 - modG_Dir) < 1E-4)) {
                    for (int j = 0; j < 3; j++) {

                        if (Double.isNaN(gDir[i][j])) {
                            throw new LoggedException("Can't normalize gradient " + i + ":" +
                                    g[i][0] + " " + g[i][1] + " " + g[i][2]);
                        }

                        gDir[i][j] /= modG_Dir;

                    }

                    neededNorm = true;
                }

            }
        }


        if (neededNorm) {
            logger.info("Some measurements had non unit gradient directions, which have been normalized.");
        }

        return gDir;

    }
}

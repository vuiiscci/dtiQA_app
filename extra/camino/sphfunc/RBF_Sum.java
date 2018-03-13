package sphfunc;

import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Function of the sphere represented using radial basis functions.
 * 
 * <dt>Description:
 * 
 * <dd>General class for spherical functions represented as the sum of radial
 * basis functions. The specific RBF is specified in subclasses.
 *
 * Note that the array of basis function centres is static to avoid 
 * excessive memory usage.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public abstract class RBF_Sum extends LinearBasisSum {


    /**
     * The default Pointset to use
     */
    private static final int DEFAULT_RBF_POINTSET = 246;


    /**
     * This is the array of sample points on the sphere.
     */
    protected static double[][] pts = SphericalPoints.getElecPointSet(DEFAULT_RBF_POINTSET);


    /**
     * The constructor takes an array containing the coefficients.
     * 
     * @param coeffs
     *            RBF coefficients.
     */
    public RBF_Sum(double[] coeffs) {

        // Check that the coefficients match the RBF pts
        if (coeffs.length != pts.length + 2) {
            throw new RuntimeException("RBF_sum was only given " + (coeffs.length-2)
                    + " coefficients for " + pts.length + " basis functions.");
        }

        // Copy the coefficients for storage.
        c = new double[coeffs.length - 2];
        for (int i = 0; i < c.length; i++) {
            c[i] = coeffs[i+2];
        }
    }


    /**
     * The radial basis function itself is abstract and defined in subclasses.
     * 
     * @param x
     *            x-coordinate.
     * 
     * @param y
     *            y-coordinate.
     * 
     * @param z
     *            z-coordinate.
     * 
     * @param p
     *            The point defining the RBF.
     * 
     * @return f(x, y, z; p).
     */
    public abstract double rbf(double x, double y, double z, double[] p);


    /**
     * Overridden to sum the weighted radial basis functions at the specified
     * point.
     * 
     * @param x
     *            x-coordinate.
     * 
     * @param y
     *            y-coordinate.
     * 
     * @param z
     *            z-coordinate.
     * 
     * @return f(x, y, z).
     */
    public double getRadius(double x, double y, double z) {

        // Compute the value of each RBF at (x, y, z) and
        // sum the values weighted by the coefficients.
        double total = 0.0;
        double rbfTotal = 0.0;
        for (int i = 0; i < pts.length; i++) {
            double r = rbf(x, y, z, pts[i]);
            rbfTotal += r;
            total += r * c[i];
        }

        // Normalize the result and return it
        return total / rbfTotal;
    }


    /**
     * Allows setting of the directions defining the radial basis.
     * 
     * @param p
     *            An array of points on the sphere.
     */
    public static void setPoints(double[][] p) {
        pts = new double[p.length][p[0].length];
        for (int i = 0; i < p.length; i++) {
            for (int j = 0; j < p[0].length; j++) {
                pts[i][j] = p[i][j];
            }
        }
    }


    /**
     * Returns the number of functions in the basis.
     * 
     * @return The number of functions in the basis.
     */
    public static int numPoints() {

        if (pts == null) {
            return 0;
        }

        return pts.length;
    }

}


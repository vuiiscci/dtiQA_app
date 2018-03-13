package sphfunc;

/**
 * <dl>
 *
 * <dt>Purpose:
 *
 * <dd> Implements the radial basis function used in Dave Tuch's
 * original q-ball implementation..
 * 
 * <dt>Description:
 *
 * <dd> The function is f(x) = exp(-acos(x.p)^2/sigma^2), where y is a
 * fixed point on the sphere and sigma is a scaling parameter.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 *
 * @version $Id$
 *  
 */
public class TuchRBF extends LinearBasisFunction {

    private double[] p;
    private double sigma;


    /**
     * Constructor requires parameters of the function.  The
     * constructor assumes normalization of newP and does not check.
     *
     * @param newP The centre of the function.
     *
     * @param newSigma The scale of the function.
     */
    public TuchRBF(double[] newP, double newSigma) {
        p = new double[3];
        for(int i=0; i<3; i++) {
            p[i] = newP[i];
        }
        sigma = newSigma;
    }


    /**
     * Evaluates the function at (x,y,z).  The function does not
     * check normalization of (x,y,z).
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

        return rbf(x, y, z, p, sigma);

    }


    /**
     * Evaluates the function at (x,y,z) assuming the centre of the
     * radial basis function is p.  The function does not check
     * normalization of (x,y,z).
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
     *            RBF centre.
     * 
     * @return f(x, y, z).
     */
    public static double rbf(double x, double y, double z, double[] p, double sig) {


        // Compute the dot product
        double xdy = x * p[0] + y * p[1] + z * p[2];

        // Check it is within the range [-1, 1].
        xdy = (xdy < 1.0) ? xdy : 1.0;
        xdy = (xdy > -1.0) ? xdy : -1.0;

        // Compute the inverse cosine
        double acosXdY = Math.acos(Math.abs(xdy));

        // Return the value
        return Math.exp(-acosXdY * acosXdY / (sig * sig));
    }

    

    
}


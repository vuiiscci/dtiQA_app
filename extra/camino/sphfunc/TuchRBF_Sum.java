package sphfunc;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Function of the sphere represented using Dave Tuch's spherical radial
 * basis function.
 * 
 * <dt>Description:
 * 
 * <dd>This class uses the radial basis function:
 * 
 * f(x) = exp(-acos(x.y)^2/sigma^2).
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class TuchRBF_Sum extends RBF_Sum {

    /**
     * This is the scaling parameter in the RBFs (Tuch calls it odfsigma in his
     * qball code).  The default value of this parameter is 15*pi/180
     */
    private static double sigma = 0.2617993877991494;


    public TuchRBF_Sum(double[] coeffs) {
        super(coeffs);
    }


    /**
     * Tuch's radial basis function.
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
    public double rbf(double x, double y, double z, double[] p) {
        return TuchRBF.rbf(x,y,z,p,sigma);
    }


    /**
     * Allows setting of the scaling parameter.
     * 
     * @param newSigma the new value of the scaling parameter
     */
    public static void setSigma(double newSigma) {
        sigma = newSigma;
    }

    /**
     * Returns the value of the scaling parameter
     *
     * @return The value of the scaling parameter
     */
    public static double getSigma() {
        return sigma;
    }

    public LinearBasisFunction basisFunction(int i) {
        return new TuchRBF(pts[i], sigma);
    }

    /**
     * Returns the parameters used in the basis sum as a string
     *
     * @return The parameters of the basis sum
     */
    public String getSettings(){
	String details = "basis type = rbf\nrbf sigma = " + sigma + "\nrbf pointset = " + pts.length + "\n\n";
	return details;
    }
}


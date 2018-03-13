package sphfunc;

/**
 * <dl>
 *
 * <dt>Purpose:
 *
 * <dd> General class for functions of the sphere represented as the
 * sum of basis functions.
 * 
 * <dt>Description:
 *
 * <dd> Contains a list of coefficients and generally applicable
 * methods.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 *
 * @version $Id$
 *  
 */
public abstract class LinearBasisSum extends SphericalFunction {

    /**
     * This is the array that holds the coefficients.
     */
    protected double[] c;


    /**
     * Default constructor
     */
    public LinearBasisSum() {
        super();
    }


    /**
     * Returns the basis function with weight c[i].
     *
     * @param i The index of the basis function.
     *
     * @return The basis function.
     */
    public abstract LinearBasisFunction basisFunction(int i);
    
    /**
     * Returns the number of coeffs (terms in the sum).
     * 
     * @return The number of coeffs
     */
    public  int numBasisFunctions(){
	return c.length;
    }

    /**
     * Returns the parameters used in the basis sum as a string
     *
     * @return The parameters of the basis sum
     */
    public abstract String getSettings();

}


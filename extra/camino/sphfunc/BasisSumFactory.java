package sphfunc;
import tools.*;
import misc.LoggedException;
import misc.SphericalPoints;
/**
 * <dl>
 *
 * <dt>Purpose:
 *
 * <dd>Creates  Basis Sum object on request
 * 
 * <dt>Description:
 *
 * <dd> Factory for creating LinearBasisSum Objects.
 * Used by linear reconstruction algorithms (such as QBallMX)
 * 
 * </dl>
 * 
 * @author Kiran Seunarine
 *
 * $Id$
 *  
 */
public class BasisSumFactory {

    /**
     * Creates a TuchRBF_Sum object
     */
    public static final int TUCH_RBF = 1;

    /**
     * Creates an EvenSHS object
     */
    public static final int SPHERICAL_HARMONICS = 2;

    /**
     * The default Pointset to use
     */
    private static final int DEFAULT_RBF_POINTSET = 246;

    /**
     * Returns a linear basis sum of type basisType.
     *
     * @param basisType the type of basis sum to create.
     *
     * @return A linear basis sum.
     */
    public static LinearBasisSum getBasisSum(int basisType){
	if(basisType == TUCH_RBF) {

	    // need to check the number of points > 0
	    int numCoeffs = RBF_Sum.numPoints();
	    if(numCoeffs > 0) {
		numCoeffs +=2;
	    }
	    else {
		RBF_Sum.setPoints(SphericalPoints.getElecPointSet(DEFAULT_RBF_POINTSET));
		numCoeffs = DEFAULT_RBF_POINTSET + 2;
	    }
	    double [] coeffs = getCoeffs(numCoeffs);

	    // use radial basis function representation
	    return new TuchRBF_Sum(coeffs);
	}
	else if (basisType == SPHERICAL_HARMONICS) {
	    int order = CL_Initializer.maxOrder;
	    int numCoeffs = ((order * order) + (3 * order) + 2)/2;
	    numCoeffs +=2;
	    double [] coeffs = getCoeffs(numCoeffs);

	    // use spherical harmonic representation
	    return new EvenSHS(coeffs, order);
	}
	else {
	    // basis code is invalid!
	    throw new LoggedException("Unrecognized basis type: " + basisType);
	}
    }

    /**
     * Generates a list of coefficients.
     *
     * @param numCoeffs the length of the coefficient array.
     *
     * @return an array of coefficients (set to 1).
     */
    private static double [] getCoeffs(int numCoeffs) {
	double [] coeffs = new double [numCoeffs];
	for(int i=0; i< numCoeffs; i++) {
		coeffs[i]=1;
	    }
	return coeffs;
    }

}


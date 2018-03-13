package numerics;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Spherical harmonic functions and associated Legendre polynomials.
 * 
 * <dt>Description:
 * 
 * <dd>Implements methods to compute the spherical harmonics and the associated
 * Legendre polynomials using the recurrence relation described in NRC.
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id: SphericalHarmonics.java,v 1.3 2005/08/18
 *         11:12:21 ucacmgh Exp $
 *  
 */
public class SphericalHarmonics {

    /**
     * Computes spherical harmonic of order l, index m at colatitude theta and
     * longitude phi.
     * 
     * @param l
     *            The spherical harmonic order.
     * 
     * @param m
     *            The spherical harmonic index.
     * 
     * @param theta
     *            The angle of colatitude.
     * 
     * @param phi
     *            The angle of longitude.
     * 
     * @return Y_lm(theta, phi)
     */
    public static Complex Y(int l, int m, double theta, double phi)
            throws SphericalHarmonicException {
        int mPos = Math.abs(m);

        //Normalisation factor for spherical harmonic function.
        double factor = getFactor(l, mPos);

        //Associated Legendre poly of cos(theta)
        double thetaPart = plgndr(l, mPos, Math.cos(theta));

        //exp(i.m.phi)
        Complex phiPart = new Complex(Math.cos((double) mPos * phi), Math
                .sin((double) mPos * phi));

        //Combine them all to get the result.
        Complex result = phiPart.times(factor * thetaPart);

        //If m is negative take the conjugate of the positive
        //result, negated if m is odd.
        if (m < 0) {
            result = result.conjugate();
            if (mPos % 2 == 1) {
                result = result.negate();
            }
        }

        return result;
    }

    /**
     * Computes the factor of the Legendre polynomial used to compute the
     * spherical harmonics.
     * 
     * @param l
     *            Spherical harmonic order
     * 
     * @param m
     *            Spherical harmonic index
     * 
     * @return The scaling factor.
     */
    private static double getFactor(int l, int m) {

        double factor = 1.0;
        for (int i = l - m + 1; i <= l + m; i++) {
            factor *= (double) i;
        }
        factor = (double) (2 * l + 1) / factor;
        factor /= 4.0 * Math.PI;
        factor = Math.sqrt(factor);

        return factor;
    }

    /**
     * NRC method for computing associated Legendre polynomial values.
     * 
     * @param l
     *            The order of the polynomial
     * 
     * @param m
     *            The index of the polynomial
     * 
     * @param x
     *            The argument.
     * 
     * @return the value of the polynomial at x.
     */
    public static double plgndr(int l, int m, double x) throws SphericalHarmonicException {
        double fact = 0.0;
        double pll = 0.0;
        double pmm = 0.0;
        double pmmp1 = 0.0;
        double somx2 = 0.0;
        int i = 0;
        int ll = 0;

        if (m < 0 || m > l || Math.abs(x) > 1.0) {
            throw new SphericalHarmonicException("Bad arguments in routine plgndr");
        }
        pmm = 1.0;
        if (m > 0) {
            somx2 = Math.sqrt((1.0 - x) * (1.0 + x));
            fact = 1.0;
            for (i = 1; i <= m; i++) {
                pmm *= -fact * somx2;
                fact += 2.0;
            }
        }
        if (l == m)
            return pmm;
        else {
            pmmp1 = x * (2 * m + 1) * pmm;
            if (l == (m + 1))
                return pmmp1;
            else {
                for (ll = m + 2; ll <= l; ll++) {
                    pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                    pmm = pmmp1;
                    pmmp1 = pll;
                }
                return pll;
            }
        }
    }

    /**
     * Computes the number of spherical harmonic functions up to a specified
     * order.
     * 
     * @param order
     *            The order up to which the number of functions is required.
     * 
     * @return The number of functions.
     */
    public static int funcsUpTo(int order) {
        return (order + 1) * (order + 1);
    }

    /**
     * Computes the number of even spherical harmonic functions up to a
     * specified order.
     * 
     * @param order
     *            The order up to which the number of functions is required.
     * 
     * @return The number of functions.
     */
    public static int evenFuncsUpTo(int order) {

        // If order is odd, change it to the even order just below.
        if (order % 2 == 1) {
            order = order - 1;
        }

        return (order + 1) * (order / 2 + 1);
    }

}
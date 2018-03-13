package sphfunc;

import java.util.logging.Logger;

import numerics.*;

/**
 * <dl>
 *
 * <dt>Purpose:
 *
 * <dd> Real part of spherical harmonic function.
 * 
 * <dt>Description:
 *
 * <dd> Encapsulates the spherical function Re(Y_lm) for specific
 * order l and index m.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 *
 * $Id$
 *  
 */
public class RealSH extends LinearBasisFunction {

    private Logger logger = Logger.getLogger(this.getClass().getName());
    
    // Index and order of spherical harmonic.
    protected int m;
    protected int l;

    // Scaling factor for function.
    protected double scalar;


    /**
     * Constructor requires order and the index.  No scaling.
     *
     * @param order The spherical harmonic order
     *
     * @param index The spherical harmonic index.
     */
    public RealSH(int order, int index) {
        l = order;
        m = index;
        scalar = 1.0;
    }


    /**
     * Constructor requires order, index and a constant scaling factor.
     *
     * @param order The spherical harmonic order
     *
     * @param index The spherical harmonic index.
     *
     * @param scalar Constant scaling factor of function.
     */
    public RealSH(int order, int index, double sc) {
        l = order;
        m = index;
        scalar = sc;
    }


    /**
     * Overridden to evaluate the function directly as it is simpler
     * in polar coordinates.
     * 
     * @param theta
     *            The angle of colatitude.
     * 
     * @param phi
     *            The angle of longitude.
     * 
     * @return f(theta, phi).
     */
    public double getRadius(double theta, double phi) {

        Complex z;
        try {
            z = SphericalHarmonics.Y(l, m, theta, phi);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        return scalar*z.real();

    }


    /**
     * Overridden to convert the Cartesian parameters to spherical polars and
     * call the getRadius(theta, phi).
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
        double[] angles = toAngles(x, y, z);
        return getRadius(angles[0], angles[1]);
    }

    
    
    /** overridden great circle integral calculation. With Spherical harmonics
     *  we can extend perform an analytical  calculation to get an exact result.
     * 
     * \[
     *   \frac{1}{2\pi}
     *    \oint_{\forall \mathbf{q} \perpto \mathbf{u}} Y_l^m(\mathbf{q})d\mathbf{q}
     * 
     *     = 2\pi P_l(0)Y_l^m(\mathbf{u})
     * \]
     * 
     * Backus, Bull Math Seismol Soc Amer 54, 571-610 (1964)
     * 
     * (in english: "the integral over the great circle on the sphere defined by the 
     *  vector u is equal to the product of P_l(0) and Y_l^m(u) where P_l(x) is the
     * l-th legendre polynomial and Y_l^m is the l,m-th spherical harmonic)
     * 
     * Which is a very handy result!
     * 
     * @param u the axis defining the plane of the circle
     */
    public double greatCircleIntegral(double[] u){
        
        double[] angles = toAngles(u[0], u[1], u[2]);
        
        double legendrePoly=0.0;
        Complex Y= new Complex();
        
        try{
        	legendrePoly = SphericalHarmonics.plgndr(l, 0, 0.0);
            Y = SphericalHarmonics.Y(l, m, angles[0], angles[1]);
        }
        catch(SphericalHarmonicException she){
            String errMess= she.getMessage();
            logger.warning(errMess);
        }
        
        double sphHarmonic = Y.real();
        
        return 2.0*Math.PI*legendrePoly*sphHarmonic*scalar;
        
    }
}


package sphfunc;

import java.util.logging.Logger;

import numerics.*;

/**
 * <dl>
 *
 * <dt>Purpose:
 *
 * <dd> Imaginary part of spherical harmonic function.
 * 
 * <dt>Description:
 *
 * <dd> Encapsulates the spherical function Im(Y_lm) for specific
 * order l and index m.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 *
 * $Id$
 *  
 */
public class ImagSH extends RealSH {

    private Logger logger = Logger.getLogger(this.getClass().getName()); 
    
    
    
    /**
     * Constructor requires order and the index.
     *
     * @param order The spherical harmonic order
     *
     * @param index The spherical harmonic index.
     */
    public ImagSH(int order, int index) {
        super(order, index);
    }


    /**
     * Constructor requires order and the index.
     *
     * @param order The spherical harmonic order
     *
     * @param index The spherical harmonic index.
     *
     * @param scalar Constant scaling factor of function.
     */
    public ImagSH(int order, int index, double sc) {
        super(order, index, sc);
    }


    /**
     * Overridden to return the imaginary part rather than the real
     * part.
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

        return scalar*z.imag();

    }
    
    /** overridden great circle integral calculation. With Spherical
     *  harmonics we can extend perform an analytical calculation to
     *  get an exact result.
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
     * (in english: "the integral over the great circle on the sphere
     * defined by the vector u is equal to the product of P_l(0) and
     * Y_l^m(u) where P_l is the l-th legendre polynomial and Y_l^m
     * is the l,m-th spherical harmonic)
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
        
        double sphHarmonic = Y.imag();
        
        return 2.0*Math.PI*legendrePoly*sphHarmonic*scalar;
        
    }


    
    
    public static void main(String[] args){

        final double Rt2= 1.414213562;
        
        RealSH Y_00= new RealSH(0, 0);
        
        RealSH Y_20= new RealSH(2, 0);
        
        RealSH Y_22r= new RealSH(2, 2);
        ImagSH Y_22i= new ImagSH(2, 2);

        RealSH Y_21r= new RealSH(2, 1);
        ImagSH Y_21i= new ImagSH(2, 1);

        RealSH Y_40r= new RealSH(4, 0);
        ImagSH Y_40i= new ImagSH(4, 0);
        
        RealSH Y_42r= new RealSH(4, 2);
        ImagSH Y_42i= new ImagSH(4, 2);

        RealSH Y_43r= new RealSH(4, 3);
        ImagSH Y_43i= new ImagSH(4, 3);

        
        double[] u= new double[] {1.0/Rt2, 1.0/Rt2, 0.0};
        double[] u2= new double[] {0.0, 0.0, 1.0};
        
        
        double numIntY_00= Y_00.greatCircleIntegral(u);
        System.err.println("Y_00: analytic = "+numIntY_00+" 2pi.r = "+2.0*Math.PI*Y_00.getRadius(0.0, 0.0, 1.0));
        
        double numIntY_20= Y_20.greatCircleIntegral(u);
        System.err.println("Y_20: analytic = "+numIntY_20);
        numIntY_20= Y_20.numGreatCircleIntegral(u);
        System.err.println("Y_20: numerical = "+numIntY_20);
       
        double numIntY_21r= Y_21r.greatCircleIntegral(u);
        System.err.println("Y_21r: analytic = "+numIntY_21r);
        double numIntY_21i= Y_21i.greatCircleIntegral(u);
        System.err.println("Y_21i: analytic = "+numIntY_21i);
        
        numIntY_21r= Y_21r.greatCircleIntegral(u2);
        System.err.println("Y_21r: analytic = "+numIntY_21r);
        numIntY_21i= Y_21i.greatCircleIntegral(u2);
        System.err.println("Y_21i: analytic = "+numIntY_21i);
        
        numIntY_21r= Y_21r.numGreatCircleIntegral(u2);
        System.err.println("Y_21r: numerical = "+numIntY_21r);
        numIntY_21i= Y_21i.numGreatCircleIntegral(u2);
        System.err.println("Y_21i: numerical = "+numIntY_21i);
        
        double numIntY_22r= Y_22r.greatCircleIntegral(u);
        System.err.println("Y_22r: analytic = "+numIntY_22r);
        double numIntY_22i= Y_22i.greatCircleIntegral(u);
        System.err.println("Y_22i: analytic = "+numIntY_22i);
        
        numIntY_22r= Y_22r.greatCircleIntegral(u2);
        System.err.println("Y_22r: analytic = "+numIntY_22r);
        numIntY_22i= Y_22i.greatCircleIntegral(u2);
        System.err.println("Y_22i: analytic = "+numIntY_22i);
        
        numIntY_22r= Y_22r.numGreatCircleIntegral(u2);
        System.err.println("Y_22r: numerical = "+numIntY_22r);
        numIntY_22i= Y_22i.numGreatCircleIntegral(u2);
        System.err.println("Y_22i: numerical = "+numIntY_22i);
        
        double numIntY_40r= Y_40r.greatCircleIntegral(u);
        System.err.println("Y_40r: analytic = "+numIntY_40r);
        double numIntY_40i= Y_40i.greatCircleIntegral(u);
        System.err.println("Y_40i: analytic = "+numIntY_40i);

        numIntY_40r= Y_40r.greatCircleIntegral(u2);
        System.err.println("Y_40r: analytic = "+numIntY_40r);
        numIntY_40i= Y_40i.greatCircleIntegral(u2);
        System.err.println("Y_40i: analytic = "+numIntY_40i);

        numIntY_40r= Y_40r.numGreatCircleIntegral(u2);
        System.err.println("Y_40r: numerical = "+numIntY_40r);
        numIntY_40i= Y_40i.numGreatCircleIntegral(u2);
        System.err.println("Y_40i: numerical = "+numIntY_40i);

        double numIntY_42r= Y_42r.greatCircleIntegral(u);
        System.err.println("Y_42r: analytic = "+numIntY_42r);
        double numIntY_42i= Y_42i.greatCircleIntegral(u);
        System.err.println("Y_42i: analytic = "+numIntY_42i);
        
        double numIntY_43r= Y_43r.greatCircleIntegral(u);
        System.err.println("Y_43r: analytic = "+numIntY_43r);
        double numIntY_43i= Y_43i.greatCircleIntegral(u);
        System.err.println("Y_43i: analytic = "+numIntY_43i);
        
        numIntY_43r= Y_43r.greatCircleIntegral(u2);
        System.err.println("Y_43r: analytic = "+numIntY_43r);
        numIntY_43i= Y_43i.greatCircleIntegral(u2);
        System.err.println("Y_43i: analytic = "+numIntY_43i);
        
        numIntY_43r= Y_43r.numGreatCircleIntegral(u2);
        System.err.println("Y_43r: numerical = "+numIntY_43r);
        numIntY_43i= Y_43i.numGreatCircleIntegral(u2);
        System.err.println("Y_43i: numerical = "+numIntY_43i);
        
        
    }

}


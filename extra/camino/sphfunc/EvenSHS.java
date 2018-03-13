package sphfunc;

import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Real antipodally symmetric function represented as a spherical harmonic
 * series.
 * 
 * <dt>Description:
 * 
 * <dd>The spherical harmonic series contains only even terms and has the
 * constraint that c_{lm} = (-1)^m c^\star_{l -m}.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class EvenSHS extends LinearBasisSum {

    /**
     * This is the maximum order of the coefficients stored in the array.
     */
    private int maxOrder;


    /**
     * The constructor takes an array containing the spherical harmonic
     * coefficients in the format output by EvenSphHarmFitter.fit.
     * 
     * @param coeffs
     *            [exitCode, log A^\star(0), c00, c20, Re(c21), Im(c21),
     *            Re(c22), Im(c22), c40, Re(c41), Im(c41), ...]
     * 
     * @param order
     *            The maximum order of the series.
     */
    public EvenSHS(double[] coeffs, int order) {

        // Initialize the instance variables.
        maxOrder = order;
        c = new double[SphericalHarmonics.evenFuncsUpTo(maxOrder)];

        // Discard the first two elements of the input array and
        // copy the rest for storage.
        for (int i = 0; i < c.length; i++) {
            c[i] = coeffs[i + 2];
        }
    }


    /**
     * Overridden to sum the spherical harmonics at the specified point.
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

        // The formula in terms of the coefficients we have is: f(t,
        // p) = \sum_{l=0}^{order} (c_{l0} Y_{l0} + \sum_{m=1}^l (2
        // Re(c_{lm}) Re(Y_{lm}) - 2 Im(c_{lm}) Im(Y_{lm}))).
        // This comes from the symmetry of the diffusion signal, ie
        // f(x) = f(-x), which implies c(lm) = (-1)^m c^\star(l -m).
        // We also have that Y(lm) = (-1)^m Y^\star(l -m).  Those two
        // identities lead to the formula.
        double total = 0.0;
        int i = 0;
        for (int l = 0; l <= maxOrder; l += 2)
            try {
                Complex z = SphericalHarmonics.Y(l, 0, theta, phi);
                total += c[i] * z.real();
                i += 1;
                for (int m = 1; m <= l; m++) {
                    z = SphericalHarmonics.Y(l, m, theta, phi);
                    total += 2 * c[i] * z.real() - 2 * c[i + 1] * z.imag();
                    i += 2;
                }
            }
            catch (Exception e) {
                throw new RuntimeException(e);
            }

        return total;

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


    /**
     * Constructs a string containing the coefficients of the series.
     */
    public String toString() {
        String s = "";
        for (int l = 0; l <= maxOrder; l += 2) {
            for (int m = -l; m <= l; m++) {
                s = s + getCoeff(l, m) + "\n";
            }
            s = s + "\n";
        }

        return s;
    }


    /**
     * Returns a particular coefficient given the spherical harmonic order and
     * index.
     * 
     * @param l
     *            order of the spherical harmonic.
     * 
     * @param m
     *            index of the spherical harmonic.
     * 
     * @return The coefficient c_lm.
     */
    public Complex getCoeff(int l, int m) {

        // Sanity check.
        if (l > maxOrder || Math.abs(m) > l || (l % 2) != 0) {
            return new Complex(0.0, 0.0);
        }

        // First work out the first index of the group of array
        // elements for order l.
        int firstIndex = (l - 1) * ((l - 2) / 2 + 1);

        // For index 0, we are done.
        if (m == 0) {
            return new Complex(c[firstIndex], 0.0);
        }

        // The offset for the specified index.
        int offset = 2 * Math.abs(m) - 1;

        // Construct the complex coefficient from the values in the
        // array.
        Complex coeff = new Complex(c[firstIndex + offset], c[firstIndex + offset + 1]);

        // Negate and conjugate if necessary.
        if (m < 0) {
            if (m % 2 == 0) {
                coeff = coeff.conjugate();
            }
            else {
                coeff = coeff.conjugate().negate();
            }
        }

        return coeff;
    }


    /**
     * We can compute the anisotropy (sqrt of normalized second
     * moment) analytically.  The normalized second moment is \int
     * \int (f(x) - \bar{x})^2 dx/f^2(x) dx.  \int f^2(x) dx is the
     * sum of squares of the coefficients of Y_l0, l = 1, ..., which
     * are real valued.
     *
     * @return The anisotropy
     */
     public double anisotropy() {
         int cind = 0;
         double c00 = c[cind];
         double intfx2 = c00*c00;
         cind += 1;
         for(int i=2; i<=maxOrder; i+=2) {
             intfx2 += c[cind]*c[cind];
             cind += 1;
             for(int j=0; j<2*i; j+=1) {
                 intfx2 += 2.0*c[cind]*c[cind];
                 cind += 1;
             }
         }

         return (intfx2==0)?0.0:Math.sqrt((intfx2 - c00*c00)/intfx2);
     }


    /**
     * Computes the order of the series from the number of values in the array
     * passed to the constructor, which is the number of free parameters in the
     * series.
     * 
     * @param numParams
     *            The number of real values defining the series.
     * 
     * @return The order.
     */
    private int getOrder(int numParams) {
        int cs = 1;
        int order = 0;
        while (cs < numParams) {
            order += 2;
            cs += 2 * order + 1;
        }

        return order;
    }


    public LinearBasisFunction basisFunction(int i) {

        int m = 0;
        int l = 0;
        boolean r = true;
        int counter = 0;

        while(counter < i) {
            counter+=1;
            if(l == 0) {
                l += 2;
            }
            else if(m == l && (!r)) {
                l += 2;
                m = 0;
                r = true;
            }
            else if(m == 0) {
                m = 1;
            }
            else if(r) {
                r = false;
            }
            else {
                r = true;
                m += 1;
            }

        }

        double scalar = 1.0;
        if(m>0)
            scalar = r?2.0:-2.0;

        if(r) {
            return new RealSH(l, m, scalar);
        }

        return new ImagSH(l, m, scalar);
    } 

     /**
     * Returns the parameters used in the basis sum as a string
     *
     * @return The parameters of the basis sum
     */
    public String getSettings(){
	String details = "basis type = sh\norder = " + maxOrder + "\n\n";
	return details;
    }
              
                
}


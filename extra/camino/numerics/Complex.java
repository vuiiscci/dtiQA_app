package numerics;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Complex number class.
 * 
 * <dt>Description:
 * 
 * <dd>Provides methods manipulate complex numbers.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @author Philip Cook
 * @version $Id$
 *  
 */
public class Complex implements Cloneable {

    /**
     * Real part.
     */
    private double real;

    /**
     * Imaginary part.
     */
    private double imaginary;

    /**
     * Constructor from real and imaginary parts.
     * 
     * @param r
     *            Real part.
     * 
     * @param i
     *            Imaginary part.
     */
    public Complex(double r, double i) {
        real = r;
        imaginary = i;
    }

    /**
     * Default constructor creates zero.
     */
    public Complex() {
        real = 0.0;
        imaginary = 0.0;
    }

    public Object clone() {
        return new Complex(real, imaginary);
    }

    /**
     * Returns the real part.
     * 
     * @return Real part.
     */
    public double real() {
        return real;
    }

    /**
     * Returns the imaginary part.
     * 
     * @return Imaginary part.
     */
    public double imag() {
        return imaginary;
    }

    /**
     * Returns the imaginary part.
     * 
     * @return Imaginary part.
     */
    public double imaginary() {
        return imaginary;
    }

    /**
     * Returns the modulus.
     * 
     * @return sqrt(re^2 + im^2).
     */
    public double mod() {
        return Math.sqrt(real * real + imaginary * imaginary);
    }

    /**
     * Returns the magnitude.
     * 
     * @return re^2 + im^2.
     */
    public double mag() {
        return real * real + imaginary * imaginary;
    }

    public boolean equals(Complex c) {
        return ((real == c.real()) && (imaginary == c.imaginary()));
    }

    /**
     * Tests for equality with zero.
     * 
     * @return True if im = re = 0.
     */
    public boolean isZero() {
        return ((real == 0.0) && (imaginary == 0.0));
    }

    public String toString() {
        return "(" + String.valueOf(real) + ") + (" + String.valueOf(imaginary) + ")i";
    }

    /**
     * Adds complex numbers.
     * 
     * @param c
     *            Complex number to add.
     * 
     * @return Result of addition.
     */
    public Complex plus(Complex c) {
        double newr = real + c.real();
        double newi = imaginary + c.imaginary();
        return new Complex(newr, newi);
    }

    /**
     * Subtracts complex numbers.
     * 
     * @param c
     *            Complex number to subtract.
     * 
     * @return Result of subtraction.
     */
    public Complex subtract(Complex c) {
        double newr = real - c.real();
        double newi = imaginary - c.imaginary();
        return new Complex(newr, newi);
    }

    /**
     * Multiplies complex numbers.
     * 
     * @param c
     *            Complex number to multiply by.
     * 
     * @return Result of multiplication.
     */
    public Complex times(Complex c) {
        double newr = real * c.real() - imaginary * c.imaginary();
        double newi = real * c.imaginary() + imaginary * c.real();
        return new Complex(newr, newi);
    }

    /**
     * Multiplies by a scalar.
     * 
     * @param d
     *            Scalar to multiply by.
     * 
     * @return Result of scalar multiplication.
     */
    public Complex times(double d) {
        double newr = real * d;
        double newi = imaginary * d;
        return new Complex(newr, newi);
    }

    /**
     * Complex division.
     * 
     * @param c
     *            Complex number to divide by.
     * 
     * @return Result of division.
     */
    public Complex divide(Complex c) throws ComplexNumberException {
        return times(c.inverse());
    }

    /**
     * Negates the complex number.
     * 
     * @return The negated complex number.
     */
    public Complex negate() {
        double newr = -real;
        double newi = -imaginary;
        return new Complex(newr, newi);
    }

    /**
     * Makes both components positive.
     * 
     * @return (abs(re), abs(im))
     */
    public Complex positive() {
        double newr = Math.abs(real);
        double newi = Math.abs(imaginary);
        return new Complex(newr, newi);
    }

    /**
     * Reciprocates the complex number.
     * 
     * @return The reciprocal of the complex number.
     */
    public Complex inverse() throws ComplexNumberException {
        double denominator = real * real + imaginary * imaginary;

        if (denominator == 0) {
            throw new ComplexNumberException("Cannot reciprocate zero");
        }

        double newr = real / denominator;
        double newi = -imaginary / denominator;

        return new Complex(newr, newi);
    }

    /**
     * Returns the complex conjugate.
     * 
     * @return The complex conjugate.
     */
    public Complex conjugate() {
        return new Complex(real, -imaginary);
    }

    /**
     * Returns the square root.
     * 
     * @return The square root.
     */
    public Complex sqrt() {

        //Compute r, theta.
        double m = mod();
        double theta = 0.0;
        if (m != 0.0) {
            theta = Math.atan(imaginary / real);
        }

        if (real < 0.0) {
            theta += Math.PI;
        }

        double newTheta = theta / 2.0;
        double newMod = Math.sqrt(m);

        double newr = newMod * Math.cos(newTheta);
        double newi = newMod * Math.sin(newTheta);

        return new Complex(newr, newi);
    }

    /**
     * Returns the sine.
     * 
     * @return The sin of the complex number.
     */
    public Complex sin() {
        double newr = Math.sin(real) * cosh(imaginary);
        double newi = sinh(imaginary) * Math.cos(real);

        return new Complex(newr, newi);
    }

    /**
     * Returns the cosine.
     * 
     * @return The cos of the complex number.
     */
    public Complex cos() {
        double newr = Math.cos(real) * cosh(imaginary);
        double newi = -Math.sin(real) * sinh(imaginary);

        return new Complex(newr, newi);
    }

    /**
     * Returns the tangent.
     * 
     * @return The tan of the complex number.
     */
    public Complex tan() throws ComplexNumberException {
        return sin().divide(cos());
    }

    /**
     * @return the hyperbolic sine of a real number.
     *  
     */
    public static double sinh(double a) {
        return 0.5 * (Math.exp(a) - Math.exp(-a));
    }

    /**
     * @return the hyperbolic cosine of a real number.
     *  
     */
    public static double cosh(double a) {
        return 0.5 * (Math.exp(a) + Math.exp(-a));
    }

    /**
     * @return the hyperbolic tangent of a real number.
     *  
     */
    public static double tanh(double a) {
        return (Math.exp(2.0 * a) - 1.0) / (Math.exp(2.0 * a) + 1);
    }

    /**
     * @return the hyperbolic sine of this.
     *  
     */
    public Complex sinh() {
        return new Complex(sinh(real) * Math.cos(imaginary), cosh(real)
                * Math.sin(imaginary));
    }

    /**
     * @return the hyperbolic cosine of this.
     *  
     */
    public Complex cosh() {
        return new Complex(cosh(real) * Math.cos(imaginary), sinh(real)
                * Math.sin(imaginary));
    }

    /**
     * @return the hyperbolic tangent of this.
     *  
     */
    public Complex tanh() {

        double denominator = cosh(2.0 * real) + Math.cos(2.0 * imaginary);

        return new Complex(sinh(2.0 * real) / denominator, Math.sin(2.0 * imaginary)
                / denominator);
    }

}
package numerics;

/**
 * 
 * Contains static method for evaluation of the error function (erf(z)) using
 * the MacClaurin series representation.
 * 
 * 
 * @author Matt Hall m.hall@cs.ucl.ac.uk
 * 
 * 
 *  
 */
public class ErrorFunction {

    private static final double THRESHOLD = 1e-15;

    private static final int ITMAX = 100;

    private static final double EPS = 3.0e-07;

    private static final double FPMIN = 1.0e-30;

    /**
     * Returns the value ln[Gamma(x)] for x > 0.
     * 
     * @param xx
     *            argument of the function
     * @return ln(\Gamma{x})
     */
    private static double gammln(double xx) {
        //Internal arithmetic will be done in double precision, a nicety that
        // you can omit if five-figure
        //accuracy is good enough.
        double x, y, tmp, ser;
        double[] cof = new double[] { 76.18009172947146, -86.50532032941677,
                24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
                -0.5395239384953e-5 };
        int j;

        y = x = xx;
        tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        ser = 1.000000000190015;
        for (j = 0; j <= 5; j++)
            ser += cof[j] / ++y;
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }

    /**
     * Returns the incomplete gamma function Q(a, x) evaluated by its continued
     * fraction representation
     * 
     * @param a
     *            power
     * @param x
     *            argument (limit)
     * 
     * @return Q(a, x) using continued fraction
     */
    private static double gcf(double a, double x) throws ErrorFunctionException {

        int i;
        double an, b, c, d, del, h;
        double gammcf, gln;

        gln = gammln(a);
        b = x + 1.0 - a;
        //Set up for evaluating continued fraction
        //by modified Lentz's method (section 5.2) with b0 = 0.
        c = 1.0 / FPMIN;
        d = 1.0 / b;
        h = d;
        for (i = 1; i <= ITMAX; i++) { //Iterate to convergence.
            an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (Math.abs(d) < FPMIN) {
                d = FPMIN;
            }

            c = b + an / c;

            if (Math.abs(c) < FPMIN) {
                c = FPMIN;
            }
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) < EPS) {
                break;
            }
        }
        if (i > ITMAX) {
            throw new ErrorFunctionException("a too large, ITMAX too small in gcf");
        }
        gammcf = Math.exp(-x + a * Math.log(x) - (gln)) * h; // Put factors in
        // front.
        return gammcf;
    }

    /**
     * evaluates the series representation of the incomplete gamma
     * 
     * 
     * @param a
     *            power
     * @param x
     *            limit
     * @return series rep of lower incomplete gamma
     */
    private static double gser(double a, double x) throws ErrorFunctionException {

        int n;
        double sum, del, ap;

        double gamser, gln;

        gln = gammln(a);

        if (x <= 0.0) {
            if (x < 0.0) {
                throw new ErrorFunctionException("x less than 0 in routine gser");
            }

            return 0.0;
        }
        else {
            ap = a;
            del = sum = 1.0 / a;
            for (n = 1; n <= ITMAX; n++) {
                ++ap;
                del *= x / ap;
                sum += del;
                if (Math.abs(del) < Math.abs(sum) * EPS) {
                    gamser = sum * Math.exp(-x + a * Math.log(x) - (gln));
                    return gamser;
                }
            }
            throw new ErrorFunctionException(
                    "a too large, ITMAX too small in routine gser");
        }
    }

    /**
     * Returns the (lower) incomplete gamma function P(a, x).
     * 
     * @param a
     *            exponent
     * @param x
     *            argument
     * @return incomplete gamma function
     */
    public static double gammp(double a, double x) throws ErrorFunctionException {

        double gamser, gammcf, gln;

        if (x < 0.0 || a <= 0.0) {
            throw new ErrorFunctionException("Invalid arguments in routine gammp");
        }

        if (x < (a + 1.0)) {
            // Use the series representation.
            try {
                gamser = gser(a, x);
            }
            catch (ErrorFunctionException e) {
                throw e;
            }
            return gamser;
        }
        else {
            //Use the continued fraction representation
            try {
                gammcf = gcf(a, x);
            }
            catch (ErrorFunctionException e) {
                throw e;
            }
            return 1.0 - gammcf; //and take its complement.
        }
    }

    /**
     * evaluates the error function for a real argumant using the MacClaurin
     * series representation. The summation continues until the ration of the
     * first and most recent term falls below a predefined threshold reflecting
     * double precision accuracy
     * 
     * @param x
     *            upper limit on the integral
     * 
     * @return erf(x) = \frac{2}{\pi}\int_0^x e^{-x^2} dx
     */
    public static double erf(double x) throws ErrorFunctionException {
        return x < 0.0 ? -gammp(0.5, x * x) : gammp(0.5, x * x);
    }

    /**
     * evaluates the error function over the interval [-4.0, 4.0] and outputs
     * the results to stderr
     * 
     * @param args
     *            command line args
     */
    public static void main(String[] args) {

        System.err.println("testing error function evaluation");

        for (int i = -20; i <= 20; i++) {
            double x = (double) i / 5.0;

            try {
                System.err.println(x + " " + erf(x));
            }
            catch (ErrorFunctionException e) {
                throw new RuntimeException(e);
            }
        }

    }
}

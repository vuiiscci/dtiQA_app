package numerics;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Contains implementation of the incomplete beta function.
 * 
 * <dt>Description:
 * 
 * <dd>Adapted from NRC code.
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id: IncompleteBeta.java,v 1.3 2005/08/18 11:12:22
 *         ucacmgh Exp $
 *  
 */
public class IncompleteBeta {

    /**
     * Computes the incomplete beta function I_x(a, b). See NRC p. 227.
     * 
     * @param a
     *            Argument a of incomplete beta.
     * 
     * @param b
     *            Argument b of incomplete beta.
     * 
     * @param x
     *            Argument x of incomplete beta.
     * 
     * @return I_x(a, b)
     */
    public static double betai(double a, double b, double x)
            throws IncompleteBetaException {
        double bt;

        if (x < 0.0 || x > 1.0)
                throw new IncompleteBetaException("Invalid x in IncompleteBeta.betai."
                        + x + " " + a + " " + b);
        if (x == 0.0 || x == 1.0)
            bt = 0.0;
        else
            bt = Math.exp(GammaFunctions.gammln(a + b) - GammaFunctions.gammln(a)
                    - GammaFunctions.gammln(b) + a * Math.log(x) + b * Math.log(1.0 - x));
        if (x < (a + 1.0) / (a + b + 2.0))
            return bt * betacf(a, b, x) / a;
        else
            return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
    }

    /**
     * Evaluates the continued fraction for the incomplete beta function by
     * modified Lentz's method. See NRC p. 227.
     * 
     * @param a
     *            Argument a of incomplete beta.
     * 
     * @param b
     *            Argument b of incomplete beta.
     * 
     * @param x
     *            Argument x of incomplete beta.
     * 
     * @return I_x(a, b)
     */
    public static double betacf(double a, double b, double x)
            throws IncompleteBetaException {
        int MAXIT = 100;
        double EPS = 3.0e-7;
        double FPMIN = 1.0e-30;

        int m, m2;
        double aa, c, d, del, h, qab, qam, qap;

        qab = a + b;
        qap = a + 1.0;
        qam = a - 1.0;
        c = 1.0;
        d = 1.0 - qab * x / qap;
        if (Math.abs(d) < FPMIN) d = FPMIN;
        d = 1.0 / d;
        h = d;
        for (m = 1; m <= MAXIT; m++) {
            m2 = 2 * m;
            aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < FPMIN) d = FPMIN;
            c = 1.0 + aa / c;
            if (Math.abs(c) < FPMIN) c = FPMIN;
            d = 1.0 / d;
            h *= d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1.0 + aa * d;
            if (Math.abs(d) < FPMIN) d = FPMIN;
            c = 1.0 + aa / c;
            if (Math.abs(c) < FPMIN) c = FPMIN;
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) < EPS) break;
        }
        if (m > MAXIT)
                throw new IncompleteBetaException(
                        "a or b too big, or MAXIT too small in betacf");
        return h;
    }
}
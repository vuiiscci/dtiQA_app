package numerics;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General Levenburg Marquardt Non-linear least squares fitter.
 * 
 * <dt>Description:
 * 
 * <dd>This is an adaptation of the Levenburg-Marquardt code in Numerical
 * Recipes in C - Press, et al. It uses Ridder's method of polynomial
 * extrapolation.
 * 
 * </dl>
 * 
 * $Id$
 * 
 * @author Danny Alexander
 *  
 */
abstract public class NumDeriv {

    private static float CON = 1.4f;

    private static float CON2 = (CON * CON);

    private static float BIG = 1.0e30f;

    private static int NTAB = 10;

    private static float SAFE = 2.0f;

    /**
     * The function to differentiate.
     * 
     * @param arg
     *            The function argument.
     * 
     * @return The value of the function at arg.
     */
    abstract protected float func(float arg);

    /**
     * Computes the derivative at x. See NRC page 188.
     * 
     * @param x
     *            The point at which to compute the derivative.
     * 
     * @param h
     *            Step size.
     * 
     * @param err
     *            Estimated error in the numerical derivative.
     * 
     * @return The value of the derivative.
     */
    public float dfridr(float x, float h, float[] err) {
        int i = 0, j = 0;
        float errt = 0.0f, fac = 0.0f, hh = 0.0f, ans = 0.0f;

        if (h == 0.0f) {
            throw new RuntimeException("h must be nonzero in dfridr.");
        }

        float[][] a = new float[NTAB + 1][NTAB + 1];
        hh = h;
        a[1][1] = (func(x + hh) - func(x - hh)) / (2.0f * hh);
        err[0] = BIG;
        for (i = 2; i <= NTAB; i++) {
            hh /= CON;
            a[1][i] = (func(x + hh) - func(x - hh)) / (2.0f * hh);
            fac = CON2;
            for (j = 2; j <= i; j++) {
                a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0f);
                fac = CON2 * fac;
                errt = Math.max(Math.abs(a[j][i] - a[j - 1][i]), Math.abs(a[j][i]
                        - a[j - 1][i - 1]));
                if (errt <= err[0]) {
                    err[0] = errt;
                    ans = a[j][i];
                }
            }
            if (Math.abs(a[i][i] - a[i - 1][i - 1]) >= SAFE * (err[0])) break;
        }

        return ans;
    }

}
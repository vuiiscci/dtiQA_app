package optimizers;

import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General Conjugate Gradients minimisation algorithm.
 * 
 * <dt>Description:
 * 
 * <dd>This is an adaptation of the conjugate gradients code in Numerical
 * Recipes in C - Press, et al. It is an abstract class that implements the
 * algorithm for minimisation of an unspecified objective function. The
 * objective function and its derivatives must be implemented in subclasses.
 * 
 * </dl>
 * 
 * @version $Date$
 * @author Danny Alexander
 *  
 */
abstract public class ConjGradMinimizer {

    /**
     * The number of parameters in the optimization.
     */
    protected int n;

    /**
     * Counter for the number of iterations in the optimization.
     */
    protected int iter;

    /**
     * Stored the current minimum value.
     */
    protected double fret;

    // Constants used in frprmn.
    private final int ITMAX = 500;

    private final double EPS = 1.0e-10;

    // Constants used in linmin.
    private final double TOL = 2.0e-4;

    // Constants used in mnbrak.
    private final double GOLD = 1.618034;

    private final double GLIMIT = 100.0;

    private final double TINY = 1.0e-20;

    // Constants used in brent.
    private final int BRENTITMAX = 100;

    private final double CGOLD = 0.3819660;

    private final double ZEPS = 1.0e-10;

    /**
     * Variables storing the output of mnbrak and brent.
     */
    private double ax, bx, cx, xmin, fa, fb, fc;

    //Methods to be implemented for particular models.

    /**
     * Returns the value of the objective function with parameters atry.
     * 
     * @param atry
     *            The point at which to evaluate the objective function
     * 
     * @return The value of the objective function.
     */
    abstract protected double fObj(double[] atry);

    /**
     * Returns the derivative of the objective function with respect to each of
     * the parameters in atry. The default implementation computes the
     * derivatives numerically. This should be replaced by an analytic
     * derivative computation if possible.
     * 
     * @param p
     *            The point at which to evaluate the objective function
     *            derivative.
     * 
     * @return The derivative of the objective function.
     */
    protected double[] dfObj(double[] p) {
        return dfObjNumerical(p);
    }

    /**
     * Compute the derivatives numerically.
     * 
     * @param p
     *            The point at which to compute the derivative.
     * 
     * @return The derivates of the objective function at p.
     */
    protected double[] dfObjNumerical(double[] p) {
        double[] derivs = new double[n + 1];

        final double[] pCopy = new double[p.length];
        for (int k = 1; k < p.length; k++) {
            pCopy[k] = p[k];
        }

        for (int j = 1; j <= n; j++) {

            final int index = j;

            NumDeriv df = new NumDeriv() {
                protected float func(float arg) {

                    pCopy[index] = arg;
                    float t = (float) fObj(pCopy);
                    return t;
                }
            };

            float[] err = new float[1];
            float h = (p[j] != 0.0) ? 0.1f * (float) p[j] : 0.001f;
            derivs[j] = df.dfridr((float) p[j], h, err);

            //Replace the original value in a[index] ready to compute
            //the next derivative.
            pCopy[index] = p[index];
        }

        return derivs;
    }

    /**
     * Runs the minimization.
     * 
     * @param p
     *            The starting point for the minimization. Also contains the
     *            parameters at the located minimum after optimization.
     */
    public void minimise(double[] p, double ftol) throws ConjGradMinimizerException {
        iter = 0;
        frprmn(p, ftol);
    }

    /**
     * Initialises working arrays.
     * 
     * @param noParams
     *            The number of parameters to fit.
     */
    protected void init(int noParams) {
        n = noParams;
    }

    /**
     * Translated from NRC routine Pg 423.
     * 
     * @param p
     *            The array containing the parameter values.
     */
    protected void frprmn(double[] p, double ftol) throws ConjGradMinimizerException {

        double[] g = new double[n + 1];
        double[] h = new double[n + 1];
        double fp = fObj(p);
        double[] xi = dfObj(p);
        for (int j = 1; j <= n; j++) {
            g[j] = -xi[j];
            xi[j] = h[j] = g[j];
        }
        for (int its = 1; its <= ITMAX; its++) {
            iter = its;
            linmin(p, xi);
            if (2.0 * Math.abs(fret - fp) <= ftol * (Math.abs(fret) + Math.abs(fp) + EPS)) {
                return;
            }
            fp = fObj(p);
            xi = dfObj(p);
            double dgg = 0.0;
            double gg = 0.0;
            for (int j = 1; j <= n; j++) {
                gg += g[j] * g[j];
                dgg += (xi[j] + g[j]) * xi[j];
            }
            if (gg == 0.0) {
                return;
            }
            double gam = dgg / gg;
            for (int j = 1; j <= n; j++) {
                g[j] = -xi[j];
                xi[j] = h[j] = g[j] + gam * h[j];
            }
        }
        throw new ConjGradMinimizerException("Too many iterations in frprmn");
    }

    /**
     * Does line minimization along xi from p.
     * 
     * @param p
     *            The starting point for the minimization.
     * 
     * @param xi
     *            The direction of the line.
     */
    void linmin(double[] p, double[] xi) throws ConjGradMinimizerException {
        ax = 0.0;
        bx = 1.0;
        mnbrak(p, xi);
        fret = brent(p, xi, TOL);
        for (int j = 1; j <= n; j++) {
            xi[j] *= xmin;
            p[j] += xi[j];
        }
    }

    /**
     * Translated from NRC page 400.
     * 
     * @param p
     *            The current parameter settings.
     * 
     * @param xi
     *            The direction of the line along which bracketing is required.
     */
    protected void mnbrak(double[] p, double[] xi) {
        double ulim, u, r, q, fu;

        fa = f1dim(ax, p, xi);
        fb = f1dim(bx, p, xi);
        if (fb > fa) {
            double dum = ax;
            ax = bx;
            bx = dum;
            dum = fb;
            fb = fa;
            fa = dum;
        }
        cx = bx + GOLD * (bx - ax);
        fc = f1dim(cx, p, xi);
        while (fb > fc) {
            r = (bx - ax) * (fb - fc);
            q = (bx - cx) * (fb - fa);
            u = (bx) - ((bx - cx) * q - (bx - ax) * r)
                    / (2.0 * SIGN(Math.max(Math.abs(q - r), TINY), q - r));
            ulim = bx + GLIMIT * (cx - bx);
            if ((bx - u) * (u - cx) > 0.0) {
                fu = f1dim(u, p, xi);
                if (fu < fc) {
                    ax = bx;
                    bx = u;
                    fa = fb;
                    fb = fu;
                    return;
                }
                else if (fu > fb) {
                    cx = u;
                    fc = fu;
                    return;
                }
                u = cx + GOLD * (cx - bx);
                fu = f1dim(u, p, xi);
            }
            else if ((cx - u) * (u - ulim) > 0.0) {
                fu = f1dim(u, p, xi);
                if (fu < fc) {
                    bx = cx;
                    cx = u;
                    u = cx + GOLD * (cx - bx);

                    fb = fc;
                    fc = fu;
                    fu = f1dim(u, p, xi);
                }
            }
            else if ((u - ulim) * (ulim - cx) >= 0.0) {
                u = ulim;
                fu = f1dim(u, p, xi);
            }
            else {
                u = cx + GOLD * (cx - bx);
                fu = f1dim(u, p, xi);
            }
            ax = bx;
            bx = cx;
            cx = u;

            fa = fb;
            fb = fc;
            fc = fu;
        }
    }

    /**
     * Translated from NRC page 404.
     * 
     * @param p
     *            The current parameter settings.
     * 
     * @param xi
     *            The direction of the line along which the minimum is required.
     */
    protected double brent(double[] p, double[] xi, double tol)
            throws ConjGradMinimizerException {
        double a, b, etemp, fu, fv, fw, fx, p1, q1, r1, tol1, tol2, u, v, w, x, xm;
        double e = 0.0;
        double d = 0.0;

        a = (ax < cx ? ax : cx);
        b = (ax > cx ? ax : cx);
        x = w = v = bx;
        fw = fv = fx = f1dim(x, p, xi);
        for (int iter = 1; iter <= BRENTITMAX; iter++) {
            xm = 0.5 * (a + b);
            tol2 = 2.0 * (tol1 = tol * Math.abs(x) + ZEPS);
            if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                xmin = x;
                return fx;
            }
            if (Math.abs(e) > tol1) {
                r1 = (x - w) * (fx - fv);
                q1 = (x - v) * (fx - fw);
                p1 = (x - v) * q1 - (x - w) * r1;
                q1 = 2.0 * (q1 - r1);
                if (q1 > 0.0) p1 = -p1;
                q1 = Math.abs(q1);
                etemp = e;
                e = d;
                if (Math.abs(p1) >= Math.abs(0.5 * q1 * etemp) || p1 <= q1 * (a - x)
                        || p1 >= q1 * (b - x))
                    d = CGOLD * (e = (x >= xm ? a - x : b - x));
                else {
                    d = p1 / q1;
                    u = x + d;
                    if (u - a < tol2 || b - u < tol2) d = SIGN(tol1, xm - x);
                }
            }
            else {
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            }
            u = (Math.abs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
            fu = f1dim(u, p, xi);
            if (fu <= fx) {
                if (u >= x)
                    a = x;
                else
                    b = x;
                v = w;
                w = x;
                x = u;

                fv = fw;
                fw = fx;
                fx = fu;
            }
            else {
                if (u < x)
                    a = u;
                else
                    b = u;
                if (fu <= fw || w == x) {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                }
                else if (fu <= fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                }
            }
        }
        throw new ConjGradMinimizerException("Too many iterations in brent");
    }

    /**
     * Computes the value of the objective function a distance x along the line
     * in direction xi from p.
     */
    private double f1dim(double x, double[] p, double[] xi) {
        double[] xt = new double[n + 1];
        for (int j = 1; j <= n; j++)
            xt[j] = p[j] + x * xi[j];
        double f = fObj(xt);
        return f;
    }

    /**
     * Translation of the nrutil function SIGN.
     */
    public static double SIGN(double a, double b) {
        return ((b) >= 0.0 ? Math.abs(a) : -Math.abs(a));
    }

}

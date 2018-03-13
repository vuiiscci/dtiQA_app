package optimizers;

import misc.LoggedException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Simple test example for conjugate gradients minimiser class.
 * 
 * <dt>Description:
 * 
 * <dd>Uses the conjugate gradients minimisation to find the minimum point of
 * f(x, y) = (3x-1)^2 + (y+1)^2 + 3.
 * 
 * </dl>
 * 
 * @version $Id: ConjGradMinimizerTest.java,v 1.5 2005/08/18 11:13:57 ucacmgh
 *          Exp $
 * @author Danny Alexander
 *  
 */
public class ConjGradMinimizerTest extends ConjGradMinimizer {

    /**
     * Basic constructor.
     */
    public ConjGradMinimizerTest() {
        init(2);
    }

    /**
     * Implements the function (3x-1)^2 + (y+1)^2 + 3.
     */
    protected double fObj(double[] params) {

        //Note that the params array starts counting from 1.

        return (3.0 * params[1] - 1.0) * (3.0 * params[1] - 1.0) + (params[2] + 1.0)
                * (params[2] + 1.0) + 3.0;
    }

    /**
     * Implements the analytic derivates.
     * 
     * @param params
     *            The values of the minimization parameters.
     * 
     * @return The derivatives.
     */
    /*
     * protected double[] dfObj(double[] params) {
     * 
     * double[] dfda = new double[3];
     * 
     * //Note that the derivatives array must start counting from 1 like //the
     * params array. dfda[1] = 6.0*(3.0*params[1]-1.0); dfda[2] =
     * 2.0*(params[2]+1.0);
     * 
     * return dfda; }
     */

    /**
     * Creates an object and calls minimisation routine.
     */
    public static void main(String[] args) {
        try {
            ConjGradMinimizerTest m = new ConjGradMinimizerTest();

            //Set the initial values for the parameters.
            double[] initPs = { 0.0, 2.0, 3.0 };

            //This is the tolerance in the minimization.
            double ftol = 1.0E-5;

            //Do the actual minimisation.
            m.minimise(initPs, ftol);

            //Get the parameter values at the minimum point and
            //print them out.
            for (int i = 0; i < initPs.length; i++) {
                System.err.print(initPs[i] + " ");
            }
            System.err.println();
        }
        catch (Exception e) {
            LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
        }
    }

}
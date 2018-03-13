package optimizers;

import misc.LoggedException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Simple test example for Marquardt minimiser class.
 * 
 * <dt>Description:
 * 
 * <dd>Uses the Marquardt minimisation algorithm to find the minimum point of
 * f(x, y) = (3x-1)^2 + (y+1)^2 + 3.
 * 
 * </dl>
 * 
 * @version $Id: MarquardtMinimiserTest.java,v 1.4 2005/08/18 11:13:57 ucacmgh
 *          Exp $
 * @author Danny Alexander
 *  
 */
public class MarquardtMinimiserTest extends MarquardtMinimiser {

    /**
     * Basic constructor.
     */
    public MarquardtMinimiserTest() {
        init(2);
    }

    /**
     * Implements the function (3x-1)^2 + (y+1)^2 + 3.
     */
    protected double fObj(double[] params, double[] dfda, double[][] d2fda2) {

        //Note that the params array starts counting from 1.

        //Note that the derivatives array must start counting from 1 like
        //the params array.
        dfda[1] = 6.0 * (3.0 * params[1] - 1.0);
        dfda[2] = 2.0 * (params[2] + 1.0);

        //IMPORTANT: Note that it is often not necessary to implement
        //the second derivatives exactly. Often the matrix of products of
        //first derivatives is used instead - d2fda2[i][j] = dfda[i]dfda[j]
        //but this must be explicitly implemented here.

        //Note that the derivatives array must start counting from 1 like
        //the params array.
        d2fda2[1][1] = 18.0;
        d2fda2[2][2] = 2.0;

        return (3.0 * params[1] - 1.0) * (3.0 * params[1] - 1.0) + (params[2] + 1.0)
                * (params[2] + 1.0) + 3.0;
    }

    /**
     * Creates an object and calls minimisation routine.
     */
    public static void main(String[] args) {
        try {
            MarquardtMinimiserTest m = new MarquardtMinimiserTest();

            //Set the initial values for the parameters.
            //Note initPs starts counting from zero, unlike arrays in the
            //Marquardt minimiser code - see fObj and dfdas above.
            double[] initPs = { 2.0, 3.0 };
            m.setInitParams(initPs);

            //Do the actual minimisation.
            m.minimise();

            //Get the parameter values at the minimum point and
            //print them out.
            double[] finalps = m.getParameters();
            for (int i = 0; i < finalps.length; i++) {
                System.err.print(finalps[i] + " ");
            }
            System.err.println();
        }
        catch (Exception e) {
            LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
        }
    }

}
package inverters;

import misc.*;
import optimizers.*;
import tools.*;

/**
 * Optimizes the Cholesky decomposition of the tensor, similar to DiffTensorFitter, but
 * uses the Minpack optimizer rather than optimizers.MarquardtMinimiser. 
 *
 * @author Philip Cook
 * 
 */
public class DiffTensorMinPackFitter extends DiffTensorFitter implements Lmdif_fcn {



    /**
     * The number of unweighted acquisitions that are made (nob0s) is required to
     * estimate the noise levels of each data item.
     * 
     * @param indepVals
     *            The matrix of non-zero gradient directions without the zero measurements.
     * 
     * @param depVals
     *            The normalized measurements.
     * 
     * @param bValues
     *            The array of b-values.
     * 
     * @param nob0s
     *            The number of b=0 acquisitions.
     */
    public DiffTensorMinPackFitter(double[][] indepVals, double[] depVals, double[] bValues,
            int nob0s) throws MarquardtMinimiserException {
        super(indepVals, depVals, bValues, nob0s);
    }


    /**
     * Objective function called by Minpack.
     *
     * @param m A positive integer set to the number of functions [number of observations]
     * @param n A positive integer set to the number of variables [number of parameters]. n must not exceed m.
     * @param x On input, it contains the initial estimate of the solution vector [the least squares parameters]. On output it contains 
     *          the final estimate of the solution vector.
     * @param fvec - An output vector that contains the m functions [residuals] evaluated at x.
     * @param iflag - Array indexed from 1. If iflag[1] == 1, then update the parameters.
     *
     */
    public void fcn(int m, int n, double[] x, double[] fvec, int[] iflag) {
        
        System.arraycopy(x, 0, a, 0, 7);

        // could put some weights in here to do weighted NLLS

        for (int i = 1; i <= ndata; i++) {
            double yf = yfit(a, i);
            fvec[i] = (yf - y[i]); 
        }
        
    }

    
    public void minimise() throws MarquardtMinimiserException {

        double[] fVec = new double[ndata+1];
        int[] info = new int[2];
     
        double[] x = new double[7];

        System.arraycopy(a, 0, x, 0, 7);

	// Use minpack minimizer with numerical derivatives
        Minpack_f77.lmdif1_f77(this, ndata, noParams, x, fVec, 1E-8, info);

        // copy final parameters back into this object
        System.arraycopy(x, 0, a, 0, 7);

        if (info[1] < 1 || info[1] > 4) {
            System.err.println("Minpack failure, code " + info[1]);
            System.err.println("Final DT: \n" + getDT_Chol(a, 1));

            // Might be better to return linear DT here?
        }
    }

}

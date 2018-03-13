package optimizers;

/**
 * Interface for objective functions passed to Minpack_f77.
 *
 * @author Steve Verrill, Philip Cook (docs)
 *
 * @see optimizers.Minpack_f77
 */
public interface  Lmdif_fcn {

    /**
     * Note: all arrays index from 1
     * 
     * @param m - A positive integer set to the number of functions [number of observations]
     * @param n - A positive integer set to the number of variables [number of parameters]. n must not exceed m.
     * @param x - On input, it contains the initial estimate of the solution vector [the least squares parameters]. On output it contains the final estimate of the solution vector.
     * @param fvec - An output vector that contains the m residuals of the function evaluated at x.
     *
     * @param iflag - iflag[1] is positive to continue, set iflag[1] to negative value to halt optimization
     */
    public void fcn(int m, int n, double[] x, double[] fvec, int[] iflag);

}

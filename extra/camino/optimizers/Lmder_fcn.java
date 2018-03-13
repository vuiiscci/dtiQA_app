package optimizers;

/**
 * Interface for Minpack LM functions. Implement this interface if you can compute both the objective function
 * and the jacobian.
 *
 * Note all arrays index from 1.
 * 
 * @author Steve Verrill, Philip Cook (docs)
 *
 * @see optimizers.Minpack_f77
 */
public interface Lmder_fcn {


    /**
     * Note: all arrays index from 1
     * 
     * @param m - A positive integer set to the number of functions [number of observations]
     * @param n - A positive integer set to the number of variables [number of parameters]. n must not exceed m.
     * @param x - On input, it contains the initial estimate of the solution vector [the least squares parameters]. 
     * On output it contains the final estimate of the solution vector.
     *
     * @param fvec - An output vector that contains the m residuals of the function evaluated at x. Set this only
     * if iflag[1] == 1. If iflag[1] == 2, then you should update fjac and not fvec.
     *
     * @param fjac - An output m by n array. The upper n by n submatrix of fjac contains an upper triangular matrix 
     * R with diagonal elements of nonincreasing magnitude such that P^T (jac^t jac)P = R^T R, where P is a permutation 
     * matrix and jac is the final calculated Jacobian. Column j of P is column ipvt[j] of the 
     * identity matrix. The lower trapezoidal part of fjac contains information generated during the computation of R.
     *
     *
     * @param iflag - A flag that is 1 if we are to compute fvec, and 2 if we are to compute fjac. Set iflag[1] to a 
     * negative value to halt optimization, otherwise don't change it.
     *
     *
     */
    public void fcn(int m, int n, double x[], double fvec[], double fjac[][], int iflag[]);

}

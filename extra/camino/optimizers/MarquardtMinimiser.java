package optimizers;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General Levenburg Marquardt minimisation algorithm.
 * 
 * <dt>Description:
 * 
 * <dd>This is an adaptation of the Levenburg-Marquardt code in Numerical
 * Recipes in C - Press, et al. It is an abstract class that implements the
 * algorithm for minimisation of an unspecified objective function. The
 * objective function and its derivatives must be implemented in subclasses.
 * 
 * </dl>
 * 
 * @version $Id$
 * @author Danny Alexander
 *  
 */
abstract public class MarquardtMinimiser {

	//Parameters of Levenberg-Marquardt algorithm
	protected double alamda;

	protected double fObjVal;

	protected double ofObjVal;

	//Threshold on fObj change at which to declare convergence.
	private static double CONVERGETHRESH = 1.0E-8;

	//Maximum number of iterations for LM fitting.
	private static int MAXITER = -1;
	
	//Parameters of the model

	//Number of parameters
	protected int ma;

	//Boolean values indicating whether each parameter is fixed or variable.
	protected boolean[] ia;

	//Parameter values.
	protected double[] a;

	//Method to get/set new convergence threshold
	public static void setCONVERGETHRESH(double convthresh) {
		CONVERGETHRESH = convthresh;                         
	}
	
	//Method to get/set new MAXITER
	public static void setMAXITER(int maxiter) {
		MAXITER = maxiter;                         
	}

	public static double getCONVERGETHRESH() {
		return CONVERGETHRESH;
	}

	//Methods to be implemented for particular models.

	/**
	 * Returns the value of the objective function with parameters atry. Array
	 * <code>dfda</code> is filled with values of the first derivative of fObj at atry wrt
	 * each of the parameters and <code>d2fda2</code> is filled with the second derivatives
	 * or an approximation thereof (such as dfda.dfda^T).
	 * 
	 * @param atry
	 *            The point at which to evaluate the objective function
	 * 
	 * @param dfda
	 *            The array to be filled with the first derivatives.
	 * 
	 * @param d2fda2
	 *            The array to be filled with the second derivatives.
	 * 
	 * @return The value of the objective function.
	 */
	abstract protected double fObj(double[] atry, double[] dfda,
			double[][] d2fda2);

	/**
	 * Initialises working arrays.
	 * 
	 * @param noParams
	 *            The number of parameters to fit.
	 */
	protected void init(int noParams) {
		ma = noParams;
		ia = new boolean[ma + 1];
		a = new double[ma + 1];
		for (int i = 1; i <= ma; i++) {
			ia[i] = true;
			a[i] = 0.0;
		}

	}

	/**
	 * Sets the initial values of the parameters. The values in aInit start
	 * counting from zero.
	 * 
	 * @param aInit
	 *            Array containing the new parameter values starting from index
	 *            0.
	 */
	public void setInitParams(double[] aInit)
	throws MarquardtMinimiserException {
		if (aInit.length != a.length - 1) {
			throw new MarquardtMinimiserException(
					"Wrong number of parameters in initializing array.  Got "
					+ aInit.length + " expected " + (a.length - 1)
					+ ".");
		}
		for (int i = 0; i < aInit.length; i++) {
			a[i + 1] = aInit[i];
		}
	}

	/**
	 * Runs the minimization.
	 */
	public void minimise() throws MarquardtMinimiserException {

		double[][] covar = new double[ma + 1][ma + 1];
		double[][] alpha = new double[ma + 1][ma + 1];
		double[] atry = new double[ma + 1];
		double[] beta = new double[ma + 1];
		double[] da = new double[ma + 1];
		int mfit = 0;
		for (int j = 1; j <= ma; j++)
			if (ia[j])
				mfit++;
		
		
		double[][] oneda = new double[mfit + 1][2];
		alamda = 0.001;
		mrqcof(a, alpha, beta);
		ofObjVal = fObjVal;
		for (int j = 1; j <= ma; j++)
			atry[j] = a[j];
		
		
		// begin test output for intermediate steps
		double[] dfdaa = new double[ma + 1];
		double[][] d2fdaa2 = new double[ma + 1][ma + 1];

		
		boolean done = false;
		int iter = 0;
		
		while (!done) {
			double diff = mrqmin(covar, alpha, atry, beta, da, oneda, mfit);

			//Check on alambda, if it blows up massively, we
			//might as well give up.
			double maxALambda = 1.0E30;
			if (alamda > maxALambda) {
				throw new MarquardtMinimiserException("alamda exceeds "
						+ maxALambda + ". Giving up.");
			}
			
			
			// begin test output for intermediate steps
			double[] dfda = new double[ma + 1];
			double[][] d2fda2 = new double[ma + 1][ma + 1];

			
			
            iter++;
			done = (diff <= 0.0) && (diff > -CONVERGETHRESH) || ((iter >= MAXITER) && (MAXITER > 0) );
			
		}

		double[] dfda = new double[ma + 1];
		double[][] d2fda2 = new double[ma + 1][ma + 1];
		double fmin = fObj(a, dfda, d2fda2);

                if((iter >= MAXITER) && (MAXITER > 0) ) {
                    throw new MarquardtMinimiserNonConvergenceException("Exceeded maximum number of iterations: " + MAXITER);
                }

	}

	/**
	 * Returns the values of the parameters.
	 * 
	 * @return The values of the parameters.
	 */
	public double[] getParameters() {
		double[] a2 = new double[a.length];
		for (int i = 0; i < a2.length; i++) {
			a2[i] = a[i];
		}
		return a2;
	}

	/**
	 * Returns the value of the objective function with the current parameter
	 * settings.
	 * 
	 * @return The value of the objective function.
	 */
	public double getFObjVal() {
		return fObjVal;
	}

	/**
	 * Sets the convergence threshold, default is 10^-8.
	 */
	public void setConvergence(double c) {
		setCONVERGETHRESH(c);
	}

	/**
	 * NRC routine for performing the Marquardt method.
	 */
	private double mrqmin(double[][] covar, double[][] alpha, double[] atry,
			double[] beta, double[] da, double[][] oneda, int mfit)
	throws MarquardtMinimiserException {
		for (int j = 1; j <= mfit; j++) {
			for (int k = 1; k <= mfit; k++)
				covar[j][k] = alpha[j][k];
			covar[j][j] = alpha[j][j] * (1.0 + alamda);
			oneda[j][1] = beta[j];
		}

		gaussj(covar, mfit, oneda, 1);

		for (int j = 1; j <= mfit; j++)
			da[j] = oneda[j][1];
		if (alamda == 0.0) {
			covsrt(covar, ma, ia, mfit);
			return 0.0;
		}
		for (int j = 0, l = 1; l <= ma; l++)
			if (ia[l])
				atry[l] = a[l] + da[++j];
		mrqcof(atry, covar, da);
		double diff = fObjVal - ofObjVal;
		if (fObjVal < ofObjVal) {
			alamda *= 0.1;
			ofObjVal = fObjVal;
			for (int j = 1; j <= mfit; j++) {
				for (int k = 1; k <= mfit; k++)
					alpha[j][k] = covar[j][k];
				beta[j] = da[j];
			}
			for (int l = 1; l <= ma; l++)
				a[l] = atry[l];
		} else {
			alamda *= 10.0;
			fObjVal = ofObjVal;
		}
		return diff;
	}

	/**
	 * See NRC
	 */
	protected void mrqcof(double[] atry, double[][] alpha, double[] beta)
	throws MarquardtMinimiserException {
		int mfit = 0;
		for (int j = 1; j <= ma; j++)
			if (ia[j])
				mfit++;

		//Get value of fObj and its first and second derivatives.
		double[] dfda = new double[ma + 1];
		double[][] d2fda2 = new double[ma + 1][ma + 1];
		fObjVal = fObj(atry, dfda, d2fda2);

		/**
		 * For testing.
		for(int i=0; i<dfda.length; i++)
		    System.out.print(dfda[i] + " ");
		System.out.println();
		for(int i=1;i<atry.length; i++) {
		    atry[i] = atry[i] + 0.000000001;
		    double test = fObj(atry, dfda, d2fda2);
		    System.out.println(fObjVal + " " + test + " " + (test-fObjVal)/0.000000001);
		    atry[i] = atry[i] - 0.000000001;
		}
		for(int i=0; i<d2fda2.length; i++) {
		    for(int j=0; j<d2fda2.length; j++)
		        System.out.print(d2fda2[i][j] + " ");
		    System.out.println();
		}
		System.exit(1);
		 */

		//If fObj becomes infinite, give up.
		if (Double.isNaN(fObjVal) || Double.isInfinite(fObjVal)) {
			throw new MarquardtMinimiserException(
			"Infinite objective function in mrqcof.");
		}

		int i = 0;
		for (int l = 1; l <= ma; l++) {
			if (ia[l]) {
				i++;
				int k = 0;
				for (int m = 1; m <= ma; m++)
					if (ia[m])
						alpha[i][++k] = d2fda2[l][m] / 2.0;
				beta[i] = -dfda[l] / 2.0;
			}
		}

	}

	/**
	 * See NRC
	 */
	private static void gaussj(double[][] aMat, int n, double[][] bMat, int m)
	throws MarquardtMinimiserException {
		int[] indxc = new int[n + 1];
		int[] indxr = new int[n + 1];
		int[] ipiv = new int[n + 1];
		int irow = 0;
		int icol = 0;
		for (int j = 1; j <= n; j++)
			ipiv[j] = 0;
		for (int i = 1; i <= n; i++) {
			double big = 0.0;
			for (int j = 1; j <= n; j++)
				if (ipiv[j] != 1)
					for (int k = 1; k <= n; k++) {
						if (ipiv[k] == 0) {
							if (Math.abs(aMat[j][k]) >= big) {
								big = Math.abs(aMat[j][k]);
								irow = j;
								icol = k;
							}
						} else if (ipiv[k] > 1)
							throw new MarquardtMinimiserException(
							"Singular Matrix in gaussj - 1");
					}
			//NRC does this:
			//++(ipiv[icol]);
			//but gcj doesn't understand so replaced here by:
			ipiv[icol] += 1;

			if (irow != icol) {
				for (int l = 1; l <= n; l++) {
					double temp = aMat[irow][l];
					aMat[irow][l] = aMat[icol][l];
					aMat[icol][l] = temp;
				}
				for (int l = 1; l <= m; l++) {
					double temp = bMat[irow][l];
					bMat[irow][l] = bMat[icol][l];
					bMat[icol][l] = temp;
				}
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if (aMat[icol][icol] == 0.0) {
				throw new MarquardtMinimiserException(
				"Singular Matrix in gaussj - 2  ");
			}
			double pivinv = 1.0 / aMat[icol][icol];
			aMat[icol][icol] = 1.0;
			for (int l = 1; l <= n; l++)
				aMat[icol][l] *= pivinv;
			for (int l = 1; l <= m; l++)
				bMat[icol][l] *= pivinv;
			for (int ll = 1; ll <= n; ll++)
				if (ll != icol) {
					double dum = aMat[ll][icol];
					aMat[ll][icol] = 0.0;
					for (int l = 1; l <= n; l++)
						aMat[ll][l] -= aMat[icol][l] * dum;
					for (int l = 1; l <= m; l++)
						bMat[ll][l] -= bMat[icol][l] * dum;
				}
		}
		for (int l = n; l >= 1; l--) {
			if (indxr[l] != indxc[l])
				for (int k = 1; k <= n; k++) {
					double temp = aMat[k][indxr[l]];
					aMat[k][indxr[l]] = aMat[k][indxc[l]];
					aMat[k][indxc[l]] = temp;
				}
		}
	}

	/**
	 * See NRC
	 */
	private static void covsrt(double[][] covar, int ma, boolean[] ia, int mfit) {
		for (int i = mfit + 1; i <= ma; i++)
			for (int j = 1; j <= i; j++)
				covar[i][j] = covar[j][i] = 0.0;
		int k = mfit;
		for (int j = ma; j >= 1; j--) {
			if (ia[j]) {
				for (int i = 1; i <= ma; i++) {
					double temp = covar[i][k];
					covar[i][k] = covar[i][j];
					covar[i][j] = temp;
				}
				for (int i = 1; i <= ma; i++) {
					double temp = covar[k][i];
					covar[k][i] = covar[j][i];
					covar[j][i] = temp;
				}
				k--;
			}
		}
	}

}

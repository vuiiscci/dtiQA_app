    package optimizers;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Interface for objects that act as minimizers for model fitting.
 * 
 * <dt>Description:
 * 
 * <dd>Must be able to accept starting parameter settings and
 * measurements to fit to, as well as have a minimise operation and a
 * function that returns a collection of sets of fitted parameter
 * estimates.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: DataSource.java 58 2006-10-26 11:38:18Z ucacmgh $
 *  
 */
public interface Minimizer {


    /**
     * Sets the initial values of the optimized parameters. The values
     * in aInit start counting from zero.
     * 
     * @param aInit
     *            Array containing the new parameter values starting from index
     *            0.
     */
    public void setInitParams(double[] aInit) throws MinimizerException;


    /**
     * Initializes the fitting procedure with a new set of
     * measurements (dependent variables).
     * 
     * @param newMeas The new set of measurements.
     */
    public void setMeasurements(double[] newMeas) throws MinimizerException;


    /**
     * Runs the minimization.
     */
    public void minimise() throws MinimizerException;


    /**
     * Returns the list of candidate solutions.
     *
     * @return Kx(N+1) matrix in which each of the K rows contains the
     * N model parameters followed by the objective function value.
     */
    public double[][] getSolutions();


    /**
     * Returns the number of candidate solutions each run returns.
     *
     * @return The number K of rows in the matrix returned by getSolutions.
     */
    public int getNumSolutions();
    
    /**
     * Returns the number of parameters the solution contains.
     *
     * @return The number N of columns in the matrix returned by getSolutions.
     */
    public int getNumParameters();
    
    /**
     * Adjusts the gradient for each voxel
     * @param gradAdj Gradient adjustment to be applied to the voxel
     */
    public void adjustScheme(double[] gradAdj);

}

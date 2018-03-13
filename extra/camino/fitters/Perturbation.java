package fitters;

import numerics.*;
import java.util.Random;


/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Base class for starting point perturbations (as in for example
 * MultiRunLM fitters) or proposal distributions required for Bayesian
 * fitting routines such as MCMC.
 * 
 * <dt>Description:
 * 
 * <dd>Requires only one function that perturbs a candidate set of
 * fitting parameters according to some proposal distribution.  Note
 * the perturbations act on the optimized parameters rather than the
 * model parameters (see Codec.java).
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: DataSource.java 58 2006-10-26 11:38:18Z ucacmgh $
 *  
 */
public abstract class Perturbation {


    /**
     * Generate new parameters from proposal distribution.
     * 
     * @param params
     *            The values of the current parameters.
     * 
     * @param paramsInitial The starting values of the parameters,
     * which are sometimes used to determine the step length.
     * 
     * @param rand The random number generator to generate the
     * perturbation.
     * 
     * @return The values of the proposed parameters.
     */				
    public abstract double[] perturb(double[] params, double[] paramsInitial, Random rand);


}




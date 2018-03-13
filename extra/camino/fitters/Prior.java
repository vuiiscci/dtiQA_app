package fitters;

import numerics.*;


/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Base class for prior probability distributions required for
 * Bayesian fitting routines such as MCMC.
 * 
 * <dt>Description:
 * 
 * <dd>Requires only one function that returns the prior probability
 * of a candidate set of fitting parameters.  Note these are the
 * optimized parameters rather than the model parameters (see Codec.java).
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: DataSource.java 58 2006-10-26 11:38:18Z ucacmgh $
 *  
 */
public abstract class Prior {


    /**
     * Returns the prior probability of the specified list of
     * optimization parameters.
     *
     * By default, all the combinations of parameter settings are
     * equally likely.
     * 
     * Note that normalization constants are not required in the prior
     * distributions, because it is only ratios of prior probabilities
     * of different parameter values that are important.  Thus for
     * uniform priors, for example, values can be one if in range and
     * zero if outside regardless of the size of the range.
     *
     * @param array of optimized parameters
     *
     * @return prior probability of parameters.
     */
    public double prior(double[] optParams) {
        return 1;
    }


}




/**
 * 
 */
package fitters;

import java.util.Random;

/**
 *  * <dd> Class for starting point perturbations using a gaussian standard distribution 
 * 
 * <dt>Description:
 *Perturbs the parameters using a gaussian standard distribution which is controlled by a single parameter: the standard
 *deviation (0.1 by default).
 * @author laura
 *
 */

public class FixedSTD_GaussianPerturbation extends Perturbation {

	private double stdGaussianDistribution;
	
	/**
	 * Creates a gaussian perturbation with std 0.1
	 */
	public FixedSTD_GaussianPerturbation()
	{
		setStd(0.1);
	}
	
	/**
	 *  Creates a gaussian perturbation with std
	 * @param std
	 */
	public FixedSTD_GaussianPerturbation(double std)
	{
		setStd(std);
	}
	
	 /**
     * Generate new parameters from gaussian distribution.
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
	public double[] perturb(double[] params, double[] paramsInitial, Random rand) {

        double[] paramsNew = new double[params.length];

        for(int i =0;i<params.length;i++) {
            double r = rand.nextGaussian();

            if (paramsInitial[i] != 0.0)
                paramsNew[i] = params[i] + r * getStd() * paramsInitial[i];
            else 
                paramsNew[i] = params[i] + r * getStd();
        }
        return paramsNew;
	
	}

	/**
	 * Sets std of gaussian distribution.
	 * @param stdGaussianDistribution the stdUniformDistribution to set
	 */
	private void setStd(double stdGaussianDistribution) {
		this.stdGaussianDistribution = stdGaussianDistribution;
	}

	/**
	 * Returns std from the uniform distribution.
	 * @return the stdUniformDistribution
	 */
	private double getStd() {
		return stdGaussianDistribution;
	}

}

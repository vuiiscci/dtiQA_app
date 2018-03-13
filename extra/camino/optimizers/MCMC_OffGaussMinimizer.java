package optimizers;

import java.util.logging.Logger;
import models.*;
import optimizers.*;
import numerics.*;
import misc.*;
import imaging.*;
import fitters.*;
import tools.CL_Initializer;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>This class is used by fitter classes.
 *
 * <dt>Description:
 * 
 * <dd>Markov Chain Monte Carlo algorithm for minimizing chi-squared
 * objective function, with offset Gaussian noise model.
 * 
 * </dl>
 *
 * @author Danny
 * @version $Id$
 *
 */
public class MCMC_OffGaussMinimizer extends MCMC_GaussianMinimizer {
	
    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.optimizers.MCMC_OffGaussMinimizer");
	
    /**
     * Default constructor.
     */
    public MCMC_OffGaussMinimizer() {
    }
	
    /**
     * Constructor needs all the following:
     *
     * @param scheme The acquisition protocol
     *
     * @param pm The diffusion model to fit
     *
     * @param cod The Codec specifying the transformation from model
     *
     * @param p Prior object specifying the prior distribution on the parameters.
     *
     * @param pb Perturbation object specifying the proposal distributions.
     *
     * @param sigma The expected standard deviation of the noise in the data
     * to optimized parameters.
     */
    public MCMC_OffGaussMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod, Prior p, Perturbation pb, double sigma) throws MarkovChainMonteCarloException {
		
	super(scheme, pm, cod, p, pb, sigma);
		
    }
	
	
    /**
     * Implements the chi-squared function.
     * 
     * @param params
     *            The current model parameter settings.
     */
    protected double fObj(double[] optParams) {
		
        //Initialise values and arrays.
        double chisq = 0.0;
		
        double[] modParams = codec.optToModel(optParams);
        RealMatrix modSignals = model.getSignals(modParams, scheme);

	boolean print = false;
        if(print) {
            System.err.println("modParams");
            for(int i=0; i<modParams.length; i++)
                System.err.print(modParams[i] + " ");
            System.err.println();
            System.err.println("measurements");
            for(int i=0; i<measurements.length; i++)
                System.err.print(measurements[i] + " ");
            System.err.println();
            System.err.println("modSignals");
            System.err.println(modSignals);
        }

        for (int i = 1; i <= scheme.numMeasurements(); i++) {
            double sig2 = sig[i-1] * sig[i-1];
            double sig2i = 1.0/sig2;
            double ymod = modSignals.entries[i-1][0];
            ymod = Math.sqrt(ymod*ymod + sig2);
            double dy = measurements[i-1] - ymod;
            chisq += dy * dy * sig2i;
        }
        if(print)
            System.err.println("chisq: " + chisq);
		
        return chisq;
    }
	
	
}


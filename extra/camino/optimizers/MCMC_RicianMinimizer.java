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
 * <dd>Markov Chain Monte Carlo algorithm for minimizing log likelihood
 * objective function, with Rician noise model.
 * 
 * </dl>
 *
 * @author Danny
 * @version $Id$
 *
 */
public class MCMC_RicianMinimizer extends MCMC_GaussianMinimizer {
	
    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.optimizers.MCMC_RicianMinimizer");
	
    /**
     * Default constructor.
     */
    public MCMC_RicianMinimizer() {
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
    public MCMC_RicianMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod, Prior p, Perturbation pb, double sigma) throws MarkovChainMonteCarloException {
		
	super(scheme, pm, cod, p, pb, sigma);
		
    }
	
	
    /**
     * Implements the log likelihood function.
     * 
     * @param params
     *            The current model parameter settings.
     */
    protected double fObj(double[] optParams) {
		
        //Initialise values and arrays.
        double loglik = 0.0;
		
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
            double z = ymod*measurements[i-1]*sig2i;
            
            // Above z threshold, the exact Bessel function
            // implementation fails so we approximate the
            // values we need.
            double logBessI0;
            double bessI1_bessI0;
            double bessI2_bessI0;
            if(z<700) {
                double bessI0 = BesselFunctions.besselI0(z);
                logBessI0 = Math.log(bessI0);
            }
            else {
                // These are log linear approximations.  See
                // ActiveImaging/models/besseli1d0.m for details.
                logBessI0 = z - Math.log(2.0*Math.PI*z)/2.0;
            }

            // Constant terms excluded for a slight speed gain.
            loglik += 0.5*sig2i*ymod*ymod -logBessI0;
        }

        if(print)
            System.err.println("loglik: " + loglik);
		
        return loglik;
    }
	
	
}


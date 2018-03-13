package fitters;

import java.util.logging.Logger;
import java.util.Random;
import java.io.*;
import models.*;
import optimizers.*;
import inverters.*;
import numerics.*;
import misc.*;
import data.*;
import tools.*;
import imaging.*;
import models.compartments.CompartmentModel;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General purpose MCMC fitter with implementations of various
 * features of all fitters that use MCMC.
 * 
 * <dt>Description:
 * 
 * <dd> .
 * 
 * </dl>
 *
 * @author Danny
 * @version $Id$
 *
 */
public abstract class MCMC_Fitter extends CompartmentFitter {

    // Prior object
    protected Prior prior;

    // Perturbation object
    protected Perturbation perturbation;

    // Output stream for acceptance rates etc
    protected DataOutputStream statsFile;


    /**
     * Default constructor.
     */
    public MCMC_Fitter() {
    }


    /**
     * Constructor sets up various required objects for MCMC.
     *
     * @param scheme The imaging protocol.
     * @param compNames the compartment names.
     * @param initialParams the initial parameters.
     */
    public MCMC_Fitter(DW_Scheme scheme, String[] compNames, double[] initialParams) {
        
        if(CL_Initializer.mcmcStatsFile != null) try {
            statsFile = new DataOutputStream(new FileOutputStream(CL_Initializer.mcmcStatsFile));
        } catch(Exception e) {
            throw new LoggedException(e);
        }
        
        this.scheme = scheme;
        makeCodec();
        makePrior();
        makePerturbation();
        cm = new CompartmentModel(compNames, initialParams);
        
        try {
	    initMCMC_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
        } catch (Exception e) {
            throw new LoggedException(e);
        }
        
    }
    

    /**
     * Initializes the minimizer object to a MCMC_Minimizer with
     * the specified noise model.
     *
     * @param nm Specifies the noise model.
     */
    protected void initMCMC_Minimizer(NoiseModel nm) throws MarkovChainMonteCarloException {
	if(CL_Initializer.sigma < 0) {
	    throw new LoggedException("No valid noise level specified; required for MCMC: use -sigma.");
	}

        if(nm == NoiseModel.GAUSSIAN) {
            minimizer = new MCMC_GaussianMinimizer(scheme, cm, cod, prior, perturbation, CL_Initializer.sigma);
        }
        else if(nm == NoiseModel.OFFGAUSS) {
            if(CL_Initializer.sigma < 0) {
                throw new LoggedException("No valid noise level specified: use -sigma.");
            }
            minimizer = new MCMC_OffGaussMinimizer(scheme, cm, cod, prior, perturbation, CL_Initializer.sigma);
        }
        else if(nm == NoiseModel.RICIAN) {
            if(CL_Initializer.sigma < 0) {
                throw new LoggedException("No valid noise level specified: use -sigma.");
            }
            minimizer = new MCMC_RicianMinimizer(scheme, cm, cod, prior, perturbation, CL_Initializer.sigma);
        }
	
        else {
            throw new LoggedException("Unknown noise model: " + nm);
        }
    }


    public double[][] fit(double[] data) throws MinimizerException {
        double[][] sol = super.fit(data);
        if(CL_Initializer.mcmcStatsFile != null) try {
            statsFile.writeDouble(((MCMC_Minimizer)minimizer).getAcceptanceRate());
        } catch(Exception e) {
            throw new LoggedException(e);
        }
        return sol;
    }
    

    /**
     * Creates the Codec object encoding the mapping
     * from model to optimized parameters.
     */
    protected abstract void makeCodec();


    /**
     * Creates the Prior object.
     */
    protected abstract void makePrior();


    /**
     * Creates the Perturbation object.
     */
    protected abstract void makePerturbation();


}


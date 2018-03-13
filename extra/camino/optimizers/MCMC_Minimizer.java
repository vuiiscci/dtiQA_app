package optimizers;

import java.util.logging.Logger;
import models.*;
import optimizers.*;
import numerics.*;
import misc.*;
import imaging.*;
import fitters.*;


/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Abstract MCMC class implementing Minimizer interfase. This extemsion implements 
 * getSolutions ans setMeasurements methods
 * 
 * <dt>Description:
 * 
 * <dd>This is extension of MarkovChainMonteCarlo
 * which implements remaining two methods of Minimizer interfase.
 * 
 * </dl>
 * 
 * @version $Id: MCMC_Minimizer.java 2010-07-20 $
 * @author Aidos Abzhanov
 *  
 */
public abstract class MCMC_Minimizer extends MarkovChainMonteCarlo implements Minimizer {
	
	/**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.optimizers.MCMC_Minimizer");
	
	/**
     * Scheme object; describes the acquisition protocol
     */
    protected DW_Scheme scheme;
	
	/**
     * Model of diffusion which is to be fit
     */
    protected ParametricModel model;
	
	/**
     * Codec object used to convert the model parameters
	 * into optimised ones and back.
     */
    protected Codec codec;
	
	/**
     * Dependent values at sampled points: the measurements.
     */
    protected double[] measurements;
	
    /**
     * Returns the list of M samples from the MCMC.
     *
     * @return Mx(N+2) matrix in which each row contains an exitcode,
     * the N model parameters and the objective function value.
     */
    public double[][] getSolutions() {
		
		double[][] solutions = new double[samplesNumber][codec.getNumModelParams() + 2];
		
		for (int i =0; i < samplesNumber; i++) {
			
		    // Compute the model parameters for return
		    double[] optParams = new double[paramsNumber];
		    for(int j=0; j<paramsNumber; j++)
			optParams[j] = samples[i][j+1];
		    double[] modelParams = codec.optToModel(optParams);
		    
		    // Copy in the exit code first.
		    solutions[i][0] = samples[i][0];

		    // Now the model parameters
		    for(int j=0; j<modelParams.length; j++)
			solutions[i][j+1] = modelParams[j];
		    
		    // Finally the objective function value.
		    solutions[i][solutions[i].length-1] = samples[i][paramsNumber+1];
		}
		
		return solutions;
    }
	
	
    public int getNumSolutions() {
        return samplesNumber;
    }
    
    public int getNumParameters()
    {
    	return model.numParams()+2; //error code and exit code 
    }


    /**
     * Initializes the fitting procedure with a new set of
     * measurements (dependent variables).
     * 
     * @param newMeas The new set of measurements.
     */
    public void setMeasurements(double[] newMeas) throws MarkovChainMonteCarloException {
        if (newMeas.length != scheme.numMeasurements()) {
            throw new MarkovChainMonteCarloException("New data contains the wrong number of values: " 
													 + newMeas.length + " instead of " + measurements.length + ".");
        }
        for (int i = 0; i < newMeas.length; i++) {
            measurements[i] = newMeas[i];
        }
		
    }
	
	
}



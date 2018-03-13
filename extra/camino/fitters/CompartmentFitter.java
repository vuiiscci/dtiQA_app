package fitters;

import misc.LoggedException;
import models.compartments.CompartmentModel;
import tools.CL_Initializer;
import optimizers.*;

/**
 * Fitter class specific to compartment models. At the moment it handles the initialisation
 * of the starting point from command line for all the compartment models.
 * @author laura
 *
 */
public abstract class CompartmentFitter extends Fitter {
	
	
	
	protected CompartmentModel cm;
	
	
	

	  /**
     * Reads the starting set of parameters from the command line 
     * if the fixed start point option was set else returns null. 
     * 
     * Compartment names have to agree in name and order with the fitted compartment model,
     *  otherwise it throws an exception.
     *
     * @param data is ignored.
     *
     * @return The starting model parameter values.
     */
	protected double[] getFixedStartPoint(double[] data) {
		
		// Set starting point from command line
		if (CL_Initializer.fixedStartPoint) {
			String[] compartmentnames = cm.getCompartmentnames();
			
			for (int i=0; i<compartmentnames.length;i++)
			{
				if (!compartmentnames[i].equalsIgnoreCase(CL_Initializer.compartmentNames[i]))
				{
					throw new LoggedException("Fixed startpoint doesn't match fitmodel." + "Given: " + CL_Initializer.compartmentNames[i] + ". Expected: " + compartmentnames[i] );
				}
			}
			return CL_Initializer.compParams;
			
		}

		return null;
	}
	

    /**
     * Initializes the minimizer object to a LM_Minimizer with
     * the specified noise model.
     *
     * @param nm Specifies the noise model.
     */
    protected void initLM_Minimizer(NoiseModel nm) throws MarquardtMinimiserException {
        if(nm == NoiseModel.GAUSSIAN) {
            minimizer = new LM_GaussianMinimizer(scheme, cm, cod);
        }
        else if(nm == NoiseModel.OFFGAUSS) {
            if(CL_Initializer.sigma < 0) {
                throw new LoggedException("No valid noise level specified: use -sigma.");
            }
            minimizer = new LM_OffGaussMinimizer(scheme, cm, cod, CL_Initializer.sigma);
        }
        else if(nm == NoiseModel.RICIAN) {
            if(CL_Initializer.sigma < 0) {
                throw new LoggedException("No valid noise level specified: use -sigma.");
            }
            minimizer = new LM_RicianMinimizer(scheme, cm, cod, CL_Initializer.sigma);
        }
        else {
            throw new LoggedException("Unknown noise model: " + nm);
        }
    }


    /**
     * Initializes the minimizer object to a MultiRunLM_Minimizer with
     * the specified noise model.
     *
     * @param nm Specifies the noise model.
     */
    protected void initMultiRunLM_Minimizer(NoiseModel nm, Perturbation p, int repeats, int seed) throws MinimizerException {
        if(nm == NoiseModel.GAUSSIAN) {
            minimizer = new MultiRunLM_GaussianMinimizer(scheme, cm, cod, p, repeats, seed);
        }
        else if(nm == NoiseModel.OFFGAUSS) {
            if(CL_Initializer.sigma < 0) {
                throw new LoggedException("No valid noise level specified: use -sigma.");
            }
            minimizer = new MultiRunLM_OffGaussMinimizer(scheme, cm, cod, p, repeats, seed, CL_Initializer.sigma);
        }
        else if(nm == NoiseModel.RICIAN) {
            if(CL_Initializer.sigma < 0) {
                throw new LoggedException("No valid noise level specified: use -sigma.");
            }
            minimizer = new MultiRunLM_RicianMinimizer(scheme, cm, cod, p, repeats, seed, CL_Initializer.sigma);
	}
        else {
            throw new LoggedException("Unknown noise model: " + nm);
        }
    }

}

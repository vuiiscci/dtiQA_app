package optimizers;

import java.util.logging.*;
import java.util.Random;
import java.io.*;
import java.net.URL;
import tools.*;
import fitters.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General Markov Chain Monte Carlo minimisation algorithm. Algorithm minimises
 * fObj function, which is implemented in MCMC_GaussianMinimizer class	
 * 
 * <dt>Description:
 * 
 * <dd>This is the Markov Chain Monte Carlo code for parameters estimation.
 * 
 * </dl>
 * 
 * @version $Id: MarkovChainMonteCarlo.java 2010-07-20 $
 * @author Aidos Abzhanov
 *  
 */
abstract public class MarkovChainMonteCarlo {
	
	
	/**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.optimizers.MarkovChainMonteCarlo");
	
	/**
     * MCMC parameter. Burn-in period before starting sample collection
     */
	protected int burninTime;
	
	/**
     * MCMC parameter. Interval period between sample collection
     */
	protected int intervalTime;
	
	/**
     * MCMC parameter. Number of samples produced by MCMC
     */
	protected int samplesNumber;
	
	/**
     * MCMC output. Array of optimized samples, each sample is a array of the model parameters
     */
	protected double[][] samples;
	
	/**
     * MCMC output. Number of acceptances in a single run.
     */
	protected int acceptances;
	
	/**
     * Random number generator
     */
	protected Random rand;

    /**
     * The prior distribution.
     */
    protected Prior prior;

    /**
     * The proposal distributions.
     */
    protected Perturbation perturbation;

	//Parameters of the model
	
	//Number of parameters
	/**
     * Nuber of the model parameters
     */
	protected int paramsNumber;
	
	/**
     * Array of current parameters values.
     */
	protected double[] paramsCurrent;
	
	/**
     * Array of initial parameters values.
     */
	protected double[] paramsInitial;
	
	/**
	 * Value of objective function for current parameters.
	 */
	protected double fObjCurrent;
	
    /**
     * Value of prior for current parameters.
     */
    protected double priorCurrent;
	
	/**
     * Serial number of voxel being processing by minimise() methos
     */
	protected int voxelIndex;
	
    /**
     * Constructor. Initialises MCMC parameters
     * given a Prior object.
     */
	public MarkovChainMonteCarlo()
	{
		//MCMC algorithm parameters initialization
		this.burninTime  = CL_Initializer.burnInTime;
                this.intervalTime = CL_Initializer.interval;
		this.samplesNumber = CL_Initializer.samples;
		
		this.rand = new Random(CL_Initializer.seed);
		
		//Initially, 0 voxels have been minimised
		voxelIndex = 0;
		
	}	
	
	/**
	 * Returns the value of the objective function with parameters
	 * atry, which provides the likelihood of atry for MCMC.  The
	 * method to be implemented for particular noise models.
	 * 
	 * @param atry
	 *            The point at which to evaluate the objective function.
	 * 
	 * @return The value of the objective function.
	 */
	abstract protected double fObj(double[] atry);
	
	
	/**
	 * Initialises working arrays.
	 * 
	 * @param noParams
	 *            The number of parameters to fit.
	 */
	protected void init(int noParams) {
		paramsNumber  = noParams;
		paramsCurrent = new double[paramsNumber];
		paramsInitial = new double[paramsNumber];
		samples = new double[samplesNumber][paramsNumber+2];
		for (int i = 0; i <paramsNumber; i++) {
			paramsCurrent[i] = 0.0;
			paramsInitial[i] = 0.0;
		}
	}
	
	/**
	 * Sets the initial values of the parameters. 
	 * 
	 * @param aInit
	 *            Array containing the new parameter values
	 */
	public void setInitParams(double[] aInit) throws MarkovChainMonteCarloException {
		
		if (aInit.length != paramsNumber) {
			throw new MarkovChainMonteCarloException(
												  "Wrong number of parameters in initializing array.  Got "
												  + aInit.length + " expected " + paramsNumber + ".");
		}
		for (int i = 0; i < aInit.length; i++) {
			paramsInitial[i] = aInit[i];
			paramsCurrent[i] = aInit[i];
		}
		fObjCurrent = fObj(paramsCurrent);
		priorCurrent = prior.prior(paramsCurrent);
	}
	
	/**
	 * Returns the values of the current parameters.
	 * 
	 * @return The values of the current parameters.
	 */
	public double[] getParameters() {
		double[] aCurr = new double[paramsCurrent.length];
		for (int i = 0; i < aCurr.length; i++) {
			aCurr[i] = paramsCurrent[i];
		}
		return aCurr;
	}
		
	/**
	 * Runs the minimization. Main MCMC routine. minimise() function
	 * is declared in the Minimizer interface. 
	 */
	public void minimise() throws MarkovChainMonteCarloException {
		
		//The index of a voxel being processed by MCMC minimise method
		voxelIndex++;

                // Initialize the acceptance count
                acceptances = 0;
		
                for (int i = 0; i < burninTime; i++) {
                    update();				
                }
                for (int i = 0; i < samplesNumber; i++) {
                    for (int j = 0; j < intervalTime; j++) {
                        update();
                    }
		    // The exit code is always zero.
		    samples[i][0] = 0;
                    for (int k = 0; k < paramsNumber; k++) {
                        samples[i][k+1] = paramsCurrent[k];
                    }
                    samples[i][paramsNumber+1] = fObjCurrent;
                }
	}


    /**
     * Returns the acceptance rate of the previous run, ie the
     * fraction of updates that resulted in a new current set
     * of parameters.
     *
     * @return number of accepts / (burn in iterations + number samples*interval)
     */
    public double getAcceptanceRate() {
        return (double)acceptances/(double)(burninTime + samplesNumber*intervalTime);
    }


	/**
	 * One step of MCMC. Proposes new parameters, and accepts them with 
	 * probabitily min(1,likelihoodRatio)
	 */
	private void update() throws MarkovChainMonteCarloException {
		
            double[] paramsNew = perturbation.perturb(paramsCurrent, paramsInitial, rand);
            double priorNew = prior.prior(paramsNew);
            if(priorNew>0) {
                
                double fObjNew = fObj(paramsNew);
                
                // Work out the likelihood ratio for accepting
                // perturbation.
                double f = fObjCurrent + Math.log(priorCurrent) - fObjNew - Math.log(priorNew);
		
                // if f less than -100, likelihoodRatio = 0.0
                // avoids numerical errors.
                double likelihoodRatio = 0.0;
                if (f > 0)
                    likelihoodRatio = 1;
                else if (f>(-100)) 
                    likelihoodRatio = Math.exp(f);
		
                //If likelihood ratio becomes infinite, give up.
                if (Double.isNaN(likelihoodRatio) || Double.isInfinite(likelihoodRatio)) {
                    System.err.println(f);
                    System.err.println(likelihoodRatio);
                    throw new MarkovChainMonteCarloException("Infinite likelihood ratio in Update function.");
                }
        
		
                // Now make the decision
                double r = rand.nextDouble();		
                if (r < likelihoodRatio) {
                    for (int i = 0; i < paramsNumber; i++) {
                        paramsCurrent[i] = paramsNew[i];
                    }
                    fObjCurrent = fObjNew;
                    priorCurrent = priorNew;
                    acceptances += 1;
                }
            }
	}

	
}



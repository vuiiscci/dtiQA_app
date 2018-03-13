package optimizers;

import tools.FileOutput;

import java.util.Random;
import java.util.Arrays;

import java.util.logging.Logger;

import java.text.*;

/**
 *
 * General superclass for simulated annealing optimization. Concrete subclasses must describe a system
 * whose parameters (either discrete or continuous) to be optimized by the minimization of an energy function. 
 * Subclasses must also provide a method to suggest changes in state to the algorithm.
 * 
 *
 * @author Philip Cook
 * @version $Id$
 */
public abstract class SimulatedAnnealingOptimizer {


    protected double temperature = 0.0;

    private final double coolingRate;

    private final int iterationsBetweenCooling;

    private final long iterationsBeforeQuit;

    private static Logger logger = Logger.getLogger("camino.optimizers.SimulatedAnnealingOptimizer");

    protected final Random ran;

    protected String saveStateRoot = "SimulatedAnnealingOptimizer";


    /**
     * @param coolRate a constant. Unless the method cool(double) is overridden, the temperature is lowered 
     * by a factor of (1 - coolRate) after <code>iterationsPerT</code> iterations.
     * @param iterationsPerT attempt this many changes to the state of the system before lowering the temperature.
     * @param iterationsBeforeQuit Optimization will end if the energy is not lowered after this many iterations.
     * state before quitting.
     * @param saveStateRoot Root of file to save the state of the system during optimization.
     * @param r a random number generator.
     *
     *
     */
    public SimulatedAnnealingOptimizer(double coolRate, int iterationsPerT, long iterationsBeforeQuit, 
				       String saveState, Random r) {

	coolingRate = coolRate;

	iterationsBetweenCooling = iterationsPerT;

	this.iterationsBeforeQuit = iterationsBeforeQuit;

	ran = r;
	saveStateRoot = saveState;

    }

    
    public void setTemperature(double t) {
	temperature = t;
    }

    /**
     * Runs the minimization.
     *
     */
    public final void minimize() {

	long start = System.currentTimeMillis();

	int hoursPassed = 0;
	
	int counter = 0;

	boolean lastRun = false;

	double stateEnergy = stateEnergy();
	
	// has best energy been saved to disk
	boolean writtenBestEnergy = false;

	long iterationsSinceLastLowered = 0;

	double bestEverStateEnergy = Double.MAX_VALUE;

	double tempLastLowered = 0.0;
	
	logger.info("Initial temperature\t" + temperature + 
		    "\nCooling rate\tT(i+1) = (1.0 - " + coolingRate + ")T" +
		    "\nIterations\t" + iterationsBetweenCooling +
		     "\nInitial state energy \t" + stateEnergy);


	while (!lastRun) {

	    if (temperature <= 1E-8) {
		lastRun = true;
		temperature = 0.0;
	    }
	    
	    boolean loweredBestEnergy = false;

	    for (int it = 0; it < iterationsBetweenCooling; it++) {

		// get suggested change
		
		double tmpStateEnergy = generateCandidateState();
		
		boolean accept = false;
		
		double deltaE = tmpStateEnergy - stateEnergy;
		
		if (deltaE < 0.0) {
		    accept = true;
		}
		else {
		    if (temperature > 0.0) {
			accept = ran.nextDouble() <= Math.exp(-1.0 * deltaE / temperature);
		    }
		}
		
		if ( accept ) {
		    // swap
		    acceptCandidateState();

		    stateEnergy = tmpStateEnergy;

		    if (stateEnergy < bestEverStateEnergy) {
			loweredBestEnergy = true;
			writtenBestEnergy = false;
			bestEverStateEnergy = stateEnergy;
			saveStateAsLowestEnergy();
		    }

		} // end if swap
		    
	    } // end for iterations
		
	    double tmpStateEnergy = stateEnergy;
	    
	    stateEnergy = stateEnergy();

	    if (Math.abs(stateEnergy - tmpStateEnergy) > 1E-8) {
		throw new misc.LoggedException("Expected energy : " + tmpStateEnergy + "\n\t actual energy : "  + stateEnergy);
	    }
	    
	    if (!loweredBestEnergy) {
		iterationsSinceLastLowered += iterationsBetweenCooling;
		if (iterationsSinceLastLowered > iterationsBeforeQuit) {
		    logger.info("Giving up (" + iterationsBeforeQuit + " iterations since energy last lowered)"); 
		    logger.info("Current temperature : " + temperature);
		    logger.info("Temperature of last improvement: " + tempLastLowered);
		    lastRun = true;
		}
	    }
	    else {
		iterationsSinceLastLowered = 0l;
		tempLastLowered = temperature;
	    }
	    
	    if (!lastRun) {
		temperature = cool(temperature);
	    }
	    
	    if ((int)(System.currentTimeMillis() - start) / (1000 * 3600) > hoursPassed) {
		hoursPassed = (int)(System.currentTimeMillis() - start) / (1000 * 3600);
			   
		// log current state 
		FileOutput out = new FileOutput(saveStateRoot + ".state");
		out.writeString(state());
		out.close();
		
		if (!writtenBestEnergy) {
		    out = new FileOutput(saveStateRoot + ".lowestEnergy");
		    out.writeString(lowestEnergyState());
		    out.close();
		    writtenBestEnergy = true;
		}
		
		logger.info(hoursPassed + " hours elapsed. Current state saved to " + saveStateRoot + 
			    ".state. Current temperature is: " + temperature);
		
	    }

	    
	} // end while
	
	    
    }


    /**
     * Lowers the temperature of the annealing process. Implements basic scheme T = (1 - e)T where e
     * is a constant.
     *
     */
    protected double cool(double temperature) {
	return (1.0 - coolingRate) * temperature;
    }


    /**
     * @return the energy of the current state.
     *
     */
    protected abstract double stateEnergy();


    /**
     * Generates a candidate state. The candidate state is not guaranteed to be randomly generated. 
     *
     * @return the energy of the candidate state.
     */
    protected abstract double generateCandidateState();

    
    /**
     * Changes the current state of the system to the last candidate.
     *
     */
    protected abstract void acceptCandidateState();
    

    /**
     * Saves the current state as the lowest energy state.
     *
     */
    protected abstract void saveStateAsLowestEnergy();


    /**
     * Sets the current state to the lowest energy state discovered so far.
     *
     */
    protected abstract void setStateToLowestEnergy();


    /**
     * Gets a String representation of the current state of the system.
     */
    public abstract String state();
    

    /**
     * Gets a String representation of the lowest energy state of the system.
     */
    public abstract String lowestEnergyState();

   

    /**
     * Calibrates the temperature by moving to the nearest local minimum and calculating the mean 
     * increase in energy of attempted state changes. The temperature is then calculated to give a
     * specified acceptance rate for the mean uphill step. 
     * <p>
     * This method assumes that the system will quickly reach a local minimum if no positive energy
     * changes are allowed (if this is not the case, perhaps simulated annealing is not the most suitable
     * optimization algorithm), and hence that the mean energy change over a large number of configuration
     * change will be positive. 
     *
     * @return the temperature that would give the specified acceptance rate for uphill jumps.
     *
     */
    public double calibratedTemp(double acceptanceRate) {
	
	// local variable, so that the calibration is not dependent on the system temperature
	// this means the method gives us the energy changes associated with trying to get
	// out of the nearest local minima

	double temp = 0.0;

	double initStateEnergy = stateEnergy();

	double meanDeltaE = 0.0;

	int runs = 1000000;

	double[] deltaE = new double[runs];
	
	double stateEnergy = initStateEnergy;

	double meanStateEnergy = stateEnergy / runs;

	logger.info("Calibrating...");

	logger.info("Input state energy (before calibration) = " + initStateEnergy);

	long start = System.currentTimeMillis();

	for (int it = 0; it < runs; it++) {

	    double tmpStateEnergy = generateCandidateState();
	    
	    boolean accept = false;
	    
	    deltaE[it] = tmpStateEnergy - stateEnergy;
	    
	    if (deltaE[it] <= 0.0) {
		accept = true;
	    }
	    else {
		if (temp > 0.0) {
		    double d = ran.nextDouble();
		    accept = d <= Math.exp(-1.0 * deltaE[it] / temp);
		}
	    }
	    
	    if ( accept ) {
		acceptCandidateState();
		stateEnergy = tmpStateEnergy;
	    } 

	    meanStateEnergy += stateEnergy / runs;

	    meanDeltaE += deltaE[it] / runs;

	}

	long end = System.currentTimeMillis();

	double stdDeltaE = std(deltaE, meanDeltaE);

	logger.info("mean state energy = " + meanStateEnergy + 
		    "\nmean (deltaE) = " + meanDeltaE +
		    "\nsigma (deltaE) = " + stdDeltaE);

	StringBuffer buffer = new StringBuffer();

	buffer.append("Acceptance rate for mean deltaE // Temperature\n");
	
	DecimalFormat df = new DecimalFormat("0.000");

	for (int i = 1; i < 10; i++) {
	    buffer.append(i / 10.0 + "\t" + df.format(-1.0 * meanDeltaE / Math.log(i / 10.0)) + "\n");
	}
	
	buffer.append("\n");

	logger.info(buffer.toString());

	int calibrationTime = (int)(end - start);

	logger.info("Calibration run time == " + calibrationTime + " msecs to complete " + runs + " runs");

	return -1.0 * meanDeltaE / Math.log(acceptanceRate);
	
     }



    /**
     * Standard deviation of a sample: \sum_{i=1}^N (x_i - \bar{x})^2 / (N-1). 
     */
    protected static double std(double[] data, double mean) {
	
	double var = 0.0;
	
	for (int i = 0; i < data.length; i++) {
	    var += (data[i] - mean) * (data[i] - mean);
	}
	
	var = var / (data.length - 1.0);
	

	return Math.sqrt(var);
    }




}

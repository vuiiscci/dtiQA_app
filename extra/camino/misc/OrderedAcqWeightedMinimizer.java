package misc;

import numerics.*;


import java.util.Random;
import java.util.Arrays;

public class OrderedAcqWeightedMinimizer extends OrderedAcqMinimizer {

    private double[][] bestEverPoints;
    
    private double[] weights;

    private final double weightPower = 2.0;

    /** Energy if each point and the points before it. */
    private double[] partialEnergy = null;
    private double[] candidatePartialEnergy = null;

    /** candidate state change swaps two points, store the index of those points here. */
    private int i = 0;
    private int j = 0;

    
    private double currentStateEnergy = 0.0;

    private double candidateStateEnergy = 0.0;

    private double bestEverStateEnergy = 0.0;


    public OrderedAcqWeightedMinimizer(int numberOfPairs, double coolRate, int iterationsPerT, 
				       String saveState, Random r) {
	
	super(numberOfPairs, coolRate, iterationsPerT, saveState, r);

	setWeights();

	partialEnergy = new double[noPairs];

	stateEnergy();
	saveStateAsLowestEnergy();

    }

    
    public OrderedAcqWeightedMinimizer(double[][] points, double coolRate, int iterationsPerT, 
				       String saveState, Random r) {

	super(points, coolRate, iterationsPerT, saveState, r);

	setWeights();

	partialEnergy = new double[noPairs];

	stateEnergy();
	saveStateAsLowestEnergy();

    }

    
    private void setWeights() {
    
	weights = new double[noPairs];
	
	for (int p = 1; p < noPairs; p++) {
	    weights[p] = 1.0 / Math.pow(p+1, weightPower);
	}
    }



    /**
     * @return the state energy of the current ordering.
     *
     */
    protected double stateEnergy() {

	for (int p = 0; p < noPairs; p++) {
	    
	    partialEnergy[p] = 0.0;
	    
	    for (int q = 0; q < p; q++) {
		// add the energy between this point and all points before it
		partialEnergy[p] += energy[p][q]; 
	    }
	    
	}


	double previousEnergy = 0.0; // last increment, used in calculation
	double stateEnergy = 0.0;
	
	for (int p = 0; p < noPairs; p++) {
	    
	    previousEnergy += partialEnergy[p];
	    
	    if (p > 4) {
		stateEnergy += previousEnergy * weights[p];
	    }
	    
	}

	currentStateEnergy = stateEnergy;

	return stateEnergy;
    }



    /**
     *
     * @return the state energy that would result from swapping points i and j.
     *
     */
    private double swappedStateEnergy() {
	
	candidatePartialEnergy = new double[noPairs];

	System.arraycopy(partialEnergy, 0, candidatePartialEnergy, 0, noPairs);

	for (int si = i+1; si < j; si++) {
	    
	    candidatePartialEnergy[si] = partialEnergy[si] + energy[si][j] - energy[si][i];
			
	    candidatePartialEnergy[i] += energy[si][i]; // point i now acquired after all si
	    candidatePartialEnergy[j] -= energy[si][j]; // point j now acquired before all si
	    

	} // end for si
	
	candidatePartialEnergy[i] += energy[i][j]; // point i now after j
	candidatePartialEnergy[j] -= energy[j][i]; // point j now before i
	
	// now swap position of i and j
	double tmp = candidatePartialEnergy[i];
	candidatePartialEnergy[i] = candidatePartialEnergy[j];
	candidatePartialEnergy[j] = tmp;

	// now ready to calculate total energy
		    
	
	double tmpStateEnergy = 0.0;
	
	double previousEnergy = 0.0;

	for (int p = 0; p < noPairs; p++) {
	    previousEnergy += candidatePartialEnergy[p];

	    if (p > 4) {
		tmpStateEnergy += previousEnergy * weights[p];
	    }
	}	
	
	candidateStateEnergy = tmpStateEnergy;

	return tmpStateEnergy;
	
    }




   
    private void calculateCandidateSwapPositions() {


	i = (int)(ran.nextFloat() * noPairs);
	j = (int)(ran.nextFloat() * noPairs);

	while ((i < 5 && j < 5) || i == j) {
	    j = (int)(ran.nextFloat() * noPairs);
	}



	// enforce i < j
	
	if (i > j) {
	    int tmp = j;
	    j = i;
	    i = tmp;
	}
	
    }


  
    /**
     * @return the change in energy between the current state and the candidate state.
     */
    protected double generateCandidateState() {

	calculateCandidateSwapPositions();

	return swappedStateEnergy();

    }

    
    /**
     * Changes the current state of the system to the last candidate.
     *
     */
    protected void acceptCandidateState() {
	partialEnergy = candidatePartialEnergy;

	// sanity check
	if (candidateStateEnergy == 0.0) {
	    throw new LoggedException("Attempting to set an invalid state.");
	}

	currentStateEnergy = candidateStateEnergy;
	
	// stops any attempt to reuse this
	candidatePartialEnergy = null;
	
	candidateStateEnergy = 0.0;
	
	swapPoints(i,j);
	
    }
    

    /**
     * Saves the current state as the lowest energy state.
     *
     */
    protected void saveStateAsLowestEnergy() {
	bestEverPoints = getPoints();
	bestEverStateEnergy = currentStateEnergy;
    }


    /**
     * Sets the current state to the lowest energy state discovered so far.
     *
     */
    protected void setStateToLowestEnergy() {
	points = bestEverPoints;

	// make a copy of the best points
	bestEverPoints = getPoints();

	buildEnergyMatrix();
	stateEnergy();
    }


    /**
     * Gets a String representation of the current state of the system.
     */
    public String state() {
	return state(points, temperature, currentStateEnergy);
    }


    /**
     * Gets a String representation of the lowest energy state of the system.
     */
    public String lowestEnergyState() {
	return state(bestEverPoints, temperature, bestEverStateEnergy);
    }


   
    private String state(double[][] points, double temp, double energy) {

	String s = noPairs + " " + temp + " " + energy + "\n";
	
	for (int n = 0; n < noPairs; n++) {
	    s += points[n][0] + " " + points[n][1] + " " + points[n][2] + "\n";
	}

	return s;

    }
    



}

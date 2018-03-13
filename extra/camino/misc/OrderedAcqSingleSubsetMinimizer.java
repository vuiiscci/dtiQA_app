package misc;

import numerics.*;


import java.util.Random;
import java.util.Arrays;

public class OrderedAcqSingleSubsetMinimizer extends OrderedAcqMinimizer {

    private double[][] bestEverPoints;
    

    /** candidate state change swaps two points, store the index of those points here. */
    private int i = 0;
    private int j = 0;

    private double currentStateEnergy = 0.0;

    private double candidateStateEnergy = 0.0;

    private double bestEverStateEnergy = 0.0;

    private final int pairsInSubset;


    public OrderedAcqSingleSubsetMinimizer(int numberOfPairs, int pairsInSubset, double coolRate, 
					   int iterationsPerT, String saveState, Random r) {
	
	super(numberOfPairs, coolRate, iterationsPerT, saveState, r);

	this.pairsInSubset = pairsInSubset;

	stateEnergy();
	saveStateAsLowestEnergy();

    }

    
    public OrderedAcqSingleSubsetMinimizer(double[][] points, int pairsInSubset, double coolRate, 
				       int iterationsPerT, String saveState, Random r) {
	
	super(points, coolRate, iterationsPerT, saveState, r);

	this.pairsInSubset = pairsInSubset;

	stateEnergy();
	saveStateAsLowestEnergy();

    }

    

    /**
     * @return the state energy of the current subset.
     *
     */
    protected double stateEnergy() {

	double stateEnergy = 0.0;

	for (int p = 0; p < pairsInSubset; p++) {
	    for (int q = 0; q < p; q++) {
		stateEnergy += energy[p][q]; 
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

	double tmpStateEnergy = currentStateEnergy;
		    
	for (int p = 0; p < pairsInSubset; p++) {
	    if (p != i) {
		tmpStateEnergy -= energy[p][i]; 
		tmpStateEnergy += energy[p][j];
	    }
	}

	candidateStateEnergy = tmpStateEnergy;

	return tmpStateEnergy;
	
    }

   
    private void calculateCandidateSwapPositions() {

	i = (int)(ran.nextFloat() * pairsInSubset);
	j = pairsInSubset + (int)(ran.nextFloat() * (noPairs - pairsInSubset));

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

	// sanity check
	if (candidateStateEnergy == 0.0) {
	    throw new LoggedException("Attempting to set an invalid state.");
	}

	currentStateEnergy = candidateStateEnergy;
	
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

	// 1 subset
	String s = noPairs + " 1 " + pairsInSubset + " " + temp + " " + energy + "\n";
	
	for (int n = 0; n < noPairs; n++) {
	    s += points[n][0] + " " + points[n][1] + " " + points[n][2] + "\n";
	}

	return s;

    }
    



}

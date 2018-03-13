package misc;


import java.util.Random;
import java.util.Arrays;

import java.util.logging.Logger;


/**
 * Divides a set of gradient directions into subsets, and searches for a configuration
 * where the sum of the electrostatic energy of each subset is minimal.
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class OrderedAcqSubsetMinimizer extends OrderedAcqMinimizer {

    // Each charge pair is assigned a subset, 0-(N-1), where there are N subsets. 
    // Pairs are repelled by other pairs that share the same subset, and have zero 
    // interaction with pairs that do not share the same subset. 

    private int noSubsets;

    // total subset energy of each subset.
    // this is the energy between pairs of the same subset.
    private double[] subsetEnergy;
    
    // number of pairs that have each subset
    private int[] pairsPerSubset;

    // subset of each pair
    private int[] subsets;


    private int[] candidateSubsets = null;
    private double[] candidateSubsetEnergy = null;


    private int i = 0;
    private int j = 0;


    private int[] bestEverSubsets = null;
    private double[] bestEverSubsetEnergy = null;


    private static Logger logger = Logger.getLogger("camino.misc.OrderedAcqSubsetMinimizer");



    /**
     * @param numPairsPerSubset contains the number of pairs in each subset.
     * 
     *
     */
    public OrderedAcqSubsetMinimizer(int numberOfPairs, int[] numPairsPerSubset, double coolRate, 
				     int iterationsPerT, String saveState, Random r) {

	super(numberOfPairs, coolRate, iterationsPerT, saveState, r);

	bestEverSubsets = new int[noPairs];
	pairsPerSubset = numPairsPerSubset;

	noSubsets = pairsPerSubset.length;

	init();
    }


    /**
     * @param points the point set, where the points in each subset are listed in sequence. Each
     * point is a unit vector that describes the axes of one pair of charged particles. The first 
     * <code>n = numPairsPerSubset[0]</code> points should be the
     * <code>n</code> points that make up the first subset. The next <code>m = numPairsPerSubset[1]</code>
     * points should be the second subset, and so on.
     *
     * @param numPairsPerSubset the number of points in each subset. 
     *
     */
    public OrderedAcqSubsetMinimizer(double[][] points, int[] numPairsPerSubset, double coolRate, int iterationsPerT, String saveState, Random r) {
	super(points, coolRate, iterationsPerT, saveState, r);

	pairsPerSubset = numPairsPerSubset;

	noSubsets = pairsPerSubset.length;

	init();
	
    }
    


    private void init() {
	
	String s = "Optimizing for " + noPairs + " pairs and " + noSubsets + " subsets.\nPairs per subset: ";
	
	for (int n = 0; n < noSubsets; n++) {
	    s += pairsPerSubset[n] + " ";
	}

	logger.info(s);

	initializeSubsets();

	stateEnergy();

	saveStateAsLowestEnergy();

    }



    /**
     * Assign subsets to points, in order. If the subsets should be assigned randomly, 
     * randomize the order of the points before calling this method.
     *
     */
    private void initializeSubsets() {

	subsets = new int[noPairs];

	int[] subsetBank = new int[noPairs]; // contains list of subsets for each pair.

	int subsetCounter = 0;
	int currentSubset = 0;
	
	for (int p = 0; p < pairsPerSubset.length; p++) {
	    for (int q = 0; q < pairsPerSubset[p]; q++) {
		subsetBank[subsetCounter++] = currentSubset;
	    }
	    currentSubset++;
	}

	for (int p = 0; p < noPairs; p++) {
	    subsets[p] = subsetBank[p];
	}

    }


    private void calculateCandidateSwapPositions() {

	i = (int)(ran.nextDouble() * noPairs);
	j = (int)(ran.nextDouble() * noPairs);

	while (i == j || subsets[i] == subsets[j]) {
	    j = (int)(ran.nextDouble() * noPairs);
	}

    }
    

    /**
     * @return the energy of the current state.
     *
     */
    protected double stateEnergy() {
	double totalEnergy = 0.0;

	subsetEnergy = new double[noSubsets];

	for (int p = 0; p < noPairs; p++) {
	    for (int q = 0; q < p; q++) {

		if (subsets[p] == subsets[q]) {
		    subsetEnergy[subsets[p]] += energy[p][q];
		}
	    }
	}
	
	for (int c = 0; c < noSubsets; c++) {
	    totalEnergy += subsetEnergy[c];
	}
	
	return totalEnergy;
    }



    /**
     * @return the change in energy between the current state and the candidate state.
     */
    protected double generateCandidateState() {

	calculateCandidateSwapPositions();
	
	return swappedStateEnergy();
	
    }

    /**
     *
     * Calculates candidate state energy.
     *
     */
    private double swappedStateEnergy() {

	int iSubsetIfSwapped = subsets[j];
	int jSubsetIfSwapped = subsets[i];
	
	candidateSubsetEnergy = new double[noSubsets];
	
	double deltaEi = 0.0;
	double deltaEj = 0.0;
	
	double totalSubsetEnergy = 0.0;

	for (int c = 0; c < noSubsets; c++) {
	    candidateSubsetEnergy[c] = subsetEnergy[c];
	}
	
	for (int si = 0; si < noPairs; si++) {
	    
	    if (si != i && si != j) {

		if (subsets[si] == iSubsetIfSwapped) {
		    
		    candidateSubsetEnergy[subsets[si]] += energy[si][i] - energy[si][j]; 
		    
		}
		else if (subsets[si] == jSubsetIfSwapped) {

		    candidateSubsetEnergy[subsets[si]] += energy[si][j] - energy[si][i]; 

		}

	    }
	    
	} // end for si

	
	double totalCandidateSubsetEnergy = 0.0;

	for (int c = 0; c < noSubsets; c++) {
	    totalCandidateSubsetEnergy += candidateSubsetEnergy[c];
	}	

	return totalCandidateSubsetEnergy; 
    }

    

    /**
     * Changes the current state of the system to the last candidate.
     *
     */
    protected void acceptCandidateState() {

	// sanity check
	if (candidateSubsetEnergy == null) {
	    throw new LoggedException("Attempting to set an invalid state." +
				      "Please report this error to the Camino team");
	}

	int tmp = subsets[i];

	subsets[i] = subsets[j];

	subsets[j] = tmp;

	subsetEnergy = candidateSubsetEnergy;

	candidateSubsetEnergy = null;

	i = 0;
	j = 0;

    }
    

    /**
     * Saves the current state as the lowest energy state.
     *
     */
    protected void saveStateAsLowestEnergy() {
	bestEverSubsets = new int[noPairs];
	System.arraycopy(subsets, 0, bestEverSubsets, 0, noPairs);
	bestEverSubsetEnergy = new double[noSubsets];
	System.arraycopy(subsetEnergy, 0, bestEverSubsetEnergy, 0, noSubsets);
    }


    /**
     * Sets the current state to the lowest energy state discovered so far.
     *
     */
    protected void setStateToLowestEnergy() {

	System.arraycopy(bestEverSubsets, 0, subsets, 0, noPairs);
	System.arraycopy(bestEverSubsetEnergy, 0, subsetEnergy, 0, noSubsets);

    }


    /**
     * Gets a String representation of the current state of the system.
     */
    public String state() {
	return state(points, pairsPerSubset, subsets, temperature, subsetEnergy);
    }


    /**
     * Gets a String representation of the lowest energy state of the system.
     */
    public String lowestEnergyState() {
	return state(points, pairsPerSubset, bestEverSubsets, temperature, bestEverSubsetEnergy);
    }


   
    private String state(double[][] points, int[] pairsPerSubset, int[] subsets, double temp, 
			 double[] subsetEnergy) {

	String s = noPairs + " " + noSubsets + " ";

	double totalSubsetEnergy = 0.0;

	for (int c = 0; c < noSubsets; c++) {
	    s += pairsPerSubset[c] + " ";
	    totalSubsetEnergy += subsetEnergy[c];
	}	

	s += temperature + " " + totalSubsetEnergy + " ";

	for (int c = 0; c < noSubsets; c++) {
	    s += subsetEnergy[c] + " ";
	}

	s += "\n";
	
	for (int c = 0; c < noSubsets; c++) {
	    for (int n = 0; n < noPairs; n++) {
		if (subsets[n] == c) {
		    s += points[n][0] + " " + points[n][1] + " " + points[n][2] + "\n";
		}
	    }
	}

	return s;

    }


}

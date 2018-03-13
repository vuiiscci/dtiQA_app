package misc;

import numerics.*;
import optimizers.*;


import java.util.Random;
import java.util.Arrays;

/**
 * Superclass for classes that order electrostatic point sets.
 *
 * @author Philip Cook
 * @version $Id$
 * @see optimizers.SimulatedAnnealingOptimizer
 */
public abstract class OrderedAcqMinimizer extends SimulatedAnnealingOptimizer {

    protected final int noPairs;

    protected double[][] points;
    protected double[][] energy;

    protected final static double root2Times2 = Math.sqrt(2.0) * 2.0;


    /**
     * Initializes points in a random order.
     *
     *
     */
    public OrderedAcqMinimizer(int numberOfPairs, double coolRate, int iterationsPerT, 
				 String saveState, Random r) {
	
	// quit after 1E11 runs without lowering energy. This takes approximately 24 hrs on iMac 2GHz
	super(coolRate, iterationsPerT, (long)1E11, saveState, r);

	noPairs = numberOfPairs;

	if (noPairs < 7) {
	    throw new LoggedException("Can't optimize less than 7 pairs");
	}

	// get points from ElecPointSet
	points = SphericalPoints.getElecPointSet(noPairs);
	
	initializeRandom();
	
    }

    
    public OrderedAcqMinimizer(double[][] points, double coolRate, int iterationsPerT, 
				 String saveState, Random r) {

	// quit after 1E11 runs without lowering energy. This takes approximately 24 hrs on iMac 2GHz
	super(coolRate, iterationsPerT, (long)1E11, saveState, r);

	noPairs = points.length;
	this.points = points;
	
	buildEnergyMatrix();
    }


    public final void buildEnergyMatrix() {

	energy = pairEnergyMatrix(points);

    }

    
    public final static double[][] pairEnergyMatrix(double[][] points) {

	int noPairs = points.length;

	double[][] energy = new double[noPairs][noPairs];

	for(int i=0; i<noPairs; i++) {
	    
	    double[] tpi = Vector3D.thetaPhi(new Vector3D(points[i][0], points[i][1], points[i][2]));
	    
	    //Convert values in parameters array to theta and phi coords.
	    double ti = tpi[0];
	    double pi = tpi[1];
	    
	    //Precompute sin's and cos's.
	    double sti = Math.sin(ti);
	    double cti = Math.cos(ti);
	    //	    double spi = Math.sin(pi);
	    // double cpi = Math.cos(pi);
	 
	    
	    for(int j=i+1; j<noPairs; j++) {
		
	
		double[] tpj = Vector3D.thetaPhi(new Vector3D(points[j][0], points[j][1], points[j][2]));
		
		double tj = tpj[0];
		double pj = tpj[1];

		double stj = Math.sin(tj);
		double ctj = Math.cos(tj);
		//		double spj = Math.sin(pj);
		//		double cpj = Math.cos(pj);
				
		double cpimpj = Math.cos(pi - pj);
		//	double spimpj = Math.sin(pi - pj);
		//		double s2tj = Math.sin(2*tj);
				
		energy[i][j] = root2Times2 * (1.0 / Math.sqrt(1 - cti*ctj - cpimpj*sti*stj) + 1.0 / Math.sqrt(1 + cti*ctj + cpimpj*sti*stj));

		energy[j][i] = energy[i][j];
	    }

	}

	return energy;

    }

    
    /**
     * Randomly order the points.
     * 
     */
    public final void initializeRandom() {

	double[][] tmpPoints = new double[noPairs][];

	boolean[] added = new boolean[noPairs];

	for (int i = 0; i < noPairs; i++) {
	    
	    int p = ran.nextInt(noPairs);

	    while(added[p]) {
		p = ran.nextInt(noPairs);
	    }

	    tmpPoints[i] = points[p];
	    added[p] = true;
    
	}
	
	points = tmpPoints;
	buildEnergyMatrix(); // need to rebuild matrix for re-ordered points


    }


    /**
     * Swaps points i and j and updates the energy matrix.
     */
    protected void swapPoints(int i, int j) {

	double[] tmpPoint = points[i];

	points[i] = points[j];

	points[j] = tmpPoint;

	for (int r = 0; r < noPairs; r++) {
	    if (r == i || r == j) {
		continue;
	    }
	    
	    double tmp = energy[r][i];
	    energy[r][i] = energy[r][j];
	    energy[r][j] = tmp;
	    
	}

	for (int c = 0; c < noPairs; c++) {

	    if (c == i || c == j) {
		continue;
	    }

	    double tmp = energy[i][c];
	    energy[i][c] = energy[j][c];
	    energy[j][c] = tmp;
	}

    }


    public double[][] getPoints() {

	double[][] copy = new double[noPairs][3];

	for (int i = 0; i < noPairs; i++) {
	    copy[i][0] = points[i][0];
	    copy[i][1] = points[i][1];
	    copy[i][2] = points[i][2];
	}
	
	return copy;
    }


}

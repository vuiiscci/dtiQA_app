package tractography;

import imaging.*;
import inverters.*;
import misc.*;
import numerics.*;
import sphfunc.*;

import java.util.Random;
import java.util.logging.*;

/**
 * Samples orientations x from a Bayesian PDF P(x|data, model). Assumes a dirac-delta prior on the
 * model parameters other than orientation. Uses a partial volume "ball and stick" model of the data.
 *
 * @see data.BallStick
 *
 * @author Philip Cook
 * @version $Id$
 *
 *
 */
public class BayesDiracBallStickRandomizer extends BayesDiracRandomizer {

   /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.tractography.BayesDiracBallStickRandomizer");


    protected final BallStickInversion inv;

    /**
     * 
     */
    public BayesDiracBallStickRandomizer(float[][][][] data, DW_Scheme scheme, Random ran) {
	
	super(data, scheme, ran);

	inv = new BallStickInversion(scheme);

    }


    /**
     * 
     */
    public BayesDiracBallStickRandomizer(float[][][][] data, DW_Scheme scheme, int pointSetInd, Random ran) {
	
	super(data, scheme, pointSetInd, ran);

	inv = new BallStickInversion(scheme);

    }


    /**
     * Compute the likelihood function P(x | data) * P(x), where P(x) is an external prior,
     * if one has been set. By default, P(x) is a uniform distribution. 
     *
     */
    protected void computeLikelihood(int i, int j, int k) {

	likelihood[i][j][k] = new float[numVectors];

	double[] voxelData = new double[numMeas];

	for (int n = 0; n < numMeas; n++) { 
	    voxelData[n] = data[i][j][k][n];
	}

	double[] comps = inv.invert(voxelData);

	double lnS0 = comps[1];

	double s0 = Math.exp(comps[1]);

	double d = comps[2];

	double f = comps[3];

	// log data 
	double[] z = new double[numMeas];
	
	for (int n = 0; n < numMeas; n++) {
	    z[n] = Math.log(data[i][j][k][n]);
	}


	// predicted values from the model
	double[] mu = new double[numMeas];
	double[] lnMu = new double[numMeas];

	double[] lf = new double[numVectors];

	double sumLikelihood = 0.0;
	
	// compute likelihood of data given v, for all vectors
	for (int v = 0; v < numVectors; v++) {
	    
	    lf[v] = 1.0;
	    
	    // sigmaSq estimated from the residuals
	    double sigmaSq = 0.0;

	    // 5 parameters in the model
	    double norm = 1.0 / (numMeas - 5.0);
	    
	    for (int n = 0; n < numMeas; n++) {

		double zModel = 0.0;

		mu[n] = s0 * ( (1.0 - f) * Math.exp(-b[n] * d) + f * 
			       Math.exp(-b[n] * d * vDotG[n][v] * vDotG[n][v]) );

		lnMu[n] = Math.log(mu[n]);

		sigmaSq += norm * (data[i][j][k][n] - mu[n]) * (data[i][j][k][n] - mu[n]);
		
	    }
	    	    


	    double c = 1.0;
	    double sumExp = 0.0;

	    for (int n = 0; n < numMeas; n++) {

		c *= (mu[n] / Math.sqrt(2.0 * Math.PI * sigmaSq));
		
		sumExp += (-1.0 * (mu[n] * mu[n]) * (z[n] - lnMu[n]) * (z[n] - lnMu[n]) / (2.0 * sigmaSq));
	    }	

	    // trying to save on calls to Math.exp
	    lf[v] = c * Math.exp(sumExp);


	    if (externalPriors != null) {
		lf[v] = lf[v] * externalPriors.pdf(i,j,k, 0, vectors[v]);
	    }
  
	    if (lf[v] > 0.0) { 
		sumLikelihood += lf[v];
	    }

	}

//  	if (!(sumLikelihood > 0.0)) {
// 	    logger.warning("Unable to determine likelihood for voxel " + i + " " + j + " " + k + ". Returning uniform distribution");
// 	}
	
	for (int v = 0; v < numVectors; v++) {
	    if (sumLikelihood > 0.0) {
		likelihood[i][j][k][v] = (float)(lf[v] / sumLikelihood);
	    }
	    else {
		likelihood[i][j][k][v] = (float)(1.0 / (4.0 * Math.PI));
	    }
	}
		 	
    }
    

 

}

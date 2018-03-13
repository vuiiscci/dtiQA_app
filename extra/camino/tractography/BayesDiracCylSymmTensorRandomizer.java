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
 * model parameters other than orientation. This class uses a Gaussian PDF with a cylindrically-symmetric 
 * diffusion tensor to model the data.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class BayesDiracCylSymmTensorRandomizer extends BayesDiracRandomizer {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.tractography.BayesDiracCylSymmTensorRandomizer");


    protected final DT_Inversion inv;


    /**
     * 
     */
    public BayesDiracCylSymmTensorRandomizer(float[][][][] data, DW_Scheme scheme, Random ran) {
	
	super(data, scheme, ran);

	inv = new LinearDT_Inversion(scheme);


    }


    /**
     * 
     */
    public BayesDiracCylSymmTensorRandomizer(float[][][][] data, DW_Scheme scheme, int pointSetInd, Random ran) {
	
	super(data, scheme, pointSetInd, ran);

	inv = new LinearDT_Inversion(scheme);


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
	
	DT d = new DT(comps[2], comps[3], comps[4], comps[5], comps[6], comps[7]);

	double[][] seig = d.sortedEigenSystem();

	double l1 = seig[0][0];

	// cylindrically symmetric tensor has l2 == l3
	double alpha = (seig[0][1] + seig[0][2]) / 2.0;

	double beta = l1 - alpha;

	// D = \beta e_1 e_1^T + \alpha I

	// log data 
	double[] z = new double[numMeas];
	
	for (int n = 0; n < numMeas; n++) {
	    z[n] = Math.log(data[i][j][k][n]);
	}


	// predicted values from the model
	double[] mu = new double[numMeas];
	double[] lnMu = new double[numMeas];

	// sigmaSq estimated from the residuals
	double sigmaSq = 0.0;

	// 5 parameters in the model: \theta, \phi (of v), \alpha, \beta, S_0
	double norm = 1.0 / (numMeas - 5.0);

	// estimate sigmaSq from variance of OLS residuals
	// could potentially use WLS
	Vector3D e1 = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);

	for (int n = 0; n < numMeas; n++) {
	    
	    double e1DotG = e1.dot(g[n]);

	    mu[n] = Math.exp(lnS0 - b[n] * (alpha + beta * e1DotG * e1DotG));
	    
	    sigmaSq += norm * (data[i][j][k][n] - mu[n]) * (data[i][j][k][n] - mu[n]);

	}

	double sqrt2PISigmaSq = Math.sqrt(2.0 * Math.PI * sigmaSq);

	// un-normalized likelihood
	double[] lf = new double[numVectors];

	// normalization constant
	double sumLikelihood = 0.0;
	
	// compute likelihood of data given v, for all vectors
	for (int v = 0; v < numVectors; v++) {
	    
	    lf[v] = 1.0;
	    
	    
	    for (int n = 0; n < numMeas; n++) {
		// Friman gives \mu = \mu_0 \exp[-\alpha b] \exp[-\beta b (g^T v)^2]
		//        so ln(\mu) = ln(\mu_0 \exp[-\alpha b -\beta b (g^T v)^2])
		//                   = ln(\mu_0) -\alpha b -\beta b (g^T v)^2
		//                   = ln(\mu_0) -b (\alpha + beta (g^T v)^2)
		
		double zModel = 0.0;
		
		lnMu[n] = lnS0 - b[n] * (alpha + beta * vDotG[n][v] * vDotG[n][v]);
		
		
		mu[n] = Math.exp(lnMu[n]);
	
	    }
	    	  
	    double c = 1.0;
	    double sumExp = 0.0;

	    for (int n = 0; n < numMeas; n++) {

		c *= (mu[n] / sqrt2PISigmaSq);
		
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

// 	if (!(sumLikelihood > 0.0)) {
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

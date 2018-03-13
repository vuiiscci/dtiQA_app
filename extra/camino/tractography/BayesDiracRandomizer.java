package tractography;

import imaging.*;
import misc.*;
import numerics.*;
import sphfunc.*;

import java.util.logging.*;
import java.util.Random;


/**
 * Samples orientations x from a Bayesian PDF P(x|data, model). Assumes a dirac-delta prior on the
 * model parameters other than orientation.
 *
 * @author Philip Cook
 * @version $Id$
 */
public abstract class BayesDiracRandomizer {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.tractography.BayesDiracRandomizer");

    protected Vector3D[] vectors;

    protected Vector3D[] negVectors;

    protected int numVectors;

    protected final float[][][][] data;

    protected final int xDataDim;
    protected final int yDataDim;   
    protected final int zDataDim;

    protected final DW_Scheme scheme;

    /** DW measurements in the imaging scheme */
    protected final int numMeas;

    /** q / |q|  */
    protected final Vector3D[] g;

    /** b-values  */
    protected final double[] b;

    /** dot product of each vector vectors[n] with g[n]  */
    protected final double[][] vDotG;

    /** curvature priors */
    protected double[][] curvePrior = null;

    // concentration of curvature prior, which is a Watson distribution
    // centred on the previous direction
    private double curvePriorKappa = 0.0;
    
    // concentration of curvature prior term, which is |v_i \dot v_{i-1}|^\gamma
    private double curvePriorGamma = 0.0;

    private boolean useCurvePrior = false;


    /** 
     * likelihood[i][j][k] = likelihood of each of the vectors given the data and any external prior
     * the curvature prior is not included in this calculation  
     */
    protected float[][][][] likelihood;
    
    // keep a counter of functions calculated. The program will certainly run out of memory if too many
    // voxel likelihood functions are calculated, since they occupy approx 2Kb each
    // if this gets too large, we can clear out the cache
    private int functionsCalculated = 0;

    // clean out the cache after this many functions are stored
    private final int FUNCTION_CACHE_LIMIT;


    /** external priors. */
    protected PICoRandomizer externalPriors = null;

  
    protected final Random ran;


    // if we know the last randomized vector, we can avoid searching for it
    private int lastRandomizedIndex = 0;



    public BayesDiracRandomizer(float[][][][] data, DW_Scheme scheme, Random ran) {
	// default point set 1 (1922 directions)
	this(data, scheme, 1, ran);
    }

    public BayesDiracRandomizer(float[][][][] data, DW_Scheme scheme, int pointSetIndex, Random ran) {
	
	this.data = data;
	this.scheme = scheme;
	this.ran = ran;

	initializeSphericalPointSet(pointSetIndex);

	// set cache size based on how much memory we have and how big the data is
	long maxMemoryLeft = Runtime.getRuntime().maxMemory() - (Runtime.getRuntime().totalMemory() - 
								 Runtime.getRuntime().freeMemory());
	// don't let the cache exceed half the free space
	FUNCTION_CACHE_LIMIT = (int)(maxMemoryLeft / (2 * 4 * vectors.length));

	logger.info("Allocating space to cache " + FUNCTION_CACHE_LIMIT + " likelihood functions");

	xDataDim = data.length;
	yDataDim = data[0].length;
	zDataDim = data[0][0].length;

	numMeas = data[0][0][0].length;

	g = new Vector3D[numMeas];

	b = new double[numMeas];

	// dot product of gradient vectors with PDF vectors
	vDotG = new double[numMeas][numVectors];

	for (int n = 0; n < numMeas; n++) {

	    b[n] = scheme.getB_Value(n);

	    if (b[n] > 0.0) {
		g[n] = new Vector3D(scheme.getG_Dir(n));
		
		for (int v = 0; v < numVectors; v++) {
		    vDotG[n][v] = g[n].dot(vectors[v]);
		}
	    }
	    else {
		g[n] = new Vector3D(1.0, 0.0, 0.0);
	    }

	}

	likelihood = new float[xDataDim][yDataDim][zDataDim][];

    }


    /**
     * Compute the likelihood function P(x | data) * P(x), where P(x) is an external prior.
     * If there is no external prior, P(x) is uniform. The curvature prior term is handled 
     * separately. The likelihood should be normalized to sum to 1 over all vectors.
     * 
     *
     */
    protected abstract void computeLikelihood(int i, int j, int k);


     
    /**
     * Samples from the likelihood function in this voxel, does not consider the prior fibre orientation.
     */
    public final Vector3D[] getRandomizedPDs(int i, int j, int k) {

	if (likelihood[i][j][k] == null) {

	    functionsCalculated++;

	    if (functionsCalculated > FUNCTION_CACHE_LIMIT) {
		clearCache();
	    }

	    computeLikelihood(i,j,k);

	}

	float r = ran.nextFloat();

	float[] cdf = new float[numVectors];
	
	// return vector n where cdf[n-1] < targetValue <= cdf[n]
	for (int n = 0; n < numVectors; n++) {

	    float previous = n > 0 ? cdf[n-1] : 0.0f;
	    
	    // add a small epsilon to the CDF, because sometimes it does not sum
	    // exactly to 1, resulting in no vector being chosen. To avoid this, 
	    // we expand the likelihood of each term by 0.002%
	    cdf[n] = previous + likelihood[i][j][k][n] * 1.00002f;

	    if (r <= cdf[n]) {
		lastRandomizedIndex = n;
		return new Vector3D[] {vectors[n]};
	    }

	}
            
 	// If we still have a small error, return most likely vector
 	if (1.0 - r < 1E-5) {

            logger.warning("Sampling failed, returning maximum likelihood vector for this iteration in voxel " + i + " " + j + " " + k);

	    // inefficiently search for ML vector but that's OK because this happens
	    // less than once in 10^5 samples
	    double ml = 0.0;
	    int mlInd = 0;
	    
	    for (int n = 0; n < numVectors; n++) {
		if (likelihood[i][j][k][n] > ml) {
		    ml = likelihood[i][j][k][n];
		    mlInd = n;
		}
	    }
	    
	    lastRandomizedIndex = mlInd;
 	    return new Vector3D[] {vectors[mlInd]};
 	}

	String cdfString = "";

	for (int n = 0; n < numVectors; n++) {
	    cdfString += n + "\t" + cdf[n] + "\n";
	}

	System.err.println("CDF: " + cdfString);

	throw new LoggedException("target value (" + r + ") is greater than CDF of all vectors in " +
				  "voxel " + i + " " + j + " " + k);

    }


  
    /**
     * Samples from the posterior distribution in this voxel. The prior distribution for
     * fibre orientations is a Watson distribution centred on <code>previousDir</code>, which
     * should be a vector returned by the getRandomizedPDs method of this class.
     *
     * 
     */
    public final Vector3D[] getRandomizedPDs(int i, int j, int k, Vector3D previousDir) {

	if (!useCurvePrior) {
	    return getRandomizedPDs(i,j,k);
	}

	if (likelihood[i][j][k] == null) {
	    functionsCalculated++;

	    if (functionsCalculated > FUNCTION_CACHE_LIMIT) {
		clearCache();
	    }

	    computeLikelihood(i,j,k);

	}

	float r = ran.nextFloat();

	float[] cdf = new float[numVectors];

	int prevIndex = lastRandomizedIndex;

	if (!previousDir.equals(vectors[prevIndex]) && !previousDir.equals(negVectors[prevIndex])) {

	    // need to search for vector
	    prevIndex = -1;
	    
	    findVector: 
	    for (int n = 0; n < numVectors; n++) { 
		if (vectors[n].equals(previousDir) || negVectors[n].equals(previousDir)) {
		    prevIndex = n;
		    break findVector;
		}
	    }
	}

	if (prevIndex > -1) {
	    
	    for (int n = 0; n < numVectors; n++) {
		float posterior = (float)(likelihood[i][j][k][n] * curvePrior[prevIndex][n]);
		
		float previous = n > 0 ? cdf[n-1] : 0.0f;
		
		cdf[n] = previous + posterior;	      
	    }
	}
	else {
	    // calculating the prior for all vectors given an arbitrary previousDir is very slow
	    // so we restrict ourselves to using directions in the vectors array

	    throw new IllegalArgumentException("Previous direction is not one of the PDF samples");
	}

	// cdf no longer normalized as we multiply (normalized) likelihood by curvature prior
	float targetValue = r * cdf[numVectors - 1];

	// return vector n where cdf[n-1] < targetValue <= cdf[n]
	// potential for a binary search to speed things up here
	for (int n = 0; n < numVectors; n++) {
	    if (targetValue <= cdf[n]) {
		lastRandomizedIndex = n;

		if (vectors[n].dot(previousDir) < 0.0) {
		    return new Vector3D[] {negVectors[n]};
		}
		else {
		    return new Vector3D[] {vectors[n]};
		}
	    }
	}


	String cdfString = "";

	for (int n = 0; n < numVectors; n++) {
	    cdfString += n + "\t" + cdf[n] + "\n";
	}
	
	
	System.err.println("CDF: " + cdfString);

	throw new LoggedException("target value (" + targetValue + ") is greater than CDF of all vectors in " +
				  "voxel " + i + " " + j + " " + k);
    }



    /**
     * Curvature prior is a Watson distribution centred on the previous direction. 
     * Prior concentration must be greater than zero. Set zero concentration to remove the curvature 
     * prior.  
     * <p>
     * Calling this method removes any previous curvature prior.
     */
    public void setCurvatureK(double k) {

	if (!(k >= 0.0)) {
	    throw new LoggedException("Can only use positive concentration for curve prior");
	}

	curvePriorKappa = k;

	if (k == 0.0) {
	    useCurvePrior = false;
	    curvePrior = null;
	    return;
	}

	useCurvePrior = true;

	curvePriorKappa = k;
	curvePriorGamma = 0.0;

	curvePrior = new double[numVectors][numVectors];

	for (int m = 0; m < numVectors; m++) {
	    for (int n = 0; n <= m; n++) {
		curvePrior[m][n] = WatsonDistribution.pdf(vectors[m], vectors[n], curvePriorKappa); 
		curvePrior[n][m] = curvePrior[m][n];
	    }
	}

    }

    
    /**
     * Sets the improper curvature prior |v_i \dot v_{i-1}|^\gamma, where v_i is the current
     * orientation and v_{i-1} is the previous orientation. The value of gamma should be positive.
     *
     * Calling this method removes any previous curvature prior. 
     */
    public void setCurvatureG(double g) {

	if (!(g >= 0.0)) {
	    throw new LoggedException("Can only use positive concentration for curve prior");
	}

	curvePriorGamma = g;

	if (g == 0.0) {
	    useCurvePrior = false;
	    curvePrior = null;
	    return;
	}

	useCurvePrior = true;

	curvePriorGamma = g;
	curvePriorKappa = 0.0;

	curvePrior = new double[numVectors][numVectors];

	for (int m = 0; m < numVectors; m++) {
	    for (int n = 0; n <= m; n++) {
		curvePrior[m][n] = Math.pow(Math.abs(vectors[m].dot(vectors[n])), g);
		curvePrior[n][m] = curvePrior[m][n];
	    }
	}

    }


    /**
     * Set external priors. Should have one PDF defined per voxel (if not, only the first
     * one is used).
     *
     */
    public void setExternalPriors(PICoRandomizer randomizer) {
	externalPriors = randomizer;
    }


    private void clearCache() {

	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    likelihood[i][j][k] = null;
		}
	    }
	}

	functionsCalculated = 0;

    }


    /**
     * Approximation of E^x. Implementation of C algorithm from 
     * Schraudolph, "A Fast, Compact Approximation of the Exponential Function",
     * Neural Computation 11:863-870 (1999). It appears to be correct to within
     * +/- 4%.
     */
    protected static double exp(double val) {
	final long tmp = (long) (1512775 * val + (1072693248 - 60801));
	return Double.longBitsToDouble(tmp << 32);
    }


    
    public void initializeSphericalPointSet(int index) {
	
	try {
	    double[][] points = ISCodes.getPointSet(index).data;
	   
	    numVectors = points.length;
	    
	    vectors = new Vector3D[numVectors];
	    negVectors = new Vector3D[numVectors];
	    
	    for (int i = 0; i < points.length; i++) {
		vectors[i] = new Vector3D(points[i][0], points[i][1], points[i][2]);
		negVectors[i] = new Vector3D(-points[i][0], -points[i][1], -points[i][2]);
	    }
	}
	catch (ISCodesException e) {
	    throw new LoggedException(e);
	}

    }


}

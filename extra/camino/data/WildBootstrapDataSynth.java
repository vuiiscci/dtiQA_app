package data;

import java.io.*;

import imaging.*;
import misc.*;
import numerics.*;
import tools.*;

import java.util.Random;


/**
 * Provides wild bootstrap data. For information on the technique, see B Whitcher, DS Tuch, L Wang, 
 * "The wild bootstrap to quantify variability in diffusion tensor MRI", Proc ISMRM 1333, (2005).
 * 
 *
 * 
 * @author Philip Cook
 * @version $Id$
 * 
 */
public abstract class WildBootstrapDataSynth implements DataSource {

    private final DataSource rawDataSource;

    // samples to generate per voxel
    private final int samples;

    // components of data in each voxel - will not match the number of measurements
    // in the scheme if data is normalized
    private int components;

    private int samplesGenerated = 0;


    // constant used in heteroscedasticity term
    // should be either HC2 = -0.5, or HC3 = -1.0
    //
    // See Flachaire, "Bootstrapping Heteroskedasticity Consistent Covariance Matrix Estimator", 
    // Computational Statistics 17:501-506 (2002)
    // See also Davidson and Flachaire, "The wild bootstrap, tamed at last", Working paper, 
    // IER#1000, Queen's University.
    private double alpha = -1.0;

    private final Random ran;


    /**
     * If true, fit the model parameters to the log of the source data. 
     *
     */
    private boolean logData = false;


    /**
     * If true, fit the model parameters to the normalized data.
     *
     */
    private boolean normalizeData = false;


    // BEGIN VOXEL-SPECIFIC DATA

    /**
     * Residuals \epsilon_i where y_i = X_i \beta_i + \epsilon_i, for measurement i out of N.
     *
     */
    private RealMatrix epsilon;


    /**
     * weightedResiduals[i] = epsilon[i] * 1.0 / (1.0 - h_i)^\alpha 
     *
     */
    private double[] weightedResiduals;

    
    /**
     * mu[i] = (X_i \beta_i)
     * <p>
     * A wild bootstrap sample is a combination of this and the residuals.
     */
    private RealMatrix mu;


    /**
     * mean unweighted signal in the input data.
     */
    private double meanS0;


    // END VOXEL-SPECIFIC DATA

    
    /**
     * a wild bootstrap sample has the form y_i = mu_i + \epsilon_i * r / (1 - h_i)^alpha,
     * where r is a random number (-1 or 1 in our implementation). 
     * <p>
     * h_i = H(i,i) = X_i (X^TX)^{-1} X_i^T, where X_i is the row of the design matrix corresponding to 
     * measurement i. H = X (X^TX)^{-1} X^T
     * <p>
     * According to Chung et al, H should have a weighting term, but Flachaire doesn't mention it. 
     * 
     */
    protected RealMatrix H;


    /**
     * y = X \beta + \epsilon, where y is the data vector, \beta is a vector containing the parameters 
     * of the model, and \epsilon is a vector of the residuals. 
     * 
     */
    protected RealMatrix X;


    /**
     * linearInv = L = Inverse or pseudo inverse of X, such that L * y = \beta. 
     * 
     * 
     * 
     */
    protected RealMatrix linearInv;



    /**
     * true when we've got the last of the source data.
     */
    private boolean reachedEndOfFile;


    /**
     * true when there are no more voxels to return, ie after (samples) calls
     * to nextVoxel using the last of the source data.
     */
    private boolean noMoreData;


    /**
     * Scheme corresponding to input data.
     *
     */
    protected DW_Scheme ip; 


    /**
     * zero[i] == true if b=0 at measurement i according to the scheme
     */
    private boolean[] zero;

    /**
     * Number of measurements in scheme
     */
    private int numMeas;


    /**
     * 
     * @param input
     *            A data source that should provide voxel-order data. If <code>null</code>, 
     *            the user must provide data with the <code>setSourceData</code> method.
     *
     * @param scheme
     *            The <code>DW_Scheme</code> object with the details of the scanner
     *            sequence of the raw data.
     * 
     * @param samples
     *            The number of bootstrap samples to generate from each voxel. Set to -1 for
     *            unlimited samples.
     * 
     */
    public WildBootstrapDataSynth(DataSource input, DW_Scheme scheme, int samples) {

	this(input, scheme, samples, new numerics.MTRandom(2350));
    }



    /**
     * 
     * @param input
     *            A data source that should provide voxel-order data. If <code>null</code>, 
     *            the user must provide data with the <code>setSourceData</code> method.
     *
     * @param scheme
     *            The <code>DW_Scheme</code> object with the details of the scanner
     *            sequence of the raw data.
     * 
     * @param samples
     *            The number of bootstrap samples to generate from each voxel. Set to -1 for
     *            unlimited samples.
     *
     * @param seed 
     *            Seed for the random number generator.
     * 
     */
    public WildBootstrapDataSynth(DataSource input, DW_Scheme scheme, int samples, int seed) {

	this(input, scheme, samples, new numerics.MTRandom(seed));
    }



    /**
     * This is the constructor that actually does everything, so ensure that this eventually gets
     * called. 
     *
     * @param input
     *            A data source that should provide voxel-order data. If <code>null</code>, 
     *            the user must provide data with the <code>setSourceData</code> method.
     *
     * @param scheme
     *            The <code>DW_Scheme</code> object with the details of the scanner
     *            sequence of the raw data.
     * 
     * @param samples
     *            The number of bootstrap samples to generate from each voxel. Set to -1 for
     *            unlimited samples.
     * 
     * @param ran 
     *          A random number generator.
     */
    public WildBootstrapDataSynth(DataSource input, DW_Scheme scheme, int samples, Random ran) {

	rawDataSource = input;
	ip = scheme;

	numMeas = ip.numMeasurements();
	
	// default, this will change if normalizeData() is called
	components = numMeas; 

	zero = new boolean[numMeas];

	meanS0 = 0.0;


	this.samples = samples;
	this.ran = ran;    
    }



    /**
     * Intializes the object. Should be called by subclass constructors. 
     * 
     */
    protected void init() {

	// Initialise the flags indicating end of input file and whether
        // more voxel data can be requested.
        reachedEndOfFile = false;

	initializeReconSpecificParams();

	for (int i = 0; i < numMeas; i++) {
	    if (ip.zero(i)) {
		zero[i] = true;
	    }
	}

        // Read in the first voxel ready for the first data request.

	if (rawDataSource != null) {
	
	    try {
		getNextSourceVoxel();
	    }
	    catch (DataSourceException e) {
		throw new LoggedException(e);
	    }
	}
	else {
	    // reset to false when data is given
	    noMoreData = true;
	}
	
    }


    /**
     * Initialize X, H, linearInv.  If the regression is to be performed on the normalized 
     * measurements, X, H and linearInv should operate on normalized data (ie, the raw data with
     * the b=0 measurements removed, and other measurements divided by the mean b=0 measurement).
     * If the regression is to be performed on log data, then X, H and linearInv should be set
     * accordingly. 
     *
     */
    protected abstract void initializeReconSpecificParams();


    /**
     * Gets the next bootstrap sample of the raw data from this data source.
     *
     * @return the next bootstrap sample from the data source. 
     * @throws DataSourceException if no more voxels are available.
     *
     */
    public final double[] nextVoxel() {

	if (noMoreData) {
	   throw new DataSourceException("No more voxels in data source.");
	}

	double[] bootstrap = getBootstrapSample();

	if (samples > -1) {
	    samplesGenerated++;
	}

        // The previous read may have been the last voxel, in which case
        // the current read failed indicating that this is the last
        // voxel's worth of data to return.
	if (samplesGenerated == samples) {

	    if (reachedEndOfFile || rawDataSource == null) {
		noMoreData = true;
	    }
	    else {
		getNextSourceVoxel();
	    }
        }

        // Return the current voxel's data.
        return bootstrap;
    }



    /**
     * Puts the data for the current voxel into an array ready for return by
     * nextVoxel.
     * 
     * @return An array containing the voxel data.
     */
    private double[] getBootstrapSample() {

	double[] bootstrap = new double[components];

	for (int i = 0; i < components; i++) {
	    
	    double t = ran.nextBoolean() ? 1.0 : -1.0;

	    bootstrap[i] = mu.entries[i][0] + t * weightedResiduals[i];

	    // if the model fits to the log data, exponentiate here
	    if (logData) {
		bootstrap[i] = Math.exp(bootstrap[i]);
	    }
	}
	
	if (normalizeData) {
	    return unNormalizeData(bootstrap, meanS0);
	}

	return bootstrap;
    }
    

    public boolean more() {
        return !noMoreData;
    }


    
    /**
     * Reads the next voxel's worth of source data into the array next. If the read
     * fails, the method sets a flag to indicate the end of the data file. If
     * the read fails part way through, the method reports a separate error
     * indicating that the file did not contain a whole number of voxels.
     *
     * @throws DataSourceException if the data cannot be read.
     */
    private void getNextSourceVoxel() {

        if (reachedEndOfFile) {
            noMoreData = true;
            throw new DataSourceException("No more voxels in data source.");
        }
	
	double[][] data = new double[1][];

	data[0] = rawDataSource.nextVoxel();

	setSourceData(data);

	if (!rawDataSource.more()) {
	    reachedEndOfFile = true;
	}

    }


    /**
     * Set the source data directly. 
     *
     */
    public void setSourceData(float[] inputData) {
	double[][] data = new double[1][inputData.length];

	// Slight waste of time copying float to double, but there is not much choice,
	// since RealMatrix needs doubles

	// can't even use System.arraycopy with arrays of different types

	// the real performance bottleneck is the creation of RealMatrices

	final int length = inputData.length;

	for (int i = 0; i < length; i++) {
	    data[0][i] = inputData[i];
	}

	setSourceData(data);

    }


    /**
     * Set the source data directly. 
     *
     */
    public void setSourceData(double[] inputData) {
	double[][] data = new double[1][inputData.length];
	
	System.arraycopy(inputData, 0, data[0], 0, inputData.length);

	setSourceData(data);

    }


    /**
     * Set the source data directly. Note that the array will be modified if this object uses
     * log or normalized data.
     *
     */
    private void setSourceData(double[][] data) {

	noMoreData = false;
	samplesGenerated = 0;

	double meanS0 = ip.geoMeanZeroMeas(data[0]);

	if (meanS0 == 0.0) {
	    meanS0 = 1.0;
	}

	if (normalizeData) {
	    data[0] = ip.normalizeData(data[0], meanS0);
	}
	if (logData) {
	    data[0] = getLogData(data[0]);
	}

	
	// calculate mu
	RealMatrix y = new RealMatrix(data).transpose();

	RealMatrix beta = linearInv.product(y); 

	mu = X.product(beta);

	epsilon = y.sub(mu);

	calculateWeightedResiduals();

    }


    /**
     * Calculate the weighted residuals using epsilon and H.
     *
     */
    private void calculateWeightedResiduals() {
	weightedResiduals = new double[components];

	for (int i = 0; i < components; i++) {
	    weightedResiduals[i] = epsilon.entries[i][0] / Math.pow(1.0 - H.entries[i][i], alpha);
	}
	
    }
    

    /**
     * @return the log of the raw data.
     */
    private double[] getLogData(double[] rawData) {

	double[] proc = new double[rawData.length];

	for (int i = 0; i < rawData.length; i++) {
	    if (rawData[i] > 0.0) {
		proc[i] = Math.log(rawData[i]);
	    }
	    else {
		proc[i] = 0.0;
	    }
	}

	return proc;
    }

    /**
     * Re-integrates the zero measurements into the data. Sets all zero measurements to the
     * estimated unweighted measurement s0, and scales other measurements by that value.
     */
    private double[] unNormalizeData(double[] data, double s0) {
	double[] unNormedData = new double[numMeas];
	
	int counter = 0;

	for (int i = 0; i < numMeas; i++) {
	    if (zero[i]) {
		unNormedData[i] = s0;
	    }
	    else {
		unNormedData[i] = s0 * data[counter++];
	    }
	} 

	return unNormedData;
    }

    /**
     * Does the model fitting and bootstrapping on the log of the raw data. 
     *
     */
    protected void useLogData() {
	logData = true;
    }


    /**
     * Does the model fitting and bootstrapping on the normalized data. b=0 data will not be bootstrapped,
     * the b=0 signal will be fixed. If this is called, then the expected number of components in the data 
     * vector y will be N - M, where N is the number of measurements in the scheme and M is the number
     * of zero measurements.
     * <p>
     * This capability is untested.
     */
    protected void useNormalizedData() {
	normalizeData = true;

	components = ip.numMeasurements() - ip.numZeroMeasurements();
    }



    /**
     * Set the constant used in heteroscedasticity term: HC2 = -0.5, or HC3 = -1.0.
     *
     * See Flachaire, "Bootstrapping Heteroskedasticity Consistent Covariance Matrix Estimator", 
     * Computational Statistics 17:501-506 (2002)
     *
     */
    public void setHC2() {
	alpha = -0.5;
	calculateWeightedResiduals();
    }


    /**
     * Set the constant used in heteroscedasticity term: HC2 = -0.5, or HC3 = -1.0. This is the default.
     *
     * See Flachaire, "Bootstrapping Heteroskedasticity Consistent Covariance Matrix Estimator", 
     * Computational Statistics 17:501-506 (2002)
     *
     */
    public void setHC3() {
	alpha = -1.0;
	calculateWeightedResiduals();
    }

}

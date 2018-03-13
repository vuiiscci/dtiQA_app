package data;

import java.io.*;

import misc.*;
import tools.*;

import java.util.Random;


/**
 * Source of bootstrap data. Given voxel-order data from a series of files, the
 * class returns a specified number of bootstrap samples of each source voxel.
 *
 *
 * 
 * @author Philip Cook
 * @version $Id$
 * 
 */
public class BootstrapDataSynthesizer implements DataSource {

    private DataSource[] rawDataSources;

    // repeats of the measurements
    protected int repeats;

    // samples to generate per voxel
    protected int samples;

    // components of data in each voxel
    protected int components;

    protected int samplesGenerated = 0;

    protected Random ran;


    /**
     * An array that always contains the data for the next voxel.
     */
    protected double[][] next;


    /**
     * true when we've got the last of the source data.
     */
    protected boolean reachedEndOfFile;


    /**
     * true when there are no more voxels to return, ie after (samples) calls
     * to nextVoxel using the last of the source data.
     */
    protected boolean noMoreData;


    /**
     * Default constructor for subclasses.
     *
     */
    protected BootstrapDataSynthesizer() {

    }


    /**
     * Constructor for subclasses just initializes the number of repeats.
     * 
     * @param repeats
     *            The number of repeats of the data.
     * 
     *
     */
    protected BootstrapDataSynthesizer(int repeats) {
	this.repeats = repeats;
    }


    /**
     * 
     * @param filenames
     *            The name of the data files containing the repeats.
     * 
     * @param components
     *            The number of values in each voxel.
     *
     * @param samples
     *            The number of bootstrap samples to generate from each voxel.
     * 
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     */
    public BootstrapDataSynthesizer(String[] filenames, int components, int samples, String type) {

	this(filenames, components, samples, type, new numerics.MTRandom(2350));
    }


    /**
     * 
     * @param filenames
     *            The name of the data files containing the repeats.
     * 
     * @param components
     *            The number of values in each voxel.
     *
     * @param samples
     *            The number of bootstrap samples to generate from each voxel.
     *
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     *
     * @param seed 
     *          Seed for the random number generator.
     */
    public BootstrapDataSynthesizer(String[] filenames, int components, int samples, String type, int seed) {
	this(filenames, components, samples, type, new numerics.MTRandom(seed));
    }


    /**
     * 
     * @param filenames
     *            The name of the data files containing the repeats.
     * 
     * @param components
     *            The number of values in each voxel.
     *
     * @param samples
     *            The number of bootstrap samples to generate from each voxel.
     *
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     *
     * @param ran 
     *          A random number generator.
     */
    public BootstrapDataSynthesizer(String[] filenames, int components, int samples,  String type, Random ran) {

	initDataSource(filenames, components, type);

	init(samples, components, ran);

    }



    // call this before init, because init gets data
    private void initDataSource(String[] filenames, int components, String type) {

	repeats = filenames.length;

	// Don't really want a bunch of 24 meg buffers
	// doing this makes the method not thread safe, because the buffer size
	// is static
	
	int tmp = ExternalDataSource.FILEBUFFERSIZE;

	ExternalDataSource.FILEBUFFERSIZE = 8*1024*1024;

	rawDataSources = new DataSource[repeats];

	for (int i = 0; i < repeats; i++) {
	    rawDataSources[i] = ExternalDataSource.getDataSource(filenames[i], components, type);
	}

	ExternalDataSource.FILEBUFFERSIZE = tmp;
    }


    /**
     * Does the work of the constructor.
     * 
     * @param samples
     *            The number of bootstrap samples to generate from each voxel.
     *
     * @param components
     *            The number of components in the data.
     *
     * @param ran 
     *          A random number generator.
     * 
     */
    protected void init(int samples, int components, Random r) {

	this.components = components;
	this.samples = samples;
	ran = r;
	
        next = new double[repeats][components];

        // Initialise the flags indicating end of input file and whether
        // more voxel data can be requested.
        reachedEndOfFile = false;
        noMoreData = false;

        // Read in the first voxel ready for the first data request.
        try {
            getNextSourceVoxel();
        }
        catch (DataSourceException e) {
            throw new LoggedException(e);
        }
    }

    /**
     * @return the next bootstrap sample from the data source. 
     *
     */
    public double[] nextVoxel() throws DataSourceException {

	double[] bootstrap = getBootstrapSample();
	samplesGenerated++;

        // The previous read may have been the last voxel, in which case
        // the current read failed indicating that this is the last
        // voxel's worth of data to return.
	if (samplesGenerated == samples) {

	    if (reachedEndOfFile) {
		noMoreData = true;
	    }
	    else {
		getNextSourceVoxel();
		samplesGenerated = 0;
	    }
        }

        // Return the current voxel's data.
        return bootstrap;
    }


    public boolean more() {
        return !noMoreData;
    }


    /**
     * Puts the data for the current voxel into an array ready for return by
     * nextVoxel.
     * 
     * @return An array containing the voxel data.
     */
    protected double[] getBootstrapSample() {

	double[] bootstrap = new double[components];

	for (int i = 0; i < components; i++) {

	    bootstrap[i] = next[ran.nextInt(repeats)][i];
	}

	return bootstrap;
    }


    /**
     * Reads the next voxel's worth of source data into the array next. If the read
     * fails, the method sets a flag to indicate the end of the data file. If
     * the read fails half way through, the method reports a separate error
     * indicating that the file did not contain a whole number of voxels.
     */
    protected void getNextSourceVoxel() throws DataSourceException {

        if (reachedEndOfFile) {
            noMoreData = true;
            throw new DataSourceException("No more voxels in data source.");
        }

        // Need this in the Exception handling.
        int indexReached = 0;
	
	try {
	    
	    for (int r = 0; r < repeats; r++) {
		
		indexReached = r;
		
		next[r] = rawDataSources[r].nextVoxel();
		
		if (!rawDataSources[r].more()) {
		    reachedEndOfFile = true;
		}
	    }
	}
	catch (DataSourceException e) {
	    
	    reachedEndOfFile = true;
	    
	    noMoreData = true;
	    throw new DataSourceException(
					  "End of file reached without completing voxel in data source " + 
					  indexReached);
	}

    }



}

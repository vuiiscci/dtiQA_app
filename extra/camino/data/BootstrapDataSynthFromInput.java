package data;

import java.io.*;

import imaging.*;
import misc.*;
import tools.*;

import java.util.Random;


/**
 * Source of bootstrap data. Synthesizes repeated measurements from each test function, 
 * adds noise and then samples with replacement from the repeats.
 *
 *
 * 
 * @author Philip Cook
 * @version $Id$
 * 
 */
public class BootstrapDataSynthFromInput extends BootstrapDataSynthesizer 
    implements DataSource {


    protected DataSynthFromInput inputSource;

    protected DW_Scheme impars;

    protected double snr;


    /**
     * 
     * @param input
     *            input data source.
     * 
     * @param repeats
     *            The number of repeats to take from each voxel of the input data.
     *
     * @param components
     *            The number of values in each voxel.
     *
     * @param samples
     *            The number of bootstrap samples to generate from each voxel.
     *
     * @param snr
     *            The signal to noise ratio with q=0.
     * 
     *
     * @param ran 
     *            A random number generator.
     */
    public BootstrapDataSynthFromInput(DataSynthFromInput input, int repeats, DW_Scheme scheme, 
				       int samples, double snr, Random ran) {

	super(repeats);

	impars = scheme;

	inputSource = input;

	this.snr = snr;

	init(samples, scheme.numMeasurements(), ran);


    }


    /**
     * 
     * @param input
     *            input data source.
     * 
     * @param repeats
     *            The number of repeats to take from each voxel of the input data.
     *
     * @param components
     *            The number of values in each voxel.
     *
     * @param samples
     *            The number of bootstrap samples to generate from each voxel.
     *
     * @param snr
     *            The signal to noise ratio with q=0.
     *
     * @param seed 
     *            Seed for the random number generator
     */
    public BootstrapDataSynthFromInput(DataSynthFromInput input, int repeats, DW_Scheme scheme, 
				       int samples, double snr, int seed) {
	this(input, repeats, scheme, samples, snr, new numerics.MTRandom(seed));
    }






    /**
     * Gets the next set of source data.
     */
    protected void getNextSourceVoxel() throws DataSourceException {

        if (reachedEndOfFile) {
            noMoreData = true;
            throw new DataSourceException("No more voxels in data source.");
        }

	// inefficiently get data from source, re-normalize it, add noise,
	// and un-normalize it again.
	
	// DataSynthesizer already has code to do this.
	// However, it uses java.util.Random, which is not very good for small samples.

	double[] noiseless = inputSource.nextVoxel();

	if (snr > 0) {
	    double aZero = impars.geoMeanZeroMeas(noiseless);

	    double[] normalized = new double[components];

	    for (int i = 0; i < components; i++) {
		normalized[i] = noiseless[i] / aZero;
	    }
			    
	    for (int r = 0; r < repeats; r++) {

		System.arraycopy(normalized, 0, next[r], 0, components);
		
		next[r] = DataSynthesizer.addNoise(next[r], 1.0 / snr, ran);
	    
		for (int i = 0; i < components; i++) {
		    next[r][i] = aZero * next[r][i];
		}
	    }
	    
	}
	else {
	    for (int r = 0; r < repeats; r++) {
		System.arraycopy(noiseless, 0, next[r], 0, components);
	    }
	}

	if (!inputSource.more()) {
	    reachedEndOfFile = true;
	}

    }

}

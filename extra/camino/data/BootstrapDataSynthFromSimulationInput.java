package data;

import java.io.*;

import misc.*;
import simulation.*;
import tools.*;

import java.util.Random;


/**
 * Source of bootstrap data. Synthesizes repeated measurements from each test function, 
 * adds noise and then samples with replacement from the repeats. Required because
 * DiffusionSimulation always adds noise.
 *
 *
 * 
 * @author Philip Cook
 * @version $Id$
 * 
 */
public class BootstrapDataSynthFromSimulationInput extends BootstrapDataSynthesizer
    implements DataSource {

    private DiffusionSimulation inputSource;



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
     * @param ran 
     *            A random number generator.
     */
    public BootstrapDataSynthFromSimulationInput(DiffusionSimulation input, int repeats, 
						 int components, int samples, Random ran) {

	super(repeats);

	inputSource = input;

	init(samples, components, ran);

    }


    /**
     * 
     * @param input
     *            input data source. Should have <code>repeats</code> voxels of data.
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
     *
     * @param seed 
     *            Seed for the random number generator
     */
    public BootstrapDataSynthFromSimulationInput(DiffusionSimulation input, int repeats, 
						 int components, int samples, int seed) {

	this(input, repeats, components, samples, new numerics.MTRandom(seed));

    }






    /**
     * Gets the source data. This method should only be called once and should
     * provide one voxel of data for each repeat.
     *
     */
    protected void getNextSourceVoxel() throws DataSourceException {

        if (reachedEndOfFile) {
            noMoreData = true;
            throw new DataSourceException("No more voxels in data source.");
        }

	for (int r = 0; r < repeats; r++) {
	    next[r] = inputSource.nextVoxel();
	}

	if (!inputSource.more()) {
	    reachedEndOfFile = true;
	}

    }

}

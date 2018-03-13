package data;

import java.util.Random;
import tools.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Source of synthetic data.
 * 
 * <dt>Description:
 * 
 * <dd>Samples the Fourier transform of a test function.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class DataSynthesizer implements DataSource {

    /**
     * Contains the imaging parameters for the sequence emulated to synthesize
     * data.
     */
    protected DW_Scheme imParams;

    /**
     * The test function.
     */
    protected ModelPDF p;

    /**
     * Signal to noise ratio of q=0 measurements.
     */
    protected double snr;

    /**
     * An array containing noise free samples of the FT of p.
     */
    protected double[] noiseFreeSamples;

    /**
     * The number of voxels the data source provides. The negative value
     * indicates infinity, i.e., the data source continues providing new data
     * forever.
     */
    protected int numVoxels = -1;

    /**
     * The number of the next voxel to give the data from.
     */
    protected int voxelNum;

    /**
     * Random number generator for adding noise.
     */
    protected Random r;

    /**
     * Constructor requires the <code>DW_Scheme</code>, the test function and the
     * signal to noise level.
     * 
     * @param pdf
     *            The test function to sample.
     * 
     * @param ip
     *            The <code>DW_Scheme</code> object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     */
    public DataSynthesizer(ModelPDF pdf, DW_Scheme ip, double s) {
        // Make an infinite data source
        init(pdf, ip, s, -1, 0);
    }

    /**
     * Constructor includes the number of voxels in the synthetic data source.
     * 
     * @param pdf
     *            The test function to sample.
     * 
     * @param ip
     *            The <code>DW_Scheme</code> object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     * 
     * @param numVox
     *            The number of voxels to provide.
     */
    public DataSynthesizer(ModelPDF pdf, DW_Scheme ip, double s, int numVox) {
        init(pdf, ip, s, numVox, 0);
    }

    /**
     * Constructor includes the number of voxels in the synthetic data source.
     * 
     * @param pdf
     *            The test function to sample.
     * 
     * @param ip
     *            The <code>DW_Scheme</code> object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     * 
     * @param numVox
     *            The number of voxels to provide.
     * 
     * @param seed
     *            Random number generator seed.
     */
    public DataSynthesizer(ModelPDF pdf, DW_Scheme ip, double s, int numVox,
            int seed) {
        init(pdf, ip, s, numVox, seed);
    }


    /**
     * Constructor includes the number of voxels in the synthetic data
     * source and the actual <code>Random</code> object to use.
     *
     * @param pdf The test function to sample.
     *
     * @param ip The <code>DW_Scheme</code> object with the details of the
     * scanner sequence to emulate.
     *
     * @param s The signal to noise ratio with q=0.
     *
     * @param numVox The number of voxels to provide.
     *
     * @param ran the Random number generator.
     */
    public DataSynthesizer(ModelPDF pdf, DW_Scheme ip, double s, int numVox, Random ran) {
		init(pdf, ip, s, numVox, 0);
		r = ran;
    }

    /**
     * Does the work of the constructors.
     * 
     * @param pdf
     *            The test function to sample.
     * 
     * @param ip
     *            The <code>DW_Scheme</code> object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     * 
     * @param numVox
     *            The number of voxels to provide.
     * 
     * @param seed
     *            Random number generator seed.
     */
    protected void init(ModelPDF pdf, DW_Scheme ip, double s, int numVox, int seed) {
        imParams = ip;
        p = pdf;
        snr = s;

        noiseFreeSamples = new double[imParams.numMeasurements()];
        for (int i = 0; i < noiseFreeSamples.length; i++) {
            noiseFreeSamples[i] = p.ftAtB_Vec(ip.getG_Dir(i), ip.getB_Value(i));
        }

        r = new Random(seed);

        numVoxels = numVox;
        voxelNum = 1;
    }


    public double[] nextVoxel() throws DataSourceException {
        if (numVoxels >= 0 && voxelNum > numVoxels) {
            throw new DataSourceException("Exceeded the specified number of voxels.");
        }

        voxelNum += 1;

        double[] noisyData;

		if(snr<=0) { 
		    
		    int numMeas = noiseFreeSamples.length;
		    
		    noisyData = new double[numMeas];
		    
		    System.arraycopy(noiseFreeSamples, 0, noisyData, 0, numMeas);
		}
		else {
	
		    if(CL_Initializer.noiseType.equals("gaussian")) {
	                noisyData = addGaussianNoise(noiseFreeSamples, (snr > 0) ? 1.0 / snr : 0.0, r);
		    }
		    else {
	                noisyData = addNoise(noiseFreeSamples, (snr > 0) ? 1.0 / snr : 0.0, r);
		    }
		}

        return noisyData;
    }

    public boolean more() {
        return (numVoxels >= 0 && voxelNum <= numVoxels);
    }


    /**
     * Adds isotropic complex Gaussian noise to each measurement and takes the
     * modulus resulting in Rician noise.
     * 
     * @param samples
     *            The measurements to add noise to.
     * 
     * @param sigma
     *            The standard deviation of the real and imaginary components of
     *            the noise.
     * 
     * @param r
     *            The random number generator to use to synthesize the noise.
     * 
     * @return An array of noisy measurements.
     */
    public static double[] addNoise(double[] samples, double sigma, Random r) {

        double[] noisySamples = new double[samples.length];

        for (int i = 0; i < samples.length; i++) {
            Complex nSamp = new Complex(samples[i] + r.nextGaussian() * sigma, r
                    .nextGaussian()
                    * sigma);
            noisySamples[i] = nSamp.mod();
        }

        return noisySamples;
    }


    /**
     * Adds Gaussian noise to each measurement.
     * 
     * @param samples
     *            The measurements to add noise to.
     * 
     * @param sigma
     *            The standard deviation of the noise.
     * 
     * @param r
     *            The random number generator to use to synthesize the noise.
     * 
     * @return An array of noisy measurements.
     */
    public static double[] addGaussianNoise(double[] samples, double sigma, Random r) {

        double[] noisySamples = new double[samples.length];

        for (int i = 0; i < samples.length; i++) {
            noisySamples[i] = samples[i]  + r.nextGaussian() * sigma;
        }

        return noisySamples;
    }

}

package data;

import java.util.Random;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Source of synthetic data bootstrapped from a small number of samples.
 * 
 * <dt>Description:
 * 
 * <dd>Uses a DataSynthesizer to generate a fixed number of synthetic
 * measurement sets, then draws combinations from those samples at random to
 * generate new measurement sets.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: SyntheticBootstrapper.java,v 1.4 2005/08/18 10:59:36 ucacmgh
 *          Exp $
 *  
 */
public class SyntheticBootstrapper implements DataSource {

    /**
     * The data synthesizer.
     */
    protected DataSynthesizer dataSynth;

    /**
     * An array containing sets of repeated measurements.
     */
    protected double[][] repeatedMeasurements;

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
     * Constructor includes the number of voxels in the synthetic data source.
     * 
     * @param pdf
     *            The test function to sample.
     * 
     * @param ip
     *            The DW_Scheme object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     * 
     * @param samples
     *            The number of samples of each measurement.
     * 
     * @param numVox
     *            The number of voxels to provide.
     */
    public SyntheticBootstrapper(ModelPDF pdf, DW_Scheme ip, double s,
            int samples, int numVox) {
        init(pdf, ip, s, samples, numVox, new MTRandom(12345));
    }

    /**
     * Constructor includes the number of voxels in the synthetic data source.
     * 
     * @param pdf
     *            The test function to sample.
     * 
     * @param ip
     *            The DW_Scheme object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     * 
     * @param samples
     *            The number of samples of each measurement.
     * 
     * @param numVox
     *            The number of voxels to provide.
     * 
     * @param seed
     *            Random number generator seed.
     */
    public SyntheticBootstrapper(ModelPDF pdf, DW_Scheme ip, double s,
            int samples, int numVox, int seed) {
        init(pdf, ip, s, samples, numVox, new MTRandom(seed));
    }

 /**
     * Constructor includes the number of voxels in the synthetic data source.
     * 
     * @param pdf
     *            The test function to sample.
     * 
     * @param ip
     *            The DW_Scheme object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     * 
     * @param samples
     *            The number of samples of each measurement.
     * 
     * @param numVox
     *            The number of voxels to provide.
     * 
     * @param ran
     *            Random number generator.
     */
    public SyntheticBootstrapper(ModelPDF pdf, DW_Scheme ip, double s,
            int samples, int numVox, Random ran) {
        init(pdf, ip, s, samples, numVox, ran);
    }


    /**
     * Does the work of the constructors.
     * 
     * @param pdf
     *            The test function to sample.
     * 
     * @param ip
     *            The DW_Scheme object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     * 
     * @param samples
     *            The number of samples of each measurement.
     * 
     * @param numVox
     *            The number of voxels to provide.
     * 
     * @param seed
     *            Random number generator seed.
     */
    protected void init(ModelPDF pdf, DW_Scheme ip, double s, int samples,
            int numVox, Random ran) {

        r = ran;

        dataSynth = new DataSynthesizer(pdf, ip, s, samples, r);
        int numMeas = ip.numMeasurements();
        repeatedMeasurements = new double[numMeas][samples];
        for (int i = 0; i < samples; i++)
            try {
                double[] v = dataSynth.nextVoxel();
                for (int j = 0; j < numMeas; j++) {
                    repeatedMeasurements[j][i] = v[j];
                }
            }
            catch (Exception e) {
                System.err.println("Error in initialization of SyntheticBootstrapper.");
                throw new RuntimeException(e);
            }



        numVoxels = numVox;
        voxelNum = 1;
    }

    public double[] nextVoxel() throws DataSourceException {
        if (numVoxels >= 0 && voxelNum > numVoxels) {
            throw new DataSourceException("Exceeded the specified number of voxels.");
        }

        voxelNum += 1;

        // Choose each measurement at random from the set of
        // possibilities.
        int numMeas = repeatedMeasurements.length;
        int samples = repeatedMeasurements[0].length;
        double[] voxelData = new double[numMeas];
        for (int i = 0; i < numMeas; i++) {
            int index = r.nextInt(samples);
            voxelData[i] = repeatedMeasurements[i][index];
        }

        return voxelData;
    }

    public boolean more() {
        return (numVoxels >= 0 && voxelNum <= numVoxels);
    }

}

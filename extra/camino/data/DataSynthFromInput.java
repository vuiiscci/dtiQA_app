package data;

import imaging.*;
import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General class for synthesizing diffusion weighted data from model
 * parameters read from an input stream.
 * 
 * <dt>Description:
 * 
 * <dd>This class is abstract and extended for each particular kind of model.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public abstract class DataSynthFromInput implements DataSource {

    /**
     * The model parameters are read from this data source.
     */
    protected DataSource dataSource;

    /**
     * Contains the imaging parameters for the sequence emulated to synthesize
     * data.
     */
    protected DW_Scheme imParams;

    /**
     * Signal to noise ratio of q=0 measurements.
     */
    protected double snr;

    /**
     * Random number generator seed.
     */
    protected int seed;

    /**
     * Specifies the number of data items in each voxel for a particular model.
     * Overridden in subclasses.
     */
    protected int DATAITEMSPERVOXEL;

    /**
     * Default constructor does nothing.
     */
    public DataSynthFromInput() {
    }

    /**
     * The constructor requires the name of the data file and the details of the
     * acquisition sequence to use for data synthesis. If the filename is null,
     * data comes from the standard input.
     * 
     * @param filename
     *            Name of the data file.
     * 
     * @param inputDataType
     *            The type of the input data
     * 
     * @param ip
     *            The DW_Scheme object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     */
    public DataSynthFromInput(String filename, String inputDataType,
            DW_Scheme ip, double s) {
        init(filename, inputDataType, ip, s, 0);
    }

    /**
     * This constructor allows the random number generator seed to be specified.
     * 
     * @param filename
     *            Name of the data file.
     * 
     * @param inputDataType
     *            The type of the input data
     * 
     * @param ip
     *            The DW_Scheme object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     * 
     * @param seed
     *            Random number generator seed.
     */
    public DataSynthFromInput(String filename, String inputDataType,
            DW_Scheme ip, double s, int seed) {
        init(filename, inputDataType, ip, s, seed);
    }

    /**
     * Called by the constructors to do the initialization.
     * 
     * @param filename
     *            Name of the data file.
     * 
     * @param inputDataType
     *            The type of the input data
     * 
     * @param ip
     *            The DW_Scheme object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     */
    protected void init(String filename, String inputDataType, DW_Scheme ip,
            double s, int initSeed) {
        imParams = ip;
        snr = s;
        dataSource = ExternalDataSource.getDataSource(filename, DATAITEMSPERVOXEL, inputDataType);
        seed = initSeed;
    }

    public double[] nextVoxel() throws DataSourceException {

        // Read in the next set of model parameters.
        double[] modelData = dataSource.nextVoxel();
        double exitCode = modelData[0];
        double aZero = Math.exp(modelData[1]);

        // Return all zeros if the voxel is in the background.
        if (exitCode < 0) {
            return new double[imParams.numMeasurements()];
        }

        // Construct the model PDF from the parameters on the input.
        ModelPDF p = getNextModel(modelData);

        // Create a data synthesizer using the model.
        DataSynthesizer ds = new DataSynthesizer(p, imParams, snr, 1, new MTRandom(seed));

        // Need to increment the seed or we get the same noise on every
        // voxel.
        seed += 1;

        // Get the synthetic data from the data synthesizer.
        double[] normData = ds.nextVoxel();

        // Unnormalize and return.
        double[] rawData = new double[normData.length];
        for (int i = 0; i < rawData.length; i++) {
            rawData[i] = aZero * normData[i];
        }

        return rawData;
    }

    /**
     * Constructs the model PDF from the next voxel's worth of input data.
     * 
     * @param modelData
     *            The parameters of the model.
     * 
     * @return The model to use for data synthesis.
     */
    protected abstract ModelPDF getNextModel(double[] modelData);

    public boolean more() {
        return dataSource.more();
    }

}

package data;

import imaging.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Class for synthesizing diffusion weighted data from diffusion tensor
 * parameters read from an input stream.
 * 
 * <dt>Description:
 * 
 * <dd>This data source provides synthetic data from each consecutive diffusion
 * tensor read from an input stream. The data is synthesized by emulating a
 * specified imaging sequence. The input should have the same format as the
 * output of <code>DT_Inversion</code>.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @see inverters.DT_Inversion
 * @version $Id: DataSynthFromDT_Input.java,v 1.4 2005/08/18 10:59:37 ucacmgh
 *          Exp $
 *  
 */
public class DataSynthFromDT_Input extends DataSynthFromInput {

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
     *            The <code>DW_Scheme</code> object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.
     */
    public DataSynthFromDT_Input(String filename, String inputDataType,
            DW_Scheme ip, double s) {
        DATAITEMSPERVOXEL = 8;
        init(filename, inputDataType, ip, s, 0);
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
     *            The <code>DW_Scheme</code> object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            The signal to noise ratio with q=0.

     * @param seed the seed for the random number generator. 
     *
     */
    public DataSynthFromDT_Input(String filename, String inputDataType,
            DW_Scheme ip, double s, int seed) {
        DATAITEMSPERVOXEL = 8;
        init(filename, inputDataType, ip, s, seed);
    }


    /**
     * Constructs the model PDF from the next voxel's worth of input data.
     * 
     * @param modelData
     *            The parameters of the model.
     * 
     * @return The model to use for data synthesis.
     */
    protected ModelPDF getNextModel(double[] modelData) {

        // Set up the required parameters to the GaussianMixture
        // constructor.
        DT[] dts = new DT[1];
        dts[0] = new DT(modelData[2], modelData[3], modelData[4], modelData[5],
                modelData[6], modelData[7]);

        double[] mixPars = { 1.0 };

        return new GaussianMixture(dts, mixPars);

    }

}

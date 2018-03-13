package data;

import imaging.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Class for synthesizing diffusion weighted data from two-tensor
 * parameters read from an input stream.
 * 
 * <dt>Description:
 * 
 * <dd>This data source provides synthetic data from each consecutive pair of
 * diffusion tensors read from an input stream. The data is synthesized by
 * emulating a specified imaging sequence. The input should have the same format
 * as the output of <code>TwoTensorInversion</code>.
 * 
 * </dl>
 *
 * @see inverters.TwoTensorInversion
 * 
 * @author Danny Alexander
 * @version $Id: DataSynthFromTwoTensorInput.java,v 1.4 2005/08/18 10:59:37
 *          ucacmgh Exp $
 *  
 */
public class DataSynthFromTwoTensorInput extends DataSynthFromInput {

    /**
     * The constructor requires the name of the data file and the details of the
     * acquisition sequence to use for data synthesis. If the filename is null,
     * data comes from the standard input.
     * 
     * @param filename
     *            name of the data file.
     * 
     * @param inputDataType
     *            the type of the input data
     * 
     * @param ip
     *            the <code>DW_Scheme</code> object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            the signal to noise ratio with q=0.
     */
    public DataSynthFromTwoTensorInput(String filename, String inputDataType,
            DW_Scheme ip, double s) {
        DATAITEMSPERVOXEL = 17;
        init(filename, inputDataType, ip, s, 0);
    }

    /**
     * The constructor requires the name of the data file and the details of the
     * acquisition sequence to use for data synthesis. If the filename is null,
     * data comes from the standard input.
     * 
     * @param filename
     *            name of the data file.
     * 
     * @param inputDataType
     *            the type of the input data
     * 
     * @param ip
     *            the <code>DW_Scheme</code> object with the details of the scanner
     *            sequence to emulate.
     * 
     * @param s
     *            the signal to noise ratio with q=0.

     * @param seed 
     *            the seed for the random number generator. 
     *
     */
    public DataSynthFromTwoTensorInput(String filename, String inputDataType,
            DW_Scheme ip, double s, int seed) {
        DATAITEMSPERVOXEL = 17;
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
        DT[] dts = new DT[2];
        dts[0] = new DT(modelData[4], modelData[5], modelData[6], modelData[7],
                modelData[8], modelData[9]);
        dts[1] = new DT(modelData[11], modelData[12], modelData[13], modelData[14],
                modelData[15], modelData[16]);

        double[] mixPars = new double[2];
        mixPars[0] = modelData[3];
        mixPars[1] = modelData[10]; 

        return new GaussianMixture(dts, mixPars);

    }

}

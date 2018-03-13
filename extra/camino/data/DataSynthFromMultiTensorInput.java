package data;

import imaging.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Class for synthesizing diffusion weighted data from multi-tensor
 * parameters read from an input stream.
 * 
 * <dt>Description:
 * 
 * <dd>This data source provides synthetic data from the diffusion tensors in each 
 * voxel. The data is synthesized by
 * emulating a specified imaging sequence. The input should have the same format
 * as the output of <code>MultiTensorInversion</code>.
 * 
 * </dl>
 * 
 * @see inverters.MultiTensorInversion
 *
 * @author Danny Alexander
 * @version $Id: DataSynthFromTwoTensorInput.java,v 1.4 2005/08/18 10:59:37
 *          ucacmgh Exp $
 *  
 */
public class DataSynthFromMultiTensorInput extends DataSynthFromInput {

    private final int maxComponents;


    /**
     * @param filename the name of the input file.
     * @param inputDataType input primitive data type.
     * @param ip the imaging parameters
     * @param maxComps the maximum number of tensors in the voxel.
     * @param s the signal to noise ratio.
     */
    public DataSynthFromMultiTensorInput(String filename, String inputDataType,
            DW_Scheme ip, int maxComps, double s) {
        DATAITEMSPERVOXEL = 3 + 7 * maxComps;
	maxComponents = maxComps;
        init(filename, inputDataType, ip, s, 0);
    }


    /**
     * @param filename the name of the input file.
     * @param inputDataType input primitive data type.
     * @param ip the imaging parameters
     * @param maxComps the maximum number of tensors in the voxel.
     * @param s the signal to noise ratio.
     * @param seed the random seed.
     */
    public DataSynthFromMultiTensorInput(String filename, String inputDataType,
            DW_Scheme ip, int maxComps, double s, int seed) {
        DATAITEMSPERVOXEL = 3 + 7 * maxComps;
        maxComponents = maxComps;
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

	int components = (int)modelData[2];

        // Set up the required parameters to the GaussianMixture
        // constructor.
        DT[] dts = new DT[components];
	double[] mixPars = new double[components];
 
	for (int i = 0; i < components; i++) {
	    int start = 4 + 7 * i;
	    
	    dts[i] = new DT(modelData[start], modelData[start+1], modelData[start+2], 
			    modelData[start+3], modelData[start+4], modelData[start+5]);
	    
	    mixPars[i] = modelData[3+7*i];
	}

   
        return new GaussianMixture(dts, mixPars);

    }

}

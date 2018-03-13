package inverters;

import imaging.*;
import misc.*;
import data.*;
import apps.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fits models selected according to an input voxel
 * classification.  Typically, one tensor in single fibre regions, two
 * or more tensors in multiple fibre regions.
 * 
 * <dt>Description:
 * 
 * <dd> Reads the classification for each voxel from an input file and
 * chooses which inversion to run on the data from that voxel.  The
 * output format is the standard multiple tensor format: {exitcode,
 * ln(A(0)), number of components, DT_1, ..., DT_N}, where each DT has
 * format: {mix, Dxx, Dxy, Dxz, Dyy, Dyz, Dzz}.
 * 
 * </dl>
 *
 * @see inverters.ModelIndex
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class MultiTensorInversion extends DiffusionInversion {


    /**
     * The inversions for each class of input voxel
     */
    protected DiffusionInversion[] inverters;


    /**
     * The number of tensors to output in each voxel.
     */
    protected int comps;


    /**
     * The data source for the voxel classification.
     */
    protected DataSource voxClassStream;


    /**
     * This constructor requires the details of the sequence, the name
     * of the file containing the voxel classification and the list of
     * indices of the inverters.  Voxel's labelled i in the
     * classification are fitted with the model with index at position
     * i in the modelIndices array.  If i is greater than the length
     * of the array, the voxel is fitted with the model with index in
     * the last position in the array.
     *
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param voxClassFile
     *            The name of the file containing the voxel classification.
     * 
     * @param modelIndices
     *            Indices of the models for each level of classification.
     * 
     * @param maxComponents
     *            The maximum number of output components in each voxel.
     */
    public MultiTensorInversion(DW_Scheme imParams, String voxClassFile, ModelIndex[][] modelIndices, 
				int maxComponents) {
        init(imParams, voxClassFile, modelIndices, maxComponents);
    }


    /**
     * Called by the constructors to create a MultiTensorInverter
     * with specified parameters.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param voxClassFile
     *            The name of the file containing the voxel classification.
     * 
     * @param modelIndices
     *            Indices of the models for each level of classification.
     * 
     * @param maxComponents
     *            The maximum number of output components in each voxel.
     */
    protected void init(DW_Scheme imParams, String voxClassFile, ModelIndex[][] modelIndices, 
			int maxComponents) {

        ip = imParams;

	comps = maxComponents;
	if(comps < 1) {
	    throw new 
		LoggedException("Must have at least one output component in MultiTensorInversion.");
	}

	// Allocate the array of inverters and create each object.
	// the fallback inverter for a two- or three-tensor inversion is the same as 
	// the first inverter in the list
	inverters = new DiffusionInversion[modelIndices.length];
	for(int i=0; i<modelIndices.length; i++) {
	    inverters[i] = 
		ModelFit.getIndexedInversion(modelIndices[i], imParams);
	}

	// Set up the voxel classification input stream.
	voxClassStream = ExternalDataSource.getDataSource(voxClassFile, 1, "int");

    }


    /**
     * Fits the model.
     * 
     * @param data
     *            The MRI data.
     * 
     * @return {exitcode, ln A^\star(0), num components, mix1, DT1,
     * mix2, DT2, ...}
     */
    public double[] invert(double[] data) {

        boolean background = false;

	// Get the voxel classfication.
	if(!voxClassStream.more()) {
	    throw new LoggedException("Voxel classification file too small.");
	}
	double[] nextVC = voxClassStream.nextVoxel();
	int voxClass = (int)nextVC[0];
	if(voxClass >= inverters.length) {
	    voxClass = inverters.length - 1;
	}
        else if (voxClass == -1) {
            // background voxels treated as order 0
            background = true;
            voxClass = 0;
        }
        // Do the inversion.
        double[] model = inverters[voxClass].invert(data);

	// Construct the output array.
	double[] params = new double[itemsPerVoxel()];
	if(inverters[voxClass] instanceof DT_Inversion) {
 
	    // Need to add in the number of components item for
	    // compatability with the multitensor format.

	    // Copy exitcode, ln(A(0)) and set components to 1 and mixing
            // parameter to 1.
	    params[0] = model[0];
	    params[1] = model[1];
	    params[2] = 1.0;
	    params[3] = 1.0;

            // if background, set exit code to -1 and number of components to zero
            if (background) {
                params[0] = -1.0;
                params[2] = 0.0;
            }

	    // Copy the single component.
	    for(int i=0; i<6; i++) {
		params[i+4] = model[i+2];
	    }

	}
	else {
 
	    // The data is already in the right format, but may have
	    // too many components.  Copy what we can.
	    for(int i=0; i<params.length && i<model.length; i++) {
		params[i] = model[i];
	    }
	}

        return params;
    }


    /**
     * Need to read off one item from the voxel classification when a 
     * background voxel is found.
     */
    public void background() {
	if(!voxClassStream.more()) {
	    throw new LoggedException("Voxel classification file too small.");
	}
	double[] nextVC = voxClassStream.nextVoxel();
    }


    public int itemsPerVoxel() {
        return 3 + 7*comps;
    }

}

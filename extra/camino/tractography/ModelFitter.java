package tractography;

import apps.ModelFit;
import data.*;
import imaging.*;
import inverters.*;
import misc.*;
import numerics.*;


/**
 * Superclass for classes that fit some model to DWI data. These are wrappers for the
 * actual model fitters in the inverters and fitters package. Extending this lets you
 * fit your choice of model to some data and use it for bootstrap tractography
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public abstract class ModelFitter {


    protected final DW_Scheme scheme;

    protected final int numMeas;

    private double[] data;
   
    protected Vector3D[] vectors;


    public ModelFitter(DW_Scheme imPars) {
        scheme = imPars;

        data = null;

        vectors = null;

        numMeas = scheme.numMeasurements();
    }



    /**
     * Fit the model to some data. A user-specified number of PDs will be computed. If the number of PDs 
     * found by the model is greater than the number specified, subclasses will use some procedure to choose
     * <code>numPDs</code> of them. If there are not enough independent PDs, then the last PD will be copied 
     * as necessary. This is done because the number of PDs per voxel is currently hard coded into the tractography
     * images.
     *
     * @param the DWI data.
     * @param numPDs the exact number of principal directions required.
     */
    public final void fitModel(double[] dat, int numPDs) {

        if (dat.length != numMeas) {
            throw new LoggedException("data does not match scheme for fitter");
        }

        data = dat;

        setPDs(data, numPDs);
    }


    protected abstract void setPDs(double[] data, int numPDs);


    /**
     * @return the PDs in this voxel. 

     *
     */
    public final Vector3D[] getPDs() {
        Vector3D[] defCopy = new Vector3D[vectors.length];

        System.arraycopy(vectors, 0, defCopy, 0, vectors.length);

        return defCopy;
    }


    public final DW_Scheme getScheme() {
        return scheme;
    }

    public final int numMeasurements() {
	return numMeas;
    }



    public static final ModelFitter getModelFitter(DW_Scheme scheme, ModelIndex[] indices) {
        if (indices[0].numDTs > 0) { 
            return new TensorModelFitter(scheme, indices);
        }
        else {
            throw new UnsupportedOperationException("No tractography interface implemented for model " + indices[0]);
        }
    }

}

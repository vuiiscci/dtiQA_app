package tractography;

import data.*;
import imaging.*;
import inverters.*;
import misc.DT;
import numerics.*;


/**
 * Fits a tensor or two-tensor model to DWI data hence the maximum number of 
 * PDs in a voxel is 2.
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class TensorModelFitter extends ModelFitter {
  
    private double[] mix;

    private DT[] tensors;

    private final DiffusionInversion linearDTInversion;

    private final DiffusionInversion oneDTInversion;
    private final DiffusionInversion twoDTInversion;


    /**
     * Construct a fitter with a specified DT or multi-tensor fitting routine.
     *
     * @param imPars the imaging scheme corresponding to the data.
     * @
     */
    public TensorModelFitter(DW_Scheme imPars, ModelIndex[] indices) {
        super(imPars);

	if (indices.length == 1) {
	    indices = new ModelIndex[] {ModelIndex.POSPOS, indices[0]};
	}

        linearDTInversion = DT_Inversion.getIndexedDT_Inversion(ModelIndex.LDT, imPars);

        oneDTInversion = DT_Inversion.getIndexedDT_Inversion(indices[1], imPars);
        twoDTInversion = new TwoTensorInversion(imPars, indices[0], indices[1]);


    }


   
    /**
     * Computes the DTs.
     *
     * @param numDTs 0, 1, or 2.
     */
    public void setPDs(double[] data, int numDTs) {

        if (numDTs == 0) {
            tensors = new DT[0];
            mix = new double[0];
            vectors = new Vector3D[0];
            return;
        }
        
	if (numDTs == 1) {
	    
	    double[] params = oneDTInversion.invert(data);
	    
            if (params[0] == 2.0) {
                params = linearDTInversion.invert(data);
            }

	    tensors = new DT[1];

            tensors[0] = new DT(params[2], params[3], params[4], params[5], params[6], params[7]);

            mix = new double[] {1.0};

            double[][] seig = tensors[0].sortedEigenSystem();
	    
	}
	else if (numDTs == 2) {
	    double[] params = twoDTInversion.invert(data);
	    
	    tensors = new DT[2];

            mix = new double[] {params[3], params[10]};

            tensors[0] = new DT(params[4], params[5], params[6], params[7], params[8], params[9]);

            tensors[1] = new DT(params[11], params[12], params[13], 
                                params[14], params[15], params[16]);
            
	}
	else {
	    throw new UnsupportedOperationException("Can't invert data higher than order 4");
	}

        vectors = new Vector3D[numDTs];
        
	for (int t = 0; t < numDTs; t++) {
            
            double[][] eig = tensors[t].sortedEigenSystem();
            
            vectors[t] = new Vector3D(eig[1][0], eig[2][0], eig[3][0]);
        }

    }


    /**
     * Get the diffusion tensors.
     *
     */
    public DT[] getDTs() {
        DT[] defCopy = new DT[tensors.length];

        System.arraycopy(tensors, 0, defCopy, 0, tensors.length);

        return defCopy;

    }


    /**
     * Get the mixing parameters, in the order corresponding to the tensors from #getDTs().
     *
     */
    public double[] getMix() {

        double[] defCopy = new double[mix.length];

        System.arraycopy(mix, 0, defCopy, 0, mix.length);

        return defCopy;

        
    }



}
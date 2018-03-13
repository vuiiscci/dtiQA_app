package fitters;

import java.util.Random;
import models.*;
import optimizers.*;
import inverters.*;
import numerics.*;
import misc.*;
import data.*;
import models.compartments.CompartmentModel;
import tools.*;
import imaging.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fits the ball and stick model using multiple runs of a Levenburg
 * Marquardt and assuming a Gaussian noise model.
 * 
 * <dt>Description:
 * 
 * <dd> Fitting is as in BallAndStickLM_GaussianFitter.  The
 * perturbations are 10% of the initial starting value.
 * 
 * </dl>
 *
 * @author  laura(panagio@cs.ucl.ac.uk)
 * @version $Id$
 *
 */
public class BallStickMultiRunLM_GaussianFitter extends BallStickLM_GaussianFitter {


    /**
     * Constructor implements the mapping between model and optimized
     * parameters in the Codec object.
     *
     * @param scheme The imaging protocol.
     */
    public BallStickMultiRunLM_GaussianFitter(DW_Scheme scheme, int repeats, int seed) {
	
	this.scheme = scheme;
	makeCodec();
	String[] compNames = new String[2];
	compNames[0] = new String("stick");
	compNames[1] = new String("ball");
	double[] initialParams = new double[7];
	initialParams[0]= 1.0; //the S0
	initialParams[1]= 0.3; //volume fraction of the first compartment
	initialParams[2]= 0.7; //volume fraction of the second compartment
	initialParams[3]= 1.7E-9;//stick diffusivity
	initialParams[4]= 1.570796326794897 ;//theta
	initialParams[5]= 0.0 ;//phi
	initialParams[6]= 1.7E-9;//ball diffusivity
	
	
	

	 cm = new CompartmentModel(compNames, initialParams);

        try {
            minimizer = new MultiRunLM_GaussianMinimizer(scheme, cm, cod, new FixedSTD_GaussianPerturbation(), repeats, seed);
        } catch(Exception e) {
            throw new LoggedException(e);
        }

	dtfitter = new LinearDT_Inversion(scheme);

    }

    public static void main(String[] args) {
        
        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
        CL_Initializer.initImagingScheme();
        CL_Initializer.initDataSynthesizer();

        OutputManager om = new OutputManager();

        BallStickMultiRunLM_GaussianFitter inv = new BallStickMultiRunLM_GaussianFitter(CL_Initializer.imPars, 3, 0);

        // Loop over the voxels.
        while (CL_Initializer.data.more()) {
            try {

                double[] nextVoxel = CL_Initializer.data.nextVoxel();
                double[][] fit = inv.fit(nextVoxel);

                for(int i=0; i<fit.length; i++)
                    om.output(fit[i]);
            } catch(Exception e) {
                System.err.println(e);
            }
        }

        om.close();
     }
 
}


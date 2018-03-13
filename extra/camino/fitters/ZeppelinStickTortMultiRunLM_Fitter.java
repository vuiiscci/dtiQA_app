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
 * <dd>Fits the zeppelin and stick model with tortuosity approximation
 * using multiple runs of a Levenburg Marquardt and a specified noise model.
 * 
 * <dt>Description:
 * 
 * <dd> Fitting is as in ZeppelinStickTortLM_Fitter.  The
 * perturbations are 10% of the initial starting value.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class ZeppelinStickTortMultiRunLM_Fitter extends ZeppelinStickTortLM_Fitter {


    /**
     * Constructor implements the mapping between model and optimized
     * parameters in the Codec object.
     *
     * @param scheme The imaging protocol.
     */
    public ZeppelinStickTortMultiRunLM_Fitter(DW_Scheme scheme, int repeats, int seed) {
	
	this.scheme = scheme;
	makeCodec();
	String[] compNames = new String[2];
	compNames[0] = new String("stick");
	compNames[1] = new String("zeppelin");
	double[] initialParams = new double[7];

	cm = new CompartmentModel(compNames, initialParams);

        try {
            initMultiRunLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel), new FixedSTD_GaussianPerturbation(), repeats, seed);
	    ((MultiRunLM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
	    ((MultiRunLM_Minimizer) minimizer).setMAXITER(5000);
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

        ZeppelinStickTortMultiRunLM_Fitter inv = new ZeppelinStickTortMultiRunLM_Fitter(CL_Initializer.imPars, 3, 0);

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


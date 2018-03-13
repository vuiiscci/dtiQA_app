package fitters;

import models.*;
import models.compartments.CompartmentModel;
import optimizers.*;
import inverters.*;
import numerics.*;
import misc.*;
import data.*;
import tools.*;
import imaging.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fits the tensor and stick model using multiple runs of a Levenburg
 * Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd> Fitting is as in TensorStickLM_Fitter.  The
 * perturbations are 10% of the initial starting value.
 * 
 * </dl>
 *
 * @author  Laura
 * @version $Id$
 *
 */
public class TensorStickMultiRunLM_Fitter extends
		TensorStickLM_Fitter {

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public TensorStickMultiRunLM_Fitter(DW_Scheme scheme, int repeats,
			int seed) {

		this.scheme = scheme;
		makeCodec();
		String[] compNames = new String[2];
		compNames[0] = new String("stick");
		compNames[1] = new String("tensor");
		double[] initialParams = new double[12];
		initialParams[0] = 1.0; //the S0
		initialParams[1] = 0.7; //volume fraction of the first compartment
		initialParams[2] = 0.3; //volume fraction of the second compartment
		initialParams[3] = 1.7E-9;//stick diffusivity
		initialParams[4] = 1.570796326794897;//theta
		initialParams[5] = 1.570796326794897;//phi
		initialParams[6] = 1.7e-9;//diff
		initialParams[7] = 1.570796326794897;//theta
		initialParams[8] = 1.570796326794897;//phi
		initialParams[9] = 1.7e-10;//diffperp1
		initialParams[10] = 1.7e-11;//diffperp2
		initialParams[11] = 0.0;//alpha

		 cm = new CompartmentModel(compNames, initialParams);
		 try {
	            initMultiRunLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel), new FixedSTD_GaussianPerturbation(), repeats, seed);
		    ((MultiRunLM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
		    ((MultiRunLM_Minimizer) minimizer).setMAXITER(5000);
	        } catch(Exception e) {
	            throw new LoggedException(e);
	        }

		dtfitter = new LinearDT_Inversion(scheme);
		zsfitter = new ZeppelinStickMultiRunLM_Fitter(scheme, 3, 0);

	}


	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		TensorStickMultiRunLM_Fitter inv = new TensorStickMultiRunLM_Fitter(
				CL_Initializer.imPars, 3, 0);

		// Loop over the voxels.
		while (CL_Initializer.data.more()) {
			try {

				double[] nextVoxel = CL_Initializer.data.nextVoxel();
				double[][] fit = inv.fit(nextVoxel);

				for (int i = 0; i < fit.length; i++)
					om.output(fit[i]);
			} catch (Exception e) {
				System.err.println(e);
			}
		}

		om.close();
	}

}

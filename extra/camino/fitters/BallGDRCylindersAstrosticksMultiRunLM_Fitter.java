package fitters;

import models.*;
import models.compartments.CompartmentModel;
import models.compartments.CompartmentType;
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
 * <dd>Fits the ballcylinderastrosticks model using multiple runs of a Levenburg
 * Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd> Fitting is as in BallCylinderAstrosticksLM_GaussianFitter.  The
 * perturbations are 10% of the initial starting value.
 * 
 * </dl>
 *
 * @author  Laura
 * @version $Id$
 *
 */
public class BallGDRCylindersAstrosticksMultiRunLM_Fitter extends
		BallGDRCylindersAstrosticksLM_Fitter {

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public BallGDRCylindersAstrosticksMultiRunLM_Fitter(DW_Scheme scheme,
			int repeats, int seed) {

		this.scheme = scheme;
		makeCodec();
		
		String[] compNames = new String[3];
		compNames[0] = new String("gammadistribradiicylinders");
		compNames[1] = new String("ball");
		compNames[2] = new String("astrosticks");
		double[] initialParams = new double[11];
		initialParams[0] = 1.0; // the S0
		initialParams[1] = 0.6; // volume fraction of the first compartment
		initialParams[2] = 0.3; // volume fraction of the second compartment
		initialParams[3] = 0.1;// volume fraction of the third compartment
		initialParams[4] = 1 ;//k
		initialParams[5] = 1 ;//beta
		initialParams[6] = 1.7e-9;//diffusivity
		initialParams[7] = 1.5;//theta
		initialParams[8] = 1.5; //phi
		initialParams[9] = 1.7e-9;//diffusivity
		initialParams[10] = 1.7e-9;//diffusivity

		 cm = new CompartmentModel(compNames, initialParams);

		 try {
	            initMultiRunLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel), new FixedSTD_GaussianPerturbation(), repeats, seed);
		    ((MultiRunLM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
		    ((MultiRunLM_Minimizer) minimizer).setMAXITER(5000);
	        } catch(Exception e) {
	            throw new LoggedException(e);
	        }

		bcfitter = new BallCylinderMultiRunLM_Fitter(scheme, 1, 0);

	}



	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		BallGDRCylindersAstrosticksMultiRunLM_Fitter inv = new BallGDRCylindersAstrosticksMultiRunLM_Fitter(
				CL_Initializer.imPars, 100, 0);

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

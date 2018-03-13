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
 * <dd>Fits the biexponential model (two cylindrically symmetric tensors) using one run of a Levenburg
 * Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each compartment.
 * Parallel diffusivity, unweighted signal and radius are constrained positive
 * by optimizing their square root. Volume fraction is constrained to [0 1] by
 * optimizing cos^{-1} sqrt(f). Perpendicular diffusivity is constrained to be
 * lower than the parallel diffusivity by optimising sin^(-1)sqrt(dperp/dparal).
 * </dl>
 * 
 * @author Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 * 
 */
public class BiZeppelinLM_Fitter extends CompartmentFitter {

	
	protected ZeppelinStickLM_Fitter zsfitter;

	private final int numModelParams = 11;
	private final int numOptParams = 7;

	/**
	 * Default constructor required for derived classes.
	 */
	public BiZeppelinLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized parameters
	 * in the Codec object.
	 * 
	 * @param scheme
	 *            The imaging protocol.
	 */
	public BiZeppelinLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = FitModel
				.getCompartmentList(CL_Initializer.fitModel);
		double[] initialParams = FitModel
				.getCompartmentModelParams(CL_Initializer.fitModel);

		 cm = new CompartmentModel(compNames, initialParams);

		 try {
             initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
             ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
             ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
		throw new LoggedException(e);
	}

		
		zsfitter = new ZeppelinStickMultiRunLM_Fitter(scheme, 1, 0);
	}

	/**
	 * Creates the Codec object.
	 */
	protected void makeCodec() {
		cod = new Codec() {

			public double[] modelToOpt(double[] modelParams) {
				double[] optParams = new double[numOptParams];
				optParams[0] = Math.sqrt(modelParams[0]); // s0
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); // f
																		// 
				optParams[2] = Math.sqrt(modelParams[3]); // diff 
				optParams[3] = modelParams[4]; // theta
				optParams[4] = modelParams[5]; // phi
				optParams[5] = (Math.asin(Math.sqrt(modelParams[6]/ modelParams[3])));// diffprp constrained
				optParams[6] = (Math.asin(Math.sqrt(modelParams[10]/ modelParams[3])));// diffprp constrained

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double zepf = Math.cos(optParams[1]);
				zepf = zepf * zepf;
				modelParams[1] = zepf; // volume fraction of zep1
											// compartment
				modelParams[2] = 1 - zepf; // volume fraction of zeppelin 2
												// compartment
				modelParams[3] = optParams[2] * optParams[2]; // diff intra
				modelParams[4] = optParams[3]; // theta
				modelParams[5] = optParams[4]; // phi
				modelParams[6] = Math.sin(optParams[5])* Math.sin(optParams[5]);
				modelParams[6] = modelParams[6]* modelParams[3];// diff perpconstrained
				modelParams[7] = modelParams[3]; // diff extra
				modelParams[8] = modelParams[4];// theta zep
				modelParams[9] = modelParams[5];// phi zep
				modelParams[10] = Math.sin(optParams[6])* Math.sin(optParams[6]);
				modelParams[10] = modelParams[10] * modelParams[3];// diff perpconstrained

				return modelParams;
			}

			public int getNumOptParams() {
				return numOptParams;
			}

			public int getNumModelParams() {
				return numModelParams; 
			}
			public int getDirectionIndex() {
				return 4;
			}
		};

	}

	/**
	 * Estimates a starting set of parameters from the 
	 * zeppelinstick fitter.The volume fraction is zeppelinstick's truncated to
	 * the range [0.1 0.9]. The parallel diffusivity and perpendicular
	 * diffusivity is the zeppelinstick diffusivity.
	 * 
	 * 
	 * @param data
	 *            The set of measurements to fit to.
	 * 
	 * @return The starting model parameter values.
	 */
	protected double[] getStartPoint(double[] data) {

		// Set starting point from command line
		if (this.getFixedStartPoint(data)!=null)
		{
			return this.getFixedStartPoint(data);
		}

		
		// Do zeppelinstick fit
		double[] zsparams;

		try {
			zsparams = MultiRunMinimizer.getBestSolution(zsfitter.fit(data));

		} catch (MinimizerException e) {
			zsparams = zsfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = zsparams[1];
		double theta = zsparams[5];
		double phi = zsparams[6];

		// set mixing parameter in range 0.1 - 0.9
		double f = zsparams[2];

		// set diffusivity
		double diff = zsparams[4];



		double eigenvalue2a = zsparams[10];
	

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;
		modelParams[1] = f; // f
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		modelParams[2] = 1 - modelParams[1]; // fzep2
		modelParams[3] = diff; // diff
		modelParams[4] = theta;
		modelParams[5] = phi;
		modelParams[6] = eigenvalue2a;// diffperp
		modelParams[6] = modelParams[6] < 1e-11 ? 1e-11 : modelParams[6];
		modelParams[7] = modelParams[3];// diffzep
		modelParams[8] = modelParams[4];// thetazep
		modelParams[9] = modelParams[5];// phizep
		modelParams[10] =modelParams[6];// diffperp

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		BiZeppelinLM_Fitter inv = new BiZeppelinLM_Fitter(
				CL_Initializer.imPars);

		// Loop over the voxels.
		while (CL_Initializer.data.more()) {
			try {

				double[] nextVoxel = CL_Initializer.data.nextVoxel();
				double[][] fit = inv.fit(nextVoxel);


				om.output(fit[0]);
			} catch (Exception e) {
				System.err.println(e);
			}
		}

		om.close();
	}

}

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
 * <dd>Fits the tensor and stick model using one run of a Levenburg
 * Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each
 * compartment. Parallel diffusivity and unweighted signal are constrained
 * positive by optimizing their square root.  Volume fraction is
 * constrained to [0 1] by optimizing cos^{-1} sqrt(f). Perpendicular 
 * diffusivity1 is constrained to be lower than the parallel diffusivity 
 * by optimising sin^(-1)sqrt(dperp/dparal) and perpendicular diffusivity2 
 * is constrained to be lower than the perpendicular diffusivity1
 * by optimising sin^(-1)sqrt(dperp1/dperp2) 
 * 
 * </dl>
 *
 * @author  Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 *
 */
public class TensorStickLM_Fitter extends CompartmentFitter {

	protected LinearDT_Inversion dtfitter;
	protected ZeppelinStickLM_Fitter zsfitter;

	private final int numModelParams = 12;
	private final int numOptParams = 8;

	/**
	 * Default constructor required for derived classes.
	 */
	public TensorStickLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public TensorStickLM_Fitter(DW_Scheme scheme) {

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
		initialParams[11] = 1.5;//alpha

		 cm = new CompartmentModel(compNames, initialParams);
		 try {
             initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
             ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
             ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
		throw new LoggedException(e);
	}

		dtfitter = new LinearDT_Inversion(scheme);
		zsfitter = new ZeppelinStickMultiRunLM_Fitter(scheme, 1, 0);
	}

	/**
	 * Creates the Codec object.
	 */
	protected void makeCodec() {
		cod = new Codec() {

			public double[] modelToOpt(double[] modelParams) {
				double[] optParams = new double[numOptParams];
				optParams[0] = Math.sqrt(modelParams[0]); //s0
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); //f stick
				optParams[2] = Math.sqrt(modelParams[3]); //diff stick
				optParams[3] = modelParams[4]; //theta
				optParams[4] = modelParams[5]; //phi
				optParams[5] = (Math.asin(Math.sqrt(modelParams[9]
						/ modelParams[3])));//diffprp constrained
				optParams[6] = (Math.asin(Math.sqrt(modelParams[10]
						/ modelParams[9])));
				optParams[7] = modelParams[11];//alpha

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double stickf = Math.cos(optParams[1]);
				stickf = stickf * stickf;

				modelParams[1] = stickf; // volume fraction of stick compartment
				modelParams[2] = 1 - stickf; //volume fraction of the extracellular compartment
				modelParams[3] = optParams[2] * optParams[2]; // diff stick
				modelParams[4] = optParams[3]; //theta
				modelParams[5] = optParams[4]; //phi
				modelParams[6] = modelParams[3]; // diff 
				modelParams[7] = modelParams[4];//theta 
				modelParams[8] = modelParams[5];//phi 
				modelParams[9] = Math.sin(optParams[5])
						* Math.sin(optParams[5]);//diff perp constrained
				modelParams[9] = modelParams[9] * modelParams[3];//diff perp constrained
				modelParams[10] = Math.sin(optParams[6])
						* Math.sin(optParams[6]);//diff perp constrained
				modelParams[10] = modelParams[10] * modelParams[9];//diff perp2
				modelParams[11] = optParams[7];//alpha 

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
	 * Estimates a starting set of parameters from the linear estimate
	 * of the diffusion tensor from the log measurements and the zeppelinstick fitter.  The stick
	 * orientation is the zeppelinstick stick's direction.  The volume fraction
	 * is the zeppelinstick's truncated to the range [0.1 0.9]. The parallel and first perpendicular diffusivity is
	 * zeppelinstick's diffusivities. The second perpendicular diffusivity and the alpha angle are from the DT fitter.
	 *
	 * @param data The set of measurements to fit to.
	 *
	 * @return The starting model parameter values.
	 */
	protected double[] getStartPoint(double[] data) {

		//Set starting point from command line
		if (this.getFixedStartPoint(data)!=null)
		{
			return this.getFixedStartPoint(data);
		}

		// Do zeppelinstick fit.
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

		// set diffusivity as mean DT diffusivity	
		double diff = zsparams[4];

		// Do a linear DT fit
		double[] dtparams = dtfitter.invert(data);

		DT dt = new DT(dtparams[2], dtparams[3], dtparams[4], dtparams[5],
				dtparams[6], dtparams[7]);

		double[][] seig = dt.sortedEigenSystem();
		
		//perpendicular diffusivities

		double eigenvalue2 = zsparams[10];
		double eigenvalue3 = seig[0][2];
		//double eigenvalue3 = zsparams[9];
		

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;
		modelParams[1] = f; //f stick
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		modelParams[2] = 1 - modelParams[1]; //ftensor
		modelParams[3] = diff; //diff stick
		modelParams[4] = theta;
		modelParams[5] = phi;
		modelParams[6] = modelParams[3];//diff
		modelParams[7] = modelParams[4];//theta
		modelParams[8] = modelParams[5];//phi
		modelParams[9] = eigenvalue2;//diffperp1 
		modelParams[9] = modelParams[9] < 1e-11 ? 1e-11 : modelParams[9];
		modelParams[10] = eigenvalue3;//diffperp2 
		modelParams[10] = modelParams[10] < 1e-11 ? 1e-11: modelParams[10];
		
		
		modelParams[11] = 0;//alpha 

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		TensorStickLM_Fitter inv = new TensorStickLM_Fitter(
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

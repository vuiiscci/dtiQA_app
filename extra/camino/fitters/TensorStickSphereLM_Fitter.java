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
 * <dd>Fits the tensorstick and sphere model using one run of a Levenburg
 * Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each
 * compartment. Parallel diffusivity and unweighted signal are constrained
 * positive by optimizing their square root.  Volume fraction is
 * constrained to [0 1] by optimizing cos^{-1} sqrt(f). Perpendicular 
 * diffusivity1 is constrained to be lower than the parallel diffusivity 
 * by optimising sin^(-1)sqrt(dperp/dparal) and perpendicular diffusivity 
 * is constrained to be lower than the perpendicular diffusivity1
 * by optimising sin^(-1)sqrt(dperp1/dperp2).The radius is constrained to 
 * 0.-20 microns by optimizing cos^{-1} sqrt(R/20e-6). 
 * 
 * </dl>
 *
 * @author  Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 *
 */
public class TensorStickSphereLM_Fitter extends CompartmentFitter {

	protected TensorStickLM_Fitter tsfitter;
	protected TensorStickAstrosticksLM_Fitter tsafitter;

	private final int numModelParams = 15;
	private final int numOptParams = 10;

	/**
	 * Default constructor required for derived classes.
	 */
	public TensorStickSphereLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public TensorStickSphereLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = new String[3];
		compNames[0] = new String("stick");
		compNames[1] = new String("tensor");
		compNames[2] = new String("spheregpd");
		double[] initialParams = new double[15];
		initialParams[0] = 1.0; //the S0
		initialParams[1] = 0.7; //volume fraction of the first compartment
		initialParams[2] = 0.3; //volume fraction of the second compartment
		initialParams[3] = 0.1; //volume fraction of the third compartment
		initialParams[4] = 1.7E-9;//stick diffusivity
		initialParams[5] = 1.570796326794897;//theta
		initialParams[6] = 1.570796326794897;//phi
		initialParams[7] = 1.7e-9;//diff
		initialParams[8] = 1.570796326794897;//theta
		initialParams[9] = 1.570796326794897;//phi
		initialParams[10] = 1.7e-10;//diffperp1
		initialParams[11] = 1.7e-11;//diffperp2
		initialParams[12] = 1.5;//alpha
		initialParams[13] = 1.7E-9;//diff astro
		initialParams[14] = 2E-6;//R astro

		 cm = new CompartmentModel(compNames, initialParams);
		 try {
             initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
             ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
             ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
		throw new LoggedException(e);
	}

		tsfitter = new TensorStickMultiRunLM_Fitter(scheme, 3, 0);
		tsafitter = new TensorStickAstrosticksMultiRunLM_Fitter(scheme, 3, 0);

	}

	/**
	 * Creates the Codec object.
	 */
	protected void makeCodec() {
		cod = new Codec() {

			public double[] modelToOpt(double[] modelParams) {
				double[] optParams = new double[numOptParams];
				optParams[0] = Math.sqrt(modelParams[0]); //s0
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); //f 
				optParams[2] = Math.acos(Math.sqrt(modelParams[2]/(1-modelParams[1]))); // f tens
				optParams[3] = Math.sqrt(modelParams[4]); //diff stick
				optParams[4] = modelParams[5]; //theta
				optParams[5] = modelParams[6]; //phi
				optParams[6] = (Math.asin(Math.sqrt(modelParams[10]
						/ modelParams[4])));//diffprp constrained
				optParams[7] = (Math.asin(Math.sqrt(modelParams[11]
						/ modelParams[10])));
				optParams[8] = modelParams[12];//alpha
				
				optParams[9] =(Math.acos(Math
						.sqrt((modelParams[14] - 1e-7) / 40e-6)));// R//R constraining the radius from 0-20 microns
				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double f = Math.cos(optParams[1]); // volume fraction ;
				f = f * f;
				modelParams[1] = f; // volume fraction
				double tensf = Math.cos(optParams[2]) * Math.cos(optParams[2])* (1 - f);
				modelParams[2] = tensf; // volume fraction of tensor
				modelParams[3] = 1 - f - tensf;//volume fraction of the third compartment
				modelParams[4] = optParams[3] * optParams[3]; // diff stick
				modelParams[5] = optParams[4]; //theta
				modelParams[6] = optParams[5]; //phi
				modelParams[7] = modelParams[4]; // diff 
				modelParams[8] = modelParams[5];//theta 
				modelParams[9] = modelParams[6];//phi 
				modelParams[10] = Math.sin(optParams[6])
						* Math.sin(optParams[6]);//diff perp constrained
				modelParams[10] = modelParams[10] * modelParams[4];//diff perp constrained
				modelParams[11] = Math.sin(optParams[7])
						* Math.sin(optParams[7]);//diff perp constrained
				modelParams[11] = modelParams[11] * modelParams[10];//diff perp2
				modelParams[12] = optParams[8];//alpha 
				modelParams[13] = modelParams[4]; // diff 
				modelParams[14] = 1e-7 + (Math.cos(optParams[9])
						* Math.cos(optParams[9]) * 40e-6);// R constrained

				return modelParams;
			}

			public int getNumOptParams() {
				return numOptParams;
			}

			public int getNumModelParams() {
				return numModelParams; 
			}
			public int getDirectionIndex() {
				return 5;
			}
		};

	}

	/**
	 * Estimates a starting set of parameters from the the tensorstick fitter.
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

		double[] tsparams;

		try {
			tsparams = MultiRunMinimizer.getBestSolution(tsfitter.fit(data));

		} catch (MinimizerException e) {

			tsparams = tsfitter.getStartPoint(data);
			e.printStackTrace();
		}
		double[] tsaparams;

		try {
			tsaparams = MultiRunMinimizer.getBestSolution(tsafitter.fit(data));

		} catch (MinimizerException e) {

			tsaparams = tsafitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = tsparams[1];
		// set mixing parameter
		double f1 = tsparams[2]-0.00005;
		double f2 = tsparams[3]-0.00005;
		double f3 = 0.0001;// f3 is not 0 because single run LM breaks and gives singular matrix -2.

		// set diffusivity 	
		double diff = tsparams[4];
		double theta = tsparams[5];
		double phi = tsparams[6];
		//perpendicular diffusivity 1
		double eigenvalue2 = tsparams[10];
		//perpendicular diffusivity 2
		double eigenvalue3 = tsparams[11];
		double alpha = tsparams[12];
		double R = 2e-6;

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;
		modelParams[1] = f1;
		modelParams[2] = f2;
		modelParams[3] = f3;//fzep
		modelParams[4] = diff; //diff stick
		modelParams[5] = theta;
		modelParams[6] = phi;
		modelParams[7] = modelParams[4];//diff
		modelParams[8] = modelParams[5];//theta
		modelParams[9] = modelParams[6];//phi
		modelParams[10] = eigenvalue2;//diffperp1 
		modelParams[11] = eigenvalue3;//diffperp2 
		modelParams[12] = alpha;//alpha 
		modelParams[13] = modelParams[4];//diff
		modelParams[14] = R;

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		TensorStickSphereLM_Fitter inv = new TensorStickSphereLM_Fitter(
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

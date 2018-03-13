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
 * <dd>Fits the ballstick and sphere model using one run of a Levenburg
 * Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each
 * compartment.  Diffusivity and unweighted signal are constrained
 * positive by optimizing their square root. Volume fractions of the intracellular
 * and extracellular compartments are constrained to [0 1] by optimizing cos^{-1} sqrt(f).
 The radius is constrained to 20 microns by optimizing cos^{-1} sqrt(R/20e-6).
 * 
 * </dl>
 *
 * @author  laura(panagio@cs.ucl.ac.uk)
 * @version $Id$
 *
 */
public class BallStickSphereLM_Fitter extends CompartmentFitter {

	protected BallStickLM_GaussianFitter bsfitter;
	protected BallStickAstrosticksLM_Fitter bsafitter;

	private final int numModelParams = 10;
	private final int numOptParams = 7;

	/**
	 * Default constructor required for derived classes.
	 */
	public BallStickSphereLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public BallStickSphereLM_Fitter(DW_Scheme scheme) {

		this(scheme, FitModel.getCompartmentList(CL_Initializer.fitModel),
				FitModel.getCompartmentModelParams(CL_Initializer.fitModel));

	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public BallStickSphereLM_Fitter(DW_Scheme scheme,
			String[] compNames, double[] initialParams) {

		this.scheme = scheme;
		makeCodec();
		 cm = new CompartmentModel(compNames, initialParams);

		 try {
             initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
             ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
             ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
		throw new LoggedException(e);
	}

		bsfitter = new BallStickMultiRunLM_GaussianFitter(scheme, 3, 0);
		bsafitter = new BallStickAstrosticksMultiRunLM_Fitter(scheme, 3, 0);

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
				optParams[2] = Math.acos(Math.sqrt(modelParams[2]/(1-modelParams[1]))); // f ball
				optParams[3] = Math.sqrt(modelParams[4]); //diff stick
				optParams[4] = modelParams[5]; //theta
				optParams[5] = modelParams[6]; //phi
				optParams[6] = (Math.acos(Math
						.sqrt((modelParams[9] - 1e-7) / 40e-6)));// R


				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double f = Math.cos(optParams[1]); // volume fraction ;
				f = f * f;
				modelParams[1] = f; // volume fraction
				double ballf = Math.cos(optParams[2]) * Math.cos(optParams[2])* (1 - f);
				modelParams[2] = ballf; // volume fraction of ball
				modelParams[3] = 1 - f - ballf; // volume fraction third compartment
				modelParams[4] = optParams[3] * optParams[3]; // diff 
				modelParams[5] = optParams[4]; //theta
				modelParams[6] = optParams[5]; //phi
				modelParams[7] = modelParams[4]; // diff ball
				modelParams[8] = modelParams[7]; //diff sphere
				modelParams[9] = 1e-7 + (Math.cos(optParams[6])
						* Math.cos(optParams[6]) * 40e-6);// R constrained


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
	 * Estimates a starting set of parameters from the ballstick model. 
	 *
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
		
		double bsparams[];
		
		try {
			bsparams = MultiRunMinimizer.getBestSolution(bsfitter.fit(data));

		} catch (MinimizerException e) {
			
			bsparams = bsfitter.getStartPoint(data);
			e.printStackTrace();
		}
double bsaparams[];
		
		try {
			bsaparams = MultiRunMinimizer.getBestSolution(bsafitter.fit(data));

		} catch (MinimizerException e) {
			
			bsaparams = bsafitter.getStartPoint(data);
			e.printStackTrace();
		}
		
		double S0 = bsparams[1];
		double f1 = bsparams[2]-0.00005;
		double f2 = bsparams[3]-0.00005;
		double f3 = 0.0001;// f3 is not 0 because single run LM breaks and gives singular matrix -2.

		
		
	    
		double diff = bsparams[4];
		double theta = bsparams[5];
		double phi = bsparams[6];
		double R =2e-6;

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;
		modelParams[1] = f1;
		modelParams[2] = f2;
		modelParams[3] = f3;
		modelParams[4] = diff;
		modelParams[5] = theta;
		modelParams[6] = phi;
		modelParams[7] = modelParams[4];
		modelParams[8] = modelParams[7];
		modelParams[9] = R;

		return modelParams;

	
	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		BallStickSphereLM_Fitter inv = new BallStickSphereLM_Fitter(
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

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
 * <dd>Fits the ball and gamma distributed radii cylinders  and sphere model using one run of a Levenburg Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each compartment.
 * Diffusivity, the unweighted signal and the radius are constrained positive by
 * optimizing their square root. Volume fraction is constrained to [0 1] by
 * optimizing cos^{-1} sqrt(f).
 * 
 * </dl>
 * 
 * @author Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 * 
 */
public class BallGDRCylindersSphereLM_Fitter extends CompartmentFitter {

	protected BallCylinderLM_Fitter bcfitter;

	private final int numModelParams = 12;
	private final int numOptParams = 9;

	/**
	 * Default constructor required for derived classes.
	 */
	public BallGDRCylindersSphereLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized parameters
	 * in the Codec object.
	 * 
	 * @param scheme
	 *            The imaging protocol.
	 */
	public BallGDRCylindersSphereLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = new String[3];
		compNames[0] = new String("gammadistribradiicylinders");
		compNames[1] = new String("ball");
		compNames[2] = new String("spheregpd");
		double[] initialParams = new double[12];
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
		initialParams[11] = 1.7e-6;//R
		

		 cm = new CompartmentModel(compNames, initialParams);

		 try {
             initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
             ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
             ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
		throw new LoggedException(e);
	}
		bcfitter = new BallCylinderMultiRunLM_Fitter(scheme, 1, 0);

	}

	/**
	 * Creates the Codec object.
	 */
	protected void makeCodec() {
		cod = new Codec() {

			public double[] modelToOpt(double[] modelParams) {
				double[] optParams = new double[numOptParams];
				optParams[0] = Math.sqrt(modelParams[0]); // s0
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); //f 
				optParams[2] = Math.acos(Math.sqrt(modelParams[2]/(1-modelParams[1]))); // f2
				
				optParams[3] = Math.acos(Math.sqrt(modelParams[4]/20)); // k
				optParams[4] = (Math.acos(Math.sqrt((modelParams[5] - 1e-7) / 20e-6)));// beta

				optParams[5] = Math.sqrt(modelParams[6]); // diff cylinder
				optParams[6] = modelParams[7]; // theta
				optParams[7] = modelParams[8]; //phi
				if (optParams[7]<1e-5)
				{
					optParams[7]=optParams[7]-Math.PI;
				}
				
				optParams[8] = (Math.acos(Math
						.sqrt((modelParams[11] - 1e-7) / 40e-6)));// Rglial constrained
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
				modelParams[4] = Math.cos(optParams[3]);//k
				modelParams[4] = modelParams[4]*modelParams[4]*20;
				modelParams[5] =  1e-7 + (Math.cos(optParams[4])
										* Math.cos(optParams[4]) * 20e-6); // beta
				modelParams[6] = optParams[5]*optParams[5];//diff
				modelParams[7] = optParams[6]; // theta
				modelParams[8] = optParams[7];//phi
				modelParams[9] = modelParams[6]; // diff ball
				modelParams[10] = modelParams[6]; // diff 
				modelParams[11] = 1e-7 + (Math.cos(optParams[8])
						* Math.cos(optParams[8]) * 40e-6);// Rglial constrained


				
				return modelParams;
			}

			public int getNumOptParams() {
				return numOptParams;
			}

			public int getNumModelParams() {
				return numModelParams; 
			}
			public int getDirectionIndex() {
				return 7;
			}
		};

	}

	/**
	 * Estimates a starting set of parameters from the ballstick fitter. The
	 * cylinder orientation is the stick's direction. The volume fraction is the
	 * ballstick's truncated to the range [0.1 0.9]. The diffusivity is the
	 * ballstick's diffusivity. The starting point for the radius is set to
	 * 2e-6.
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

		// Do ballstick fit, if fails do dtfit.
		double[] bcparams;

		try {
			bcparams = MultiRunMinimizer.getBestSolution(bcfitter.fit(data));

		} catch (MinimizerException e) {
			// dtfitter.invert(data);
			bcparams = bcfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = bcparams[1];
		double theta = bcparams[5];
		double phi = bcparams[6];

		// set mixing parameter in range 0.1 - 0.9
		double f = bcparams[2]-0.00005;
		double f2 = bcparams[3]-0.00005;
		double f3 = 0.0001;// f3 is not 0 because single run LM breaks and gives singular matrix -2.
		// set diffusivity as mean DT diffusivity
		double diff = bcparams[4];

		// set cylinder radius
		double k= 1.8;
		double beta =  bcparams[7];
		double R = 2e-6;

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;// s0
		modelParams[1] = f;// intracellular volume fraction
		modelParams[2] = f2;
		modelParams[3] = f3;
		modelParams[4] = k;// 
		modelParams[5] = beta;
		modelParams[6] = diff;
		modelParams[7] = theta;
		modelParams[8] = phi;// extracellular diffusivity
		modelParams[9] = diff;
		modelParams[10] = diff;
		modelParams[11] = R;
		
		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		BallGDRCylindersSphereLM_Fitter inv = new BallGDRCylindersSphereLM_Fitter(
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

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
 * <dd>Fits the ball and gamma distributed radii cylinders model using one run
 * of a Levenburg Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd>The diffusivity parameter is assumed the same in each compartment.
 * Diffusivity, the unweighted signal are constrained positive by optimizing
 * their square root. Volume fraction is constrained to [0 1] by optimizing
 * cos^{-1} sqrt(f). The shape parameter k is constrained to 0-20 microns by
 * optimizing cos^{-1}sqrt(k/20e-6), and the scale parameter beta, is
 * constrained to 0.1-20 microns by optimizing cos^{-1}sqrt(beta-1e-7/20e-6).
 * 
 * </dl>
 * 
 * @author Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 * 
 */
public class BallGDRCylindersLM_Fitter extends CompartmentFitter {

	protected BallCylinderLM_Fitter bcfitter;

	private final int numModelParams = 9;
	private final int numOptParams = 7;

	/**
	 * Default constructor required for derived classes.
	 */
	public BallGDRCylindersLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized parameters
	 * in the Codec object.
	 * 
	 * @param scheme
	 *            The imaging protocol.
	 */
	public BallGDRCylindersLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = new String[2];
		compNames[0] = new String("gammadistribradiicylinders");
		compNames[1] = new String("ball");
		double[] initialParams = new double[9];
		initialParams[0] = 1.0; // the S0
		initialParams[1] = 0.7; // volume fraction of the first compartment
		initialParams[2] = 0.3; // volume fraction of the second compartment
		initialParams[3] = 1;// k
		initialParams[4] = 1;// beta
		initialParams[5] = 1.7e-9;// diffusivity
		initialParams[6] = 1.5;// theta
		initialParams[7] = 1.5; // phi
		initialParams[8] = 1.7e-9;// diffusivity

		cm = new CompartmentModel(compNames, initialParams);

		try {
            initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
            ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
            ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
	throw new LoggedException(e);
}

		bcfitter = new BallCylinderMultiRunLM_Fitter(scheme, 3, 0);

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
				optParams[2] = Math.acos(Math.sqrt(modelParams[3] / 20)); // k
				optParams[3] = (Math.acos(Math
						.sqrt((modelParams[4] - 1e-7) / 20e-6)));// beta
				optParams[4] = Math.sqrt(modelParams[5]); // diff cylinder
				optParams[5] = modelParams[6]; // theta
				optParams[6] = modelParams[7]; // phi

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double cylindf = Math.cos(optParams[1]);
				cylindf = cylindf * cylindf;
				modelParams[1] = cylindf; // volume fraction of cylinder
				modelParams[2] = 1 - cylindf; // volume fraction of ball
				modelParams[3] = Math.cos(optParams[2]);// k
				modelParams[3] = modelParams[3] * modelParams[3] * 20;
				modelParams[4] = 1e-7 + (Math.cos(optParams[3])
						* Math.cos(optParams[3]) * 20e-6); // beta
				modelParams[5] = optParams[4] * optParams[4];// diff
				modelParams[6] = optParams[5]; // theta
				modelParams[7] = optParams[6];// phi
				modelParams[8] = modelParams[5]; // diff ball

				return modelParams;
			}

			public int getNumOptParams() {
				return numOptParams;
			}

			public int getNumModelParams() {
				return numModelParams; 
			}
			public int getDirectionIndex() {
				return 6;
			}
		};

	}

	/**
	 * Estimates a starting set of parameters from the ballcylinder fitter. The
	 * GDRCylinders orientation is the cylinders's direction. The volume fraction is the
	 * ballcylinders's truncated to the range [0.1 0.9]. The diffusivity is the
	 * ballcylinder's diffusivity. The starting point for k= 1.8 and beta is the radius of the cylinder.
	 * @param data
	 *            The set of measurements to fit to.
	 * 
	 * @return The starting model parameter values.
	 */
	protected double[] getStartPoint(double[] data) {

		// Set starting point from command line
		if (this.getFixedStartPoint(data) != null) {
			return this.getFixedStartPoint(data);
		}

		// Do ballcylinder fit, if fails do dtfit.
		double[] bcparams;

		try {
			bcparams = MultiRunMinimizer.getBestSolution(bcfitter.fit(data));

		} catch (MinimizerException e) {
			bcparams = bcfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = bcparams[1];
		double theta = bcparams[5];
		double phi = bcparams[6];

		// set mixing parameter in range 0.1 - 0.9
		double f = bcparams[2];

		// set diffusivity as ballcylinder diffusivity
		double diff = bcparams[4];

		// set k,beta from cylinders radius
		double k = 1.8;
		double beta = bcparams[7];

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;// s0
		modelParams[1] = f;// intracellular volume fraction
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		modelParams[2] = 1 - modelParams[1]; // extracellular volume fraction
		modelParams[3] = k;//
		modelParams[4] = beta;
		modelParams[5] = diff;
		modelParams[6] = theta;
		modelParams[7] = phi;
		modelParams[8] = diff;
		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		BallGDRCylindersLM_Fitter inv = new BallGDRCylindersLM_Fitter(
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

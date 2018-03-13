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
 * <dd>Fits the ball and cylinder model using one run of a Levenburg Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd>The diffusivity parameter is assumed the same in each compartment.
 * Diffusivity, the unweighted signal are constrained positive by optimizing
 * their square root. Volume fraction is constrained to [0 1] by optimizing
 * cos^{-1} sqrt(f). The radius is constrained to 0.1-20 microns by optimizing
 * cos^{-1} sqrt(R-1e-7/20e-6).
 * 
 * </dl>
 * 
 * @author Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 * 
 */
public class BallCylinderLM_Fitter extends CompartmentFitter {

	protected BallStickLM_GaussianFitter bsfitter;

	private final int numModelParams = 8;
	private final int numOptParams = 6;

	/**
	 * Default constructor required for derived classes.
	 */
	public BallCylinderLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized parameters
	 * in the Codec object.
	 * 
	 * @param scheme
	 *            The imaging protocol.
	 */
	public BallCylinderLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = new String[2];
		compNames[0] = new String("cylindergpd");
		compNames[1] = new String("ball");
		double[] initialParams = new double[8];
		initialParams[0] = 1.0; // the S0
		initialParams[1] = 0.7; // volume fraction of the first compartment
		initialParams[2] = 0.3; // volume fraction of the second compartment
		initialParams[3] = 1.7E-9;// diffusivity
		initialParams[4] = 1.570796326794897;// theta
		initialParams[5] = 0.0;// phi
		initialParams[6] = 2.0e-6;// R
		initialParams[7] = 1.7E-9;// ball diffusivity

		cm = new CompartmentModel(compNames, initialParams);

		try {
            initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
            ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
            ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
	throw new LoggedException(e);
}

		bsfitter = new BallStickMultiRunLM_GaussianFitter(scheme, 3, 0);

	}

	/**
	 * Creates the Codec object.
	 */
	protected void makeCodec() {
		cod = new Codec() {

			public double[] modelToOpt(double[] modelParams) {
				double[] optParams = new double[numOptParams];
				optParams[0] = Math.sqrt(modelParams[0]); // s0
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); // f cylinder
				optParams[2] = Math.sqrt(modelParams[3]); // diff cylinder
				optParams[3] = modelParams[4]; // theta
				optParams[4] = modelParams[5]; // phi
				optParams[5] = (Math.acos(Math
						.sqrt((modelParams[6] - 1e-7) / 20e-6)));// R

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double cylindf = Math.cos(optParams[1]);
				cylindf = cylindf * cylindf;
				modelParams[1] = cylindf; // volume fraction of cylinder
				modelParams[2] = 1 - cylindf; // volume fraction of ball
				modelParams[3] = optParams[2] * optParams[2]; // diff cylinder
				modelParams[4] = optParams[3]; // theta
				modelParams[5] = optParams[4]; // phi
				modelParams[6] = 1e-7 + (Math.cos(optParams[5])
						* Math.cos(optParams[5]) * 20e-6);// R 
				modelParams[7] = modelParams[3]; // diff ball

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
		if (this.getFixedStartPoint(data) != null) {
			return this.getFixedStartPoint(data);
		}

		// Do ballstick fit, if fails do dtfit.
		double[] bsparams;

		try {
			bsparams = MultiRunMinimizer.getBestSolution(bsfitter.fit(data));

		}catch (MinimizerException e) {

			bsparams = bsfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = bsparams[1];
		double theta = bsparams[5];
		double phi = bsparams[6];

		// set mixing parameter in range 0.1 - 0.9
		double f = bsparams[2];

		// set diffusivity as mean DT diffusivity
		double diff = bsparams[4];

		// set cylinder radius
		double R = 2.0e-6;

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;// s0
		modelParams[1] = f;// intracellular volume fraction
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		modelParams[2] = 1 - modelParams[1]; // extracellular volume fraction
		modelParams[3] = diff;// intracellular diffusivity
		modelParams[4] = theta;
		modelParams[5] = phi;
		modelParams[6] = R;
		modelParams[7] = modelParams[3];// extracellular diffusivity
		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		BallCylinderLM_Fitter inv = new BallCylinderLM_Fitter(
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

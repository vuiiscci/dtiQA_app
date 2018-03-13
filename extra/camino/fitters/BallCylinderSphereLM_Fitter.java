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
 * <dd>Fits the ballcylinderspherer model using one run of a Levenburg
 * Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each compartment.
 * Diffusivity and the unweighted signal are constrained positive by optimizing
 * their square root. Volume fractions of the intracellular and extracellular
 * compartments are constrained to [0 1] by optimizing cos^{-1} sqrt(f). The
 * radius for the axons is constrained to 0.1-20 microns by optimizing cos^{-1}
 * sqrt(R-1e-7/20e-6). The radius for the glial is constrained to 0.1-20 microns by optimizing cos^{-1}
 * sqrt(R-1e-7/40e-6).
 * </dl>
 * 
 * @author Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 * 
 */
public class BallCylinderSphereLM_Fitter extends CompartmentFitter {

	protected BallCylinderLM_Fitter bcfitter;

	private final int numModelParams = 11;
	private final int numOptParams = 8;

	/**
	 * Default constructor required for derived classes.
	 */
	public BallCylinderSphereLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized parameters
	 * in the Codec object.
	 * 
	 * @param scheme
	 *            The imaging protocol.
	 */
	public BallCylinderSphereLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = new String[3];
		compNames[0] = new String("cylindergpd");
		compNames[1] = new String("ball");
		compNames[2] = new String("spheregpd");
		double[] initialParams = new double[11];
		initialParams[0] = 1.0; // the S0
		initialParams[1] = 0.6; // volume fraction of the first compartment
		initialParams[2] = 0.3; // volume fraction of the second compartment
		initialParams[3] = 0.1; // volume fraction of the third compartment
		initialParams[4] = 1.7E-9;// diffusivity
		initialParams[5] = 1.570796326794897;// theta
		initialParams[6] = 0.0;// phi
		initialParams[7] = 2.0e-6;// R
		initialParams[8] = 1.7E-9;// ball diffusivity
		initialParams[9] = 1.7E-9;// Astro diffusivity
		initialParams[10] = 4.0e-6;// Rglial
	

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
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); //f 
				optParams[2] = Math.acos(Math.sqrt(modelParams[2]/(1-modelParams[1]))); // f
                optParams[3] = Math.sqrt(modelParams[4]); // diff cylinder
				optParams[4] = modelParams[5]; // theta
				optParams[5] = modelParams[6]; // phi
				optParams[6] = (Math.acos(Math
						.sqrt((modelParams[7] - 1e-7) / 20e-6)));// R constrained
				optParams[7] = (Math.acos(Math
						.sqrt((modelParams[10] - 1e-7) / 20e-6)));// Rglial constrained


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
				modelParams[5] = optParams[4]; // theta
				modelParams[6] = optParams[5]; // phi
				modelParams[7] = 1e-7 + (Math.cos(optParams[6])
						* Math.cos(optParams[6]) * 80e-6);// R constraining
				// the radius from
				// 0-20 microns
				modelParams[8] = modelParams[4]; // diff ball
				modelParams[9] = modelParams[4]; // diff astro
				modelParams[10]  = 1e-7 + (Math.cos(optParams[7])
						* Math.cos(optParams[7]) * 20e-6);// Rglial constraining

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
	 * Estimates a starting set of parameters from the ballcylinder fitter.
	 * 
	 * @param data
	 *            The set of measurements to fit to.
	 * 
	 * @return The starting model parameter values.
	 */
	protected double[] getStartPoint(double[] data) {

		if (this.getFixedStartPoint(data)!=null)
		{
			return this.getFixedStartPoint(data);
		}

		// Do ballcylinder fit
		double[] bcparams;

		try {
			bcparams = MultiRunMinimizer.getBestSolution(bcfitter.fit(data));

		} catch (MinimizerException e) {

			bcparams = bcfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = bcparams[1];

		double f1 = bcparams[2]-0.00005;
		double f2 = bcparams[3]-0.00005;
		double f3 = 0.0001;// f3 is not 0 because single run LM breaks and gives singular matrix -2.
		double diff = bcparams[4];
		double theta = bcparams[5];
		double phi = bcparams[6];
		double R = bcparams[7];
		double Rg= 2e-6;

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];

		modelParams[0] = S0;
		modelParams[1] = f1;
		modelParams[2] = f2;
		modelParams[3] = f3;
		modelParams[4] = diff; // diff
		modelParams[5] = theta;
		modelParams[6] = phi;
		modelParams[7] = R;
		modelParams[8] = modelParams[4];// extracellular diffusivity
		modelParams[9] = modelParams[4];
		modelParams[10] = Rg;
		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		BallCylinderSphereLM_Fitter inv = new BallCylinderSphereLM_Fitter(
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

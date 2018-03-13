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
 * <dd>Fits the zeppelincylindersphere model using one run of a Levenburg
 * Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd>The diffusivity parameter is assumed the same in each compartment.
 * Parallel diffusivity and unweighted signal are constrained positive by
 * optimizing their square root. Volume fraction is constrained to [0 1] by
 * optimizing cos^{-1} sqrt(f). Perpendicular diffusivity is constrained to be
 * lower than the parallel diffusivity by optimising sin^(-1)sqrt(dperp/dparal).
 * </dl>
 * 
 * @author Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 * 
 */
public class ZeppelinCylinderSphereLM_Fitter extends CompartmentFitter {

	protected ZeppelinCylinderLM_Fitter zcfitter;

	private final int numModelParams = 14;
	private final int numOptParams = 9;

	/**
	 * Default constructor required for derived classes.
	 */
	public ZeppelinCylinderSphereLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized parameters
	 * in the Codec object.
	 * 
	 * @param scheme
	 *            The imaging protocol.
	 */
	public ZeppelinCylinderSphereLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = new String[3];
		compNames[0] = new String("cylindergpd");
		compNames[1] = new String("zeppelin");
		compNames[2] = new String("spheregpd");
		double[] initialParams = new double[14];
		initialParams[0] = 1.0; // the S0
		initialParams[1] = 0.7; // volume fraction of the first compartment
		initialParams[2] = 0.2; // volume fraction of the second compartment
		initialParams[3] = 0.1; // volume fraction of the third compartment
		initialParams[4] = 1.7E-9;// diffusivity
		initialParams[5] = 1.570796326794897;// theta
		initialParams[6] = 1.570796326794897;// phi
		initialParams[7] = 2e-6;// R
		initialParams[8] = 1.7e-9;// diff
		initialParams[9] = 1.570796326794897;// theta
		initialParams[10] = 1.570796326794897;// phi
		initialParams[11] = 1.7e-10;// diffperp
		initialParams[12] = 1.7e-9;// diff
		initialParams[13] = 2e-6;// Rglial

		cm = new CompartmentModel(compNames, initialParams);

		try {
			initLM_Minimizer(NoiseModel
					.getNoiseModel(CL_Initializer.noiseModel));
			((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
			((LM_Minimizer) minimizer).setMAXITER(5000);
		} catch (Exception e) {
			throw new LoggedException(e);
		}

		zcfitter = new ZeppelinCylinderMultiRunLM_Fitter(scheme, 3, 0);

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
				optParams[2] = Math.acos(Math.sqrt(modelParams[2]/(1-modelParams[1]))); // f zep
				optParams[3] = Math.sqrt(modelParams[4]); // diff
				optParams[4] = modelParams[5]; // theta
				optParams[5] = modelParams[6]; // phi
				optParams[6] = (Math.acos(Math
						.sqrt((modelParams[7] - 1e-7) / 20e-6)));// R
				optParams[7] = (Math.asin(Math.sqrt(modelParams[11]
						/ modelParams[4])));// diffprp constrained
				optParams[8] = (Math.acos(Math
						.sqrt((modelParams[13] - 1e-7) / 40e-6)));// Rglial
																	// constrained

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double f = Math.cos(optParams[1]); // volume fraction ;
				f = f * f;
				modelParams[1] = f; // volume fraction
				double zepf = Math.cos(optParams[2]) * Math.cos(optParams[2])* (1 - f);
				modelParams[2] = zepf; // volume fraction of zep
				modelParams[3] = 1 - f - zepf;//volume fraction of the third compartment
				modelParams[4] = optParams[3] * optParams[3]; // diff
				modelParams[5] = optParams[4]; // theta
				modelParams[6] = optParams[5]; // phi
				modelParams[7] = 1e-7 + (Math.cos(optParams[6])
						* Math.cos(optParams[6]) * 20e-6);// R constraining
															// the radius
				modelParams[8] = modelParams[4]; // diff zeppelin
				modelParams[9] = modelParams[5];// theta zep
				modelParams[10] = modelParams[6];// phi zep
				modelParams[11] = Math.sin(optParams[7])
						* Math.sin(optParams[7]);// diff perp constrained
				modelParams[11] = modelParams[11] * modelParams[4];
				modelParams[12] = modelParams[4];
				modelParams[13] = 1e-7 + (Math.cos(optParams[8])
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
				return 5;
			}
		};

	}

	/**
	 * Estimates a starting set of parameters from the zeppelincylinder fitter.
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

		double[] zcparams;

		try {
			zcparams = MultiRunMinimizer.getBestSolution(zcfitter.fit(data));

		} catch (MinimizerException e) {

			zcparams = zcfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = zcparams[1];
		// set mixing parameter
		double f1 = zcparams[2]-0.00005;
		double f2 = zcparams[3]-0.00005;
		double f3 = 0.0001;// f3 is not 0 because single run LM breaks and gives singular matrix -2.
		// set diffusivity
		double diff = zcparams[4];
		double theta = zcparams[5];
		double phi = zcparams[6];
		// R
		double R = zcparams[7];
		// perpendicular diffusivity
		double eigenvalue2 = zcparams[11];
		double Rg = 2e-6;

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;
		modelParams[1] = f1;
		modelParams[2] = f2;
		modelParams[3] = f3;
		modelParams[4] = diff;
		modelParams[5] = theta;
		modelParams[6] = phi;
		modelParams[7] = R;
		modelParams[8] = modelParams[4];// diffzep
		modelParams[9] = modelParams[5];// thetazep
		modelParams[10] = modelParams[6];// phizep
		modelParams[11] = eigenvalue2;// diffperp
		modelParams[12] = modelParams[4];
		modelParams[13] = Rg;

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		ZeppelinCylinderSphereLM_Fitter inv = new ZeppelinCylinderSphereLM_Fitter(
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

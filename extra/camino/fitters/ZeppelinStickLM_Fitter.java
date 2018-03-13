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
 * <dd>Fits the zeppelin and stick model using one run of a Levenburg Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd>The diffusivity parameter is assumed the same in each compartment.
 * Parallel diffusivity and unweighted signal are constrained positive by
 * optimizing their square root. Volume fraction is constrained to [0 1] by
 * optimizing cos^{-1} sqrt(f). Perpendicular diffusivity is constrained to be
 * lower than the parallel diffusivity by optimising sin^(-1)sqrt(dperp/dparal).
 * 
 * </dl>
 * 
 * @author Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 * 
 */
public class ZeppelinStickLM_Fitter extends CompartmentFitter {

	protected LinearDT_Inversion dtfitter;
	protected BallStickLM_GaussianFitter bsfitter;

	private final int numModelParams = 10;
	private final int numOptParams = 6;

	/**
	 * Default constructor required for derived classes.
	 */
	public ZeppelinStickLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized parameters
	 * in the Codec object.
	 * 
	 * @param scheme
	 *            The imaging protocol.
	 */
	public ZeppelinStickLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = new String[2];
		compNames[0] = new String("stick");
		compNames[1] = new String("zeppelin");
		double[] initialParams = new double[10];
		initialParams[0] = 1.0; // the S0
		initialParams[1] = 0.7; // volume fraction of the first compartment
		initialParams[2] = 0.3; // volume fraction of the second compartment
		initialParams[3] = 1.7E-9;// stick diffusivity
		initialParams[4] = 1.570796326794897;// theta
		initialParams[5] = 1.570796326794897;// phi
		initialParams[6] = 1.7e-9;// diff
		initialParams[7] = 1.570796326794897;// theta
		initialParams[8] = 1.570796326794897;// phi
		initialParams[9] = 1.7e-10;// diffperp

		cm = new CompartmentModel(compNames, initialParams);

		try {
            initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
            ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
            ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
	throw new LoggedException(e);
}


		dtfitter = new LinearDT_Inversion(scheme);
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
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); // f stick
				optParams[2] = Math.sqrt(modelParams[3]); // diff stick
				optParams[3] = modelParams[4]; // theta
				optParams[4] = modelParams[5]; // phi
				optParams[5] = (Math.asin(Math.sqrt(modelParams[9]
						/ modelParams[3])));// diffperp constrained

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double stickf = Math.cos(optParams[1]); // volume fraction of
														// stick compartment;
				stickf = stickf * stickf;
				modelParams[1] = stickf; // volume fraction of stick
											// compartment
				modelParams[2] = 1 - stickf; // volume fraction of the
												// extracellular compartment
				modelParams[3] = optParams[2] * optParams[2]; // diff stick
				modelParams[4] = optParams[3]; // theta
				modelParams[5] = optParams[4]; // phi
				modelParams[6] = modelParams[3]; // diff zeppelin
				modelParams[7] = modelParams[4];// theta zep
				modelParams[8] = modelParams[5];// phi zep
				modelParams[9] = Math.sin(optParams[5])
						* Math.sin(optParams[5]);// diff perp constrained
				modelParams[9] = modelParams[9] * modelParams[3];

				return modelParams;
			}

			public RealMatrix getJacobian(double[] optParams) {
			    RealMatrix J = new RealMatrix(numModelParams, numOptParams);
			    J.entries[0][0] = 2.0*optParams[0];
			    J.entries[1][1] = -2.0*Math.sin(optParams[1])*Math.cos(optParams[1]);
			    J.entries[2][1] = 2.0*Math.sin(optParams[1])*Math.cos(optParams[1]);
			    J.entries[3][2] = 2.0*optParams[2];
			    J.entries[4][3] = 1.0;
			    J.entries[5][4] = 1.0;
			    J.entries[6][2] = 2.0*optParams[2];
			    J.entries[7][3] = 1.0;
			    J.entries[8][4] = 1.0;
			    J.entries[9][2] = 2.0*optParams[2]*Math.sin(optParams[5])*Math.sin(optParams[5]);
			    J.entries[9][5] = 2.0*optParams[2]*optParams[2]*Math.sin(optParams[5])*Math.cos(optParams[5]);

			    return J;
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
	 * Estimates a starting set of parameters from the ballstick  and DT fitter. The
	 * stick orientation is the ballstick stick's direction. The volume fraction
	 * is ballstick's truncated to the range [0.1 0.9]. The parallel diffusivity
	 * is the ballstick's diffusivity and perpendicular diffusivity is the
	 * second eiganvalue from the DT fit.
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

		} catch (MinimizerException e) {
			bsparams = bsfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = bsparams[1];
		double theta = bsparams[5];
		double phi = bsparams[6];

		// set mixing parameter in range 0.1 - 0.9
		double f = bsparams[2];

		// set diffusivity 
		double diff = bsparams[4];

		// Do a linear DT fit
		double[] dtparams = dtfitter.invert(data);

		DT dt = new DT(dtparams[2], dtparams[3], dtparams[4], dtparams[5],
				dtparams[6], dtparams[7]);

		double[][] seig = dt.sortedEigenSystem();
		Vector3D v = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);

		// perpendicular diffusivity
		double eigenvalue2 = seig[0][1];

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;
		modelParams[1] = f; // f stick
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		modelParams[2] = 1 - modelParams[1]; // fzep
		modelParams[3] = diff; // diff stick
		modelParams[4] = theta;
		modelParams[5] = phi;
		modelParams[6] = modelParams[3];// diffzep
		modelParams[7] = modelParams[4];// thetazep
		modelParams[8] = modelParams[5];// phizep
		modelParams[9] = eigenvalue2;// diffperp

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		ZeppelinStickLM_Fitter inv = new ZeppelinStickLM_Fitter(
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

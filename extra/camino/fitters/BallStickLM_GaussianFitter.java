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
 * <dd>Fits the ball and stick model using one run of a Levenburg
 * Marquardt and assuming a Gaussian noise model.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each
 * compartment.  Diffusivity and unweighted signal are constrained
 * positive by optimizing their square root.  Volume fraction is
 * constrained to [0 1] by optimizing cos^{-1} sqrt(f).
 * 
 * </dl>
 *
 * @author  laura(panagio@cs.ucl.ac.uk)
 * @version $Id$
 *
 */
public class BallStickLM_GaussianFitter extends CompartmentFitter {

	protected LinearDT_Inversion dtfitter;

	private final int numModelParams = 7;
	private final int numOptParams = 5;

	/**
	 * Default constructor required for derived classes.
	 */
	public BallStickLM_GaussianFitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public BallStickLM_GaussianFitter(DW_Scheme scheme) {

		this(scheme, FitModel.getCompartmentList(CL_Initializer.fitModel),
				FitModel.getCompartmentModelParams(CL_Initializer.fitModel));

	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 * @param compNames the compartment names.
	 * @param initialParams the initial parameters.
	 */
	public BallStickLM_GaussianFitter(DW_Scheme scheme, String[] compNames,
			double[] initialParams) {

		this.scheme = scheme;
		makeCodec();
                cm = new CompartmentModel(compNames, initialParams);

		try {
			minimizer = new LM_GaussianMinimizer(scheme, cm, cod);
			((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
                        // This does not seem necessary, but set high in case
			((LM_Minimizer) minimizer).setMAXITER(5000);
 		} catch (Exception e) {
			throw new LoggedException(e);
		}

		dtfitter = new LinearDT_Inversion(scheme);

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

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double stickf = Math.cos(optParams[1]); //volume fraction of stick compartment;
				stickf = stickf * stickf;

				modelParams[1] = stickf; //  volume fraction of stick compartment
				modelParams[2] = 1 - stickf; // volume fraction of ball compartment
				modelParams[3] = optParams[2] * optParams[2]; // diff stick
				modelParams[4] = optParams[3]; //theta
				modelParams[5] = optParams[4]; //phi
				modelParams[6] = modelParams[3]; // diff ball

				return modelParams;
			}

                public RealMatrix getJacobian(double[] optParams) {
                    RealMatrix J = new RealMatrix(numModelParams, numOptParams);
                    J.entries[0][0] = 2.0*optParams[0];
                    J.entries[1][1] = -2.0*Math.sin(optParams[1])*Math.cos(optParams[1]);
                    J.entries[2][1] = 2.0*Math.sin(optParams[1])*Math.cos(optParams[1]);
                    J.entries[3][2] = 2.0*optParams[2];
                    J.entries[6][2] = 2.0*optParams[2];
                    J.entries[4][3] = 1.0;
                    J.entries[5][4] = 1.0;

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
	 * Estimates a starting set of parameters from the linear estimate
	 * of the diffusion tensor from the log measurements.  The stick
	 * orientation is the DT principal direction.  The volume fraction
	 * is the FA truncated to the range [0.1 0.9].  The diffusivity is
	 * the mean diffusivity of the DT.
	 *
	 * @param data The set of measurements to fit to.
	 *
	 * @return The starting model parameter values.
	 */
	protected double[] getStartPoint(double[] data) {
		
		if (this.getFixedStartPoint(data)!=null)
		{
			return this.getFixedStartPoint(data);
		}

                // Get B0 directly
                double B0 = scheme.geoMeanZeroMeas(data);

		// Do a linear DT fit
		double[] dtparams = dtfitter.invert(data);

		DT dt = new DT(dtparams[2], dtparams[3], dtparams[4], dtparams[5],
				dtparams[6], dtparams[7]);

		double logS0 = dtparams[1];

		double[][] seig = dt.sortedEigenSystem();
		Vector3D v = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);
		double[] tp = Vector3D.thetaPhi(v);

		double theta = tp[0];
		double phi = -tp[1];
		
		// set mixing parameter in range 0.1 - 0.9	
		double f = dt.fa();

		// set diffusivity as mean DT diffusivity
                // Take abs, as can pathologically be negative.
		double diff = Math.abs(dt.trace()) / 3.0;

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = Math.exp(logS0);//B0
		modelParams[1] = f;
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		modelParams[2] = 1 - modelParams[1]; //ADDED
		modelParams[3] = diff;
		modelParams[4] = theta;
		modelParams[5] = phi;
		modelParams[6] = modelParams[3];

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		BallStickLM_GaussianFitter inv = new BallStickLM_GaussianFitter(
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

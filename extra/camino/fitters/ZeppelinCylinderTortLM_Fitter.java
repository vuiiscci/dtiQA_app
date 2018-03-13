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
 * <dd>Fits the zeppelin and cylinder model with tortuosity
 * approximation using one run of a Levenburg
 * Marquardt.  Accommodates any noise model.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each
 * compartment.  Diffusivity and unweighted signal are constrained
 * positive by optimizing their square root.  Volume fraction is
 * constrained to [0 1] by optimizing cos^{-1} sqrt(f).
 * The perpendicular zeppelin diffusivity is (1-f) times parallel
 * diffusivity.
 *
 * The starting point comes from a linear DT fit.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class ZeppelinCylinderTortLM_Fitter extends CompartmentFitter {

	protected LinearDT_Inversion dtfitter;

	private final int numModelParams = 11;
	private final int numOptParams = 6;

	/**
	 * Default constructor required for derived classes.
	 */
	public ZeppelinCylinderTortLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public ZeppelinCylinderTortLM_Fitter(DW_Scheme scheme) {

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
	public ZeppelinCylinderTortLM_Fitter(DW_Scheme scheme, String[] compNames,
			double[] initialParams) {

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

		dtfitter = new LinearDT_Inversion(scheme);

	}

	/**
	 * Creates the Codec object.
	 */
	protected void makeCodec() {
		cod = new Codec() {

			// Bounds on R
			private double minR = 1e-7;
			private double maxR = 20e-6;

			public double[] modelToOpt(double[] modelParams) {
				double[] optParams = new double[numOptParams];
				optParams[0] = Math.sqrt(modelParams[0]); //s0
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); //f
				optParams[2] = Math.sqrt(modelParams[3]); //diff cyl
				optParams[3] = modelParams[4]; //theta
				optParams[4] = modelParams[5]; //phi
                                optParams[5] =(Math.acos(Math.sqrt((modelParams[6] - minR) / maxR)));// R

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double f = Math.cos(optParams[1]); //volume fraction
				f = f*f;
				modelParams[1] = f; //  volume fraction
				modelParams[2] = 1 - f; // volume fraction of zep
				modelParams[3] = optParams[2] * optParams[2]; // diff cyl
				modelParams[4] = optParams[3]; //theta
				modelParams[5] = optParams[4]; //phi
                                modelParams[6] = minR + (Math.cos(optParams[5])* Math.cos(optParams[5]) * maxR);// R 
                                modelParams[7] = modelParams[3]; // diff zeppelin
                                modelParams[8] = modelParams[4];// theta zep
                                modelParams[9] = modelParams[5];// phi zep
                                modelParams[10] = modelParams[2]*modelParams[3];// diff perp zep

				return modelParams;
			}

			
                public RealMatrix getJacobian(double[] optParams) {
                    RealMatrix J = new RealMatrix(numModelParams, numOptParams);
                    J.entries[0][0] = 2.0*optParams[0];
		    double sc1 = Math.sin(optParams[1])*Math.cos(optParams[1]);
                    J.entries[1][1] = -2.0*sc1;
                    J.entries[2][1] = 2.0*sc1;
                    J.entries[3][2] = 2.0*optParams[2];
                    J.entries[4][3] = 1.0;
                    J.entries[5][4] = 1.0;
		    J.entries[6][5] = -2.0*maxR*Math.sin(optParams[5])*Math.cos(optParams[5]);
                    J.entries[7][2] = 2.0*optParams[2];
                    J.entries[8][3] = 1.0;
                    J.entries[9][4] = 1.0;
                    J.entries[10][1] = J.entries[2][1]*optParams[2]*optParams[2];
                    double f = Math.cos(optParams[1]); //volume fraction of cyl compartment;                                                                    
                    J.entries[10][2] = f*f*J.entries[3][2];

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
	 * of the diffusion tensor from the log measurements.  The cylinder
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

		// Use the first a second diffusivities (abs values to
		// avoid pathological settings) as parallel and 
		// perpedicular diffusivity starting points.
		double diff1 = Math.abs(seig[0][0]);
		double diff2 = Math.abs(seig[0][1]);

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = Math.exp(logS0);//B0
		modelParams[1] = f;
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		modelParams[2] = 1 - modelParams[1]; //ADDED
		modelParams[3] = diff1;
		modelParams[4] = theta;
		modelParams[5] = phi;
		modelParams[6] = 1e-6;
		modelParams[7] = modelParams[3];
		modelParams[8] = theta;
		modelParams[9] = phi;
		modelParams[10] = modelParams[2]*modelParams[3];

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		ZeppelinCylinderTortLM_Fitter inv = new ZeppelinCylinderTortLM_Fitter(CL_Initializer.imPars);

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

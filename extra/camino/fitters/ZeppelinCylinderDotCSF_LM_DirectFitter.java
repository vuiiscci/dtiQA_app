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
 * <dd>Fits the zeppelin-cylinder-dot model, with additional
 * compartment for CSF contamination, using one run of a Levenburg
 * Marquardt.  Accommodates any noise model.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each tissue
 * compartment, but different in CSF and fitted here rather than
 * fixed.  Diffusivity and unweighted signal are constrained positive
 * by optimizing their square root.  Volume fractions are constrained
 * to [0 1] by optimizing cos^{-1} sqrt(f). Volume fraction sum
 * constrained to one.
 *
 * The starting point comes from a linear DT fit.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class ZeppelinCylinderDotCSF_LM_DirectFitter extends CompartmentFitter {

	protected LinearDT_Inversion dtfitter;

	private final int numModelParams = 14;
	private final int numOptParams = 10;

	/**
	 * Default constructor required for derived classes.
	 */
	public ZeppelinCylinderDotCSF_LM_DirectFitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public ZeppelinCylinderDotCSF_LM_DirectFitter(DW_Scheme scheme) {

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
	public ZeppelinCylinderDotCSF_LM_DirectFitter(DW_Scheme scheme, String[] compNames,
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
                                optParams[2] = Math.acos(Math.sqrt(modelParams[2]/(1-modelParams[1]))); // f zep
                                optParams[3] = Math.acos(Math.sqrt(modelParams[3]/(1-modelParams[1]-modelParams[2]))); // f dot
				optParams[4] = Math.sqrt(modelParams[5]); //diff 
				optParams[5] = modelParams[6]; //theta
				optParams[6] = modelParams[7]; //phi
                                optParams[7] =(Math.acos(Math.sqrt((modelParams[8] - minR) / maxR)));// R
				// perpendicular zep diff constrained
				// smaller than parallel diffusivity.
                                optParams[8] = (Math.asin(Math.sqrt(modelParams[12] / modelParams[5])));
                                // CSF diffusivity constrained only positive
				optParams[9] = Math.sqrt(modelParams[13]); //diff CSF

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double f = Math.cos(optParams[1]); //volume fraction ;
				f = f*f;
				modelParams[1] = f; //  volume fraction 
                                double zepf = Math.cos(optParams[2]) * Math.cos(optParams[2])*(1-f);
				modelParams[2] = zepf; // volume fraction of zep
                                double dotf = Math.cos(optParams[3]) * Math.cos(optParams[3])*(1-f-zepf);
				modelParams[3] = dotf; // volume fraction of dot
				modelParams[4] = 1 - f - zepf - dotf; // volume fraction of CSF
				modelParams[5] = optParams[4] * optParams[4]; // cyl diff
				modelParams[6] = optParams[5]; //theta
				modelParams[7] = optParams[6]; //phi
                                modelParams[8] = minR + (Math.cos(optParams[7])* Math.cos(optParams[7]) * maxR);// R 
                                modelParams[9] = modelParams[5]; // diff zeppelin
                                modelParams[10] = modelParams[6];// theta zep
                                modelParams[11] = modelParams[7];// phi zep
                                modelParams[12] = Math.sin(optParams[8])
				    * Math.sin(optParams[8]);// diff perp fraction of diff par
                                modelParams[12] = modelParams[12] * modelParams[5];// diff perp zep
				modelParams[13] = optParams[9] * optParams[9]; // CSF diff

				return modelParams;
			}

                        	
                public RealMatrix getJacobian(double[] optParams) {
                    RealMatrix J = new RealMatrix(numModelParams, numOptParams);
                    J.entries[0][0] = 2.0*optParams[0];
		    double sc1 = Math.sin(optParams[1])*Math.cos(optParams[1]);
		    double cc1 = Math.cos(optParams[1])*Math.cos(optParams[1]);
		    double sc2 = Math.sin(optParams[2])*Math.cos(optParams[2]);
		    double cc2 = Math.cos(optParams[2])*Math.cos(optParams[2]);
		    double sc3 = Math.sin(optParams[3])*Math.cos(optParams[3]);
		    double cc3 = Math.cos(optParams[3])*Math.cos(optParams[3]);
                    J.entries[1][1] = -2.0*sc1;
                    J.entries[2][1] = 2.0*sc1*cc2;
                    J.entries[3][1] = 2.0*sc1*cc3*(1.0 - cc2);
                    J.entries[4][1] = 2.0*sc1*(1.0 - cc2 - cc3 + cc2*cc3);
                    J.entries[2][2] = -2.0*(1.0 - cc1)*sc2;
                    J.entries[3][2] = 2.0*cc3*sc2*(1.0 - cc1);
                    J.entries[4][2] = 2.0*sc2*(1.0 - cc1)*(1.0 - cc3);
                    J.entries[3][3] = -2.0*sc3*(1 - cc1)*(1 - cc2);
                    J.entries[4][3] = -J.entries[3][3];
                    J.entries[5][4] = 2.0*optParams[4];
                    J.entries[6][5] = 1.0;
                    J.entries[7][6] = 1.0;
		    J.entries[8][7] = -2.0*maxR*Math.sin(optParams[7])*Math.cos(optParams[7]);
                    J.entries[9][4] = 2.0*optParams[4];
                    J.entries[10][5] = 1.0;
                    J.entries[11][6] = 1.0;
		    J.entries[12][4] = 2.0*optParams[4]*Math.sin(optParams[8])*Math.sin(optParams[8]);
		    J.entries[12][8] = 2.0*optParams[4]*optParams[4]*Math.sin(optParams[8])*Math.cos(optParams[8]);
                    J.entries[13][9] = 2.0*optParams[9];

                    return J;
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
		// start dot fraction close to zero.
		double initDotFraction = 0.01;
		double initCSF_Fraction = 0.01;
		modelParams[2] = 1 - modelParams[1] - initDotFraction - initCSF_Fraction;
		modelParams[3] = initDotFraction; 
		modelParams[4] = initCSF_Fraction; 
		modelParams[5] = diff1;
		modelParams[6] = theta;
		modelParams[7] = phi;
		modelParams[8] = 1e-6;
		modelParams[9] = modelParams[4];
		modelParams[10] = theta;
		modelParams[11] = phi;
		modelParams[12] = diff2;
                // Initialize CSF diffusivity at approximate in-vivo value.
		modelParams[13] = 3e-9;

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		ZeppelinCylinderDotCSF_LM_DirectFitter inv = new ZeppelinCylinderDotCSF_LM_DirectFitter(CL_Initializer.imPars);

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

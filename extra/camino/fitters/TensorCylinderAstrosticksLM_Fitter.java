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
 * <dd>Fits the tensorcylinderastrosticks model using one run of a Levenburg
 * Marquardt.
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
 * @author  Laura
 * @version $Id$
 *
 */
public class TensorCylinderAstrosticksLM_Fitter extends CompartmentFitter {

	protected TensorCylinderLM_Fitter tcfitter;

	private final int numModelParams = 15;
	private final int numOptParams = 10;

	/**
	 * Default constructor required for derived classes.
	 */
	public TensorCylinderAstrosticksLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public TensorCylinderAstrosticksLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = new String[3];
		compNames[0] = new String("cylindergpd");
		compNames[1] = new String("tensor");
		compNames[2] = new String("astrosticks");
		double[] initialParams = new double[15];
		initialParams[0] = 1.0; //the S0
		initialParams[1] = 0.6; //volume fraction of the first compartment
		initialParams[2] = 0.3; //volume fraction of the second compartment
		initialParams[3] = 0.1; //volume fraction of the third compartment
		initialParams[4] = 1.7E-9;//diffusivity
		initialParams[5] = 1.570796326794897;//theta
		initialParams[6] = 1.570796326794897;//phi
		initialParams[7] = 2e-6;//R
		initialParams[8] = 1.7e-9;//diff
		initialParams[9] = 1.570796326794897;//theta
		initialParams[10] = 1.570796326794897;//phi
		initialParams[11] = 1.7e-10;//diffperp1
		initialParams[12] = 1.7e-11;//diffperp2
		initialParams[13] = 1.5;//alpha
		initialParams[14] = 1.7e-9;//diff astro

		 cm = new CompartmentModel(compNames, initialParams);

		 try {
             initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
             ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
             ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
		throw new LoggedException(e);
	}

		tcfitter = new TensorCylinderMultiRunLM_Fitter(scheme, 3, 0);

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
				optParams[2] = Math.acos(Math.sqrt(modelParams[2]/(1-modelParams[1]))); // f tens
				optParams[3] = Math.sqrt(modelParams[4]); //diff cyl
				optParams[4] = modelParams[5]; //theta
				optParams[5] = modelParams[6]; //phi
				optParams[6] = (Math.acos(Math
						.sqrt((modelParams[7] - 1e-7) / 20e-6)));// R
				optParams[7] = (Math.asin(Math.sqrt(modelParams[11]
						/ modelParams[4])));//diffprp constrained
				optParams[8] = (Math.asin(Math.sqrt(modelParams[12]
						/ modelParams[11])));
				optParams[9] = modelParams[13];//alpha

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double f = Math.cos(optParams[1]); // volume fraction ;
				f = f * f;
				modelParams[1] = f; // volume fraction
				double tensf = Math.cos(optParams[2]) * Math.cos(optParams[2])* (1 - f);
				modelParams[2] = tensf; // volume fraction of tensor
				modelParams[3] = 1 - f - tensf;//volume fraction of the third compartment
				modelParams[4] = optParams[3] * optParams[3]; // diff stick
				modelParams[5] = optParams[4]; //theta
				modelParams[6] = optParams[5]; //phi
				modelParams[7] = 1e-7 + (Math.cos(optParams[6])
						* Math.cos(optParams[6]) * 20e-6);// R constraining
				modelParams[8] = modelParams[4]; // diff 
				modelParams[9] = modelParams[5];//theta 
				modelParams[10] = modelParams[6];//phi 
				modelParams[11] = Math.sin(optParams[7])
						* Math.sin(optParams[7]);//diff perp constrained
				modelParams[11] = modelParams[11] * modelParams[4];//diff perp constrained
				modelParams[12] = Math.sin(optParams[8])
						* Math.sin(optParams[8]);//diff perp constrained
				modelParams[12] = modelParams[12] * modelParams[11];//diff perp2
				modelParams[13] = optParams[9];//alpha 
				modelParams[14] = modelParams[4];

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
	 * Estimates a starting set of parameters from the tensorcylinder model.
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

		double[] tcparams;

		try {
			tcparams = MultiRunMinimizer.getBestSolution(tcfitter.fit(data));

		} catch (MinimizerException e) {
			tcparams = tcfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = tcparams[1];
		double f1 = tcparams[2]-0.00005;
		double f2 = tcparams[3]-0.00005;
		double f3 = 0.0001;// f3 is not 0 because single run LM breaks and gives singular matrix -2.
		double diff = tcparams[4];
		double theta = tcparams[5];
		double phi = tcparams[6];
		double R = tcparams[7];
		double diffperp1 = tcparams[11];
		double diffperp2 = tcparams[12];
		double alpha = tcparams[13];

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;
		modelParams[1] = f1;
		modelParams[2] = f2;
		modelParams[3] = f3;
		modelParams[4] = diff; //diff stick
		modelParams[5] = theta;
		modelParams[6] = phi;
		modelParams[7] = R;
		modelParams[8] = modelParams[4];//diff
		modelParams[9] = modelParams[5];//theta
		modelParams[10] = modelParams[6];//phi
		modelParams[11] = diffperp1;//diffperp1 
		modelParams[12] = diffperp2;//diffperp2 
		modelParams[13] = alpha;//alpha 
		modelParams[14] = modelParams[4];

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		TensorCylinderAstrosticksLM_Fitter inv = new TensorCylinderAstrosticksLM_Fitter(
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

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
 * <dd>Fits the tensor and cylinder model using one run of a Levenburg
 * Marquardt.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each
 * compartment.  Parallel diffusivity and the unweighted signal  are constrained
 * positive by optimizing their square root.  Volume fraction is
 * constrained to [0 1] by optimizing cos^{-1} sqrt(f). The radius is constrained to 
 * 0.1-20 microns by optimizing cos^{-1} sqrt(R-1e-7/20e-6). Perpendicular diffusivity1 
 * is constrained  by optimizing sin^{-1}sqrt(dperp1/dpar) and perpendicular diffusivity2
 * is constrained  by optimizing sin^{-1}sqrt(dperp2/dperp1).
 * 
 * </dl>
 *
 * @author  Laura (panagio@cs.ucl.ac.uk)
 * @version $Id$
 *
 */
public class TensorCylinderLM_Fitter extends CompartmentFitter {

	protected BallCylinderLM_Fitter bcfitter;
	protected TensorStickLM_Fitter tsfitter;

	private final int numModelParams = 13;
	private final int numOptParams = 9;

	/**
	 * Default constructor required for derived classes.
	 */
	public TensorCylinderLM_Fitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public TensorCylinderLM_Fitter(DW_Scheme scheme) {

		this.scheme = scheme;
		makeCodec();

		String[] compNames = new String[2];
		compNames[0] = new String("cylindergpd");
		compNames[1] = new String("tensor");
		double[] initialParams = new double[13];
		initialParams[0] = 1.0; //the S0
		initialParams[1] = 0.7; //volume fraction of the first compartment
		initialParams[2] = 0.3; //volume fraction of the second compartment
		initialParams[3] = 1.7E-9;//diffusivity
		initialParams[4] = 1.570796326794897;//theta
		initialParams[5] = 1.570796326794897;//phi
		initialParams[6] = 2e-6;//R
		initialParams[7] = 1.7e-9;//diff
		initialParams[8] = 1.570796326794897;//theta
		initialParams[9] = 1.570796326794897;//phi
		initialParams[10] = 1.7e-10;//diffperp1
		initialParams[11] = 1.7e-11;//diffperp2
		initialParams[12] = 1.5;//alpha

		 cm = new CompartmentModel(compNames, initialParams);

		 try {
             initLM_Minimizer(NoiseModel.getNoiseModel(CL_Initializer.noiseModel));
             ((LM_Minimizer) minimizer).setCONVERGETHRESH(1e-8);
             ((LM_Minimizer) minimizer).setMAXITER(5000);
	} catch (Exception e) {
		throw new LoggedException(e);
	}

		bcfitter = new BallCylinderMultiRunLM_Fitter(scheme, 1, 0);
		tsfitter = new TensorStickMultiRunLM_Fitter(scheme, 1, 0);
	}

	/**
	 * Creates the Codec object.
	 */
	protected void makeCodec() {
		cod = new Codec() {

			public double[] modelToOpt(double[] modelParams) {
				double[] optParams = new double[numOptParams];
				optParams[0] = Math.sqrt(modelParams[0]); //s0
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); //f 
				optParams[2] = Math.sqrt(modelParams[3]); //diff 
				optParams[3] = modelParams[4]; //theta
				optParams[4] = modelParams[5]; //phi
				optParams[5] = (Math.acos(Math.sqrt((modelParams[6] - 1e-7) / 20e-6)));// R
				optParams[6] = (Math.asin(Math.sqrt(modelParams[10]
						/ modelParams[3])));//diffprp1 constrained
				optParams[7] = (Math.asin(Math.sqrt(modelParams[11]
						/ modelParams[10])));//diffprp2 constrained
				optParams[8] = modelParams[12];//alpha

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0] * optParams[0];
				double cylindf = Math.cos(optParams[1]);
				cylindf = cylindf * cylindf;

				modelParams[1] = cylindf; // volume fraction of intracellular compartment
				modelParams[2] = 1 - cylindf; //volume fraction of the extracellular compartment
				modelParams[3] = optParams[2] * optParams[2]; // diff 
				modelParams[4] = optParams[3]; //theta
				modelParams[5] = optParams[4]; //phi
				modelParams[6] =  1e-7 + (Math.cos(optParams[5])* Math.cos(optParams[5]) * 20e-6);// R constraining//R
				modelParams[7] = modelParams[3]; // diff 
				modelParams[8] = modelParams[4];//theta 
				modelParams[9] = modelParams[5];//phi 
				modelParams[10] = Math.sin(optParams[6])
						* Math.sin(optParams[6]);
				modelParams[10] = modelParams[10] * modelParams[3];//diff perp1 constrained
				modelParams[11] = Math.sin(optParams[7])
						* Math.sin(optParams[7]);
				modelParams[11] = modelParams[11] * modelParams[10];//diff perp2 constrained
				modelParams[12] = optParams[8];//alpha 

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
	 * Estimates a starting set of all parameters from the tensor stick model 
	 * except the radius which is estimated from the ballcylinder model. 
	 *
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

		// Do tensorstick fit.
		double[] tsparams;

		try {
			tsparams = MultiRunMinimizer.getBestSolution(tsfitter.fit(data));

		} catch (MinimizerException e) {
			tsparams = tsfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double S0 = tsparams[1];
		// set mixing parameter in range 0.1 - 0.9	
		double f = tsparams[2];
		double diff = tsparams[4];
		double theta = tsparams[5];
		double phi = tsparams[6];
		double diffperp1 = tsparams[10];
		double diffperp2 = tsparams[11];
		double alpha = tsparams[12];

		// Do a  BallCylinder fit
		double[] bcparams;

		try {
			bcparams = MultiRunMinimizer.getBestSolution(bcfitter.fit(data));

		} catch (MinimizerException e) {
			bcparams = bcfitter.getStartPoint(data);
			e.printStackTrace();
		}

		double R = bcparams[7];

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = S0;
		modelParams[1] = f; //f stick
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		modelParams[2] = 1 - modelParams[1]; //ftensor
		modelParams[3] = diff; //diff stick
		modelParams[4] = theta;
		modelParams[5] = phi;
		modelParams[6] = R;
		modelParams[7] = modelParams[3];//diff
		modelParams[8] = modelParams[4];//theta
		modelParams[9] = modelParams[5];//phi
		modelParams[10] = diffperp1;//diffperp1 
		modelParams[11] = diffperp2;//diffperp2 
		modelParams[12] = alpha;//alpha 

		return modelParams;

	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		TensorCylinderLM_Fitter inv = new TensorCylinderLM_Fitter(
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

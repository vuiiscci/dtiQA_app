package fitters;

import models.*;
import optimizers.*;
import inverters.*;
import numerics.*;
import misc.*;
import data.*;
import tools.*;
import imaging.*;
import models.compartments.CompartmentModel;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fits the minimal model of white matter diffusion, as described
 * in Alexander et al NIMG 2010, with in-vivo settings using one run
 * of a Levenburg Marquardt with any specified noise model.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each
 * compartment.  Radius and unweighted signal are constrained positive
 * by optimizing their square root.  White matter intrinsic
 * diffusivity is fixed at 1.7E-9 m^2/s and CSF diffusivity to 3.0E-9
 * m^2/s.  Volume fraction is constrained to [0 1] by optimizing
 * cos^{-1} sqrt(f).
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class MMWMD_InVivoLM_DirectFitter extends CompartmentFitter {

	protected LinearDT_Inversion dtfitter;

	private final int numModelParams = 13;
	private final int numOptParams = 6;

        // Fixed diffusivities of two compartments.
        protected double WMDIFF=1.7E-9;
        protected double CSFDIFF=3.0E-9;

	/**
	 * Default constructor required for derived classes.
	 */
	public MMWMD_InVivoLM_DirectFitter() {
	}

	/**
	 * Constructor implements the mapping between model and optimized
	 * parameters in the Codec object.
	 *
	 * @param scheme The imaging protocol.
	 */
	public MMWMD_InVivoLM_DirectFitter(DW_Scheme scheme) {

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
	public MMWMD_InVivoLM_DirectFitter(DW_Scheme scheme, String[] compNames,
			double[] initialParams) {

            if(CL_Initializer.mmwmddiff>0)
		WMDIFF = CL_Initializer.mmwmddiff;

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
				optParams[1] = Math.acos(Math.sqrt(modelParams[1])); //f cylinder
				optParams[2] = Math.acos(Math.sqrt(modelParams[2]/(1-modelParams[1]))); //f zeppelin
				optParams[3] = modelParams[5]; //theta
				optParams[4] = modelParams[6]; //phi
                                optParams[5] =(Math.acos(Math.sqrt((modelParams[7] - minR) / maxR)));// R

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];

                                // S0
				modelParams[0] = optParams[0] * optParams[0];

                                // Volume fractions.
				double sqrtf1 = Math.cos(optParams[1]);
                                //  Volume fraction of cylinder compartment
				modelParams[1] = sqrtf1*sqrtf1;
                                //  Volume fraction of zeppelin compartment
				double sqrtf2 = Math.cos(optParams[2]); 
				modelParams[2] = (1-modelParams[1])*sqrtf2*sqrtf2;
                                //  Volume fraction of CSF compartment
				modelParams[3] = 1 - modelParams[1] - modelParams[2];

                                // Cylinder compartment
                                // Intrinsic ICS diff
				modelParams[4] = WMDIFF;
                                // Cylinder orientation
				modelParams[5] = optParams[3]; //theta
				modelParams[6] = optParams[4]; //phi
                                // Cylinder radius
                                modelParams[7] = minR + (Math.cos(optParams[5])* Math.cos(optParams[5]) * maxR);// R 

                                // Zeppelin compartment
                                // Parallel diff
				modelParams[8] = WMDIFF;
                                // Parallel orientation
				modelParams[9] = optParams[3]; //theta
				modelParams[10] = optParams[4]; //phi
                                // Perpendicular diff
				modelParams[11] = WMDIFF*modelParams[2]/(modelParams[1]+modelParams[2]);

                                // CSF Compartment (ball)
				modelParams[12] = CSFDIFF;

				return modelParams;
			}

                                                
                public RealMatrix getJacobian(double[] optParams) {
                    RealMatrix J = new RealMatrix(numModelParams, numOptParams);
                    J.entries[0][0] = 2.0*optParams[0];
		    double sc1 = Math.sin(optParams[1])*Math.cos(optParams[1]);
		    double cc1 = Math.cos(optParams[1])*Math.cos(optParams[1]);
		    double sc2 = Math.sin(optParams[2])*Math.cos(optParams[2]);
		    double cc2 = Math.cos(optParams[2])*Math.cos(optParams[2]);
                    J.entries[1][1] = -2.0*sc1;
                    J.entries[2][1] = 2.0*sc1*cc2;
                    J.entries[3][1] = 2.0*sc1*(1.0 - cc2);
                    J.entries[2][2] = -2.0*(1.0 - cc1)*sc2;
                    J.entries[3][2] = 2.0*sc2*(1.0 - cc1);
                    J.entries[5][3] = 1.0;
                    J.entries[6][4] = 1.0;
		    J.entries[7][5] = -2.0*maxR*Math.sin(optParams[5])*Math.cos(optParams[5]);
                    J.entries[9][3] = 1.0;
                    J.entries[10][4] = 1.0;
                    double f1 = cc1;
                    double f2 = (1.0 - cc1)*cc2;
                    double f1pf2 = f1 + f2;
		    J.entries[11][1] = -WMDIFF*(f1pf2*J.entries[1][1] - f1*(J.entries[1][1] + J.entries[2][1]))/(f1pf2*f1pf2);
		    J.entries[11][2] = WMDIFF*J.entries[2][2]*f1/(f1pf2*f1pf2);

                    return J;
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

 		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = Math.exp(logS0);
		modelParams[1] = f;
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		double initCSF_Fraction = 0.01;
		modelParams[2] = 1 - modelParams[1] - initCSF_Fraction;
		modelParams[3] = initCSF_Fraction; 
		modelParams[4] = WMDIFF;
		modelParams[5] = theta;
		modelParams[6] = phi;
		modelParams[7] = 1E-6;
		modelParams[8] = WMDIFF;
		modelParams[9] = theta;
		modelParams[10] = phi;
		modelParams[11] = modelParams[2]*WMDIFF;
		modelParams[12] = CSFDIFF;

                return modelParams;
	}

	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		MMWMD_InVivoLM_DirectFitter inv = new MMWMD_InVivoLM_DirectFitter(
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

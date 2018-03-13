package fitters;

import java.util.logging.Logger;
import java.util.Random;
import java.io.*;
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
 * <dd>Fits the ball and stick model using MCMC and with a specified
 * noise model.
 * 
 * <dt>Description:
 * 
 * <dd> The diffusivity parameter is assumed the same in each
 * compartment.  Diffusivity and unweighted signal are constrained
 * positive by the prior.  Volume fraction is constrained to [0 1] by
 * the prior.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class BallStickMCMC_Fitter extends MCMC_Fitter {

	/**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.fitters.BallStickMCMC_Fitter");
	
	protected LinearDT_Inversion dtfitter;

	private final int numModelParams = 7;
	private final int numOptParams = 5;


	/**
	 * Default constructor required for derived classes.
	 */
	public BallStickMCMC_Fitter() {
	}


	/**
	 * Basic constructor extracts model particulars and calls
         * full constructor.
         *
         * @param scheme The imaging protocol.
	 */
	public BallStickMCMC_Fitter(DW_Scheme scheme) {
            this(scheme, FitModel.getCompartmentList(CL_Initializer.fitModel), FitModel.getCompartmentModelParams(CL_Initializer.fitModel));
	}


	/**
	 * Constructor sets up various required objects for MCMC.
	 *
	 * @param scheme The imaging protocol.
	 * @param compNames the compartment names.
	 * @param initialParams the initial parameters.

	 */
        public BallStickMCMC_Fitter(DW_Scheme scheme, String[] compNames, double[] initialParams) {
            super(scheme, compNames, initialParams);

            // Initialize the dtfitter for starting point
            // estimation.
            dtfitter = new LinearDT_Inversion(scheme);
	 }


	/**
	 * Creates the Codec object. The prior encodes bounds, so
         * we use a much simpler mapping than for other fitters. 
         * This keeps the trial distributions simple.
	 */
	protected void makeCodec() {
		cod = new Codec() {

			public double[] modelToOpt(double[] modelParams) {
				double[] optParams = new double[numOptParams];
				optParams[0] = modelParams[0]; //s0
				optParams[1] = modelParams[1]; //f stick
				optParams[2] = modelParams[3]; //diff stick
				optParams[3] = modelParams[4]; //theta
				optParams[4] = modelParams[5]; //phi

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0]; //s0
				double stickf = optParams[1]; //volume fraction of stick compartment;

				modelParams[1] = stickf; //  volume fraction of stick compartment
				modelParams[2] = 1 - stickf; // volume fraction of ball compartment
				modelParams[3] = optParams[2]; // diff stick
				modelParams[4] = optParams[3]; //theta
				modelParams[5] = optParams[4]; //phi
				modelParams[6] = modelParams[3]; // diff ball

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
	 * Creates the Prior object.
	 */
	protected void makePrior() {
            prior = new Prior() {

                    public double prior(double[] optParams) {

                        // s0 must be positive.
                        double p1 = optParams[0]<0?0:1;
                        
                        // f must be in [0 1].
                        double p2 = (optParams[1]<0 || optParams[1]>1)?0:1;
                        
                        // diffusivity must be positive.
                        double p3 = optParams[2]<0?0:1;
                        
                        // the prior on the direction is sin(theta)
                        // to ensure even sampling over the whole
                        // surface.
                        double p4 = Math.abs(Math.sin(optParams[3]));
                        
                        return p1*p2*p3*p4;
                    }
                    
		};

	}


	/**
	 * Creates the Perturbation object.
	 */
	protected void makePerturbation() {
            perturbation = new Perturbation() {

                    public double[] perturb(double[] params, double[] paramsInitial, Random rand) {
		
                        double[] paramsNew = new double[params.length];
		
			// Scale of perturbation compared to parameter scale
			// For non-directional parameters
			double nonDirFraction = 0.05;

			// Perturbation scale for directional parameters
			double dirFraction = 0.001;

			// S0
			paramsNew[0] = params[0] + paramsInitial[0]*nonDirFraction*rand.nextGaussian();

			// volume fraction
			paramsNew[1] = params[1] + nonDirFraction*rand.nextGaussian();

			// diffusivity
			paramsNew[2] = params[2] + paramsInitial[2]*nonDirFraction*rand.nextGaussian();

			// orientation
			paramsNew[3] = params[3] + dirFraction*rand.nextGaussian();
			paramsNew[4] = params[4] + dirFraction*rand.nextGaussian();


                        return paramsNew;
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

		// find angle that rotates this to the XZ plane (phi), then
		// onto x axis (theta)
		//double theta = Math.PI / 2.0 - tp[0];
		double theta = tp[0];
		double phi = -tp[1];

		// set mixing parameter in range 0.1 - 0.9	
		double f = dt.fa();

		// set diffusivity as mean DT diffusivity	
		double diff = dt.trace() / 3.0;

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = Math.exp(logS0);
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

		BallStickMCMC_Fitter inv = new BallStickMCMC_Fitter(CL_Initializer.imPars);


		// Loop over the voxels.
		int step =1;
		while (CL_Initializer.data.more()) {
			try {
				
				long start = System.currentTimeMillis();
				double[] nextVoxel = CL_Initializer.data.nextVoxel();
				double[][] fit = inv.fit(nextVoxel);

                                for(int i=0; i<fit.length; i++)
                                    om.output(fit[i]);
				long end = System.currentTimeMillis();
				
				logger.info(" current step is "+step +", time of minimization is "+(end-start)+" ms.");
				step++;
				
			} catch (Exception e) {
				System.err.println(e);
			}
		}

		om.close();
	}

}

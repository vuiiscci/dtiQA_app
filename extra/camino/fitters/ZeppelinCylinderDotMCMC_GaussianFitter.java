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
 * <dd>Fits the ball and stick model using MCMC and assuming a
 * Gaussian noise model.
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
public class ZeppelinCylinderDotMCMC_GaussianFitter extends MCMC_Fitter {

        /**
         * Logging object
         */
        private static Logger logger = Logger.getLogger("camino.fitters.ZeppelinCylinderDotMCMC_GaussianFitter");
	
	protected LinearDT_Inversion dtfitter;

	private final int numModelParams = 12;
	private final int numOptParams = 8;


	/**
	 * Default constructor required for derived classes.
	 */
	public ZeppelinCylinderDotMCMC_GaussianFitter() {
	}


	/**
	 * Basic constructor extracts model particulars and calls
         * full constructor.
         *
         * @param scheme The imaging protocol.
	 */
	public ZeppelinCylinderDotMCMC_GaussianFitter(DW_Scheme scheme) {
            this(scheme, FitModel.getCompartmentList(CL_Initializer.fitModel), FitModel.getCompartmentModelParams(CL_Initializer.fitModel));
	}


	/**
	 * Constructor sets up various required objects for MCMC.
	 *
	 * @param scheme The imaging protocol.
	 * @param compNames the compartment names.
	 * @param initialParams the initial parameters.

	 */
        public ZeppelinCylinderDotMCMC_GaussianFitter(DW_Scheme scheme, String[] compNames, double[] initialParams) {
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
				optParams[1] = modelParams[1]; //f cyl
				optParams[2] = modelParams[2]; //f zep
				optParams[3] = modelParams[4]; //par diff 
				optParams[4] = modelParams[5]; //theta
				optParams[5] = modelParams[6]; //phi
				optParams[6] = modelParams[7]; //radius
				optParams[7] = modelParams[11]; //perp diff

				return optParams;
			}

			public double[] optToModel(double[] optParams) {
				double[] modelParams = new double[numModelParams];
				modelParams[0] = optParams[0]; //s0
				modelParams[1] = optParams[1]; //f cyl
				modelParams[2] = optParams[2]; //f zep
				modelParams[3] = 1 - optParams[1] - optParams[2]; //f dot
				modelParams[4] = optParams[3]; // diff cyl
				modelParams[5] = optParams[4]; //theta cyl
				modelParams[6] = optParams[5]; //phi cyl
				modelParams[7] = optParams[6]; //radius
				modelParams[8] = optParams[3]; //diff zep par
				modelParams[9] = optParams[4]; //theta zep
				modelParams[10] = optParams[5]; //phi zep
				modelParams[11] = optParams[7]; //diff zep perp

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
	 * Creates the Prior object.
	 */
	protected void makePrior() {
            prior = new Prior() {

                    public double prior(double[] optParams) {

                        // s0 must be positive.
                        double p1 = optParams[0]<0?0:1;
                        
                        // sum of fs must be in [0 1].
                        double p2 = (optParams[1]<0 || optParams[1]>1)?0:1;
                        double p3 = (optParams[2]<0 || optParams[2]>1-optParams[1])?0:1;
                       
                        // diffusivity must be positive.
                        double p4 = optParams[3]<0?0:1;
                        
                        // perp diffusivity must be positive and less
                        // than parallel.
                        double p5 = (optParams[7]<0 || optParams[7]>1-optParams[3])?0:1;
                        
                        // radius must be positive.
                        double p6 = optParams[6]<0?0:1;
                        
                        // the prior on the direction is sin(theta)
                        // to ensure even sampling over the whole
                        // surface.
                        double p7 = Math.abs(Math.sin(optParams[4]));
                        
                        return p1*p2*p3*p4*p5*p6*p7;
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

			// volume fractions
			paramsNew[1] = params[1] + nonDirFraction*rand.nextGaussian();
			paramsNew[2] = params[2] + nonDirFraction*rand.nextGaussian();

			// diffusivities
			paramsNew[3] = params[3] + paramsInitial[3]*nonDirFraction*rand.nextGaussian();
			paramsNew[7] = params[7] + paramsInitial[7]*nonDirFraction*rand.nextGaussian();

			// radius
			paramsNew[6] = params[6] + paramsInitial[6]*nonDirFraction*rand.nextGaussian();
			// orientation
			paramsNew[4] = params[4] + dirFraction*rand.nextGaussian();
			paramsNew[5] = params[5] + dirFraction*rand.nextGaussian();


                        return paramsNew;
                    }
		};

	}


	/**
	 * Estimates a starting set of parameters from the linear
	 * estimate of the diffusion tensor from the log measurements.
	 * The cylinder orientation is the DT principal direction.
	 * The cylinder volume fraction is the FA truncated to the
	 * range [0.1 0.9].  The dot fraction is zero.  The radius is
	 * 1 micron.  The parallel diffusivity is the largest
	 * eigenvalue of the DT; perpendicular is the smallest.
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

		DT dt = new DT(dtparams[2], dtparams[3], dtparams[4], dtparams[5], dtparams[6], dtparams[7]);

		double logS0 = dtparams[1];

		double[][] seig = dt.sortedEigenSystem();

                // eigenvalues
                double eigmax = seig[0][0];
                double eigmin = seig[0][2];

                // principal eigenvector
		Vector3D v = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);
		double[] tp = Vector3D.thetaPhi(v);
		double theta = tp[0];
		double phi = -tp[1];

		// set mixing parameter in range 0.1 - 0.9	
		double f = dt.fa();

		// Construct the array of model parameters.
		double[] modelParams = new double[numModelParams];
		modelParams[0] = Math.exp(logS0); //s0
		modelParams[1] = f; //f cyl
		modelParams[1] = modelParams[1] < 0.1 ? 0.1 : modelParams[1];
		modelParams[1] = modelParams[1] > 0.9 ? 0.9 : modelParams[1];
		modelParams[2] = 1-f; //f zep
		modelParams[3] = 0; //f dot
		modelParams[4] = eigmax; //diff cyl
		modelParams[5] = theta; //theta cyl
		modelParams[6] = phi; //phi cyl
		modelParams[7] = 1E-6; //radius
		modelParams[8] = eigmax; //diff zep par
		modelParams[9] = theta; //theta cyl
		modelParams[10] = phi; //phi cyl
		modelParams[11] = eigmin; //diff zep perp

		return modelParams;

	}

	
	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		ZeppelinCylinderDotMCMC_GaussianFitter inv = new ZeppelinCylinderDotMCMC_GaussianFitter(CL_Initializer.imPars);


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

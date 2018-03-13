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
 * <dd>Fits the MMWMD fixed tissue model using MCMC and with a
 * specified noise model.
 * 
 * <dt>Description:
 * 
 * <dd> Compartments are:
 * Cylinder
 * Zeppelin
 * Dot
 * Ball (for CSF)
 *
 * Parameters are:
 * S0 - constrained positive by prior
 * f1 - cylinder volume fraction; constrained to [0 1].
 * f2 - zeppelin volume fraction; constrained to [0 1-f1].
 * f3 - dot volume fraction; constrained to [1 1-f1-f2].
 * f4 - CSF volume fraction; equals 1-f1-f2-f3.
 * d - intrinsic diffusivity within cylinder; fixed at 1.7E-9m^2/s.
 * theta - colatitude of cylinder orientation.
 * phi - longitude of cylinder orientation.
 * R - cylinder radius; constrained to [1E-7 2E-5] microns by prior.
 * d - zeppelin parallel diffusivity; fixed at 1.7E-9m^2/s.
 * theta - colatitude of zeppelin orientation; same as cylinder.
 * phi - longitude of zeppelin orientation; same as cylinder.
 * d_perp - perpendicular zeppelin diffusivity; equals (1-f1/(f1+f2))d.
 * d_csf - CSF diffusivity; fixed at 3.0E-9m^2/s.
 * 
 * The MCMC varies only the proportions of f1 and f2 with f3 and f4
 * held fixed and R.  All other parameters are fixed at their starting
 * values, which come from a maximum likelihood gradient descent with
 * multiple starting points.
 *
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class MMWMD_FixedMCMC_Fitter extends MCMC_Fitter {

	/**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.fitters.MMWMD_InVivoMCMC_Fitter");
	
    protected MMWMD_FixedMultiRunLM_DirectFitter gdfitter;
    protected int GDRUNS = 1;

    private final int numModelParams = 14;
    private final int numOptParams = 7;

    // Fixed diffusivities of two compartments.
    protected double WMDIFF=0.6E-9;
    protected double CSFDIFF=2.0E-9;


	/**
	 * Default constructor required for derived classes.
	 */
	public MMWMD_FixedMCMC_Fitter() {
	}


	/**
	 * Basic constructor extracts model particulars and calls
         * full constructor.
         *
         * @param scheme The imaging protocol.
	 */
	public MMWMD_FixedMCMC_Fitter(DW_Scheme scheme) {
            this(scheme, FitModel.getCompartmentList(CL_Initializer.fitModel), FitModel.getCompartmentModelParams(CL_Initializer.fitModel));
	}


	/**
	 * Constructor sets up various required objects for MCMC.
	 *
	 * @param scheme The imaging protocol.
	 * @param compNames the compartment names.
	 * @param initialParams the initial parameters.

	 */
        public MMWMD_FixedMCMC_Fitter(DW_Scheme scheme, String[] compNames, double[] initialParams) {
            super(scheme, compNames, initialParams);

            if(CL_Initializer.mmwmddiff>0)
		WMDIFF = CL_Initializer.mmwmddiff;

	    // Initialize the gradient descent fitter for starting point
            // estimation.
            gdfitter = new MMWMD_FixedMultiRunLM_DirectFitter(scheme, GDRUNS, CL_Initializer.seed);
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
			    optParams[1] = modelParams[1]; //f cylinder
			    optParams[2] = modelParams[2]; //f zeppelin
			    optParams[3] = modelParams[3]; //f dot
			    optParams[4] = modelParams[6]; //theta
			    optParams[5] = modelParams[7]; //phi
			    optParams[6] = modelParams[8];// R

			    return optParams;
                        }

                        public double[] optToModel(double[] optParams) {
			    double[] modelParams = new double[numModelParams];

			    // S0
			    modelParams[0] = optParams[0];

			    // Volume fractions.
			    modelParams[1] = optParams[1];
			    modelParams[2] = optParams[2];
			    modelParams[3] = optParams[3];
			    modelParams[4] = 1 - modelParams[1] - modelParams[2] - modelParams[3];

			    // Cylinder compartment
			    // Intrinsic ICS diff
			    modelParams[5] = WMDIFF;
			    // Cylinder orientation
			    modelParams[6] = optParams[4]; //theta
			    modelParams[7] = optParams[5]; //phi
			    // Cylinder radius
			    modelParams[8] = optParams[6];// R 

			    // Zeppelin compartment
			    // Parallel diff
			    modelParams[9] = WMDIFF;
			    // Parallel orientation
			    modelParams[10] = optParams[4]; //theta
			    modelParams[11] = optParams[5]; //phi
			    // Perpendicular diff
			    modelParams[12] = WMDIFF*modelParams[2]/(modelParams[1]+modelParams[2]);

			    // CSF Compartment (ball)
			    modelParams[13] = CSFDIFF;

			    return modelParams;
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
	 * Creates the Prior object.
	 */
	protected void makePrior() {
            prior = new Prior() {

		    // Bounds on R
		    private double minR = 1e-7;
		    private double maxR = 20e-6;

                    public double prior(double[] optParams) {

                        // s0 must be positive.
                        double p1 = optParams[0]<0?0:1;
                        
                        // f1 must be in [0 1].
                        double p2 = (optParams[1]<0 || optParams[1]>1)?0:1;
                        
                        // f2 must be in [0 1-f1].
                        double p3 = (optParams[2]<0 || optParams[2]>(1-optParams[1]))?0:1;
                        
                        // f3 must be in [0 1-f1-f2].
                        double p4 = (optParams[3]<0 || optParams[3]>(1-optParams[1]-optParams[2]))?0:1;
                        
                        // the prior on the direction is sin(theta)
                        // to ensure even sampling over the whole
                        // surface.
                        double p5 = Math.abs(Math.sin(optParams[4]));

			// R constrained to realistic range
                        double p6 = (optParams[6]<minR || optParams[6]>maxR)?0:1;
                                                
                        return p1*p2*p3*p4*p5*p6;
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
                        // In the original matlab implementation, this
                        // was lower at 0.05, but convergence is better
                        // at this value.
			double nonDirFraction = 0.2;

			// Perturbation scale for directional parameters
			double dirFraction = 0.001;

			// Fix S0 to starting value
			paramsNew[0] = params[0];

			// Fix the white matter volume fraction but
			// perturb the intra versus extracellular proportion.
			// f1.
                        // Slight difference here to original matlab code,
                        // which scales perturbation by paramsInitial[1].
			paramsNew[1] = params[1] + nonDirFraction*rand.nextGaussian();
			// f2
			paramsNew[2] = params[1] + params[2] - paramsNew[1];
			paramsNew[3] = params[3];

			// Keep orientation fixed
			paramsNew[4] = params[4];
			paramsNew[5] = params[5];

			// R
			paramsNew[6] = params[6] + paramsInitial[6]*nonDirFraction*rand.nextGaussian();

                        return paramsNew;
                    }
		};

	}


	/**
	 * Estimates a starting set of parameters from multirun
	 * gradient descent fit of the same model.
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

		// Use the gradient descent fitter.
		double[] modelParams = new double[numModelParams];
                try {
		    double[] gdparams = MultiRunMinimizer.getBestSolution(gdfitter.fit(data));
		    for(int i=0; i<modelParams.length; i++) {
			modelParams[i] = gdparams[i+1];
		    }

                } catch (MinimizerException e) {
		    modelParams = gdfitter.getStartPoint(data);
		    logger.warning(e.toString());
                }

		return modelParams;

	}

	
	public static void main(String[] args) {

		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();

		OutputManager om = new OutputManager();

		MMWMD_FixedMCMC_Fitter inv = new MMWMD_FixedMCMC_Fitter(CL_Initializer.imPars);


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

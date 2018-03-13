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
 * <dd>Fits the basic MMWMD model using MCMC and with a specified
 * noise model.
 * 
 * <dt>Description:
 * 
 * <dd> Compartments are:
 * Cylinder
 * Zeppelin
 *
 * Parameters are:
 * S0 - constrained positive by prior
 * f1 - cylinder volume fraction; constrained to [0 1].
 * f2 - zeppelin volume fraction; equals 1-f1.
 * d - intrinsic diffusivity within cylinder; specified on command line
 * theta - colatitude of cylinder orientation.
 * phi - longitude of cylinder orientation.
 * R - cylinder radius; constrained to [1E-7 2E-5] microns by prior.
 * d - zeppelin parallel diffusivity; specified on command line.
 * theta - colatitude of zeppelin orientation; same as cylinder.
 * phi - longitude of zeppelin orientation; same as cylinder.
 * 
 * The MCMC varies only the proportions of f1, f2 and R.  All
 * other parameters are fixed at their starting
 * values, which come from a maximum likelihood gradient descent with
 * multiple starting points.
 *
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class MMWMD_BasicMCMC_Fitter extends MCMC_Fitter {

	/**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.fitters.MMWMD_InVivoMCMC_Fitter");
	
    protected MMWMD_BasicMultiRunLM_DirectFitter gdfitter;
    protected int GDRUNS = 1;

    private final int numModelParams = 11;
    private final int numOptParams = 5;

    // Fixed diffusivity.
    protected double WMDIFF=1.7E-9;


	/**
	 * Default constructor required for derived classes.
	 */
	public MMWMD_BasicMCMC_Fitter() {
	}


	/**
	 * Basic constructor extracts model particulars and calls
         * full constructor.
         *
         * @param scheme The imaging protocol.
	 */
	public MMWMD_BasicMCMC_Fitter(DW_Scheme scheme) {
            this(scheme, FitModel.getCompartmentList(CL_Initializer.fitModel), FitModel.getCompartmentModelParams(CL_Initializer.fitModel));
	}


	/**
	 * Constructor sets up various required objects for MCMC.
	 *
	 * @param scheme The imaging protocol.
	 * @param compNames the compartment names.
	 * @param initialParams the initial parameters.

	 */
        public MMWMD_BasicMCMC_Fitter(DW_Scheme scheme, String[] compNames, double[] initialParams) {
            super(scheme, compNames, initialParams);

            if(CL_Initializer.mmwmddiff>0)
		WMDIFF = CL_Initializer.mmwmddiff;

	    // Initialize the gradient descent fitter for starting point
            // estimation.
            gdfitter = new MMWMD_BasicMultiRunLM_DirectFitter(scheme, GDRUNS, CL_Initializer.seed);
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
			    optParams[2] = modelParams[4]; //theta
			    optParams[3] = modelParams[5]; //phi
			    optParams[4] = modelParams[6];// R

			    return optParams;
                        }

                        public double[] optToModel(double[] optParams) {
			    double[] modelParams = new double[numModelParams];

			    // S0
			    modelParams[0] = optParams[0];

			    // Volume fractions.
			    modelParams[1] = optParams[1];
			    modelParams[2] = 1-optParams[1];

			    // Cylinder compartment
			    // Intrinsic ICS diff
			    modelParams[3] = WMDIFF;
			    // Cylinder orientation
			    modelParams[4] = optParams[2]; //theta
			    modelParams[5] = optParams[3]; //phi
			    // Cylinder radius
			    modelParams[6] = optParams[4];// R 

			    // Zeppelin compartment
			    // Parallel diff
			    modelParams[7] = WMDIFF;
			    // Parallel orientation
			    modelParams[8] = optParams[2]; //theta
			    modelParams[9] = optParams[3]; //phi
			    // Perpendicular diff
			    modelParams[10] = WMDIFF*modelParams[2];

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

		    // Bounds on R
		    private double minR = 1e-7;
		    private double maxR = 20e-6;

                    public double prior(double[] optParams) {

                        // s0 must be positive.
                        double p1 = optParams[0]<0?0:1;
                        
                        // f1 must be in [0 1].
                        double p2 = (optParams[1]<0 || optParams[1]>1)?0:1;
                        
                        // the prior on the direction is sin(theta)
                        // to ensure even sampling over the whole
                        // surface.
                        double p3 = Math.abs(Math.sin(optParams[2]));

			// R constrained to realistic range
                        double p4 = (optParams[4]<minR || optParams[4]>maxR)?0:1;
                                                
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

                        // Fix S0 to starting value
                        paramsNew[0] = params[0];

                        // Fix the white matter volume fraction but
                        // perturb the intra versus extracellular proportion.
                        // f1
                        paramsNew[1] = params[1] + nonDirFraction*rand.nextGaussian();

                        // Keep orientation fixed
                        paramsNew[2] = params[2];
                        paramsNew[3] = params[3];

			// R
			paramsNew[4] = params[4] + paramsInitial[4]*nonDirFraction*rand.nextGaussian();

                        return paramsNew;
                    }
		};

	}


	/**
	 * Estimates a starting set of parameters from fitting the
	 * same model by multirun gradient descent.
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

		MMWMD_BasicMCMC_Fitter inv = new MMWMD_BasicMCMC_Fitter(CL_Initializer.imPars);


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

package optimizers;

import models.*;
import numerics.*;
import misc.*;
import imaging.*;
import fitters.*;
import java.util.Random;
import java.util.logging.Logger;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Runs another minimizer multiple times from different starting
 * points and concatenates the results.
 * 
 * <dt>Description:
 * 
 * <dd> An argument to the constructor specifies the perturbations of
 * the starting point for each run.  The perturbations are zero-mean
 * Gaussian distributed and the constructor argument gives a separate
 * standard deviation for each optimization parameter.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public abstract class MultiRunMinimizer implements Minimizer {
	
	private static Logger logger = Logger.getLogger("optimizers.MultiRunMinimizer");

    // The fixed starting point that is perturbed repeatedly.
    private double[] initOptParams;

    // Specifies the start point perturbations.
    protected Perturbation perturbation;

    // Random number generator for the perturbations
    private Random r;

    // The number of repeated runs.
    private int numRepeats;

    // The minimizer to run repeatedly.
    protected Minimizer minimizer;

    // Data structure to store the results
    private double[][] concatRes;

    protected DW_Scheme scheme;
    protected HCPScheme masterScheme;
    ParametricModel pm;
    Codec cod;
    
    //fail rate
    private static double FAILRATE =0.3;


    /**
     * Basic constructor required for inheritance.
     */
    public MultiRunMinimizer() {
    }


    /**
     * Constructor needs all the following:
     *
     * @param scheme The acquisition protocol
     *
     * @param pm The model to fit
     *
     * @param cod The Codec specifying the transformation from model
     * to optimized parameters.
     *
     * @param pertStds array of standard deviations of the
     * perturbations on each starting parameter.
     *
     * @param repeats The number of times to repeat the optimization.
     *
     * @param seed Seed for random number generator.
     */
    public MultiRunMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod, Perturbation p, int repeats, int seed) throws MinimizerException {
	init(cod, p, repeats, seed);
	makeMinimizer(scheme, pm, cod);
    }


    /**
     * Initializes various class variables
     *
     * @param cod The Codec specifying the transformation from model
     * to optimized parameters.
     *
     * @param pertStds array of standard deviations of the
     * perturbations on each starting parameter.
     *
     * @param repeats The number of times to repeat the optimization.
     *
     * @param seed Seed for random number generator.
     */
    protected void init(Codec cod, Perturbation p, int repeats, int seed) {
        initOptParams = new double[cod.getNumOptParams()];
        perturbation = p;
        r = new Random(seed);
        numRepeats = repeats;
    }


    /**
     * Creates the Minimizer object to run repeatedly.
     *
     * @param scheme The acquisition protocol
     *
     * @param pm The model to fit
     *
     * @param cod The Codec specifying the transformation from model
     * to optimized parameters.
     */
    protected abstract void makeMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod) throws MinimizerException;


    /**
     * Sets the perturbation object for getting starting points.
     */
    public void setPerturbation(Perturbation p) {
        perturbation = p;
    }


    public void setInitParams(double[] aInit) throws MinimizerException {
        for(int i=0; i<initOptParams.length; i++) {
            initOptParams[i] = aInit[i];
        }
    }


    public void setMeasurements(double[] newMeas) throws MinimizerException {
        minimizer.setMeasurements(newMeas);
    }


    /**
     * Run the minimizer repeatedly from perturbed starting points and
     * concatenate the result.
     */
    public void minimise() throws MinimizerException {
    	int i=0;
    	int fail=0;
    	while (i<numRepeats)
    	    {
            double[] pertStart = perturbation.perturb(initOptParams, initOptParams, r);
            minimizer.setInitParams(pertStart);
            
            try
            {
            	minimizer.minimise();
            	double[][] res = minimizer.getSolutions();
            	if(i==0)
                    concatRes = new double[res.length*numRepeats][res[0].length];
                for(int j=0; j<res.length; j++)
                    for(int k=0; k<res[0].length; k++)
                        concatRes[i*res.length+j][k] = res[j][k];
                i++;
            }
            catch (Exception e)
            {
            	//System.err.println("Multirun caught exception "+e.getMessage());
            	//System.err.println("" + fail + "of the fittings failed");
            	fail++;
            }
            
            if (fail > numRepeats*FAILRATE)
            {
            	logger.warning("Aborting fitting, more than " + (int)(FAILRATE * 100) + " percent of the fittings failed");
            	
            	// fill rest of concatRes with errorcode 4 = "multirun failrate exceeded"
            	// only works when res.length==1 atm 
            	// which is the case for all single run fitters atm
            	//  --- Torben ---            	
            	int numRes=minimizer.getNumParameters();            	

            	// allocate if all failed 
            	if(i==0)
            	{
                    concatRes = new double[numRepeats][numRes];
            	}
            	
            	//fill rest of result array with error code and zeros 
            	for (int errorIndex=i; errorIndex < numRepeats; errorIndex++ )
            	{            		
            		concatRes[errorIndex][0] = 4; //error code for failrate exceeded
            		for (int k=1; k<numRes;k++)
            		{
            			concatRes[errorIndex][k] = 0;
            		}
            	}
            	 
            	//end loop
            	return;
            	
            	
            }
        }
    }


    public double[][] getSolutions() {
        return concatRes;
    }


    public int getNumSolutions() {
        return numRepeats;
    }
    
    public int getNumParameters() {
        return minimizer.getNumParameters();
    }


    /**
     * Returns the solution with minimum MSE from an array of multirun
     * solutions
     */
	public static double[] getBestSolution(double[][] solutions) {
		int solution_index = 0;
		double minMSE = Double.MAX_VALUE;
		double maxMSE = Double.MIN_VALUE;
		int fobj_index = solutions[0].length - 1;

		for (int m = 0; m < solutions.length; m++) {
			if (minMSE > solutions[m][fobj_index]) {
				minMSE = solutions[m][fobj_index];
				solution_index = m;
			}

			if (maxMSE < solutions[m][fobj_index]) {
				maxMSE = solutions[m][fobj_index];
			}
		}

		return solutions[solution_index];

	}

}









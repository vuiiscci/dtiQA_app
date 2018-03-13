package optimizers;

import models.*;
import optimizers.*;
import numerics.*;
import misc.*;
import imaging.*;
import fitters.*;
import java.util.logging.Logger;

public abstract class LM_Minimizer extends MarquardtMinimiser implements Minimizer {

    private static Logger logger = Logger.getLogger("optimizers.LM_Minimizer");
    
    
    protected DW_Scheme scheme;
    protected ParametricModel model;
    protected Codec codec;
    protected HCPScheme masterScheme;

    // Standard deviations of the measurements.
    protected double[] sig;

    // Dependent values at sampled points: the measurements.
    protected double[] measurements;

    // Stores the exitcode after minimization
    protected int exitcode;

    /**
     * Returns the single minimum found by one run of the gradient descent.
     *
     * @return 1x(N+1) matrix in which the single row contains the N
     * model parameters followed by the objective function value.
     */
    public double[][] getSolutions() {
        double[] optParams = new double[a.length-1];
        for(int i=0; i<a.length-1; i++)
            optParams[i] = a[i+1];
        double f = fObj(a, new double[a.length], new double[a.length][a.length]);
        double[][] solutions = new double[1][codec.getNumModelParams() + 2];
        double[] modelParams = codec.optToModel(optParams);

	solutions[0][0] = exitcode;
        for(int i=0; i<modelParams.length; i++)
            solutions[0][i+1] = modelParams[i];
        solutions[0][solutions[0].length-1] = f;
        return solutions;
    }


    public int getNumSolutions() {
        return 1;
    }
    
    public int getNumParameters() {
        return model.numParams()+2; //model parameters + error code and objective function
    }


    public void minimise() throws MarquardtMinimiserException {
	exitcode = 0;
	try {
	    super.minimise();
	}
	catch(MarquardtMinimiserNonConvergenceException e) {
	    logger.warning(e.toString());
	    exitcode = 3;
	}
	catch(Exception e) {
	    logger.warning(e.toString());
	    exitcode = 2;
	}
    }


    /**
     * Initializes the fitting procedure with a new set of
     * measurements (dependent variables).
     * 
     * @param newMeas The new set of measurements.
     */
    public void setMeasurements(double[] newMeas) throws MarquardtMinimiserException {
        if (newMeas.length != scheme.numMeasurements()) {
            throw new MarquardtMinimiserException("New data contains the wrong number of values: " + newMeas.length + " instead of " + measurements.length + ".");
        }
        for (int i = 0; i < newMeas.length; i++) {
            measurements[i] = newMeas[i];
        }

    }


    /**
     * Computes derivatives of the signals wrt to the optimization
     * parameters directly.  For debugging.
     *
     * @param optParams Is the vector of optimization parameters.
     * Similar to that used within the LM algorithm but starts at
     * index zero rather than index one.
     */
    protected RealMatrix numericalOptJac(double[] optParams) {
        final RealMatrix jac = new RealMatrix(scheme.numMeasurements(), optParams.length);

        final double[] a = new double[optParams.length];
        for (int k = 0; k < optParams.length; k++) {
            a[k] = optParams[k];
        }

        for(int i=0; i < scheme.numMeasurements(); i++) {
            for(int j = 0; j < optParams.length; j++) {

                final int index1 = i;
                final int index2 = j;

                NumDeriv df = new NumDeriv() {
                        protected float func(float arg) {

                            a[index2] = arg;
                            float t = (float) model.getSignal(codec.optToModel(a), scheme, index1);
                            return t;
                        }
                    };

                float[] err = new float[1];
                
                jac.entries[index1][index2] = df.dfridr((float) optParams[index2], (optParams[index2] != 0.0) ? 0.1f * (float) optParams[index2] : 0.1f, err);
                

                //Replace the original value in a[index2] ready to compute
                //the next derivative.
                a[index2] = optParams[index2];
            }
        }

        return jac;

    }
    

}



package optimizers;

import models.*;
import numerics.*;
import misc.*;
import imaging.*;
import fitters.*;
import java.util.Random;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd> Multi-run minimizer specifically for LM minimization.
 * 
 * <dt>Description:
 * 
 * <dd> Adds a few methods to set parameters of the LM minimization.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public abstract class MultiRunLM_Minimizer extends MultiRunMinimizer {


    public MultiRunLM_Minimizer() {
    }

    public MultiRunLM_Minimizer(DW_Scheme scheme, ParametricModel pm, Codec cod, Perturbation p, int repeats, int seed)  throws MinimizerException {
	super(scheme, pm, cod, p, repeats, seed);
    }

    /**
     * Sets the convergence threshold in the LM algorithm.
     */
    public void setCONVERGETHRESH(double convthresh) {
	((LM_Minimizer)minimizer).setCONVERGETHRESH(convthresh);
    }
        
    /**
     * Sets the maximum number of iterations in each run.
     */
    public void setMAXITER(int maxiter) {
	((LM_Minimizer)minimizer).setMAXITER(maxiter);
    }

    /**
     * Returns the convergence threshold in the LM algorithm.
     */
    public double getCONVERGETHRESH() {
	return ((LM_Minimizer)minimizer).getCONVERGETHRESH();
    }

}









package optimizers;

import models.*;
import numerics.*;
import misc.*;
import imaging.*;
import fitters.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Runs the LM minimizer with offset Gaussian noise model multiple
 * times from different starting points and concatenates the results.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class MultiRunLM_OffGaussMinimizer extends MultiRunLM_Minimizer {

    // Noise std.
    private double sigma;

    public MultiRunLM_OffGaussMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod, Perturbation p, int repeats, int seed, double sigma) throws MinimizerException {
	this.sigma = sigma;
	init(cod, p, repeats, seed);
        makeMinimizer(scheme, pm, cod);
    }


    protected void makeMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod) throws MarquardtMinimiserException {

        minimizer = new LM_OffGaussMinimizer(scheme, pm, cod, sigma);

    }
    

    public void adjustScheme(double[] gradAdj) {
        
    }
}




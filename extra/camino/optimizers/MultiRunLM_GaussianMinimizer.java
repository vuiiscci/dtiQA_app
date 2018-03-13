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
 * <dd>Runs the LM_Gaussian minimizer multiple times from different starting
 * points and concatenates the results.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class MultiRunLM_GaussianMinimizer extends MultiRunLM_Minimizer {


    public MultiRunLM_GaussianMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod, Perturbation p, int repeats, int seed) throws MinimizerException {
        super(scheme, pm, cod, p, repeats, seed);
        this.scheme = scheme;
        this.pm = pm;
        this.cod = cod;

        if (scheme instanceof HCPScheme) {
            this.masterScheme = ((HCPScheme)scheme).copyScheme();    
        }
    }


    protected void makeMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod) throws MarquardtMinimiserException {

        minimizer = new LM_GaussianMinimizer(scheme, pm, cod);

    }

   public void adjustScheme(double[] gradAdj) {
   		if (scheme instanceof HCPScheme) {
            ((HCPScheme)scheme).resetScheme(masterScheme);
            ((HCPScheme)scheme).modifyScheme(gradAdj, masterScheme);
            try {
            	makeMinimizer(scheme, pm, cod);
            }
            catch (MarquardtMinimiserException e) {
            	System.out.println("[MultiRunLM_GaussianMinimizer caught" + e);
            }
        }
        else {
            throw new LoggedException("[MultiRunLM_GaussianMinimizer adjustScheme given a non-HCPScheme scheme");
        }
    	// System.out.println("asdfkh");
   }

}




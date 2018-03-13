package optimizers;

import models.*;
import optimizers.*;
import numerics.*;
import misc.*;
import imaging.*;
import fitters.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Levenberg-Marquardt algorithm for minimizing chi-squared
 * objective function, ie Gaussian noise model, with equal variance on
 * each measurement.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class LM_GaussianMinimizer extends LM_Minimizer {



    /**
     * Default constructor.
     */
    public LM_GaussianMinimizer() {
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
     */
    public LM_GaussianMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod) throws MarquardtMinimiserException {
	
	    this.scheme = scheme;
        this.codec = cod;
        this.model = pm;
        if (scheme instanceof HCPScheme) {
            this.masterScheme = ((HCPScheme)scheme).copyScheme();    
        }
        // Simple, constant variance, Gaussian noise model
        sig = new double[scheme.numMeasurements()];
        for(int i=0; i<sig.length; i++) {
            sig[i] = 1.0;
        }

        measurements = new double[scheme.numMeasurements()];

        init(codec.getNumOptParams());
    }


    /**
     * Implements the chi-squared function and its derivatives. The
     * second derivative matrix is replaced by values proportional to
     * the products of first derivatives: <code>d2fda2[i][j] =
     * dfda[i]*dfda[j]*sig2i</code>, which avoids computation and can
     * improve stability.
     * 
     * @param params
     *            The current model parameter settings.
     * 
     * @param dfda
     *            The first derivatives of the model with respect to each
     *            parameter.
     * 
     * @param d2fda2
     *            The second derivatives of the model with respect to each pair
     *            of parameters.
     */
    protected double fObj(double[] a, double[] dfda, double[][] d2fda2) {

        //Initialise values and arrays.
        double chisq = 0.0;
        for (int i = 1; i <= ma; i++) {
            for (int j = 1; j <= ma; j++) {
                d2fda2[i][j] = 0.0;
            }
            dfda[i] = 0.0;
        }

        double[] optParams = new double[a.length-1];
        for(int i=0; i<a.length-1; i++)
            optParams[i] = a[i+1];

        double[] modParams = codec.optToModel(optParams);

        RealMatrix modSignals = model.getSignals(modParams, scheme);
        RealMatrix modJac = model.getJacobian(modParams, scheme);
        RealMatrix codJac = codec.getJacobian(optParams);

        // Multiply to get required Jacobian wrt optimized parameters.
        RealMatrix optJac = modJac.product(codJac);
	// Useful for debugging purposes.
        //optJac = numericalOptJac(optParams);
        boolean print = false;
        if(print) {
            System.err.println("modParams");
            for(int i=0; i<modParams.length; i++)
                System.err.print(modParams[i] + " ");
            System.err.println();
            System.err.println("measurements");
            for(int i=0; i<measurements.length; i++)
                System.err.print(measurements[i] + " ");
            System.err.println();
            System.err.println("modSignals");
            System.err.println(modSignals);
            System.err.println("modJac");
            System.err.println(modJac);
            System.err.println("codJac");
            System.err.println(codJac);

            System.err.println("optJac");
            System.err.println(optJac);

	    RealMatrix nOptJac = numericalOptJac(optParams);
	    System.err.println("numerical optJac");
	    System.err.println(nOptJac);
        }

        for (int i = 1; i <= scheme.numMeasurements(); i++) {
            double ymod = modSignals.entries[i-1][0];

            double[] dyda = optJac.entries[i-1];
            double sig2i = 1.0 / (sig[i-1] * sig[i-1]);
            double dy = measurements[i-1] - ymod;
            for (int l = 1; l <= ma; l++) {
                double wt = 2.0 * dyda[l-1] * sig2i;
                for (int m = 1; m <= l; m++) {
                    d2fda2[m][l] += wt * dyda[m-1];
                }
                dfda[l] -= dy * wt;
            }
            chisq += dy * dy * sig2i;
        }

        //Fill in other off-diagonal of second derivatives assuming
        //symmetry.
        for (int i = 2; i <= ma; i++) {
            for (int j = 1; j < i; j++) {
                d2fda2[i][j] = d2fda2[j][i];
            }
        }
        if(print)
            System.err.println("chisq: " + chisq);

        return chisq;
    }
    
    /**
     * Adjusts the gradient values of the scheme per voxel for HCP data.
     * 
     * @param gradAdj speficies the value of the gradient adjustment.
     */
    public void adjustScheme(double[] gradAdj){
    
        if (scheme instanceof HCPScheme) {
            ((HCPScheme)scheme).resetScheme(masterScheme);
            ((HCPScheme)scheme).modifyScheme(gradAdj, masterScheme);
        }
        else {
            throw new LoggedException("[LM_GaussianMinimizer adjustScheme given a non-HCPScheme scheme");
        }
    }

}


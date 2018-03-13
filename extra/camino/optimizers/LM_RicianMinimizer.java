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
 * <dd>Levenberg-Marquardt algorithm for minimizing the log-likelihood
 * objective function using a Rician noise model assuming equal
 * underlying variances on the real and imaginary Gaussian models on
 * each measurement.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class LM_RicianMinimizer extends LM_OffGaussMinimizer {



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
     * @param sigma Standard deviations of the Gaussian noise
     * distribution in each channel.
     */
    public LM_RicianMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod, double sigma) throws MarquardtMinimiserException {
	
	super(scheme, pm, cod, sigma);
    }


    /**
     * Implements the Rician noise model that accounts for the
     * bias introduced by Rician noise. The second derivative matrix
     * ignores the terms in the second derivatives of the model.
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
        double loglik = 0.0;
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
        }

        for (int i = 1; i <= scheme.numMeasurements(); i++) {
            double sig2 = sig[i-1] * sig[i-1];
            double sig2i = 1.0/sig2;
            double ymod = modSignals.entries[i-1][0];
            double z = ymod*measurements[i-1]*sig2i;
            
            // Above z threshold, the exact Bessel function
            // implementation fails so we approximate the
            // values we need.
            double logBessI0;
            double bessI1_bessI0;
            double bessI2_bessI0;
            if(z<700) {
                double bessI0 = BesselFunctions.besselI0(z);
                logBessI0 = Math.log(bessI0);
                double bessI1 = BesselFunctions.besselI1(z);
                bessI1_bessI0 = bessI1 / bessI0;
                double bessI2 = BesselFunctions.besselIn(2, z);
                bessI2_bessI0 = bessI2 / bessI0;
            }
            else {
                // These are log linear approximations.  See
                // ActiveImaging/models/besseli1d0.m for details.
                logBessI0 = z - Math.log(2.0*Math.PI*z)/2.0;
                bessI1_bessI0 = 1 - Math.exp((z-700)*(-0.00143010502604)-7.24386992534278);
                bessI2_bessI0 = 1 - Math.exp((z-700)*(-0.001428569477897)-5.858647951071322);
            }
            double[] dyda = optJac.entries[i-1];
            for (int l = 1; l <= ma; l++) {
                for (int m = 1; m <= l; m++) {
                    d2fda2[m][l] += dyda[l-1]*dyda[m-1]*sig2i*(1 - measurements[i-1]*measurements[i-1]*sig2i*(0.5 + 0.5*bessI2_bessI0 - bessI1_bessI0*bessI1_bessI0));
                }
                dfda[l] += dyda[l-1]*sig2i*(ymod - measurements[i-1]*bessI1_bessI0);
            }

            // Constant terms excluded for a slight speed gain.
            loglik += 0.5*sig2i*ymod*ymod -logBessI0;
            // Here is an additional step that includes the
            // constant terms if desired.
            //loglik += 0.5*sig2i*measurements[i-1]*measurements[i-1] - Math.log(measurements[i-1]) + Math.log(sig2);
        }

        //Fill in other off-diagonal of second derivatives assuming
        //symmetry.
        for (int i = 2; i <= ma; i++) {
            for (int j = 1; j < i; j++) {
                d2fda2[i][j] = d2fda2[j][i];
            }
        }
        if(print)
            System.err.println("chisq: " + loglik);

        return loglik;
    }


}


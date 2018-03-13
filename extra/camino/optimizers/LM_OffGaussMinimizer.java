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
 * objective function with the offset Gaussian noise model, which adds
 * an expected offset to each measurement from Rician noise bias; see
 * Jones and Basser MRM 2004; Alexander "Modelling fitting and
 * sampling..." 2009.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public class LM_OffGaussMinimizer extends LM_GaussianMinimizer {



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
    public LM_OffGaussMinimizer(DW_Scheme scheme, ParametricModel pm, Codec cod, double sigma) throws MarquardtMinimiserException {
	
	this.scheme = scheme;
        this.codec = cod;
        this.model = pm;
	
        // Underlying variance the same for each measurement
        sig = new double[scheme.numMeasurements()];
        for(int i=0; i<sig.length; i++) {
            sig[i] = sigma;
        }

        measurements = new double[scheme.numMeasurements()];

        init(codec.getNumOptParams());
    }


    /**
     * Implements the chi-squared function and its derivatives but
     * using the offset Gaussian noise model that accounts for the
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

        for (int i = 1; i <= scheme.numMeasurements(); i++) {
            double sig2 = sig[i-1] * sig[i-1];
            double sig2i = 1.0/sig2;
            double ymod = modSignals.entries[i-1][0];
            ymod = Math.sqrt(ymod*ymod + sig2);

            double[] dyda = optJac.entries[i-1];
            double dy = measurements[i-1] - ymod;
            for (int l = 1; l <= ma; l++) {
                double doffdy = modSignals.entries[i-1][0] / ymod;
                double wt = 2.0 * dyda[l-1] * sig2i * doffdy;
                for (int m = 1; m <= l; m++) {
                    d2fda2[m][l] += wt * dyda[m-1] * doffdy;
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

        return chisq;
    }


}


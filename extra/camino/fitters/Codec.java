package fitters;

import numerics.*;


/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Base class for codec that maps model parameters to optimization
 * parameters during model fitting in order to enforce constraints on
 * parameter values.
 * 
 * <dt>Description:
 * 
 * <dd>Requires functions that implement both the mapping and the
 * inverse mapping and also the derivative of the forward mapping (ie
 * from model parameters to optimized parameters).
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: DataSource.java 58 2006-10-26 11:38:18Z ucacmgh $
 *  
 */
public abstract class Codec {


    /**
     * @return The number of optimized parameters.
     */
    public abstract int getNumOptParams();


    /**
     * @return The number of model parameters.
     */
    public abstract int getNumModelParams();


    /**
     * Returns the first index of the theta-phi pair of parameters
     * encoding a single fiber orientation.  It is the index in the
     * model parameters list.  Theta and phi are assumed consecutive
     * and in that order.  If the model contains no orientation
     * parameter, return -1, which is the default.
     *
     * @return The index of the direction model parameter.
     */
    public int getDirectionIndex() {
        return -1;
    }


    /**
     * Maps the model parameters to the optimized parameters.
     *
     * @param array of model parameters
     *
     * @return array of optimized parameters.
     */
    public abstract double[] modelToOpt(double[] modelParams);


    /**
     * Maps the optimized parameters to the model parameters.
     *
     * @param array of optimized parameters
     *
     * @return array of model parameters.
     */
    public abstract double[] optToModel(double[] optParams);


    // Compute it numerically by default and make non abstract.
    /**
     * Computes the jacobian of the encoding, specifically d(model
     * parameters)/d(optimization parameters).
     *
     * @return Jacobian matrix with size NxM where N is the number of
     * model parameters and M is the number of optimization
     * parameters.
     */
    public RealMatrix getJacobian(double[] optParams) {

        final RealMatrix jac = new RealMatrix(getNumModelParams(), getNumOptParams());
        final double[] modParams = optToModel(optParams);
        final double[] a = new double[optParams.length];
        for (int k = 0; k < optParams.length; k++) {
            a[k] = optParams[k];
        }

        for(int i=0; i < modParams.length; i++) {
            for(int j = 0; j < optParams.length; j++) {

                final int index1 = i;
                final int index2 = j;

                NumDeriv df = new NumDeriv() {
                        protected float func(float arg) {

                            a[index2] = arg;
                            double[] pert = optToModel(a);
                            float t = (float) pert[index1];
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


    // Could add this if necessary.
    // Compute it numerically by default and make non abstract.
    //public abstract RealMatrix[] getHessian(double[] optParams);

}




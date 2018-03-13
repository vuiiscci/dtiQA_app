package models;

import imaging.DW_Scheme;
import numerics.NumDeriv;
import numerics.RealMatrix;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Abstract parametric model class for data synthesis or model fitting.
 * 
 * <dt>Description:
 * 
 * <dd>Subclasses must implement the getSignal and getSignals method.
 * This class implements a default getJacobian method, which estimates
 * the jacobian of the signal with respect to the model parameters
 * numerically.  For greater efficiency during model fitting,
 * subclasses should override the numerical estimation with analytic
 * expressions.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 * 
 */
public abstract class ParametricModel {   
    
    
    /** how many parameters in this model? */
    private final int numParams;
    
    /** constructor. takes number of parameters */
    public ParametricModel(int numParams){
        this.numParams= numParams;
    }
    
    
    /**
     * Returns a list of signals for the specified scheme from the
     * model using the specified model parameters.
     *
     * @param modParams The model parameters
     *
     * @param scheme The acquisition scheme
     *
     * @return A 1xK RealMatrix containing the K signals corresponding
     * to the K elements of scheme.
     */
    public abstract RealMatrix getSignals(double[] modParams, DW_Scheme scheme);


    /**
     * Returns the signals for the i-th measurement of specified
     * scheme from the model using the specified model parameters.
     *
     * @param modParams The model parameters
     *
     * @param scheme The acquisition scheme
     *
     * @param i Index of the scheme entry to compute the signal for
     *
     * @return The signal from the model.
     */
    public abstract double getSignal(double[] modParams, DW_Scheme scheme, int i);


    /**
     * Computes the Jacobian of the model numerically.  Specifically,
     * d(measurements)/d(model parameters) for the specified scheme.
     *
     * @param modParams The model parameters
     *
     * @param scheme The acquisition scheme
     *
     * @return An KxN matrix where N is the number of model parameters
     * and K is the number of measurements in the scheme.
     */
    public RealMatrix getJacobian_old(final double[] modParams, final DW_Scheme scheme) {

        final RealMatrix jac = new RealMatrix(scheme.numMeasurements(), modParams.length);

        final double[] a = new double[modParams.length];
        for (int k = 0; k < modParams.length; k++) {
            a[k] = modParams[k];
        }

        for(int i=0; i < scheme.numMeasurements(); i++) {
            for(int j = 0; j < modParams.length; j++) {

                final int index1 = i;
                final int index2 = j;

                NumDeriv df = new NumDeriv() {
                        protected float func(float arg) {

                            a[index2] = arg;
                            float t = (float) getSignal(a, scheme, index1);
                            return t;
                        }
                    };

                float[] err = new float[1];
                
                jac.entries[index1][index2] = df.dfridr((float) modParams[index2], (modParams[index2] != 0.0) ? 0.1f * (float) modParams[index2] : 0.1f, err);
                

                //Replace the original value in a[index2] ready to compute
                //the next derivative.
                a[index2] = modParams[index2];
            }
        }

        return jac;

    }
    
    public RealMatrix getJacobian(final double[] modParams, final DW_Scheme scheme) {

    	final RealMatrix jac = new RealMatrix(scheme.numMeasurements(), modParams.length);

    	for(int j = 0; j < modParams.length; j++) {
    		RealMatrix cJ = getJacobian(modParams, scheme, j);
    		
    		for(int r=0; r < scheme.numMeasurements();r++){
            	jac.setEntry(r, j, cJ.entry(r, j));
            } 

    	}
    	
    	return jac;
    }
    
    /**
     * Computes the Jacobian of the model numerically.  Specifically,
     * d(measurements)/d(model parameters) for the specified scheme.
     *
     * @param modParams The model parameters
     *
     * @param scheme The acquisition scheme
     *
     * @return An KxN matrix where N is the number of model parameters
     * and K is the number of measurements in the scheme.
     */
    public RealMatrix getJacobian(final double[] modParams, final DW_Scheme scheme, int index) {

        final RealMatrix jac = new RealMatrix(scheme.numMeasurements(), modParams.length);

        final double[] a = new double[modParams.length];
        for (int k = 0; k < modParams.length; k++) {
            a[k] = modParams[k];
        }

        for(int i=0; i < scheme.numMeasurements(); i++) {
            

                final int index1 = i;
                final int index2 = index;

                NumDeriv df = new NumDeriv() {
                        protected float func(float arg) {

                            a[index2] = arg;
                            float t = (float) getSignal(a, scheme, index1);
                            return t;
                        }
                    };

                float[] err = new float[1];
                jac.entries[index1][index2] = df.dfridr((float) modParams[index2], (modParams[index2] != 0.0) ? 0.1f * (float) modParams[index2] : 0.1f, err);

                //Replace the original value in a[index2] ready to compute
                //the next derivative.
                a[index2] = modParams[index2];
            
        }

        return jac;

    }
    
    

    /** @return number of parameters in model */
    public final int numParams(){
        return numParams;
    }
}



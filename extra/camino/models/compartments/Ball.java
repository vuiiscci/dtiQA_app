package models.compartments;

import models.ParametricModel;
import numerics.RealMatrix;
import imaging.DW_Scheme;
import models.compartments.CompartmentType;

/**
 * implents the compartment interface for an isotropic 
 * Ball compartment
 * 
 * 
 * @author laura(panagio@cs.ucl.ac.uk)
 *
 */
public class Ball extends ParametricModel {
    
    
    /** constructor. needs array of params. 
     * in this case, an array of one parameter, the diffusivity.
     * 
     * @param params params array
     */
    public Ball(){

        super(CompartmentType.BALL.numParams);        
    }
    
    
    /**
     * generates signals from this compartment for each line
     * in the scheme given.
     * 
     * @param scheme scan specifics
     */
    public RealMatrix getSignals(double[] params, DW_Scheme scheme) {
        
        RealMatrix signals= new RealMatrix(scheme.numMeasurements(), 1);
        
        for(int i=0; i<scheme.numMeasurements(); i++){
            
            
           // double b= scheme.getB_Value(i);
            
            
            signals.setEntry(i, 0, getSignal( params, scheme,i));
            
        }
        
        
        return signals;
    }

    /**
     * generates signals from this compartment for each line
     * in the scheme given.
     * 
     * @param scheme scan specifics
     */
    public double getSignal(double[] params, DW_Scheme scheme, int i) {
        
            double b= scheme.getB_Value(i);
            // Commented out the next line added by Matt.  No idea
            // why this was here, but it introduces numerical problems.
            //b= b * Math.sqrt(scheme.getG_Dir(i)[0]*scheme.getG_Dir(i)[0]+scheme.getG_Dir(i)[1]*scheme.getG_Dir(i)[1]+scheme.getG_Dir(i)[2]*scheme.getG_Dir(i)[2]);
            
            
            return Math.exp(-b*params[0]);
            
    }

    /**
     * Compute the jacobian of the signal.
     * 
     * @param params the model parameter settings
     * @param scheme scan specifics
     */
    public RealMatrix getJacobian(double[] modParams, DW_Scheme scheme) {
        
        RealMatrix J = new RealMatrix(scheme.numMeasurements(), 1);
        RealMatrix S = getSignals(modParams, scheme);
        for(int i=0; i<scheme.numMeasurements(); i++){
            double b= scheme.getB_Value(i);
            J.setEntry(i, 0, -b*S.entries[i][0]);
        }
        
        return J;
    }

}

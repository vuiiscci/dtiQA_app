package models.compartments;

import numerics.RealMatrix;
import misc.LoggedException;
import models.ParametricModel;
import imaging.*;

/**
 * implements the compartment interface for a restricted compartment (according
 * to Behren's ball and stick model) the Stick compartment
 * 
 * 
 * @author laura(panagio@cs.ucl.ac.uk)
 *
 */
public class Stick extends ParametricModel {

    
    /** constructor. needs array of params. 
     * in this case, an array of 3 parameters: the diffusivity and the angles theta and phi
     * that determine the fibre orientation.
     * 
     * @param params params array
     */
    public Stick(){
        
        super(CompartmentType.STICK.numParams);        
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
            signals.setEntry(i, 0, getSignal(params, scheme, i));
        }
        
        return signals;
    }

    public double getSignal(double[] params, DW_Scheme  scheme, int i){

        double stickDiff = params[0];
        double theta = params[1];
        double phi = params[2];
        
        double sinT = Math.sin(theta);
        double cosT = Math.cos(theta);
        double sinP = Math.sin(phi);
        double cosP = Math.cos(phi);
        
        // fibre orientation
        double[] n = new double[3];
        n[0] = cosP * sinT;
        n[1] = sinP * sinT;
        n[2] = cosT;
        
        double b = scheme.getB_Value(i);
        double[] gd = scheme.getG_Dir(i);
            
        //dot product of the fibre orientation n and the wavenumber q
        double dotqn= n[0]*gd[0] + n[1]*gd[1] + n[2]*gd[2];
        
        double stickSignal = Math.exp(-b*stickDiff*dotqn*dotqn);
    
        return stickSignal;        
        
    }
    
    public RealMatrix getJacobian(final double[] modParams, final DW_Scheme scheme) {
    	RealMatrix jac = new RealMatrix(scheme.numMeasurements(), modParams.length);
    	
        double stickDiff = modParams[0];
        double theta = modParams[1];
        double phi = modParams[2];
        
        double sinT = Math.sin(theta);
        double cosT = Math.cos(theta);
        double sinP = Math.sin(phi);
        double cosP = Math.cos(phi);
        
        // fibre orientation
        double[] n = new double[3];
        n[0] = cosP * sinT;
        n[1] = sinP * sinT;
        n[2] = cosT;
        
        for(int i=0; i < scheme.numMeasurements();i++) {

            double b = scheme.getB_Value(i);
            double[] gd = scheme.getG_Dir(i);
            
            //dot product of the fibre orientation n and the wavenumber q
            double dotqn= n[0]*gd[0] + n[1]*gd[1] + n[2]*gd[2];
        
            double stickSignal = Math.exp(-b*stickDiff*dotqn*dotqn);

            // deriv wrt diffusivity
            jac.setEntry(i, 0, -b*dotqn*dotqn*stickSignal);

            // deriv wrt theta
            double cfac = -2*b*stickDiff*dotqn*stickSignal;
            jac.setEntry(i, 1, cfac*(gd[0]*cosT*cosP + gd[1]*cosT*sinP - gd[2]*sinT));

            // deriv wrt phi
            jac.setEntry(i, 2, cfac*(gd[1]*sinT*cosP - gd[0]*sinT*sinP));
        } 

        return jac;
    }


    /*public RealMatrix getJacobian(final double[] modParams, final DW_Scheme scheme) {
    	RealMatrix jac = new RealMatrix(scheme.numMeasurements(), modParams.length);
    	
    	//diff
    	RealMatrix diffJ = super.getJacobian(modParams, scheme, 0);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 0, diffJ.entry(r, 0));
        } 
        
        //theta
        double theta = modParams[1];
        
        if (theta < 1e-6)
        {
        	modParams[1] = 1e-6;
        }
        
    	RealMatrix thetaJ = super.getJacobian(modParams, scheme, 1);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 1, thetaJ.entry(r, 1));
        } 
        
        //phi
        
        
        RealMatrix phiJ = super.getJacobian(modParams, scheme, 2);
        
        for(int r=0; r < scheme.numMeasurements();r++){
            System.err.println("This is wrong in Stick.java");
        	jac.setEntry(r, 2, theta > 1e-6 ? phiJ.entry(r, 2) : 1);
        } 
    	
    	return jac;
    	
        }*/
    
    
}


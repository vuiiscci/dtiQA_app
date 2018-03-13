package models;

import optimizers.*;
import fitters.*;
import numerics.*;
import misc.*;
import models.compartments.CompartmentType;
import imaging.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Implements Behrens' ball-and-stick model
 * 
 * <dt>Description:
 * 
 * <dd>Model parameters are: {unweighted signal S0, intra-stick volume
 * fraction f, ball diffusivity, stick diffusivity, theta, phi}.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 * 
 */
public class BallAndStickModel extends ParametricModel {

    public BallAndStickModel() {
        
        super(CompartmentType.BALL.numParams + CompartmentType.STICK.numParams);
    }

    public RealMatrix getSignals(double[] modParams, DW_Scheme scheme) {

        // Compartment model parameters
        double S0 = modParams[0];
        double f = modParams[1];

        // Ball parameters
        double ballDiff = modParams[6];

        // Stick parameters
        double stickDiff = modParams[3];
        double theta = modParams[4];
        double phi = modParams[5];

        double sinT = Math.sin(theta);
        double cosT = Math.cos(theta);
        double sinP = Math.sin(phi);
        double cosP = Math.cos(phi);
        
        // fibre orientation
        double[] n = new double[3];
        n[0] = cosP * sinT;
        n[1] = sinP * sinT;
        n[2] = cosT;
        
        RealMatrix signals = new RealMatrix(scheme.numMeasurements(),1);
        
        for(int i=0; i<scheme.numMeasurements(); i++){
            double b = scheme.getB_Value(i);
            double ballSignal = Math.exp(-b*ballDiff);
                        
            double[] gd = scheme.getG_Dir(i);
            
            //dot product of the fibre orientation n and the wavenumber q
            double dotqn= n[0]*gd[0] + n[1]*gd[1] + n[2]*gd[2];
        
            double stickSignal = Math.exp(-b*stickDiff*dotqn*dotqn);

            signals.entries[i][0] = S0*(f*stickSignal + (1-f)*ballSignal);
        }
        
        return signals;
    }


    public double getSignal(double[] modParams, DW_Scheme scheme, int i) {

        // Compartment model parameters
        double S0 = modParams[0];
        double f = modParams[1];

        // Ball parameters
        double ballDiff = modParams[6];

        // Stick parameters
        double stickDiff = modParams[3];
        double theta = modParams[4];
        double phi = modParams[5];

        double sinT = Math.sin(theta);
        double cosT = Math.cos(theta);
        double sinP = Math.sin(phi);
        double cosP = Math.cos(phi);
        
        // fibre orientation
        double[] n = new double[3];
        n[0] = cosP * sinT;
        n[1] = sinP * sinT;
        n[2] = cosT;
        
        double b= scheme.getB_Value(i);
        double ballSignal = Math.exp(-b*ballDiff);    
            
        double[] gd = scheme.getG_Dir(i);
            
        //dot product of the fibre orientation n and the wavenumber q
        double dotqn= n[0]*gd[0] + n[1]*gd[1] + n[2]*gd[2];
        
        double stickSignal = Math.exp(-b*stickDiff*dotqn*dotqn);

        return S0*(f*stickSignal + (1-f)*ballSignal);
    }


}


package models.compartments;

import java.util.logging.Logger;

import numerics.RealMatrix;

import misc.LoggedException;
import models.ParametricModel;

import imaging.DW_Scheme;
import imaging.StejskalTannerScheme;

public class ListDistributionCompartment extends ParametricModel{

    /** logging object */
    private final Logger logger= Logger.getLogger(this.getClass().getName());
    
    /** array of compartment parameters we're summing over */
    private final double[][] param;
    
    /** array of weights for each parameter combination */
    private final double[] f;
    
    /** compartment model to use for evaluation */
    protected final ParametricModel model;
    
    
    
    /**
     * 
     */
    public ListDistributionCompartment(int numParams, ParametricModel model){
    	super(numParams);
    	
    	this.param= null;
    	
    	this.f= null;
    	
    	this.model=model;
    	
    }
    
    /**  
     *  constructor. takes a list of parameter combinations to sum over,
     *  and a set of weights (one for each parameter combination). In 
     *  addition to this, we specify a string that specifies
     */
    public ListDistributionCompartment(double[][] param, String compName, int numParams){
        
        super(numParams);
        
        CompartmentType type= CompartmentType.getCompartmentType(compName);
        
        this.f= new double[param.length];
                
        this.param= new double[param.length][];
        
        for(int i=0; i<param.length; i++){
            if(param[i].length-1!= type.numParams){
                throw new LoggedException("number of parameters at line "+i+" passed to list distribution("+param[i].length+" does not match number required by model");
            }
            
            f[i]= param[i][0];
            
            this.param[i]= new double[param[i].length-1];
            
            for(int j=1; j<param[i].length; j++){
                this.param[j-1]=param[j];
            }
            
        }
        
        this.model= CompartmentFactory.getCompartment(compName);
    }



	public RealMatrix getSignals(double[] params, DW_Scheme scheme){
        
        RealMatrix signals= new RealMatrix(scheme.numMeasurements(), 1);
        
        for(int i=0; i<scheme.numMeasurements(); i++){
            signals.setEntry(i, 0, getSignal(params, scheme, i));
        }
        
        return signals;
    }

    
    public double getSignal(double[] params, DW_Scheme scheme, int line){
    
        double signal=0.0;
        
        for(int i=0; i<this.param.length; i++){
            signal+=f[i]*model.getSignal(this.param[i], scheme, line);
        }
        
        return signal;
        
    }
    
    
}

package models.compartments;

import misc.LoggedException;
import models.ParametricModel;

/**
 * enumerates all types of single compartments currently implemented.
 * New compartment types should be added into here.
 * 
 * @author laura(panagio@cs.ucl.ac.uk)
 *
 */
public enum CompartmentType {

    /** isotropic diffusion */
    BALL(1),   
    
    /** delta-peaked anisotropy */
    STICK(3),
    
    /** long-time limit of diffusion constrained by impermeable cylinder */
    CYLINDERGPD(4),
    
    /** cylindrically symmetric tensor */
    ZEPPELIN(4),
    
    /** full 3d diffusion tensor */
    TENSOR(6),
    
    /** DIRECTIONALLY INTEGRATED STICKS */
    ASTROSTICKS(1),
    
    /** DIRECTIONALLY INTEGRATED CYLINDERS */
    ASTROCYLINDERS(2),
    
    /** SPHERE OF ZERO RADIUS */
    DOT(0),
    
    /** SPHERE */
    SPHEREGPD(2),
    
    /** CYLINDERS WITH GAMMA DISTRIBUTED RADII */
    GAMMADISTRIBRADIICYLINDERS(5);
    

    
    public final int numParams;
    
    CompartmentType(int numParams){
        
        this.numParams= numParams;
    }
    
    
    /**
     * returns the type enumeration based on the string handed over.
     * this does a direct comparison with the names of the enumeration 
     * types, so the string given has to match (case ignored). throws
     * an (uncaught) exception if the string is unmatched. 
     * 
     * @param compartmentName name of the compartment
     * 
     * @return matching type
     */
    public static CompartmentType getCompartmentType(String compartmentName){
        
        
        for(CompartmentType ct : CompartmentType.values()){
            if(ct.toString().equalsIgnoreCase(compartmentName)){
                return ct;
            }
        }
        
        throw new LoggedException("unrecognised compartment type name "+compartmentName);
    }
    
    
}
package simulation.geometry.elements;

import misc.LoggedException;
import numerics.Vector3D;
import simulation.SimulationParams;

/**
 * cylinder factory class. produces a cylinder of the desired 
 * type, size and orientation.
 * 
 * TODO: There are too many methods in here. This is clumsy and 
 * needs to change.
 *
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class CylinderFactory {

    public enum CylType{
        BASIC, 
        FACET,
        SQUASHY,
        MYELINATED,
        NESTED
    }

    
    /**
     * returns a cylinder of the specified type, radius, permeability, position and orientation.
     * 
     * TODO: this currently assumes constant radius and permeability. Needs fixing. Should specify
     * P, p, v, outer radius and gt whatever else from CL_Initialiser or SimulationParams.
     * 
     * @param P
     * @param simParams
     * @return
     */
    public static Cylinder getCylinder(double[] P, SimulationParams simParams){

        CylType type= SimulationParams.cylinderType;
        double r=SimulationParams.sim_r;
        double p= SimulationParams.sim_p;
        
        if(type==CylType.BASIC){
            return new BasicCylinder(P, r, p);
        }
        else if(type==CylType.FACET){
            
            Vector3D Pvec= new Vector3D(P[0], P[1], P[2]); 
            
            return new FacetCylinder(Pvec, r, SimulationParams.sim_num_facets, p);
        }
        else if(type==CylType.SQUASHY){
            
            return new SquashyCylinder(P, r, p);
        }
        else if(type==CylType.MYELINATED){
            
            return new MyelinatedCylinder(P, SimulationParams.sim_cyl_r1, 
                    SimulationParams.sim_r, SimulationParams.sim_cyl_D1, 
                    SimulationParams.sim_cyl_D2, p);
        }
        else if (type==CylType.NESTED){
        	
        	return new NestedCylinder(P, SimulationParams.sim_r, p, SimulationParams.cyl_nest_depth_max);
        
        }
        else{
            throw new LoggedException("unrecognised cylinder type "+type);
        }

        
    }    

    
    /**
     * returns a cylinder of the specified type, radius, permeability, position and orientation.
     * 
     * @param P poition of cylinder
     * @param r (outer) radius of cylinder
     * @param simParams simulation paramters object -- contains everything else!
     * 
     * @return a freshly-minted cylinder of the specified type.
     */
    public static Cylinder getCylinder(double[] P, double r, SimulationParams simParams){

        CylType type= SimulationParams.cylinderType;
        double p= SimulationParams.sim_p;
        
        if(type==CylType.BASIC){
            return new BasicCylinder(P, r, p);
        }
        else if(type==CylType.FACET){
            
            Vector3D Pvec= new Vector3D(P[0], P[1], P[2]); 
            
            return new FacetCylinder(Pvec, r, SimulationParams.sim_num_facets, p);
        }
        else if(type==CylType.SQUASHY){
            
            return new SquashyCylinder(P, r, p);
        }
        else if(type==CylType.MYELINATED){
            
            return new MyelinatedCylinder(P, SimulationParams.sim_gRatio*r, 
                    r, SimulationParams.sim_cyl_D1, 
                    SimulationParams.sim_cyl_D2, p);
        }
        else{
            throw new LoggedException("unrecognised cylinder type "+type);
        }

        
    }    

    
    
    
    public static Cylinder getCylinder(double[] V, double[] P, SimulationParams simParams){
        
        CylType type= SimulationParams.cylinderType;
        double r=SimulationParams.sim_r;
        double p= SimulationParams.sim_p;
        
        
        if(type==CylType.BASIC){
            return new BasicCylinder(V, P, r, p);
        }
        else if(type==CylType.FACET){
            
            Vector3D Pvec= new Vector3D(P[0], P[1], P[2]); 
            Vector3D Vvec= new Vector3D(V[0], V[1], V[2]);
            
            return new FacetCylinder(Vvec, Pvec, r, SimulationParams.sim_num_facets, p);
        }
        else if(type==CylType.SQUASHY){
            
            return new SquashyCylinder(V, P, r, p);
        }
        else if(type==CylType.MYELINATED){
            
            return new MyelinatedCylinder(V, P, SimulationParams.sim_cyl_r1, 
                    SimulationParams.sim_r, SimulationParams.sim_cyl_D1, 
                    SimulationParams.sim_cyl_D2, p);
        }
        else{
            throw new RuntimeException("unrecognised cylinder type "+type);
        }
        
        
    }


    /** 
     * get a new cylinder with directly specified type
     * 
     * @param P
     * @param type
     * @return
     */
    public static Cylinder getCylinder(double[] P, double r, double p, CylType type){

        if(type==CylType.BASIC){
            return new BasicCylinder(P, r, p);
        }
        else if(type==CylType.FACET){
            
            Vector3D Pvec= new Vector3D(P[0], P[1], P[2]); 
            
            return new FacetCylinder(Pvec, r, SimulationParams.sim_num_facets, p);
        }
        else if(type==CylType.SQUASHY){
            
            return new SquashyCylinder(P, r, p);
        }
        else if(type==CylType.MYELINATED){
            
            return new MyelinatedCylinder(P, SimulationParams.sim_cyl_r1, 
                    SimulationParams.sim_r, SimulationParams.sim_cyl_D1, 
                    SimulationParams.sim_cyl_D2, p);
        }
        else if (type==CylType.NESTED){
        	
        	return new NestedCylinder(P, SimulationParams.sim_r, p, SimulationParams.cyl_nest_depth_max);
        
        }
        else{
            throw new LoggedException("unrecognised cylinder type "+type);
        }

        
    }    

    
}

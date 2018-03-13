/* SubstrateFactory.java created on 28-Nov-2005
 * (simulation)
 * 
 * author: Matt Hall (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation.geometry.substrates;

import java.util.logging.Logger;

import simulation.SimulationParams;
import simulation.geometry.elements.CylinderFactory.CylType;


import misc.LoggedException;

/**
 *  Camino fibre reconstruction and tracking toolkit
 * 
 * SubstrateFactory (simulation)
 * 
 * 
 * 
 *
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public class SubstrateFactory {

    private static Logger logger=Logger.getLogger("simulation.GeometryFactory");
    
    
    public enum SubstrateType {    // cell-based geometries
        CELL_STRIPED,         /** striped lattice */
        CELL_PERC,            /** percolation cluster */
        
        // parallel cylinder geometries 
        CYL_1_FIXED,          /** parallel fixed radius cylinders in 1 dir */
        CYL_1_DISTRIB,        /** parallel distributed radius cylinders in 1 dir */
        CYL_1_INFLAM,         /** parallel cylinders on inflaming substrate 
                               * (randomly positioned, incrementing constant radii) */
        CYL_1_FACET,          /** facetted cylinders, rather than smooth surfaces */
        CYL_1_PERC,           /** randomly placed, parallel, non-overlapping cylinders */
        CYL_1_MYELIN,         /** sen and basser's myelinated cylinder substrates */
        CYL_1_STICKY,         /** sticky cylinder, square packed, parallel to z-axis */
        
        // crossing cylinder geometries
        CYL_2_FIXED,         /** basic crossing fibres */
        
        // sphere geometries
        SPH_1_pC,			 /** primitive cubic-packed spheres with fixed radius */
        
        // misc substrates
        TRI_PLY_MESH,         /** triangular mesh from PLY file */
        EMPTY                 /** substrate containing nothing at all */
    }
    
    /** factory method for making geometries.
     * 
     * @param substrateType geometry type code
     * @param simParams simulation parameters
     * 
     * @return new Geometry instance
     */
    public static Substrate getSubstrate(SubstrateType substrateType, SimulationParams simParams){
    
        if(substrateType==SubstrateType.CELL_STRIPED){
            return new StripedCellularLattice(simParams);
        }
        else if(substrateType==SubstrateType.CELL_PERC){
            return new PercCellularLattice(simParams);
        }
        else if(substrateType==SubstrateType.CYL_1_FIXED){
        	return new ParallelCylinderSubstrate(simParams);
        }
        else if(substrateType==SubstrateType.CYL_1_DISTRIB){
        	return new DistributedRadiusCylinderSubstrate(simParams);
        }
        else if(substrateType==SubstrateType.CYL_1_INFLAM){
        	return new SquashyInflammationSubstrate(simParams);
        	
        }
        else if(substrateType==SubstrateType.CYL_1_FACET){
        	//return new FacetCylinderSubstrate(simParams);
            return new ParallelCylinderSubstrate(CylType.FACET, simParams);

        }
        else if(substrateType==SubstrateType.TRI_PLY_MESH){
        	return new MeshSubstrate(simParams);
        }
        else if(substrateType==SubstrateType.CYL_1_PERC){
        	return new PercolationSubstrate(simParams);
        	
        }
        else if(substrateType==SubstrateType.CYL_1_MYELIN){
        	return new ParallelCylinderSubstrate(CylType.MYELINATED, simParams);

        }
        else if(substrateType==SubstrateType.CYL_2_FIXED){
        	return new CrossingCylinderSubstrate(simParams);
        }
        else if(substrateType==SubstrateType.SPH_1_pC){
        	return new RegularSphereSubstrate(simParams);
        }
        else if(substrateType==SubstrateType.EMPTY){
            return new EmptySubstrate(simParams);
        }
        else if(substrateType==SubstrateType.CYL_1_STICKY){
            return new StickyCylinderSubstrate(simParams);
        }
        else{
            String errMess=new String("unrecognised substrate type "+substrateType);
            logger.severe(errMess);
            throw new RuntimeException(errMess);
        }    
        
    }

 
    
    
    
    public static void main(String[] args) {
    }
}

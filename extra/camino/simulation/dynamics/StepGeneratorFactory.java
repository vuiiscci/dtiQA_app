/* StepGeneratorFactory.java created on 28-Nov-2005
 * (simulation)
 * 
 * author: Matt Hall (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation.dynamics;

import imaging.DW_Scheme;

import java.util.logging.Logger;

import misc.LoggedException;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;

/**
 * Camino fibre reconstruction and tracking toolkit
 * 
 * StepGeneratorFactory (simulation)
 * build step generators from given list of types.
 * 
 * the simulation step generator is a singleton. if 
 *
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public class StepGeneratorFactory {

    /** enumerator for step generator types */
    public enum StepType{ 
                FIXEDLENGTH,         /** fixed length steep generator */
                SPATIALLYVARYING,    /** steps with spatially varying length */
                CYLINDRICAL          /** steps on the surface of a cylinder */
    }
        
    /** logging object */
    private static Logger logger=Logger.getLogger("simulation.StepGeneratorFactory");
    
    /** current active step generator (singleton for each generator type) */
    private static StepGenerator flStepGen= null;
    
    /** singleton instance for cylindrical surface step generator */
    private static StepGenerator csStepGen= null;
    
    /** 
     * factory method
     * 
     * @param simParams simulation parameters object
     * @return a new step generator of the specified type
     */
    public static final StepGenerator getStepGenerator(SimulationParams simParams){
    
        StepType type= simParams.getStepType();
    	
        if(type==StepType.FIXEDLENGTH){
        	if(StepGeneratorFactory.flStepGen==null){
        		logger.info("instantiating fixed length step generator");
        		
        		StepGenerator stepGen= new FixedLengthStepGenerator(simParams);
        		StepGeneratorFactory.flStepGen=stepGen;
        		
        		return stepGen;
        	}
        	else{
        			return flStepGen;
        	}
        }
        else if(type==StepType.CYLINDRICAL){
            if(StepGeneratorFactory.csStepGen==null){
                logger.info("instantiating cylinderical surface step generator");
            
                StepGenerator stepGen= new CylindricalSurfaceStepGenerator(simParams);
            
                csStepGen= stepGen;
                
                return stepGen;
            }
            else{
                return csStepGen;
            }
        }
        else{
            String errMess=new String("unknown diffusion simulation step generator type code "+ type);
            
            logger.severe(errMess);
            throw new RuntimeException(errMess);
        }
    }
        
    
    /** 
     * parameters array generator
     *  
     * @param type the type of step generator being manufactured
     * @param simParams the simulation parameters
     * @param imParams the imaging parameters
     * 
     * @return the appropriate array of parameters for the desired
     *         step generator given the simulation and imaging
     *         parameters for the current simualtion
     */
    public static final double[] getStepParamsArray(StepType type, SimulationParams simParams){
        
        if(type==StepType.FIXEDLENGTH){
        	// step length is sqrt((2d)D dt) where d is dimensionality of system, 
        	//D diffusivity and dt timestep duration 
            double stepLength= Math.sqrt(6.0*DiffusionSimulation.DIFF_CONST * simParams.getDt());
            
            return new double[] {stepLength};
        }
        else if(type==StepType.CYLINDRICAL){
            double freeStepLength= Math.sqrt(6.0*DiffusionSimulation.DIFF_CONST * simParams.getDt());
            double surfaceStepLength= Math.sqrt(4.0*SimulationParams.sim_surfaceDiffusivity * simParams.getDt());
            
            return new double[] {freeStepLength, surfaceStepLength};
        }
        else{
            logger.warning("unknown step generator type "+type+" paramters array " +
            		"for fixed length step generator will be returned instead.");
            
            return getStepParamsArray(StepType.FIXEDLENGTH, simParams);
            
        }
        
    }
    
    
    
    }

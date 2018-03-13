package simulation.dynamics;

import java.util.logging.Logger;

import simulation.SimulationParams;
import simulation.geometry.substrates.ParallelCylinderSubstrate;
import simulation.geometry.substrates.StickyCylinderSubstrate;
import simulation.geometry.substrates.Substrate;

/**
 * factory class for producing step amender objects
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class StepAmenderFactory {

    public enum AmenderType{
        ELESTIC_REFLECTOR,
        STICK_AND_DIFFUSE
    }
    
	private static final Logger logger= Logger.getLogger("simulation.dynamics.StepAmenderFactory");
	
	public static final StepAmender getStepAmender(AmenderType type, Substrate substrate){
		
		if(type==AmenderType.ELESTIC_REFLECTOR){
			return new ElasticReflector(substrate);
		}
		else if(type==AmenderType.STICK_AND_DIFFUSE){
		    if(substrate instanceof StickyCylinderSubstrate){
		        return new StickAndDiffuse(substrate);
		    }
		    else{
		        logger.warning("Sticky and diffuse step amender is not implemented for non-cylinders. returning elastic reflector.");
		        return new ElasticReflector(substrate);
		    }
		}
		else{
			logger.warning("Unknown step amender type "+type+" returning elastic reflector");
			return new ElasticReflector(substrate);
		}
		
	}
	
}

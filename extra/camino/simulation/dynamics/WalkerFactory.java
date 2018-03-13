package simulation.dynamics;

import java.io.DataOutputStream;
import java.util.logging.Logger;

import simulation.SimulationParams;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.geometry.substrates.StickyCylinderSubstrate;
import simulation.geometry.substrates.Substrate;
import simulation.geometry.substrates.SubstrateFactory;
import simulation.geometry.substrates.SubstrateFactory.SubstrateType;
import simulation.measurement.SyntheticScan;
import simulation.dynamics.StepGeneratorFactory.StepType;
import tools.CL_Initializer;

/**
 * generate walkers of the chosen type
 *
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class WalkerFactory {

    /** logging object */
    private static final Logger logger= Logger.getLogger("simulation.dynamics.WalkerFactory");
    
    /** 
     * constructs the correct type of walker for the desired simulation 
     */
    public static final Walker getWalker(double[] r0, StepGenerator stepGen, Substrate substrate, 
                                    SyntheticScan synthScan, DataOutputStream trajWriter, SimulationParams simParams){
        
        if(SimulationParams.sim_geomType==SubstrateFactory.SubstrateType.CYL_1_STICKY){
            
            // construct new step params array. 2nd element is cylindrical step length
            double[] stepParams= StepGeneratorFactory.getStepParamsArray(StepType.CYLINDRICAL, simParams);
            
            // set the step params object
            simParams.setStepParams(stepParams);
            
            SimulationParams.sim_stepType= StepGeneratorFactory.StepType.CYLINDRICAL;
            
            // get surface step generator (singleton)
            StepGenerator csStepGen= new CylindricalSurfaceStepGenerator(simParams, (StickyCylinderSubstrate)substrate);
            
            return new StickyWalker(r0, stepGen, csStepGen, substrate, synthScan, trajWriter, 
                    true, SimulationParams.sim_p_unstick);
        }
        else{
            return new Walker(r0, stepGen, substrate, synthScan, trajWriter);
        }
    }
}

package simulation.measurement;

import java.util.logging.Logger;

import simulation.SimulationParams;
import simulation.dynamics.Walker;

/**
 * constructs stats modules of a specified type from simulation 
 * params and walkers array. also specifies the stats module type
 * enum that is used to specify which module we're after.
 *
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class StatisticsModuleFactory {

    /** logging object */
    private static final Logger logger = Logger.getLogger("simulation.measurement.StatisticsModuleFactory");
    
    /** enum for stats module types */
    public static enum StatsModuleType{ 
                                        MS_DISP,        /** mean-squared displacements only */
                                        STICKY_STATS    /** equilibrium times, occupation time stats, occupation levels */
                                      }
    
    /** factory method */
    public static final StatisticsModule getStatsModule(Walker[] walker){
        
        if(SimulationParams.sim_StatsModType== StatsModuleType.MS_DISP){
            
            return new MSdispStatsModule(walker);
        }
        else if(SimulationParams.sim_StatsModType== StatsModuleType.STICKY_STATS){
            
            return new StickyStatsModule(walker);
        }
        else{
            
            logger.warning("unrecognised statistics module type. returning mean-squared displacement type");
            return new MSdispStatsModule(walker);
        }
    
    
    }
    
    
}

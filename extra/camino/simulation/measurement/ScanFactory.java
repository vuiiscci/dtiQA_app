package simulation.measurement;

import imaging.RectGradSteTanScheme;
import imaging.RectQuadraticGradSteTanScheme;
import imaging.SimulableScheme;

import java.util.logging.Logger;

import misc.LoggedException;
import simulation.SimulationParams;
import simulation.dynamics.Walker;
import simulation.geometry.substrates.Substrate;

/**
 * factory class for creating scans (and other measurement modules).
 * it also contains the constant identifiers for various scan types.
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class ScanFactory {

    /** enumeration of synthetic scan types */
    public enum ScanType{
        /** Simple idealised two-lobed sequence */
        PGSE_SCAN,          
        /** Four-pulse sequence with two 180 pulses */
        TWICE_REFOCUSED,   
        /** independent of explicit form -- uses new scheme framework */
        AGNOSTIC
    }
    
    
	/** logging object */
	private static Logger logger= Logger.getLogger("simulation.measurement.ScanFactory");
	
	/** factory method. scan objects are not singletons 
	 * 
	 * @param simParams simulation parameters object
	 * @param scheme acquisition scheme
	 * @param substrate the substrate
	 * @param walker array of walkers
	 * @param scanType what type of scan do we want?
	 * 
	 * @return a new scan object based on the above
	 */
	public static final SyntheticScan getMeasurementModule(SimulationParams simParams, SimulableScheme scheme, 
																Substrate substrate, Walker[] walker){
		

		// everything with linear gradients is an agnostic scan, otherwise Quadratic.
		if(scheme instanceof RectQuadraticGradSteTanScheme){
			return new QuadraticGradientScan((RectQuadraticGradSteTanScheme)scheme, walker, substrate);
		}
		
		
	    return new AgnosticScan(scheme, walker, substrate);
		
	}
	
	
}

package simulation.geometry.substrates;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Logger;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.Walker;
import simulation.geometry.elements.SubstrateObject;
import simulation.geometry.elements.Triangle;


/**
 * this is for debugging purposes. A substrate with no structure in it
 * so that diffusion is isotropic in all directions. Substrate is cubic
 * with a given cell size.
 *
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class EmptySubstrate extends Substrate {

    /** logging object */
    private final Logger logger = Logger.getLogger(this.getClass().getName());
    
    /** dimensionality of space */
    private final int D= DiffusionSimulation.D;
    
    public EmptySubstrate(SimulationParams simParams){
        super(simParams, new double[]{SimulationParams.sim_L, SimulationParams.sim_L, SimulationParams.sim_L});
        super.subsObj = new SubstrateObject[0];
    }
    
    /**
     * test constructor, size only
     * 
     * @param substrateDims
     */
    public EmptySubstrate(double[] substrateDims){
        super(substrateDims);
    }
    
    
    /**
     * @return always false.
     */
    public boolean crossesMembrane(Walker walker, double[] offset,
            double[] step, double[] normal, double[] d, boolean skipCurrent,
            double origLength, boolean[] in, double[] p, boolean report, FileWriter debugWriter) {
        
        return false;
    }

    /**
     * @return half substrate size
     */
    public double getPeakCoord() {
        return super.L[0]/2;
        
        //return 5e-6;
    }

    /**
     * @return substrate dims from superclass
     */
    public double[] getSubstrateSize() {
        double[] subsSize= new double[D];
        
        for(int i=0; i<D; i++){
            subsSize[i]=super.L[i];
        }
        
        return subsSize;
    }

    /**
     * overrides Substrate's position checking code so that the 
     * generated position is always ok. No barriers mean no need to check!
     * 
     * @param position to check
     * @param R walker radius
     * 
     * @return always true
     */
    public boolean positionOk(double[] r0, double R){
    	
    	return true;
    }

    public Collection<Triangle> getTriangles(){
    	
    	return new ArrayList<Triangle>();
    	
    }
    

    /**
     * does nothing
     */
    public void init() {

    }

    /**
     * @return always false
     */
    public boolean intracellular(Walker walker) {
        
        return false;
    }

}

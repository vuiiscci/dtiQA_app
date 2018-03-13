package simulation.dynamics;

import java.io.DataOutputStream;
import java.util.logging.Logger;

import misc.LoggedException;
import numerics.MTRandom;

import simulation.DiffusionSimulation;
import simulation.geometry.substrates.StickyCylinderSubstrate;
import simulation.geometry.substrates.Substrate;
import simulation.measurement.SyntheticScan;
import tools.CL_Initializer;

/**
 * This is a kind of Walker that stick to barriers and
 * diffuse around on the surface for a while. This is
 * different to regular walkers, which have a simple
 * interaction with surfaces and only one form of 
 * behaviour.
 * 
 * Because these walkers have two modes of behaviour 
 * (stuck and free) they need an extra step generator,
 * an extra flag to indicate if it's stuck or not and
 * and new version of the update method.
 *
 * note that this currently supports walkers sticking to
 * only one type of object or freely diffusing in the 
 * bulk.
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class StickyWalker extends Walker {

    /** logging object */
    private final Logger logger= Logger.getLogger(this.getClass().getName());
    
    /** dimensionality of space */
    private final int D= DiffusionSimulation.D;
    
    /** are we in the bulk or stuck to a the surface? */
    public boolean free;
    
    /** step generator for when bound to a surface */
    private final StepGenerator surfaceStepGen;
    
    /** probability of unsticking from the surface */
    private final double p_unstick;

    /** random number generator for unsticking */
    private static MTRandom surfaceTwister= new MTRandom(CL_Initializer.seed+27834);
    
    /** flag to say whether the spin was initially extracellular or intracellular */
    public boolean wasExtracellular;
    
    /** all substrates are sticky for these walkers */
    private final StickyCylinderSubstrate substrate;
    
    /** counts the number of times we get stuck to the surface */
    public int numSticks = 0;
    
    /** the last time we got stuck */
    public double lastStickTime = 0.0;
    
    /** total time stuck */
    public double totalTimeStuck = 0.0;
    
    /**
     * main constructor. this should be used if a full simulation
     * is being run.
     * 
     * @param r0 initial position
     * @param freeStepGen step generator for when diffusing in the bulk
     * @param surfaceStepGen step gen for when bound on the surface
     * @param substrate substrate we're diffusing in
     * @param scan measurement module
     * @param trajWriter trajectories writer if present
     * @param free flag saying if we're stuck to a surface or not
     */
    public StickyWalker(double[] r0, StepGenerator freeStepGen, StepGenerator surfaceStepGen,
            Substrate substrate, SyntheticScan scan, DataOutputStream trajWriter, 
            boolean free, double p_unstick) {
        super(r0, freeStepGen, substrate, scan, trajWriter);

        this.surfaceStepGen= surfaceStepGen;
        
        this.free=free;
        
        this.p_unstick= p_unstick;
        
        this.substrate= (StickyCylinderSubstrate)substrate;
    }

    /**
     * test constructor with specified stick state. do not
     * use for full simulation.
     * 
     * @param r0 initial position
     * @param free are we stuck?
     */
    public StickyWalker(double[] r0, StepGenerator stepGen, Substrate substrate, boolean free) {
        super(r0, substrate);
        
        this.free=free;
        
        this.surfaceStepGen= stepGen;
        
        this.p_unstick=0.0;
        
        this.substrate=(StickyCylinderSubstrate)substrate;
        
    }

    /**
     * test constructor with unspecified stick state, assumes 
     * unstuck. do not use for full simulation
     * 
     * @param r0 initial position
     */
    public StickyWalker(double[] r0) {
        super(r0);
        
        this.free=false;

        this.surfaceStepGen= null;

        this.p_unstick=0.0;

        this.substrate= null;
    }

    
    
    /**
     * update walkers position.
     *
     *  @param t current time in seconds
     *  @param ti current time in timesteps
     *  @param i index of walker
     */
    public void update(double t, int ti, int n){
        
        if(free){
            super.update(t, ti, n, false);
        }
        else{
            // move on the surface
            step= surfaceStepGen.getStep(this);

            // amend the step
            substrate.amend(this, step, t, n, false, null);
            
            makeStep(step);
        }        
    }
    
    
    public void getSubsCoords(double[] subsCoords){
        
        final double[] zeroOffset= new double[]{0.0, 0.0, 0.0};
                
        substrate.getSubstrateCoords(r, zeroOffset, subsCoords);
        
        
    }
    
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        // TODO Auto-generated method stub

    }

}

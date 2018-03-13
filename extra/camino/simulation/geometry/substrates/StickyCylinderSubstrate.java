/**
 * 
 */
package simulation.geometry.substrates;

import java.io.FileWriter;
import java.util.logging.Logger;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StickyWalker;
import simulation.dynamics.Walker;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.geometry.elements.Cylinder;
import simulation.geometry.elements.CylinderFactory;
import simulation.geometry.elements.CylinderFactory.CylType;

/**
 *
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class StickyCylinderSubstrate extends ParallelCylinderSubstrate {

    /** logging object */
    private final Logger logger= Logger.getLogger(this.getClass().getName());
    
    /** dimensionality of space */
    private final int D= DiffusionSimulation.D;

    /** permeability of cylinder */
    private final double p;
    
    /** sticking probability */
    private final double p_stick;
    
    /** array of T2s */
    private final double[] T2;
    
    /** constructor from superclass */
    public StickyCylinderSubstrate(SimulationParams simParams) {
        super(simParams);

        logger.info("replacing regular reflecting cylinders with sticky surfaces"); 
        
        // this substrate has only a single cylinder
        Cylinder[] cylinder= new Cylinder[1];
    
        // position the cylinder in the middle of the substrate
        double pos= SimulationParams.sim_R;
        double[] P= new double[]{pos, pos, pos};
        
        // instantiate cylinder 
        cylinder[0]= CylinderFactory.getCylinder(P, simParams);
    
        // set the cylinder arrays and related quantities in superclass
        setCylinders(cylinder);
        
        //set cylinder permeability here
        this.p = SimulationParams.sim_p;
        
        // set the sticking probability
        this.p_stick= SimulationParams.sim_p_stick;
        
        logger.info("permeability p= "+p+", p_stick= "+p_stick+", p_unstick= "+SimulationParams.sim_p_unstick);        
        
        // set the relaxivities
        this.T2= new double[2];        
        for(int i=0; i<T2.length; i++){
            T2[i]= SimulationParams.sim_T2[i];
        }
        
        logger.info("bulk T2= "+T2[0]+"s, surface T2= "+T2[1]+"s.");
        
    }

    
    
    
    
    
    
    /**
     * wrapper for crosses membrane method. if walker
     * is stuck to the surface, it returns true once and then
     * false thereafter. this allows the step amender to handle
     * unsticking as well as sticking and also unifies sticking
     * and unsticking with permeability checks in a neat way
     * 
     * @param walker the (stick) walker making the step
     * @param offset net step made so far
     * @param stepVector the step to check
     * @param normal space to store normal at intersection point
     * @param d space to store distance of int point to origin dot normal
     * @param skipCurrent skip the membrane we're sitting on?
     * @param origLength length fo original step
     * @param in space to store initially in flag
     * @param p space to store permeability of interacting membrane
     * 
     * @return true if an intersection, false otherwise
     */
    public boolean crossesMembrane(Walker walker, double[] offset,
            double[] stepVector, double[] normal, double[] d,
            boolean skipCurrent, double origLength, boolean[] in, double[] p, boolean report, FileWriter debugWriter) {
       
        if(((StickyWalker)walker).free){     
            // if we're free, default to superclass behaviour
           boolean crosses= false;
           
           try{
        	   super.crossesMembrane(walker, offset, stepVector, normal, d,
                skipCurrent, origLength, in, p, false, null);
           }
           catch(StepRejectedException sre){
        	   sre.printStackTrace();
           }
           
           return crosses;
        }
        else{
            // set permeability (short cut!)
            p[0]=this.p;            
            
            // return true the first time and false thereafter
            return !skipCurrent;
        }
    }


    /**
     * returns the change in magnetisation at a given location and time.
     * overrides superclass method to include T2 effects. separate T2s are
     * specified on the surface and in the bulk (see constructor) and are
     * used to obtain an additive change in log-magenetization d(logM)
     * 
     * d(\log M) = -\frac{dt}{T_2(r)}
     * 
     * where dt is the time increment (t-tLast) and T_2(r) = T2[0] in the
     * bulk and T2[1] on the surface.
     * 
     * This assumes (without checking) that the walker passed to it will be a
     * StickyWalker. Seeing as the sticky substrate won't run without this
     * being true, it isn't a strong assumption but this will throw a
     * class-cast exception if passed another kind of walker.
     * 
     * @param walker the walker
     * @param t current time (end of timestep)
     * @param tLast the last time we were here (beginning of timestep)
     * 
     * @return change in magnetization according to compartment
     */
    public double getLogMagnetisationChange(Walker walker, double t, double tLast) {
        
        double dt= t-tLast;
        double T2;
        
        if(((StickyWalker)walker).free){
            T2=this.T2[0];
        }
        else{
            T2=this.T2[1];
        }
        
        // negative T2s mean that this effect is ignored
        if(T2<0.0){
            return 0.0;
        }
        
        // otherwise return the decay increment
        return -dt/T2;
    }




    /**
     * override intracellular check to account for being stuck to the surface
     * which would otherwise cause annoying floating-point precision related 
     * difficulties and spurious crossings
     * 
     * @param walker the walker to check. this will be explicitly cast into a StickyWalker
     * 
     * @return true if intracellular, otherwise false
     */
    public boolean intracellular(Walker walker){
        
        // cast walker as StickyWalker (should always work)
        StickyWalker stWalker= (StickyWalker)walker;
        
        if(stWalker.free){
            // if free, use the usual check
            return super.intracellular(walker);
        }
        
        // otherwise return the walker's original compartment
        return !stWalker.wasExtracellular;
    }
    
}

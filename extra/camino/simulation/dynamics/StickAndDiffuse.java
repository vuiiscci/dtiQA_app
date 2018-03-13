package simulation.dynamics;

import java.util.logging.Logger;

import numerics.MTRandom;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.geometry.substrates.StickyCylinderSubstrate;
import simulation.geometry.substrates.Substrate;
import tools.CL_Initializer;

public class StickAndDiffuse extends ElasticReflector {

    /** logging object */
    private final Logger logger= Logger.getLogger(this.getClass().getName());
    
    /** dimensionality of space */
    private final int D= DiffusionSimulation.D;
    
    /** sticking probability */
    private final double p_stick;
    
    /** unsticking probability */
    private final double p_unstick;
    
    /** random number generator for sticking and unsticking */
    private final MTRandom surfaceTwister;
    
    /** copy of the substrate the amender's attached to */
    private final StickyCylinderSubstrate substrate;
    
    /** constructor. takes substrate 
     * 
     * @param substrate the substrate
     */
    public StickAndDiffuse(Substrate substrate){
        super(substrate);
        
        
        this.p_stick= SimulationParams.sim_p_stick;
        this.p_unstick= SimulationParams.sim_p_unstick;
        
        this.surfaceTwister= new MTRandom(CL_Initializer.seed+475620);
        
        this.substrate=(StickyCylinderSubstrate)substrate;
        
        logger.info("cylindrical surface stick and diffuse step amender constructed");
        logger.info("p_stick = "+p_stick+" p_unstick = "+p_unstick);        

    }
    
    
    /**
     * if a walker is free, it sticks it with a fixed prob. if
     * stuck, it returns a null amended step and stick the walker,
     * otherwise it
     */
    public void amendStep(Walker walker, double[] subsCoords, double[] offset,
            double[] step, double[] normal, double[] d, double origLength,
            double[] toBarrier, double[] amended, double[] unamended,
            double[] L, double t, int n) {
        
        if(((StickyWalker)walker).free){
            // if the walker is free, check if it gets stuck
            
            // amend using elastic reflector
            super.amendStep(walker, subsCoords, offset, step, normal, d, 
                    origLength, toBarrier, amended, unamended, L, t, n);
            
            // only perform the sticking check if we've hit the cylinder
            // not if we've just hit the cell boundary
            if(!substrate.intersectsBoundary){

                double p= surfaceTwister.nextDouble();
                if(p<p_stick){
                    // if in here, we're sticking. mark our current compartment
                    ((StickyWalker)walker).wasExtracellular= !substrate.intracellular(walker);
    
                    // set flag to say whether walker was extracellular
                    ((StickyWalker)walker).free= false;
                    
                    ((StickyWalker)walker).numSticks++;
                    
                    ((StickyWalker)walker).lastStickTime= t;
                    
                    // set amended and unamended steps to null vectors
                    for(int i=0; i<D; i++){
                        amended[i]=0.0;
                        unamended[i]=0.0;
                    }    
                }
            }
        }
        else{
            // if the walker is already stuck to the surface then
            // make the step no matter what and check if we're 
            // going to unstick the walker
            for(int i=0; i<D; i++){
                toBarrier[i]= 0.0;
                amended[i]= step[i];
                unamended[i]= step[i];
            }
            
            // check if the walker will leave the surface
            double p= surfaceTwister.nextDouble();
        
            if(p<p_unstick){
                // free the walker
                ((StickyWalker)walker).free= true;
                
                ((StickyWalker)walker).totalTimeStuck += (t-((StickyWalker)walker).lastStickTime);
                
                // displace the spin by a small amount away from surface
                // in amended and unamended step.
                final double[] towardsCentre= new double[D];
                final double[] awayFromCentre= new double[D];
                
                double[] cylPos= substrate.cylinder[0].getPosition();
                double cylRad= substrate.cylinder[0].getRadius();
                
                final double[] posInCylCoords= new double[D];                    
                
                for(int i=0; i<D; i++){
                    posInCylCoords[i]=subsCoords[i]-cylPos[i];
                }
                
                // convert to polar coords to construct 
                final double[] polarPos= new double[2];
                
                polarPos[0]= Math.sqrt(posInCylCoords[0]*posInCylCoords[0]+posInCylCoords[1]*posInCylCoords[1]);
                polarPos[1]= Math.atan2(posInCylCoords[1], posInCylCoords[0]);
                
                //double dr= cylRad/100;
                double dr= origLength/10;
                
                // step towards centre as radius increment
                towardsCentre[0]=-dr*Math.cos(polarPos[1]);
                towardsCentre[1]=-dr*Math.sin(polarPos[1]);
                towardsCentre[2]=0.0;
                
                // step away from centre is the same negated
                for(int i=0; i<D; i++){
                    awayFromCentre[i]=-towardsCentre[i];
                }
                
                // add a small radial offset to step to insure the 
                // spin is in a chosen compartment (same if amended, 
                // other if unamended).
                if(((StickyWalker)walker).wasExtracellular){
                    
                    for(int i=0; i<D; i++){
                        amended[i]=awayFromCentre[i];
                        unamended[i]=towardsCentre[i];
                    }
                }
                else{
                    for(int i=0; i<D; i++){
                        amended[i]=towardsCentre[i];
                        unamended[i]=awayFromCentre[i];
                    }
                }
            }
        }
    }
}

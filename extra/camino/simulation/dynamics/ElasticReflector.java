package simulation.dynamics;

import java.util.logging.Logger;

import numerics.Rotations;

import simulation.DiffusionSimulation;
import simulation.geometry.substrates.Substrate;

/** 
 * reflects a step off a barrier with a given normal.
 * sum of lengths of components is equal to length of
 * original step.
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 * 
 */ 
public class ElasticReflector implements StepAmender {

	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;
		
	/** the substrate the amender belongs to */
	private final Substrate substrate;
	
	
	/** constructor
	 * 
	 * @param substrate the substrate the amender belongs to
	 */
	public ElasticReflector(Substrate substrate){
	    
	    this.substrate= substrate;
	}
	
	
	
    /** 
     * amends the step for the walker by relecting it elastically off
     * the appropriate membrane. The geometry is represented as a ray starting 
     * at the walker's current position, projecting in the direction of the step.
     * 
     * Firstly the intersection is calculated via the parametric representation
     * R_0 + tR, where R_0 is the walkers position, R is it's direction and t
     * is a real number. The distance from the current position to the plane is
     * gien by the dot product R_0.n where n is the normal.
     * 
     * The new step is the sum of the incident and reflected rays with lengths determined
     * by t and step length. The reflected step is given by
     * 
     * w=(v - (v.n)n(|step|-t) where v is the incident part of the step. The amended step is 
     * therefore S'=v+w = (t/|step|)step + ((t/|step|)S - 2(t/|step|)(S.n)n)(|S|-t)
     * 
     * which is much easier to explain with a diagram!
     * 
     * 
     * @param walker the walker
     * @param step the step the walker would take
     * @param normal the normal of the membrane
     * @param d distance of plane from origin (array of length 1)
     * @param origLength the length of the original step before ANY amendment
     * @param toBarrier space to store step that will take the walker to the barrier but not through it
     * @param amended space to store amended component of step.
     * @param unamended space to store the step through the barrier.
     * 
     */
    public void amendStep(Walker walker, double[] subsCoords, double[] rawOffset, double[] rawStep, double[] normal, double[] d, 
    						double origLength, double[] toBarrier, double[] amended, 
    						double[] unamended, double[] L, double time, int n){
        
        
        double normalDotPos=0.0;    // the point-plane distance
        double normalDotStep=0.0;   // gives the angle between normal and step vector 
        
        double modStep=0.0;

        double[] walkerPos= new double[D];

        double[] step= substrate.mapStepIntoSubstrate(walker, rawOffset, rawStep);
        double[] offset;
       
        if(step!=rawStep){
            offset= substrate.mapStepIntoSubstrate(walker, rawOffset, rawOffset);
        }
        else{
            offset=rawOffset;
        }
        
        for(int i=0;i<D;i++){
            modStep+=step[i]*step[i];
            walkerPos[i]=subsCoords[i];
            normalDotPos+=normal[i]*walkerPos[i];
            normalDotStep+=normal[i]*step[i];            
        }
        
        
        
        
        modStep=Math.sqrt(modStep);
        
        // angle between step and normal
        double cosTheta= normalDotStep/modStep;
               
        // point-plane distance - distance of walker to plane,
        double distToPlane=Math.abs(normalDotPos - d[0]);

        // catch the (extremely rare) case where a walker is actually sitting on the membrane
        if(distToPlane<=modStep*1E-12){
            distToPlane=modStep*1E-12;
        }
        
        
        // divide by cos theta (angle between step and normal to plane) to give 
        // distance to intersection of step and plane
        double t=1.0;
        
        if(cosTheta!=0.0){
            t=Math.abs(distToPlane/cosTheta);
        }
        else{
            logger.warning("WARNING: attempt to amend a step that is parallel to membrane!");
            logger.warning("walker.r=("+walker.r[0]+","+walker.r[1]+","+walker.r[2]+"), " +
            		       "step = ("+ step[0]+","+step[1]+","+step[2]+") normal=("+normal[0]+
            		       ","+normal[1]+","+normal[2]+")");
        }
       
        // check that distance is not longer than step
        if(t>modStep){
            String errMess= new String("erroneously detected membrane crossing in diffusion sim!\n" +
            		"distance to membrane "+t+" is greater than step length "+modStep+" d[0]= "+d[0] +
            		" distToPlane= "+distToPlane+ " cosTheta= "+cosTheta+" t="+time+" walker "+n);
            logger.severe(errMess);
            throw new RuntimeException(errMess);
        }
        
        // now costruct incoming and outgoing paths
        double[] v=new double[D];
        double vDotNormal= 0.0;
        
        // incoming step is original step scaled to length t
        for(int i=0; i<D; i++){
            v[i]=t*(step[i]/modStep);      				// step to barrier (internal)
            vDotNormal+=v[i]*normal[i];					// dot prod of above with barrier normal
            //System.err.println("Factor: "+Substrate.factor);
            //toBarrier[i]=Substrate.factor*v[i];			// step to barrier (external)
            toBarrier[i]=v[i];
            unamended[i]=(modStep-t)*(step[i]/modStep);	// unreflected step through barrier
        }
        
        
        for(int i=0; i<D; i++){
            double w_i;
            
            // reflected incoming path
            w_i=(v[i] - 2.0*vDotNormal*normal[i]);

            // rescale length
            w_i/=t;
            w_i*=(modStep-t);

            amended[i]= w_i;
        }
        
        
        double newLength=0.0;
        double amendLength=0.0;
        double unamendLength=0.0;
        for(int i=0; i<D; i++){
            newLength+=step[i]*step[i];
            amendLength+=amended[i]*amended[i];
            unamendLength+=unamended[i]*unamended[i];
        }
        
        newLength=Math.sqrt(newLength);
        amendLength=Math.sqrt(amendLength);
        unamendLength=Math.sqrt(unamendLength);
        
        if(newLength>modStep){
            logger.warning("amended step longer than original!");
        }
        
        if(amendLength/origLength<=1E-12){
        	for(int i=0; i<D; i++){
        		amended[i]=0.0;
        	}
        }

        if(unamendLength/origLength<=1E-12){
        	for(int i=0; i<D; i++){
        		unamended[i]=0.0;
        	}
        }
        
        // finally, check if we need to un-map the step
        if(step!=rawStep){
            double[] unmappedBarrier= substrate.unmapStepFromSubstrate(subsCoords, offset, toBarrier);
            double[] unmappedAmended= substrate.unmapStepFromSubstrate(subsCoords, offset, amended);
            double[] unmappedUnamended= substrate.unmapStepFromSubstrate(subsCoords, offset, unamended);
            
            for(int i=0; i<D; i++){
                toBarrier[i]=unmappedBarrier[i];
                amended[i]=unmappedAmended[i];
                unamended[i]=unmappedUnamended[i];
            }
        }
    }

}

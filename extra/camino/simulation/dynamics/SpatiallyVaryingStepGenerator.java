package simulation.dynamics;

import java.util.logging.Logger;

import numerics.MTRandom;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.geometry.substrates.Substrate;
import tools.CL_Initializer;

/**
 * step generator for substrates where the diffusivity can vary
 * spatially. This could mean, for example, a substrate with a
 * different diffusivity in the intracellular and extracellular
 * compartments.
 * 
 * Steps are of fixed length for each compartment
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class SpatiallyVaryingStepGenerator implements StepGenerator {

	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;
	
	/** simulation parameters object */
	private final SimulationParams simParams;
	
	/** substrate object */
	private final Substrate substrate;
	
	/** random number generator */
	private final MTRandom stepTwister=new MTRandom(CL_Initializer.seed+17);
	
	/** space to store next step */
	private final double[] step= new double[D];
	
	
	/** 
	 * constructor. needs simulation params and substrate object
	 */
	public SpatiallyVaryingStepGenerator(SimulationParams simParams, Substrate substrate){
		
		this.simParams=simParams;
		
		this.substrate= substrate;
	}
	
	/* (non-Javadoc)
	 * @see simulation.dynamics.StepGenerator#getBorder()
	 */
	public double getBorder() {
		// TODO Auto-generated method stub
		return 0;
	}

	/**
	 * queries the substrate for the diffusivity at the current location,
	 * the constructs a step of the appropriate length.
	 * 
	 * @return step vector
	 * 
	 * @see simulation.dynamics.StepGenerator#getStep(simulation.dynamics.Walker)
	 */
	public double[] getStep(Walker walker) {

		double diff=substrate.getDiffusivityAt(walker.r);
		
		double length=Math.sqrt(6.0*diff/simParams.getDt());
		
        if(D==1){
            if(stepTwister.nextDouble()<0.5){
                step[0]=-length;
            }
            else{
                step[0]= length;
            }
        }
        else if(D==2){
            double theta= 2.0*Math.PI*stepTwister.nextDouble();
            
            step[0]=length*Math.cos(theta);
            step[1]=length*Math.sin(theta);
        }
        else if(D==3){
            double theta= 2.0*Math.PI*stepTwister.nextDouble();
            double cosPhi = 2.0*stepTwister.nextDouble()-1.0;
            
            double cosTh= Math.cos(theta);
            double sinTh= Math.sin(theta);
            
            double sinPhi= Math.sqrt(1.0-cosPhi*cosPhi);
                        
            step[0]= length*cosTh*sinPhi;
            step[1]= length*sinTh*sinPhi;
            step[2]= length*cosPhi;
            
        }
        else{
            String errMess= new String("steps of dimension "+D+" are not yet implemented. sorry.");
            logger.severe(errMess);
            throw new RuntimeException(errMess);
        }
        
        return step;
	}

	/**
	 * returns the type code of the generator
	 * 
	 * @return StepGeneratorFactory.SPATIALLYVARYING
	 * 
	 * @see simulation.dynamics.StepGenerator#getType()
	 */
	public StepType getType() {
		return StepType.SPATIALLYVARYING;
	}

	
	
	/**
     * @return finite size for a walker based on step length
     */
     public final double getWalkerRadius(){

    	 double len=simParams.getStepParams()[0];

    	 return len/LengthStepRatio;
     }
}

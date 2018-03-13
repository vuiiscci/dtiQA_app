/* Step.java created on 25-Nov-2005
 * (simulation)
 * 
 * author: ucacmgh
 * 
 */
package simulation.dynamics;

import simulation.dynamics.StepGeneratorFactory.StepType;

/**
 *  Camino fibre reconstruction and tracking toolkit
 * 
 * Step (simulation)
 * 
 *  Contains a method prototype to get a step vector for 
 *  a Walker ro use during walker update
 * 
 *
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public interface StepGenerator {
    
	public static final double LengthStepRatio=1E6;
	
    /** returns a step vector as array of doubles 
     * @return step vector 
     */
    public double[] getStep(Walker walker);

    /**
     * get the border width for cloning 
     */
    public double getBorder();
    
    /**
     * get the step generator type
     */
    public StepType getType();
    
    /**
     * get the sze of a walker for steps of this kind 
     */
     public double getWalkerRadius();
}

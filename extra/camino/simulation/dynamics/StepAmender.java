package simulation.dynamics;


/**
 * interface describing a step amendment object.
 * a StepAmender is an object that handles the
 * interaction between walkers and barriers.
 * 
 * This is implemented in a separate object so that
 * substrates with different interaction mechanisms 
 * don't need to be re-implemented.
 * 
 * examples of implementations are ElasticReflector
 * or Scatterers
 * 
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public interface StepAmender {

	/**
	 * amend a step given necessary geometric information
	 * and places to store results etc
	 * 
	 * @param walker the walker making the step
	 * @param the position of the walker on the substrate
	 * @param offset displacement from walker location to start of step
	 * @param step the step it wants to make
	 * @param normal surface normal at interaction point
	 * @param d distance to interaction point
	 * @param origLength original length of step
	 * @param toBarrier storage for step to barrier
	 * @param amended storage for amended portion of step
	 * @param unamended storage for unamended
	 * @param L substrate size
	 * @param t current time
	 * @param i index of walker
	 */
	public void amendStep(Walker walker, double[] subsCoords, double[] offset, double[] step, double[] normal, double[] d, 
			double origLength, double[] toBarrier, double[] amended, 
			double[] unamended, double[] L, double t, int i);
	
	
}

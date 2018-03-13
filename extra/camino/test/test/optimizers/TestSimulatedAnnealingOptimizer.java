package optimizers;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;

import java.util.Random;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>OrderedAcqWeightedMinimizer.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestSimulatedAnnealingOptimizer.java,v 1.1 2006/07/26 21:28:19 ucacpco Exp $
 * @author  Philip Cook
 * @see misc.DT
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public abstract class TestSimulatedAnnealingOptimizer extends TestCase {


    protected SimulatedAnnealingOptimizer minimizer = null;

    public TestSimulatedAnnealingOptimizer(String name) {
	super(name);
    }

   
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	// any class variables you declare should be initialized here. This is called before each test
	minimizer = getMinimizer();
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
	minimizer = null;
    }
    

    /**
     * Tests that candidate states get generated (ie the states are different) and that the energy 
     * agrees with the stateEnergy() method. 
     *
     */
    public void testGenerateCandidateState() {
	double energy = minimizer.stateEnergy();

	double candidateEnergy = minimizer.generateCandidateState();
	
	assertTrue(energy != candidateEnergy);

	minimizer.acceptCandidateState();

	energy = minimizer.stateEnergy();

	assertEquals(energy, candidateEnergy, 1E-10);

	for (int i = 0; i < 1000; i++) {
	    candidateEnergy = minimizer.generateCandidateState();
	    minimizer.acceptCandidateState();
	}
	
	energy = minimizer.stateEnergy();
	    
	assertEquals(energy, candidateEnergy, 1E-10);

    }


    /**
     * Make sure that the best state gets saved and then stays saved.
     *
     */
    public void testBestStateSave() {

	String s = minimizer.lowestEnergyState();
	
	minimizer.generateCandidateState();

	minimizer.acceptCandidateState();
	
	assertEquals(s, minimizer.lowestEnergyState());
	
    }


    /**
     * Check that temperature is selected correctly.
     *
     */
    public void testCalibration() {
	
	double rate = 0.6;

	double temp = minimizer.calibratedTemp(rate);

	double meanDeltaE = 0.0;

	int runs = 1000;

	double stateEnergy = minimizer.stateEnergy();

	minimizer.setTemperature(0.0);

	Random ran = new Random(248);

	for (int i = 0; i < runs; i++) {
	    double newEnergy = minimizer.generateCandidateState();

	    double deltaE = (newEnergy - stateEnergy);

	    meanDeltaE += deltaE / runs;
	    
	    if (deltaE < 0.0) {
		minimizer.acceptCandidateState();
	    }
	    
	}

	assertEquals(-1.0 * meanDeltaE / Math.log(rate), temp, temp / 20.0);
	

    }


    protected abstract SimulatedAnnealingOptimizer getMinimizer();

    
}

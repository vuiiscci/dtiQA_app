package misc;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;
import optimizers.*;

import java.util.Random;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>OrderedAcqWeightedMinimizer.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestOrderedAcqWeightedMinimizer.java,v 1.2 2006/08/15 20:53:54 ucacpco Exp $
 * @author  Philip Cook
 * @see misc.DT
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestOrderedAcqWeightedMinimizer extends TestSimulatedAnnealingOptimizer {

    public TestOrderedAcqWeightedMinimizer(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    public static Test suite() {
	return new TestSuite(TestOrderedAcqWeightedMinimizer.class);
    }

    protected SimulatedAnnealingOptimizer getMinimizer() {
	return new OrderedAcqWeightedMinimizer(60, 1E-7, 1000, "", new Random(1394));
    }


    public void testSwapPoints() {
	OrderedAcqWeightedMinimizer m = (OrderedAcqWeightedMinimizer)getMinimizer();

	m.generateCandidateState();

	m.acceptCandidateState();
	
	double energy = m.stateEnergy();

	// rebuild energy matrix and make sure that it is the same
	m.buildEnergyMatrix();

	assertEquals(energy, m.stateEnergy(), 1E-10);

	double[] points2 = new double[3];

	System.arraycopy(m.points[2], 0, points2, 0, 3);

	double[] points10 = new double[3];

	System.arraycopy(m.points[10], 0, points10, 0, 3);

	// now check points get swapped
	m.swapPoints(2, 10);

	assertEquals(points2[0], m.points[10][0], 1E-10);
	assertEquals(points2[1], m.points[10][1], 1E-10);
	assertEquals(points2[2], m.points[10][2], 1E-10);

	assertEquals(points10[0], m.points[2][0], 1E-10);
	assertEquals(points10[1], m.points[2][1], 1E-10);
	assertEquals(points10[2], m.points[2][2], 1E-10);

	

    }

    
}

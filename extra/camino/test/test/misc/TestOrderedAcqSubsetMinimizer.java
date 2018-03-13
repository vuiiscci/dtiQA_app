package misc;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;
import optimizers.*;

import java.util.Random;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>OrderedAcqSubsetMinimizer.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestOrderedAcqSubsetMinimizer.java,v 1.1 2006/07/26 21:28:19 ucacpco Exp $
 * @author  Philip Cook
 * @see misc.DT
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestOrderedAcqSubsetMinimizer extends TestSimulatedAnnealingOptimizer {

    public TestOrderedAcqSubsetMinimizer(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    public static Test suite() {
	return new TestSuite(TestOrderedAcqSubsetMinimizer.class);
    }

    protected SimulatedAnnealingOptimizer getMinimizer() {
	return new OrderedAcqSubsetMinimizer(60, new int[] {15, 15, 15, 15}, 1E-7, 1000, "", new Random(1394));
    }


    
}

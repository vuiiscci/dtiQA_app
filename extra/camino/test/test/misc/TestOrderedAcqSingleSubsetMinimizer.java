package misc;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;
import optimizers.*;

import java.util.Random;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>OrderedAcqSingleSubsetMinimizer.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestOrderedAcqSingleSubsetMinimizer.java,v 1.1 2006/08/15 20:53:54 ucacpco Exp $
 * @author  Philip Cook
 * @see misc.DT
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestOrderedAcqSingleSubsetMinimizer extends TestSimulatedAnnealingOptimizer {

    public TestOrderedAcqSingleSubsetMinimizer(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    public static Test suite() {
	return new TestSuite(TestOrderedAcqSingleSubsetMinimizer.class);
    }

    protected SimulatedAnnealingOptimizer getMinimizer() {
	return new OrderedAcqSingleSubsetMinimizer(60, 15, 1E-7, 1000, "", new Random(1394));
    }


    
}

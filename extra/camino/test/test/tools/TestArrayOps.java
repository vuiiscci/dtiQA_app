package tools;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for all the array operations.
 * <BR>
 * </dl>
 *
 * @version $Id: TestArrayOps.java,v 1.2 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestArrayOps extends TestCase {


    /**
     * Test double array.
     */
    private double[] testDoubleArray = {1.01, 1.03, 1.02, 1.01, 1.05, 1.00, 1.09};

    /**
     * Test weights.
     */
    private double[] testWeights = {1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0};


    public TestArrayOps(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestArrayOps.class);
    }


    public void testMin() {
	assertEquals(1.00, ArrayOps.min(testDoubleArray), 1E-6);
    }

    public void testMax() {
	assertEquals(1.09, ArrayOps.max(testDoubleArray), 1E-6);
    }

    public void testSum() {
	assertEquals(7.21, ArrayOps.sum(testDoubleArray), 1E-6);
    }


    public void testMedian() {
	assertEquals(1.02, ArrayOps.median(testDoubleArray), 1E-6);
    }

    public void testMean() {
	assertEquals(1.03, ArrayOps.mean(testDoubleArray), 1E-6);
    }

    public void testVar() {
	assertEquals(0.00096667, ArrayOps.var(testDoubleArray, ArrayOps.mean(testDoubleArray)), 1E-6);
    }

    public void testWeightedMean() {
	assertEquals(1.032, ArrayOps.weightedMean(testDoubleArray, testWeights), 1E-6);
    }

    public void testWeightedVar() {
	assertEquals(0.0011853333333,
		     ArrayOps.weightedVar(testDoubleArray, testWeights, 
					  ArrayOps.weightedMean(testDoubleArray, testWeights)), 1E-6);
    }



}

package numerics;

import junit.framework.*;
import junit.extensions.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>SphericalHarmonics.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestSphericalHarmonics.java,v 1.1 2005/06/15 09:12:29 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see data.SphericalHarmonics
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestSphericalHarmonics extends TestCase {


    public TestSphericalHarmonics(String name) {
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
	return new TestSuite(TestSphericalHarmonics.class);
    }


    public void testEval() {

	assertEquals(SphericalHarmonics.funcsUpTo(0), 1);
	assertEquals(SphericalHarmonics.funcsUpTo(1), 4);
	assertEquals(SphericalHarmonics.funcsUpTo(2), 9);
	assertEquals(SphericalHarmonics.funcsUpTo(3), 16);
	assertEquals(SphericalHarmonics.funcsUpTo(4), 25);
	assertEquals(SphericalHarmonics.funcsUpTo(5), 36);
	assertEquals(SphericalHarmonics.funcsUpTo(6), 49);

	assertEquals(SphericalHarmonics.evenFuncsUpTo(0), 1);
	assertEquals(SphericalHarmonics.evenFuncsUpTo(1), 1);
	assertEquals(SphericalHarmonics.evenFuncsUpTo(2), 6);
	assertEquals(SphericalHarmonics.evenFuncsUpTo(3), 6);
	assertEquals(SphericalHarmonics.evenFuncsUpTo(4), 15);
	assertEquals(SphericalHarmonics.evenFuncsUpTo(5), 15);
	assertEquals(SphericalHarmonics.evenFuncsUpTo(6), 28);

	try {
	    assertEquals(SphericalHarmonics.Y(0, 0, 0.3, 0.5).real(), 0.282095, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(0, 0, 0.3, 0.5).imag(), 0.0, 1.0E-6);

	    assertEquals(SphericalHarmonics.Y(1, -1, 0.3, 0.5).real(), 0.089602, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(1, -1, 0.3, 0.5).imag(), -0.048950, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(1, 0, 0.3, 0.5).real(), 0.46678, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(1, 0, 0.3, 0.5).imag(), 0.0, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(1, 1, 0.3, 0.5).real(), -0.089602, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(1, 1, 0.3, 0.5).imag(), -0.04895, 1.0E-6);

	    assertEquals(SphericalHarmonics.Y(2, -2, 0.3, 0.5).real(), 0.018227, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(2, -2, 0.3, 0.5).imag(), -0.028386, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(2, -1, 0.3, 0.5).real(), 0.191407, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(2, -1, 0.3, 0.5).imag(), -0.104566, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(2, 0, 0.3, 0.5).real(), 0.548152, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(2, 0, 0.3, 0.5).imag(), 0.0, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(2, 1, 0.3, 0.5).real(), -0.191407, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(2, 1, 0.3, 0.5).imag(), -0.104566, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(2, 2, 0.3, 0.5).real(), 0.018227, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(2, 2, 0.3, 0.5).imag(), 0.028386, 1.0E-6);

	    assertEquals(SphericalHarmonics.Y(3, -2, 0.3, 0.5).real(), 0.046069, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(3, -2, 0.3, 0.5).imag(), -0.071749, 1.0E-6);

	    assertEquals(SphericalHarmonics.Y(4, 1, 0.3, 0.5).real(), -0.397194, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(4, 1, 0.3, 0.5).imag(), -0.216988, 1.0E-6);

	    assertEquals(SphericalHarmonics.Y(5, 4, 0.3, 0.5).real(), -0.00445, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(5, 4, 0.3, 0.5).imag(), 0.009724, 1.0E-6);

	    assertEquals(SphericalHarmonics.Y(6, -6, 0.3, 0.5).real(), -3.186E-4, 1.0E-6);
	    assertEquals(SphericalHarmonics.Y(6, -6, 0.3, 0.5).imag(), -4.5408E-5, 1.0E-8);

	} catch(Exception e) {
	    fail("Failed to compute spherical harmonic function.");
	}

    }


    
}

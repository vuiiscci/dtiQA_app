package sphfunc;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for RBF_Sum
 * <BR>
 * </dl>
 *
 * @version $Id: TestRBF_Sum.java,v 1.4 2006/02/01 18:49:07 ucacpco Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestRBF_Sum extends TestCase {


    /**
     * Test two objects with different maximum order
     */
    private TuchRBF_Sum trs1;
    private TuchRBF_Sum trs2;

    public TestRBF_Sum(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

        // Coefficients for test order-four object.
        double[] coeffs1 = {0, 0, 1.0, 0.5, 0.3, 0.3, 0.2, 0.1, 0.1, 0.05, 0.01, -0.05, -0.03, 0.02, -0.02, 0.03, 0.01};

        // Coefficients for order-six object.
        double[] coeffs2 = {0, 0, 1.0, 0.05, 0.005, 0.001, -0.001, 0.002, 0.002, -0.002, 0.001, 0.005, 0.001, 0.002, -0.002, 0.001, 0.002};

        RBF_Sum.setPoints(SphericalPoints.getElecPointSet(15));
        TuchRBF_Sum.setSigma(1.0);

        // Create the objects.
        trs1 = new TuchRBF_Sum(coeffs1);
        trs2 = new TuchRBF_Sum(coeffs2);
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestRBF_Sum.class);
    }


    public void testRadius() {

        // Test getRadius at various different points on the sphere.
	assertEquals(0.16102364, trs1.getRadius(0, 0), 0.00001);
	assertEquals(0.16909367, trs1.getRadius(0.1, 0.2), 0.00001);
	assertEquals(0.185840385, trs1.getRadius(0.5, 1.0), 0.00001);
	assertEquals(0.1646620659, trs1.getRadius(3.2, 1.3), 0.00001);

	assertEquals(0.032595587, trs2.getRadius(0, 0), 0.00001);
	assertEquals(0.04022372598, trs2.getRadius(0.1, 0.2), 0.00001);
	assertEquals(0.04828796, trs2.getRadius(0.5, 1.0), 0.00001);
	assertEquals(0.03344865429, trs2.getRadius(3.2, 1.3), 0.00001);


    }


    public void testStats() {

        double[] trs1stats = trs1.getStats();
        assertEquals(0.1680012566, trs1.moment(1), 0.00001);
        assertEquals(0.029118511, trs1.moment(2), 0.00001);
        assertEquals(-2.367092775230205E-4, trs1.normMoment(3), 0.0000001);
        assertEquals(0.001880063678, trs1.normMoment(4), 0.00000001);
        assertEquals(0.1680012566, trs1.mean(), 0.00001);
        assertEquals(0.1752289078, trs1.anisotropy(), 0.00001);
        assertEquals(-0.0618593131, trs1.skewness(), 0.00001);
        assertEquals(1.1883308877, trs1.kurtosis(), 0.00001);
        assertEquals(0.1680012566, trs1stats[0], 0.00001);
        assertEquals(8.940888370269975E-4, trs1stats[1], 0.00001);
        assertEquals(0.2273049228, trs1stats[2], 0.00001);
        assertEquals(0.10560063, trs1stats[3], 0.00001);

        double[] trs2stats = trs2.getStats();
        assertEquals(0.0711381301886, trs2.moment(1), 0.00001);
        assertEquals(0.0068364293777, trs2.moment(2), 0.00001);
        assertEquals(0.06119496211, trs2.normMoment(3), 0.0000001);
        assertEquals(0.0703306646769, trs2.normMoment(4), 0.00000001);
        assertEquals(0.07113813, trs2.mean(), 0.00001);
        assertEquals(0.509661523, trs2.anisotropy(), 0.00001);
        assertEquals(0.3940686536, trs2.skewness(), 0.00001);
        assertEquals(1.010425581028, trs2.kurtosis(), 0.00001);
        assertEquals(0.07113813, trs2stats[0], 0.00001);
        assertEquals(0.001775795811, trs2stats[1], 0.00001);
        assertEquals(0.16646703095, trs2stats[2], 0.00001);
        assertEquals(0.0158835436, trs2stats[3], 0.00001);

    }


    public void testBasisFunction() {

        double[][] sps = SphericalPoints.getElecPointSet(15);
        TuchRBF t = new TuchRBF(sps[4], 1.0);
        assertEquals(t.getRadius(0.3, 0.4), trs1.basisFunction(4).getRadius(0.3, 0.4), 1E-5);

    }
}

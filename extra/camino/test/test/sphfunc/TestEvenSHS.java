package sphfunc;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for EvenSHS
 * <BR>
 * </dl>
 *
 * @version $Id: TestEvenSHS.java,v 1.7 2006/07/24 14:45:35 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestEvenSHS extends TestCase {


    /**
     * Test two objects with different maximum order
     */
    private EvenSHS e4;
    private EvenSHS e6;

    public TestEvenSHS(String name) {
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
        double[] coeffs4 = {0, 0, 1.0, 0.5, 0.3, 0.3, 0.2, 0.1, 0.1, 0.05, 0.01, -0.05, -0.03, 0.02, -0.02, 0.03, 0.01};

        // Coefficients for order-six object.
        double[] coeffs6 = {0, 0, 1.0, 0.5, 0.3, 0.3, 0.2, 0.1, 0.1, 0.05, 0.01, -0.05, -0.03, 0.02, -0.02, 0.03, 0.01, 0.005, 0.001, -0.001, 0.002, 0.002, -0.002, 0.001, 0.005, 0.001, 0.002, -0.002, 0.001, 0.002};

        // Create the objects.
        e4 = new EvenSHS(coeffs4, 4);
        e6 = new EvenSHS(coeffs6, 6);
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestEvenSHS.class);
    }


    public void testRadius() {

        // Test getRadius at various different points on the sphere.
	assertEquals(0.68211479, e4.getRadius(0, 0), 0.00001);
	assertEquals(0.619564088, e4.getRadius(0.1, 0.2), 0.00001);
	assertEquals(0.5418, e4.getRadius(0.5, 1.0), 0.00001);
	assertEquals(0.66127678, e4.getRadius(3.2, 1.3), 0.00001);

	assertEquals(0.68720033, e6.getRadius(0, 0), 0.00001);
	assertEquals(0.62352333, e6.getRadius(0.1, 0.2), 0.00001);
	assertEquals(0.53484415, e6.getRadius(0.5, 1.0), 0.00001);
	assertEquals(0.6665555993, e6.getRadius(3.2, 1.3), 0.00001);


    }


    public void testStats() {

        double[] e4stats = e4.getStats();

        assertEquals(0.28209479177, e4.moment(1), 0.00001);
        assertEquals(0.1381215487, e4.moment(2), 0.00001);
        assertEquals(0.0777415385, e4.normMoment(3), 0.00001);
        assertEquals(0.15445904, e4.normMoment(4), 0.00001);
        assertEquals(0.28209479177, e4.mean(), 0.00001);
        assertEquals(0.65102256, e4.anisotropy(), 0.00001);
        assertEquals(0.426793415, e4.skewness(), 0.00001);
        assertEquals(0.96295789, e4.kurtosis(), 0.00001);
        assertEquals(0.28209479177, e4stats[0], 0.00001);
        assertEquals(0.05854407719, e4stats[1], 0.00001);
        assertEquals(0.82909753086, e4stats[2], 0.00001);
        assertEquals(-0.09995229345, e4stats[3], 0.00001);

        double[] e6stats = e6.getStats();
        assertEquals(0.2820988318, e6.moment(1), 0.00001);
        assertEquals(0.138133576, e6.moment(2), 0.00001);
        assertEquals(0.0781308575, e6.normMoment(3), 0.00001);
        assertEquals(0.154795919, e6.normMoment(4), 0.00001);
        assertEquals(0.2820988318, e6.mean(), 0.00001);
        assertEquals(0.651056468, e6.anisotropy(), 0.00001);
        assertEquals(0.42750467, e6.skewness(), 0.00001);
        assertEquals(0.963432337, e6.kurtosis(), 0.00001);
        assertEquals(0.2820988318, e6stats[0], 0.00001);
        assertEquals(0.0585538251, e6stats[1], 0.00001);
        assertEquals(0.83187064879, e6stats[2], 0.00001);
        assertEquals(-0.0986043003597, e6stats[3], 0.00001);

    }


    public void testBasisFunction() {

        RealSH rsh = new RealSH(2, 1);
        ImagSH ish = new ImagSH(2, 1);
        assertEquals(rsh.getRadius(0.3, 0.4)*2.0, e6.basisFunction(2).getRadius(0.3, 0.4), 1E-5);
        assertEquals(rsh.getRadius(0.1, 0.5)*2.0, e6.basisFunction(2).getRadius(0.1, 0.5), 1E-5);
        assertEquals(ish.getRadius(0.3, 0.4)*(-2.0), e6.basisFunction(3).getRadius(0.3, 0.4), 1E-5);
        assertEquals(ish.getRadius(0.1, 0.5)*(-2.0), e6.basisFunction(3).getRadius(0.1, 0.5), 1E-5);


        rsh = new RealSH(4, 0);
        ish = new ImagSH(4, 1);
        assertEquals(rsh.getRadius(0.3, 0.4), e6.basisFunction(6).getRadius(0.3, 0.4), 1E-5);
        assertEquals(rsh.getRadius(0.1, 0.5), e6.basisFunction(6).getRadius(0.1, 0.5), 1E-5);
        assertEquals(ish.getRadius(0.3, 0.4)*(-2.0), e6.basisFunction(8).getRadius(0.3, 0.4), 1E-5);
        assertEquals(ish.getRadius(0.1, 0.5)*(-2.0), e6.basisFunction(8).getRadius(0.1, 0.5), 1E-5);

        rsh = new RealSH(6, 5);
        ish = new ImagSH(6, 6);
        assertEquals(rsh.getRadius(0.3, 0.4)*2.0, e6.basisFunction(24).getRadius(0.3, 0.4), 1E-5);
        assertEquals(rsh.getRadius(0.1, 0.5)*2.0, e6.basisFunction(24).getRadius(0.1, 0.5), 1E-5);
        assertEquals(ish.getRadius(0.3, 0.4)*(-2.0), e6.basisFunction(27).getRadius(0.3, 0.4), 1E-5);
        assertEquals(ish.getRadius(0.1, 0.5)*(-2.0), e6.basisFunction(27).getRadius(0.1, 0.5), 1E-5);
    }

}

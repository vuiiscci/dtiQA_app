package sphfunc;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for MaxEntProfile
 * <BR>
 * </dl>
 *
 * @version $Id: TestMaxEntProfile.java,v 1.2 2005/10/28 14:19:44 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestMaxEntProfile extends TestCase {


    /**
     * Test two objects with different maximum order
     */
    private MaxEntProfile mep1;
    private MaxEntProfile mep2;

    public TestMaxEntProfile(String name) {
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
        double[] coeffs1 = {0, 0, 1.0, 0.5, 0.3, 0.3, 0.2, 0.1, 0.1, 0.05, 0.01, -0.05, -0.03, 0.02, -0.02, 0.03, 0.01, 0.0};

        // Coefficients for order-six object.
        double[] coeffs2 = {0, 0, 1.0, 0.05, 0.005, 0.001, -0.001, 0.002, 0.002, -0.002, 0.001, 0.005, 0.001, 0.002, -0.002, 0.001, 0.002, 0.0};

        double[] kernelParams = {0.0, 1.0};

        MaxEntProfile.setReconDirs(SphericalPoints.getElecPointSet(15));

        // Create the objects.
        mep1 = new MaxEntProfile(coeffs1, kernelParams);
        mep2 = new MaxEntProfile(coeffs2, kernelParams);
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestMaxEntProfile.class);
    }


    public void testRadius() {

        // Test getRadius at various different points on the sphere.
	assertEquals(8.113887156, mep1.getRadius(0, 0), 0.00001);
	assertEquals(7.96516277, mep1.getRadius(0.1, 0.2), 0.00001);
	assertEquals(7.771591518, mep1.getRadius(0.5, 1.0), 0.00001);
	assertEquals(8.021361473, mep1.getRadius(3.2, 1.3), 0.00001);

	assertEquals(2.88571664, mep2.getRadius(0, 0), 0.00001);
	assertEquals(2.88028099, mep2.getRadius(0.1, 0.2), 0.00001);
	assertEquals(2.87510045765, mep2.getRadius(0.5, 1.0), 0.00001);
	assertEquals(2.8853679, mep2.getRadius(3.2, 1.3), 0.00001);


    }


    public void testStats() {

        double[] mep1stats = mep1.getStats();
        assertEquals(8.4830378867, mep1.moment(1), 0.00001);
        assertEquals(72.39524886, mep1.moment(2), 0.00001);
        assertEquals(3.974646756025354E-4, mep1.normMoment(3), 0.0000001);
        assertEquals(9.645836659468326E-5, mep1.normMoment(4), 0.00000001);
        assertEquals(8.4830378867, mep1.mean(), 0.00001);
        assertEquals(0.07736559454, mep1.anisotropy(), 0.00001);
        assertEquals(0.073524629665, mep1.skewness(), 0.00001);
        assertEquals(1.2809645424, mep1.kurtosis(), 0.00001);
        assertEquals(8.4830378867, mep1stats[0], 0.00001);
        assertEquals(0.4333170722, mep1stats[1], 0.00001);
        assertEquals(10.223160384, mep1stats[2], 0.00001);
        assertEquals(7.639042068, mep1stats[3], 0.00001);

        double[] mep2stats = mep2.getStats();
        assertEquals(2.85790141, mep2.moment(1), 0.00001);
        assertEquals(8.168420998, mep2.moment(2), 0.00001);
        assertEquals(-3.250605785188354E-7, mep2.normMoment(3), 0.000000001);
        assertEquals(1.7695500063320353E-8, mep2.normMoment(4), 0.0000000001);
        assertEquals(2.85790141, mep2.mean(), 0.00001);
        assertEquals(0.01002248266, mep2.anisotropy(), 0.00001);
        assertEquals(-0.006875771486, mep2.skewness(), 0.00001);
        assertEquals(1.150774941, mep2.kurtosis(), 0.00001);
        assertEquals(2.8579014117, mep2stats[0], 0.00001);
        assertEquals(8.205191857157956E-4, mep2stats[1], 0.0000001);
        assertEquals(2.89756025, mep2stats[2], 0.00001);
        assertEquals(2.8039257928, mep2stats[3], 0.00001);

    }

}

package sphfunc;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for RealSH
 * <BR>
 * </dl>
 *
 * @version $Id: TestRealSH.java,v 1.2 2006/06/07 10:20:20 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestRealSH extends TestCase {


    /**
     * Test two objects with different maximum order
     */
    private RealSH rsh1;
    private RealSH rsh2;
    private RealSH rsh3;

    private RealSH rsh4;
    private RealSH rsh4a;
    private RealSH rsh5;
    private RealSH rsh5a;

    public TestRealSH(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

        // Create the objects.
        rsh1 = new RealSH(2,1);
        rsh2 = new RealSH(4,-3);
        rsh3 = new RealSH(4,-3,3);

        rsh4 = new RealSH(4,0);
        rsh4a = new RealSH(4,0,3);
        rsh5 = new RealSH(4,2);
        rsh5a = new RealSH(4,2,3);

    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestRealSH.class);
    }


    public void testRadius() {

        // Test getRadius at various different points on the sphere.
	assertEquals(0.0, rsh1.getRadius(0, 0), 0.0001);
	assertEquals(-0.0752111297, rsh1.getRadius(0.1, 0.2), 0.0001);
	assertEquals(-0.17561906897, rsh1.getRadius(0.5, 1.0), 0.0001);
	assertEquals(0.01204278418, rsh1.getRadius(3.2, 1.3), 0.0001);

	assertEquals(0.0, rsh2.getRadius(0, 0), 0.0001);
	assertEquals(0.00102275978, rsh2.getRadius(0.1, 0.2), 0.0001);
	assertEquals(-0.119832256, rsh2.getRadius(0.5, 1.0), 0.0001);
	assertEquals(1.8042920584531209E-4, rsh2.getRadius(3.2, 1.3), 0.0000001);

	assertEquals(0.0, rsh3.getRadius(0, 0), 0.0001);
	assertEquals(0.00102275978*3.0, rsh3.getRadius(0.1, 0.2), 0.0001);
	assertEquals(-0.119832256*3.0, rsh3.getRadius(0.5, 1.0), 0.0001);
	assertEquals(1.8042920584531209E-4*3.0, rsh3.getRadius(3.2, 1.3), 0.0000001);


    }


    public void testGreatCircleIntegral() {

        final double Rt2= 1.414213562;
        double[] u= new double[] {1.0/Rt2, 1.0/Rt2, 0.0};

        assertEquals(0.0, rsh1.greatCircleIntegral(u), 1E-4);
        assertEquals(0.0, rsh2.greatCircleIntegral(u), 1E-4);
        assertEquals(0.0, rsh3.greatCircleIntegral(u), 1E-4);

        assertEquals(0.747754, rsh4.greatCircleIntegral(u), 1E-4);
        assertEquals(0.747754*3.0, rsh4a.greatCircleIntegral(u), 1E-4);
        assertEquals(0.0, rsh5.greatCircleIntegral(u), 1E-4);
        assertEquals(0.0, rsh5a.greatCircleIntegral(u), 1E-4);
    }


    public void testStats() {

        double[] rsh1stats = rsh1.getStats();
        double[] rsh2stats = rsh2.getStats();
        double[] rsh3stats = rsh3.getStats();

        assertEquals(0.0, rsh1.moment(1), 0.0001);
        assertEquals(0.03978873577, rsh1.moment(2), 0.0001); 
        assertEquals(1.015801887, rsh1.normMoment(3), 0.0001); // large delta because of Java 1.7 vs 1.6 differences
        assertEquals(1.0, rsh1.normMoment(4), 0.00000001);
        assertEquals(0.0, rsh1.mean(), 0.0001);
        assertEquals(1.0, rsh1.anisotropy(), 0.0001);
        assertEquals(1.00523979, rsh1.skewness(), 0.0001);
        assertEquals(1.0, rsh1.kurtosis(), 0.0001);
        assertEquals(0.0, rsh1stats[0], 0.0001);
        assertEquals(0.03978873577, rsh1stats[1], 0.0001);
        assertEquals(0.3861942014458, rsh1stats[2], 0.0001);
        assertEquals(-0.3861942014458, rsh1stats[3], 0.0001);

        assertEquals(0.0, rsh2.moment(1), 0.0001);
        assertEquals(0.03984512759, rsh2.moment(2), 0.0001);
        assertEquals(0.991614, rsh2.normMoment(3), 0.0001);
        assertEquals(1.0, rsh2.normMoment(4), 0.0001);
        assertEquals(0.0, rsh2.mean(), 0.0001);
        assertEquals(1.0, rsh2.anisotropy(), 0.0001);
        assertEquals(0.997197, rsh2.skewness(), 0.0001);
        assertEquals(1.0, rsh2.kurtosis(), 0.0001);
        assertEquals(0.0, rsh2stats[0], 0.0001);
        assertEquals(0.03984512759, rsh2stats[1], 0.0001);
        assertEquals(0.40631396184, rsh2stats[2], 0.0001);
        assertEquals(-0.4063139618, rsh2stats[3], 0.0001);

        assertEquals(0.0, rsh3.moment(1), 0.0001);
        assertEquals(0.03984512759*9.0, rsh3.moment(2), 0.0001);
        assertEquals(0.9761014, rsh3.normMoment(3), 0.0001);
        assertEquals(1.0, rsh3.normMoment(4), 0.0001);
        assertEquals(0.0, rsh3.mean(), 0.0001);
        assertEquals(1.0, rsh3.anisotropy(), 0.0001);
        assertEquals(0.991969, rsh3.skewness(), 0.0001);
        assertEquals(1.0, rsh3.kurtosis(), 0.0001);
        assertEquals(0.0, rsh3stats[0], 0.0001);
        assertEquals(0.03984512759*9.0, rsh3stats[1], 0.0001);
        assertEquals(0.40631396184*3.0, rsh3stats[2], 0.0001);
        assertEquals(-0.4063139618*3.0, rsh3stats[3], 0.0001);

    }

}

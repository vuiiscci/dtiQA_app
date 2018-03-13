package sphfunc;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for ImagSH
 * <BR>
 * </dl>
 *
 * @version $Id: TestImagSH.java,v 1.2 2006/06/07 10:20:07 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestImagSH extends TestCase {


    /**
     * Test two objects with different maximum order
     */
    private ImagSH ish1;
    private ImagSH ish2;
    private ImagSH ish3;

    private ImagSH ish4;
    private ImagSH ish4a;
    private ImagSH ish5;
    private ImagSH ish5a;

    public TestImagSH(String name) {
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
        ish1 = new ImagSH(2,1);
        ish2 = new ImagSH(4,-3);
        ish3 = new ImagSH(4,-3,3);

        ish4 = new ImagSH(4,0);
        ish4a = new ImagSH(4,0,3);
        ish5 = new ImagSH(4,2);
        ish5a = new ImagSH(4,2,3);
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestImagSH.class);
    }


    public void testRadius() {

        // Test getRadius at various different points on the sphere.
	assertEquals(0.0, ish1.getRadius(0, 0), 0.0001);
	assertEquals(-0.015246050775, ish1.getRadius(0.1, 0.2), 0.0001);
	assertEquals(-0.2735104946, ish1.getRadius(0.5, 1.0), 0.0001);
	assertEquals(0.0433793423895, ish1.getRadius(3.2, 1.3), 0.0001);

	assertEquals(0.0, ish2.getRadius(0, 0), 0.0001);
	assertEquals(-6.997076128518016E-4, ish2.getRadius(0.1, 0.2), 0.0001);
	assertEquals(-0.017081673857, ish2.getRadius(0.5, 1.0), 0.0001);
	assertEquals(-1.7094307718620608E-4, ish2.getRadius(3.2, 1.3), 0.0001);
	assertEquals(0.0, ish3.getRadius(0, 0), 0.0001);
	assertEquals(-6.997076128518016E-4*3.0, ish3.getRadius(0.1, 0.2), 0.0001);
	assertEquals(-0.017081673857*3.0, ish3.getRadius(0.5, 1.0), 0.0001);
	assertEquals(-1.7094307718620608E-4*3.0, ish3.getRadius(3.2, 1.3), 0.0001);


    }


    public void testGreatCircleIntegral() {

        final double Rt2= 1.414213562;
        double[] u= new double[] {1.0/Rt2, 1.0/Rt2, 0.0};

        assertEquals(1.0508516987E-16, ish1.greatCircleIntegral(u), 1E-20);
        assertEquals(-1.2769307572E-16, ish2.greatCircleIntegral(u), 1E-20);
        assertEquals(-1.2769307572E-16*3, ish3.greatCircleIntegral(u), 1E-20);

        assertEquals(0.0, ish4.greatCircleIntegral(u), 1E-4);
        assertEquals(0.0, ish4a.greatCircleIntegral(u), 1E-4);
        assertEquals(-0.788202, ish5.greatCircleIntegral(u), 1E-4);
        assertEquals(-0.788202*3.0, ish5a.greatCircleIntegral(u), 1E-4);
    }


    public void testStats() {

        double[] ish1stats = ish1.getStats();
        double[] ish2stats = ish2.getStats();
        double[] ish3stats = ish3.getStats();

        assertEquals(0.0, ish1.moment(1), 0.0001);
        assertEquals(0.03978873577, ish1.moment(2), 0.0001);
        assertEquals(-286765.0582534, ish1.normMoment(3), 100.0); // large delta because of Java 1.7 vs 1.6 differences
        assertEquals(1.0, ish1.normMoment(4), 0.000001);
        assertEquals(0.0, ish1.mean(), 0.0001);
        assertEquals(1.0, ish1.anisotropy(), 0.0001);
        assertEquals(-65.94712, ish1.skewness(), 0.0001);
        assertEquals(1.0, ish1.kurtosis(), 0.0001);
        assertEquals(0.0, ish1stats[0], 0.0001);
        assertEquals(0.0397887357, ish1stats[1], 0.0001);
        assertEquals(0.3861942014458, ish1stats[2], 0.0001);
        assertEquals(-0.3861942014, ish1stats[3], 0.0001);

        assertEquals(0.0, ish2.moment(1), 0.0001);
        assertEquals(0.0396683424, ish2.moment(2), 0.0001);
        assertEquals(-393475.3156, ish2.normMoment(3), 0.0001);
        assertEquals(1.0, ish2.normMoment(4), 0.00000001);
        assertEquals(0.0, ish2.mean(), 0.0001);
        assertEquals(1.0, ish2.anisotropy(), 0.0001);
        assertEquals(-73.27781265, ish2.skewness(), 0.0001);
        assertEquals(1.0, ish2.kurtosis(), 0.0001);
        assertEquals(0.0, ish2stats[0], 0.0001);
        assertEquals(0.0396683424, ish2stats[1], 0.0001);
        assertEquals(0.406476197026, ish2stats[2], 0.0001);
        assertEquals(-0.40647619702, ish2stats[3], 0.0001);

        assertEquals(0.0, ish3.moment(1), 0.0001);
        assertEquals(0.0396683424*3.0*3.0, ish3.moment(2), 0.0001);
        assertEquals(-436257.1115, ish3.normMoment(3), 0.0001);
        assertEquals(1.0, ish3.normMoment(4), 0.00000001);
        assertEquals(0.0, ish3.mean(), 0.0001);
        assertEquals(1.0, ish3.anisotropy(), 0.0001);
        assertEquals(-75.84277422, ish3.skewness(), 0.0001);
        assertEquals(1.0, ish3.kurtosis(), 0.0001);
        assertEquals(0.0, ish3stats[0], 0.0001);
        assertEquals(0.0396683424*3.0*3.0, ish3stats[1], 0.0001);
        assertEquals(0.406476197026*3.0, ish3stats[2], 0.0001);
        assertEquals(-0.40647619702*3.0, ish3stats[3], 0.0001);

    }

}

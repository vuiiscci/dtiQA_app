package sphfunc;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for TuchRBF
 * <BR>
 * </dl>
 *
 * @version $Id: TestTuchRBF.java,v 1.1 2005/12/14 10:11:20 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTuchRBF extends TestCase {


    /**
     * Test two objects with different maximum order
     */
    private TuchRBF trbf1;
    private TuchRBF trbf2;

    public TestTuchRBF(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

        double[] p1 = {1, 0, 0};
        double[] p2 = {Math.cos(Math.PI/3)*Math.sin(Math.PI/7), Math.sin(Math.PI/3)*Math.sin(Math.PI/7), Math.cos(Math.PI/7)};

        double sigma1 = 1.0;
        double sigma2 = 3.0;

        // Create the objects.
        trbf1 = new TuchRBF(p1, sigma1);
        trbf2 = new TuchRBF(p2, sigma2);
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestTuchRBF.class);
    }


    public void testRadius() {

        // Test getRadius at various different points on the sphere.
	assertEquals(0.08480497247, trbf1.getRadius(0, 0), 0.00001);
	assertEquals(0.114277169, trbf1.getRadius(0.1, 0.2), 0.00001);
	assertEquals(0.1803441238, trbf1.getRadius(0.5, 1.0), 0.00001);
	assertEquals(0.08904735797898, trbf1.getRadius(3.2, 1.3), 0.00001);

	assertEquals(0.97786852, trbf2.getRadius(0, 0), 0.00001);
	assertEquals(0.9832952797, trbf2.getRadius(0.1, 0.2), 0.00001);
	assertEquals(0.9996572745, trbf2.getRadius(0.5, 1.0), 0.00001);
	assertEquals(0.9830281686, trbf2.getRadius(3.2, 1.3), 0.00001);


    }


    public void testStats() {

        double[] trbf1stats = trbf1.getStats();
        double[] trbf2stats = trbf2.getStats();

        assertEquals(0.4023256838, trbf1.moment(1), 0.00001);
        assertEquals(0.2291218512, trbf1.moment(2), 0.00001);
        assertEquals(0.0717784573613, trbf1.normMoment(3), 0.0000001);
        assertEquals(0.084507377328, trbf1.normMoment(4), 0.00000001);
        assertEquals(0.4023256838, trbf1.mean(), 0.00001);
        assertEquals(0.5417911723, trbf1.anisotropy(), 0.00001);
        assertEquals(0.4155896343, trbf1.skewness(), 0.00001);
        assertEquals(0.995157295866, trbf1.kurtosis(), 0.00001);
        assertEquals(0.4023256838, trbf1stats[0], 0.00001);
        assertEquals(0.067255895364, trbf1stats[1], 0.00001);
        assertEquals(0.9995143759, trbf1stats[2], 0.00001);
        assertEquals(0.08480497247, trbf1stats[3], 0.00001);

        assertEquals(0.8836040379, trbf2.moment(1), 0.00001);
        assertEquals(0.78551464, trbf2.moment(2), 0.00001);
        assertEquals(-2.7605081603800954E-5, trbf2.normMoment(3), 0.0000001);
        assertEquals(6.490984900523166E-5, trbf2.normMoment(4), 0.00000001);
        assertEquals(0.8836040379, trbf2.mean(), 0.00001);
        assertEquals(0.077832348, trbf2.anisotropy(), 0.00001);
        assertEquals(-0.0302224507, trbf2.skewness(), 0.00001);
        assertEquals(1.153234205857, trbf2.kurtosis(), 0.00001);
        assertEquals(0.8836040379, trbf2stats[0], 0.00001);
        assertEquals(0.0047585491, trbf2stats[1], 0.00001);
        assertEquals(0.9999914938, trbf2stats[2], 0.00001);
        assertEquals(0.760310887, trbf2stats[3], 0.00001);

    }

}

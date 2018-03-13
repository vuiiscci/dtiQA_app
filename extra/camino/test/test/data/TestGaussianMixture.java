package data;

import junit.framework.*;
import junit.extensions.*;
import misc.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>GaussianMixture.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestGaussianMixture.java,v 1.1 2005/06/15 09:12:29 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see numerics.GaussianMixture
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestGaussianMixture extends TestCase {


    public TestGaussianMixture(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	// any class variables you declare should be initialized here. This is called before each test
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestGaussianMixture.class);
    }


    public void testOneComponent() {

	double dxx = 1.7E-9;
	double dxy = 0.0;
	double dxz = 0.0;
	double dyy = 2.0E-10;
	double dyz = 0.0;
	double dzz = 2.0E-10;

	DT dt1 = new DT(dxx, dxy, dxz, dyy, dyz, dzz);

	DT[] dts = new DT[1];
	dts[0] = dt1;
	double[] mixs = {1.0};
	double tau = 0.04;

	GaussianMixture gm = new GaussianMixture(dts, mixs);

	double[] q = {2.0E5, 0.0, 0.0};
	assertEquals(gm.ftAt(q, tau), 0.06587475442640295, 0.000001);
	q[0] = 0.0;
	q[1] = 2.0E5;
	assertEquals(gm.ftAt(q, tau), 0.7261490370736909, 0.000001);
	q[1] = 0.0;
	q[2] = 2.0E5;
	assertEquals(gm.ftAt(q, tau), 0.7261490370736909, 0.000001);
	q[0] = 2.0E5/1.73205;
	q[1] = 2.0E5/1.73205;
	q[2] = 2.0E5/1.73205;
	assertEquals(gm.ftAt(q, tau), 0.32627945385628226, 0.000001);

	// Use some typical displacements and find the pdf value.
	double[] x = {2.0E-5, 0, 0};
	assertEquals(gm.at(x, tau), 7.819381543762102E13, 1.0E7);
	x[0] = 0.0;
	x[1] = 2.0E-5;
	assertEquals(gm.at(x, tau), 1.2681182036628723E9, 1.0E3);
	x[1] = 0.0;
	x[2] = 2.0E-5;
	assertEquals(gm.at(x, tau), 1.2681182036628723E9, 1.0E3);
	x[0] = 2.0E-5/1.73205;
	x[1] = 2.0E-5/1.73205;
	x[2] = 2.0E-5/1.73205;
	assertEquals(gm.at(x, tau), 5.0098768778581665E10, 1.0E4);

    }


    public void testTwoComponent() {

	double dxx = 1.7E-9;
	double dxy = 0.0;
	double dxz = 0.0;
	double dyy = 2.0E-10;
	double dyz = 0.0;
	double dzz = 2.0E-10;

	DT dt1 = new DT(dxx, dxy, dxz, dyy, dyz, dzz);

	dxx = 7.0E-10;
	dxy = 5.0E-11;
	dxz = -2.0E-11;
	dyy = 2.0E-10;
	dyz = 4.5E-11;
	dzz = 2.0E-10;

	DT dt2 = new DT(dxx, dxy, dxz, dyy, dyz, dzz);

	DT[] dts = new DT[2];
	dts[0] = dt1;
	dts[1] = dt2;
	double[] mixs = {0.5, 0.5};
	double tau = 0.04;

	GaussianMixture gm = new GaussianMixture(dts, mixs);

	double[] q = {2.0E5, 0.0, 0.0};
	assertEquals(gm.ftAt(q, tau), 0.19607727452472123, 0.000001);
	q[0] = 0.0;
	q[1] = 2.0E5;
	assertEquals(gm.ftAt(q, tau), 0.7261490370736909, 0.000001);
	q[1] = 0.0;
	q[2] = 2.0E5;
	assertEquals(gm.ftAt(q, tau), 0.7261490370736909, 0.000001);
	q[0] = 2.0E5/1.73205;
	q[1] = 2.0E5/1.73205;
	q[2] = 2.0E5/1.73205;
	assertEquals(gm.ftAt(q, tau), 0.41984812685718353, 0.000001);

	// Use some typical displacements and find the pdf value.
	double[] x = {2.0E-5, 0, 0};
	assertEquals(gm.at(x, tau), 4.616259588044246E13, 1.0E7);
	x[0] = 0.0;
	x[1] = 2.0E-5;
	assertEquals(gm.at(x, tau), 1.024052575220468E9, 1.0E3);
	x[1] = 0.0;
	x[2] = 2.0E-5;
	assertEquals(gm.at(x, tau), 1.11164064887521E9, 1.0E3);
	x[0] = 2.0E-5/1.73205;
	x[1] = 2.0E-5/1.73205;
	x[2] = 2.0E-5/1.73205;
	assertEquals(gm.at(x, tau), 1.445790372353563E11, 1.0E4);

    }



    
}

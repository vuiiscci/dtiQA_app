package inverters;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>EvenSphHarmFitter.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestEvenSphHarmFitter.java,v 1.1 2005/06/15 09:12:29 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestEvenSphHarmFitter extends TestCase {


    // Test function used for data synthesis.
    private GaussianMixture gm2;

    // Set of imaging parameters.
    private DW_Scheme ip;

    // Source of synthetic data for testing
    private DataSynthesizer ds1, ds2, ds3;


    public TestEvenSphHarmFitter(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	// Set up the test functions.
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

        gm2 = new GaussianMixture(dts, mixs);


	StandardTestFunctions.reset();
	StandardTestFunctions.setMix3(0.6);
	RealMatrix rot = Rotations.randomRotMat(0);
	StandardTestFunctions.setTransformation(rot);
	ModelPDF p3 = StandardTestFunctions.getFunction(3);
	ModelPDF p4 = StandardTestFunctions.getFunction(4);

	ip = DW_Scheme.readScheme("bmx7.scheme");

	// Test data source for two tensor fitter tests.
	ds1 = new DataSynthesizer(p3, ip, 100000);

	// Test data source for one tensor fitter tests.
	ds2 = new DataSynthesizer(gm2, ip, 16);
	
	// Test data source for three tensor fitter tests.
	ds3 = new DataSynthesizer(p4, ip, 100000);

    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestEvenSphHarmFitter.class);
    }


    public void testFittingOrder8() {

	EvenSphHarmFitter eshf = new EvenSphHarmFitter(ip, 8);

	try {
	    double[] fitted = eshf.fit(ds1.nextVoxel());
	    assertEquals(fitted.length, 47);

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.83731E-6, fitted[1], 1.0E-12);
	    assertEquals(2.11342E-9, fitted[2], 1.0E-14);
	    assertEquals(2.27277E-10, fitted[3], 1.0E-14);
	    assertEquals(-1.84868E-13, fitted[44], 1.0E-14);
	    assertEquals(1.89379E-13, fitted[45], 1.0E-14);
	    assertEquals(1.08121E-13, fitted[46], 1.0E-14);

	    fitted = eshf.fit(ds2.nextVoxel());
	    assertEquals(fitted.length, 47);

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(0.036071, fitted[1], 0.000001);
	    assertEquals(1.81271E-9, fitted[2], 1.0E-14);
	    assertEquals(-4.23579E-10, fitted[3], 1.0E-14);
	    assertEquals(1.60159E-11, fitted[44], 1.0E-14);
	    assertEquals(7.20163E-11, fitted[45], 1.0E-14);
	    assertEquals(-4.28744E-11, fitted[46], 1.0E-14);

	    fitted = eshf.fit(ds3.nextVoxel());
	    assertEquals(fitted.length, 47);

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.838183E-6, fitted[1], 1.0E-12);
	    assertEquals(2.02392E-9, fitted[2], 1.0E-14);
	    assertEquals(4.15213E-14, fitted[3], 1.0E-16);
	    assertEquals(2.24024E-13, fitted[44], 1.0E-16);
	    assertEquals(4.93959E-13, fitted[45], 1.0E-16);
	    assertEquals(-1.0088E-14, fitted[46], 1.0E-16);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


 
}

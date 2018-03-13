package data;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import imaging.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>DataSynthesizer.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestDataSynthesizer.java,v 1.1 2005/06/15 09:12:29 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see data.DataSynthesizer
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestDataSynthesizer extends TestCase {


    // Test functions used for data synthesis.
    private GaussianMixture gm1, gm2;

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestDataSynthesizer(String name) {
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

	DT[] dts = new DT[1];
	dts[0] = dt1;
	double[] mixs = {1.0};
	double tau = 0.04;

	gm1 = new GaussianMixture(dts, mixs);

	dxx = 7.0E-10;
	dxy = 5.0E-11;
	dxz = -2.0E-11;
	dyy = 2.0E-10;
	dyz = 4.5E-11;
	dzz = 2.0E-10;

	DT dt2 = new DT(dxx, dxy, dxz, dyy, dyz, dzz);

	dts = new DT[2];
	dts[0] = dt1;
	dts[1] = dt2;
	double[] mixs2 = {0.5, 0.5};

        gm2 = new GaussianMixture(dts, mixs2);


	ip = DW_Scheme.readScheme("bmx7.scheme");

    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestDataSynthesizer.class);
    }


    public void testDataSynth() {

	DataSynthesizer ds1 = new DataSynthesizer(gm1, ip, 100000);
	
	// Check parts of the first voxels worth of data.
	try {
	    double[] data = ds1.nextVoxel();
	    assertEquals(data[0], 1.0, 0.00001);
	    assertEquals(data[6], 0.216314, 0.000001);

	    double[] norm = ip.normalizeData(data);
	    assertEquals(norm[0], 0.216314, 0.000001);
	} catch(Exception e) {
	    fail("Failed to get next voxel from ds1.");
	}


	DataSynthesizer ds2 = new DataSynthesizer(gm2, ip, 16);
	
	// Check parts of the first voxels worth of data.
	try {
	    double[] data = ds2.nextVoxel();
	    assertEquals(data[0], 1.051669, 0.000001);
	    assertEquals(data[6], 0.294656, 0.000001);

	    double[] norm = ip.normalizeData(data);
	    assertEquals(norm[0], 0.284217, 0.000001);
	} catch(Exception e) {
	    fail("Failed to get next voxel from ds2.");
	}
    }


    public void testSynthBootstrapper() {

	SyntheticBootstrapper sbs = new SyntheticBootstrapper(gm1, ip, 100000, 5, 1, 1234);
	
	// Check parts of the first voxels worth of data.
	try {
	    double[] data = sbs.nextVoxel();
	    assertEquals(1.0, data[0], 0.0001);
	    assertEquals(0.2163, data[6], 0.0001);

	    double[] norm = ip.normalizeData(data);
	    assertEquals(0.2163, norm[0], 0.0001);
	} catch(Exception e) {
	    fail("Failed to get next voxel from sbs.");
	}

    }



    
}

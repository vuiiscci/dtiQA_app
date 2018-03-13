package inverters;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import fitters.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for all the diffusion inverters.
 * <BR>
 * </dl>
 *
 * @version $Id: TestInverters.java,v 1.4 2005/10/20 08:21:41 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestInverters extends TestCase {


    // Test function used for data synthesis.
    private GaussianMixture gm2;

    // Set of imaging parameters.
    private DW_Scheme ip;

    // Source of synthetic data for testing
    private DataSynthesizer ds1, ds2, ds3;


    public TestInverters(String name) {
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
	return new TestSuite(TestInverters.class);
    }


 
    public void testLinearDT_Inversion() {

	DiffusionInversion difInv = new LinearDT_Inversion(ip);

	try {

            double[] data = ds2.nextVoxel();

	    double[] fitted = difInv.invert(data);
	    assertEquals(fitted[0], 0.0, 0.000001);
	    assertEquals(fitted[1], 0.036071, 0.000001);
	    assertEquals(fitted[2], 1.020834E-9, 1.0E-13);
	    assertEquals(fitted[3], 6.046258E-11, 1.0E-13);

            // check bad data correctly handled

            data[30] = 0.0;

	    fitted = difInv.invert(data);
	    assertEquals(6.0, fitted[0], 0.000001);
	    assertEquals(0.036071, fitted[1], 0.000001);
	    assertEquals(1.021E-9, fitted[2], 1.0E-11);
	    assertEquals(6.05E-11, fitted[3], 1.0E-11);

            

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testWeightedLinearDT_Inversion() {

	WeightedLinearDT_Inversion difInv = new WeightedLinearDT_Inversion(ip);

	try {

	    double[] data = ds2.nextVoxel();

	    double[] fitted = difInv.invert(data);

	    assertEquals(0.0, fitted[0], 0.00001);
	    assertEquals(0.036071, fitted[1], 0.00001);
	    assertEquals(1.0676E-9, fitted[2], 1.0E-12);
	    assertEquals(5.8051E-11, fitted[3], 1.0E-14);

	    // test variance
	    DiffusionInversion olsInv = new LinearDT_Inversion(ip);
	    
	    double[] olsSolution = olsInv.invert(data);

	    RealMatrix[] wlsMats = difInv.computeWeightedFit(data, olsSolution);

	    double sigmaSq = difInv.computeNoiseVariance(wlsMats);
	    
	    assertEquals(0.0039417, sigmaSq, 1E-6);


	} catch(Exception e) {
	    fail("Data source failed.");
	}


    }


    public void testNonLinearDT_Inversion() {

	DiffusionInversion difInv = new NonLinearDT_Inversion(ip);

	try {
	    double[] fitted = difInv.invert(ds2.nextVoxel());
	    assertEquals(fitted[0], 0.0, 0.000001);
	    assertEquals(fitted[1], 0.036071, 0.000001);
	    assertEquals(fitted[2], 1.045218E-9, 1.0E-15);
	    assertEquals(fitted[3], 5.661728E-11, 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }



    public void testNonLinearDT_InversionUnCon() {

	DiffusionInversion difInv = new NonLinearDT_Inversion(ip, ModelIndex.NLDT);

	try {
	    double[] fitted = difInv.invert(ds2.nextVoxel());
	    assertEquals(fitted[0], 0.0, 0.000001);
	    assertEquals(fitted[1], 0.036071, 0.000001);
	    assertEquals(fitted[2], 1.045217E-9, 1.0E-15);
	    assertEquals(fitted[3], 5.6618006E-11, 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testNonLinearDT_InversionPos() {

	DiffusionInversion difInv = new NonLinearDT_Inversion(ip, ModelIndex.NLDT_POS);

	try {
	    double[] fitted = difInv.invert(ds2.nextVoxel());
	    assertEquals(fitted[0], 0.0, 0.000001);
	    assertEquals(fitted[1], 0.036071, 0.000001);
	    assertEquals(fitted[2], 1.045217E-9, 1.0E-15);
	    assertEquals(fitted[3], 5.6618006E-11, 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testRestoreDT_Inversion() {

	DiffusionInversion difInv = new RestoreDT_Inversion(ip, 0.01);

	try {
	    double[] fitted = difInv.invert(ds2.nextVoxel());
	    assertEquals(fitted[0], 1036.0, 0.000001);
	    assertEquals(fitted[1], 0.036071, 0.000001);
	    assertEquals(fitted[2], 1.043643E-9, 1.0E-15);
	    assertEquals(fitted[3], 5.058489E-11, 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testTwoTensorFitter() {

	DiffusionInversion difInv = new TwoTensorInversion(ip, ModelIndex.DTDT, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds1.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.999218E-6, fitted[1], 1.0E-12);
	    assertEquals(2.0, fitted[2], 1.0E-6);
	    assertEquals(0.621453, fitted[3], 1.0E-6);
	    assertEquals(8.4648873E-10, fitted[4], 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


 
    public void testTwoTensorFitterNonLinDT() {

	DiffusionInversion difInv = 
	    new TwoTensorInversion(ip, ModelIndex.DTDT, ModelIndex.NLDT_POS);

	try {
	    double[] fitted = difInv.invert(ds1.nextVoxel());

	    assertEquals(fitted[0], 0.0, 0.000001);
	    assertEquals(fitted[1], 5.999218E-6, 1.0E-10);
	    assertEquals(fitted[2], 2.0, 1.0E-6);
	    assertEquals(fitted[3], 0.631276, 1.0E-6);
	    assertEquals(fitted[4], 8.56280E-10, 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testTwoTensorCholFitter() {

	DiffusionInversion difInv = 
	    new TwoTensorInversion(ip, ModelIndex.POSPOS, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds1.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.999218E-6, fitted[1], 1.0E-10);
	    assertEquals(2.0, fitted[2], 1.0E-6);
	    assertEquals(0.581594, fitted[3], 1.0E-6);
	    assertEquals(8.051024E-10, fitted[4], 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testTwoTensorCholFixMP_Fitter() {

	DiffusionInversion difInv =
	    new TwoTensorInversion(ip, ModelIndex.POSPOS_EQ, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds1.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.999218E-6, fitted[1], 1.0E-10);
	    assertEquals(2.0, fitted[2], 1.0E-6);
	    assertEquals(0.5, fitted[3], 1.0E-6);
	    assertEquals(7.107233E-10, fitted[4], 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testTwoTensorOneAxiSymFitter() {

	DiffusionInversion difInv =
	    new TwoTensorInversion(ip, ModelIndex.POSCYL, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds1.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.999218E-6, fitted[1], 1.0E-10);
	    assertEquals(2.0, fitted[2], 1.0E-6);
	    assertEquals(0.589616, fitted[3], 1.0E-6);
	    assertEquals(8.136713E-10, fitted[4], 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testTwoTensorOneAxiSymFixMP_Fitter() {

	DiffusionInversion difInv =
	    new TwoTensorInversion(ip, ModelIndex.POSCYL_EQ, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds1.nextVoxel());

	    assertEquals(fitted[0], 0.0, 0.000001);
	    assertEquals(fitted[1], 5.999218E-6, 1.0E-10);
	    assertEquals(fitted[2], 2.0, 1.0E-6);
	    assertEquals(fitted[3], 0.5, 1.0E-6);
	    assertEquals(fitted[4], 7.107396E-10, 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testTwoTensorAxiSymFitter() {

	DiffusionInversion difInv = 
	    new TwoTensorInversion(ip, ModelIndex.CYLCYL, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds1.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.999218E-6, fitted[1], 1.0E-10);
	    assertEquals(2.0, fitted[2], 1.0E-6);
	    assertEquals(0.581342, fitted[3], 1.0E-6);
	    assertEquals(8.048490E-10, fitted[4], 1.0E-15);
	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testTwoTensorAxiSymFixMP_Fitter() {

	DiffusionInversion difInv = 
	    new TwoTensorInversion(ip, ModelIndex.CYLCYL_EQ, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds1.nextVoxel());

	    assertEquals(fitted[0], 0.0, 0.000001);
	    assertEquals(fitted[1], 5.999218E-6, 1.0E-10);
	    assertEquals(fitted[2], 2.0, 1.0E-6);
	    assertEquals(fitted[3], 0.5, 1.0E-6);
	    assertEquals(fitted[4], 7.107404E-10, 1.0E-15);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testThreeTensorFitter() {

	DiffusionInversion difInv = 
	    new ThreeTensorInversion(ip, ModelIndex.DTDTDT, ModelIndex.LDT);
      
	try {

	    double[] fitted = difInv.invert(ds3.nextVoxel());

            // require no error
	    assertEquals(0.0, fitted[0], 0.000001);

            // true value is 0.0, but cannot be less than zero
	    assertTrue(fitted[1] < 1E-5);
            assertTrue(fitted[1] > 0.0);
            
            // require right number of DTs
	    assertEquals(3.0, fitted[2], 1.0E-6);
            
            
	    assertEquals(0.36368, fitted[3], 1.0E-5);
	    assertEquals(8.78890E-10, fitted[4], 1.0E-15);

	    assertEquals(0.32431, fitted[10], 1.0E-5);


	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testThreeTensorCholFitter() {

	DiffusionInversion difInv =
	    new ThreeTensorInversion(ip, ModelIndex.POSPOSPOS, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds3.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.96879E-6, fitted[1], 1.0E-10);
	    assertEquals(3.0, fitted[2], 1.0E-6);
	    assertEquals(0.33629, fitted[3], 1.0E-5);
	    assertEquals(8.06209E-10, fitted[4], 1.0E-15);
	    assertEquals(0.34400, fitted[10], 1.0E-5);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


     public void testThreeTensorOneAxiSymFitter() {

	DiffusionInversion difInv = 
	    new ThreeTensorInversion(ip, ModelIndex.POSPOSCYL, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds3.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.968792E-6, fitted[1], 1.0E-10);
	    assertEquals(3.0, fitted[2], 1.0E-6);
	    assertEquals(0.33949, fitted[3], 1.0E-5);
            assertEquals(8.12149E-10, fitted[4], 1.0E-15);
	    assertEquals(0.32801, fitted[10], 1.0E-5);


	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


     public void testThreeTensorTwoAxiSymFitter() {

	DiffusionInversion difInv = 
	    new ThreeTensorInversion(ip, ModelIndex.POSCYLCYL, ModelIndex.LDT);

	try {
	    double[] fitted = difInv.invert(ds3.nextVoxel());

            assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.968792E-6, fitted[1], 1.0E-10);
	    assertEquals(3.0, fitted[2], 1.0E-6);
	    assertEquals(0.26976, fitted[3], 1.0E-5);
	    assertEquals(6.68574E-10, fitted[4], 1.0E-15);
	    assertEquals(0.36254, fitted[10], 1.0E-5);


	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testThreeTensorAxiSymFitter() {
	
	DiffusionInversion difInv = 
	    new ThreeTensorInversion(ip, ModelIndex.CYLCYLCYL, ModelIndex.LDT);
	
	try {
	    double[] fitted = difInv.invert(ds3.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.986792E-6, fitted[1], 1.0E-6);
	    assertEquals(3.0, fitted[2], 1.0E-6);
	    assertEquals(0.33086, fitted[3], 1.0E-5);
	    assertEquals(7.95983E-10, fitted[4], 1.0E-15);
            assertEquals(0.250551, fitted[10], 1.0E-5);


	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testBallStickInversion() {

        String dataFile = "./test/inverters/ballStickTest.Bfloat";

        
        // BALL: vol frac= 0.2 params =  7.0E-10
        // STICK: vol frac= 0.8 params =  7.0E-10 1.0472 0.7854


        VoxelOrderDataSource data = new VoxelOrderDataSource(dataFile, 60, "float");

        double[] voxel = data.nextVoxel();

        DW_Scheme scheme = DW_Scheme.readScheme("M6_N54_b1500.scheme");

        // For no good reason, fitters decide to look to CL_Initializer to see what model they should fit
        tools.CL_Initializer.fitModel = "BALLSTICK";

        BallStickLM_Fitter fitter = new BallStickLM_Fitter(scheme);

        Vector3D stick = Vector3D.vectorFromSPC(1.0, 1.0472, 0.7854);

        try { 

            double[] solution = fitter.fit(voxel)[0];
            
            assertEquals(0.8, solution[2], 1E-6);
            assertEquals(0.2, solution[3], 1E-6);
            assertEquals(7E-10, solution[4], 1E-12);
            
            assertEquals(1.0, Math.abs(stick.dot(Vector3D.vectorFromSPC(1.0, solution[5], solution[6]))), 1E-6);
      
        }
        catch (optimizers.MinimizerException e) {
            fail(e.toString());
        }


        // check that we can get the same thing with BallStickInversion

        BallStickInversion inv = new BallStickInversion(scheme);

        double[] fit = inv.invert(voxel);
        
        assertEquals(7E-10, fit[2], 1E-12);
        assertEquals(0.8, fit[3], 1E-6);

        assertEquals(1.0, Math.abs(stick.dot(new Vector3D(fit[4], fit[5], fit[6]))), 1E-6);


    }

 
}

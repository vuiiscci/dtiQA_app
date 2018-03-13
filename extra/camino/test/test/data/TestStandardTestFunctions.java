package data;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import inverters.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>StandardTestFunctions.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestStandardTestFunctions.java,v 1.1 2005/06/15 09:12:29 ucacdxa Exp $
 * @author  Daniel Alexander
 * @see data.StandardTestFunctions
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestStandardTestFunctions extends TestCase {


    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestStandardTestFunctions(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

	ip = DW_Scheme.readScheme("bmx7.scheme");

    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestStandardTestFunctions.class);
    }


    public void testP0() {

	StandardTestFunctions.reset();
	ModelPDF p = StandardTestFunctions.getFunction(0);
	DataSynthesizer ds = new DataSynthesizer(p, ip, 100000);
	DiffusionInversion difInv = new LinearDT_Inversion(ip);

	try {
	    double[] fitted = difInv.invert(ds.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.838377E-6, fitted[1], 1.0E-9);
	    assertEquals(6.999986E-10, fitted[2], 1.0E-15);
	    assertEquals(6.307581E-15, fitted[3], 1.0E-18);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testP1() {

	StandardTestFunctions.reset();
	ModelPDF p = StandardTestFunctions.getFunction(1);
	DataSynthesizer ds = new DataSynthesizer(p, ip, 100000);
	DiffusionInversion difInv = new LinearDT_Inversion(ip);

	try {
	    double[] fitted = difInv.invert(ds.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.838377E-6, fitted[1], 1.0E-10);
	    assertEquals(1.699985E-9, fitted[2], 1.0E-15);
	    assertEquals(0.0, fitted[3], 1.0E-14);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testP2() {

	StandardTestFunctions.reset();
	ModelPDF p = StandardTestFunctions.getFunction(2);
	DataSynthesizer ds = new DataSynthesizer(p, ip, 100000);
	DiffusionInversion difInv = new LinearDT_Inversion(ip);

	try {
	    double[] fitted = difInv.invert(ds.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(5.838377E-6, fitted[1], 1.0E-10);
	    assertEquals(9.499958E-10, fitted[2], 1.0E-15);
	    assertEquals(0, fitted[3], 1.0E-14);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testP3() {

	StandardTestFunctions.reset();
	ModelPDF p = StandardTestFunctions.getFunction(3);
	DataSynthesizer ds = new DataSynthesizer(p, ip, 100000);
	DiffusionInversion difInv = new LinearDT_Inversion(ip);

	try {
	    double[] fitted = difInv.invert(ds.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(0.0, fitted[1], 1.0E-5);
	    assertEquals(7.683466E-10, fitted[2], 1.0E-15);
	    assertEquals(0, fitted[3], 1.0E-12);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testP4() {

	StandardTestFunctions.reset();
	ModelPDF p = StandardTestFunctions.getFunction(4);
	DataSynthesizer ds = new DataSynthesizer(p, ip, 100000);
	DiffusionInversion difInv = new LinearDT_Inversion(ip);

	try {
	    double[] fitted = difInv.invert(ds.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(0.0, fitted[1], 1.0E-5);
	    assertEquals(5.704022E-10, fitted[2], 1.0E-15);
	    assertEquals(0.0, fitted[3], 1.0E-5);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testP1_ScaleL1() {

	StandardTestFunctions.reset();
	StandardTestFunctions.setLambda1(19.0E-10);
	ModelPDF p = StandardTestFunctions.getFunction(1);
	DataSynthesizer ds = new DataSynthesizer(p, ip, 100000);
	DiffusionInversion difInv = new LinearDT_Inversion(ip);

	try {
	    double[] fitted = difInv.invert(ds.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(0.0, fitted[1], 1.0E-5);
	    assertEquals(1.899978E-9, fitted[2], 1.0E-15);
	    assertEquals(0.0, fitted[3], 1.0E-12);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testP1_ScaleTrD() {

	StandardTestFunctions.reset();
	StandardTestFunctions.setTraceD(23.0E-10);
	ModelPDF p = StandardTestFunctions.getFunction(1);
	DataSynthesizer ds = new DataSynthesizer(p, ip, 100000);
	DiffusionInversion difInv = new LinearDT_Inversion(ip);

	try {
	    double[] fitted = difInv.invert(ds.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(0, fitted[1], 1.0E-5);
	    assertEquals(1.699984E-9, fitted[2], 1.0E-15);
	    assertEquals(0, fitted[3], 1.0E-12);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    public void testP1_Rotated() {

	StandardTestFunctions.reset();
	RealMatrix rot = Rotations.randomRotMat(0);
	StandardTestFunctions.setTransformation(rot);
	ModelPDF p = StandardTestFunctions.getFunction(1);
	DataSynthesizer ds = new DataSynthesizer(p, ip, 100000);
	DiffusionInversion difInv = new LinearDT_Inversion(ip);

	try {
	    double[] fitted = difInv.invert(ds.nextVoxel());

	    assertEquals(0.0, fitted[0], 0.000001);
	    assertEquals(0.0, fitted[1], 1.0E-5);
	    assertEquals(8.245454E-10, fitted[2], 1.0E-15);
	    assertEquals(0.0, fitted[3], 1.0E-9);

	} catch(Exception e) {
	    fail("Data source failed.");
	}
    }


    
}

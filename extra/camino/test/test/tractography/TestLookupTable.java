package tractography;

import junit.framework.*;
import junit.extensions.*;



/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>LookupTable</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>LookupTable</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestLookupTable.java,v 1.2 2005/11/24 12:07:35 ucacpco Exp $
 * @author  Philip Cook
 * @see numerics.WatsonDistribution
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestLookupTable extends TestCase {

    private double[][][] lutArray2D1Val = null;
    private double[][][] lutArray2D2Vals = null;
    
    private double[][][] lutArray1D1Val = null;
    private double[][][] lutArray1D2Vals = null;

    private double[][][] lutArray3D2Vals = null;

    
    private final double step = 2.5;
    private final double xMin = 1.0;
    private final double xMax = 6.0;
    private final double yMin = 1.0;
    private final double yMax = 6.0;
    private final double zMin = 1.0;
    private final double zMax = 3.0;

    private final double zStep = 2.0;

    private LookupTable lut2D1Val = null;
    private LookupTable lut2D2Vals = null;
    
    private LookupTable lut1D1Val = null;
    private LookupTable lut1D2Vals = null;
    private LookupTable lut3D2Vals = null;

    
   
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	lutArray2D1Val = new double[1][3][3];
	lutArray1D1Val = new double[1][1][3];

	lutArray1D2Vals = new double[1][1][6];
	lutArray2D2Vals = new double[1][3][6];
	
	lutArray1D1Val[0][0][0] = 0.0;
	lutArray1D1Val[0][0][1] = 2.0;
	lutArray1D1Val[0][0][2] = 4.0;

	lutArray1D2Vals[0][0][0] = 0.0;
	lutArray1D2Vals[0][0][1] = 8.0;
	lutArray1D2Vals[0][0][2] = 2.0;
	lutArray1D2Vals[0][0][3] = 16.0;
	lutArray1D2Vals[0][0][4] = 4.0;
	lutArray1D2Vals[0][0][5] = 32.0;
	

	lutArray3D2Vals = new double[3][3][4];

	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 2; k++) {
		    lutArray3D2Vals[i][j][k*2] = i + i * j + (i + j) * k;
		    lutArray3D2Vals[i][j][k*2+1] = 2.0 * lutArray3D2Vals[i][j][k*2];
		}
	    }
	}

	

	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {
		lutArray2D1Val[0][i][j] = i + i * j;
	    }
	}


	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {
		for (int v = 0; v < 2; v++) {
		    lutArray2D2Vals[0][i][2 *j + v] = i + i * j * (v+1);
		}
	    }
	}

        lut1D1Val = new LookupTable(lutArray1D1Val, new double[] {0.0, 0.0, 0.0, 0.0, yMin, yMax}, new double[] {0.0, 0.0, step}, 1); 

	lut1D2Vals = new LookupTable(lutArray1D2Vals, new double[] {0.0, 0.0, 0.0, 0.0, yMin, yMax}, new double[] {0.0, 0.0, step}, 2); 

	lut2D1Val = new LookupTable(lutArray2D1Val, new double[] {0.0, 0.0, xMin, xMax, yMin, yMax}, new double[] {0.0, step, step}, 1);
 
	lut2D2Vals = new LookupTable(lutArray2D2Vals, new double[] {0.0, 0.0, xMin, xMax, yMin, yMax}, new double[] {0.0, step, step}, 2); 

	lut3D2Vals = new LookupTable(lutArray3D2Vals, new double[] {xMin, xMax, yMin, yMax, zMin, zMax}, new double[] {step, step, zStep}, 2); 


    }

   
    protected void tearDown() {
    }


    public TestLookupTable(String name) {
	super(name);
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestLookupTable.class);
    }


    public void testBoundsChecking() {

	testExceptionThrown(lut3D2Vals, 1.0, 1.0, 4.0); // z > zMax

	testExceptionThrown(lut2D1Val, 1.0, 7.0); // y > yMax
	testExceptionThrown(lut2D1Val, 1.0, -17.0); // y < yMin
	testExceptionThrown(lut2D1Val, 8.0, 1.0); // x > xMax
	testExceptionThrown(lut2D1Val, 0.0, 2.0); // x < xMin

	// test 1D case
	testExceptionThrown(lut1D1Val, 0.0, 0.0); // y < yMin
	testExceptionThrown(lut1D1Val, 0.0, 10.0); // y > yMax
    }

    public void testClamping() {
	
	try {
	    // i + i * j + (i + j) * k

	    // test z clamped
	    assertEquals(10.0, lut3D2Vals.getValues(xMax, yMax, zMax + 1, false, true)[0], 0.000001); 
	    assertEquals(20.0, lut3D2Vals.getValues(xMax, yMax, zMax + 1, false, true)[1], 0.000001); 
    
	    assertEquals(6.0, lut3D2Vals.getValues(xMax, yMax, zMin - 1.0, false, true)[0], 0.000001); 

	    // test x and y clamped
	    assertEquals(2.0, lut3D2Vals.getValues(xMin - 1.0, yMax + 1.0, zMax, false, true)[0], 0.000001); 
	    assertEquals(4.0, lut3D2Vals.getValues(xMax + 1.0, yMin - 1, zMax, false, true)[0], 0.000001); 


	}
	catch (OutsideLUTRangeException e) {
	    fail("Unexpected exception " + e);
	}

	// 2D case
	try {
	    // i + i * j
	    assertEquals(6.0, lut2D1Val.getValues(xMax + 1, yMax + 1, false, true)[0], 0.000001); 
	    assertEquals(0.0, lut2D1Val.getValues(xMin - 1.0, yMin - 1.0, false, true)[0], 0.000001); 
	    
	    // text x clamped
	    assertEquals(4.0, lut2D1Val.getValues(xMax + 2.0, 3.5, false, true)[0], 0.000001); 
	}
	catch (OutsideLUTRangeException e) {
	    fail("Unexpected exception " + e);
	}

	// 1D case

	try {
	    
	    assertEquals(4.0, lut1D1Val.getValues(7.0, false, true)[0], 0.000001); 
	    assertEquals(0.0, lut1D1Val.getValues(0.0, false, true)[0], 0.000001); 
	    
	}
	catch (OutsideLUTRangeException e) {
	    fail("Unexpected exception " + e);
	}

    }

    public void testInterpolation() {
	try {
	    // if interp == false, we should get nearest neighbour

	    // 3D

	    // i + i * j + (i + j) * k
	    assertEquals(6.0, lut3D2Vals.getValues(5.5, 5.5, 1.0, false, true)[0], 0.000001);
	    assertEquals(12.0, lut3D2Vals.getValues(5.5, 5.5, 1.0, false, true)[1], 0.000001);

	    
	    assertEquals(3.0, lut3D2Vals.getValues(3.5, 5.5, 1.0, false, true)[0], 0.000001);
	    assertEquals(6.0, lut3D2Vals.getValues(3.5, 5.5, 1.0, false, true)[1], 0.000001);


	    // i + i * j + (i + j) * k
	    assertEquals(6.84, lut3D2Vals.getValues(5.5, 5.5, 2.0,  true, true)[0], 0.000001);
	    assertEquals(13.68, lut3D2Vals.getValues(5.5, 5.5, 2.0, true, true)[1], 0.000001);

	    // make sure that the last value in the LUT is accessible (checks for rounding errors)
	    assertEquals(6.0, lut3D2Vals.getValues(6.0, 6.0, 1.0, true, true)[0], 0.000001);
	    assertEquals(12.0, lut3D2Vals.getValues(6.0, 6.0, 1.0, true , true)[1], 0.000001);


	    // 2D

	    // i + i * j
	    assertEquals(3.96, lut2D2Vals.getValues(5.5, 4.0, true, true)[0], 0.000001);
	    
	    // i + 2(i*j)
	    assertEquals(6.12, lut2D2Vals.getValues(5.5, 4.0, true, true)[1], 0.000001);
	    
	    // i + i * j
	    assertEquals(6.0, lut2D2Vals.getValues(5.0, 5.0, false, true)[0], 0.000001);
	    
	    // i + 2(i*j)
	    assertEquals(10.0, lut2D2Vals.getValues(5.0, 5.0, false, true)[1], 0.000001);


	    assertEquals(3.0, lut2D2Vals.getValues(3.5, 5.0, false, true)[0], 0.000001);
	    assertEquals(6.0, lut2D2Vals.getValues(5.0, 3.5, false, true)[1], 0.000001);


	    assertEquals(6.0, lut2D2Vals.getValues(6.0, 6.0, false, true)[0], 0.000001);
	    assertEquals(10.0, lut2D2Vals.getValues(6.0, 6.0, false, true)[1], 0.000001);

	
	    // 1D
	    assertEquals(2.0, lut1D2Vals.getValues(3.25, false, true)[0], 0.000001);
	    assertEquals(8.8, lut1D2Vals.getValues(1.25, true, true)[1], 0.000001);

	    
	}
	catch (OutsideLUTRangeException e) {
	    fail("Unexpected exception " + e);
	}
    }

    public void testLimits() {
        try {
	    assertEquals(10.0, lut3D2Vals.getValues(xMax, yMax, zMax, true, true)[0], 0.000001);
        }
        catch(OutsideLUTRangeException e) {
            fail(e.toString());
        }
    }

    
    private void testExceptionThrown(LookupTable lut, double x, double y) {
	try {

	    lut.getValues(x,y,false, false);
	    fail("Should have thrown an OutsideLUTRangeException");
	}
	catch (OutsideLUTRangeException e) {
	}
    }


    private void testExceptionThrown(LookupTable lut, double x, double y, double z) {
	try {

	    lut.getValues(x,y,z,false, false);
	    fail("Should have thrown an OutsideLUTRangeException");
	}
	catch (OutsideLUTRangeException e) {
	}
    }


}

package imaging;

import apps.*;
import misc.*;
import numerics.*;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 *
 * Automated tests for schemes that extend DW_Scheme. 
 *
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see imaging.DW_Scheme
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public abstract class TestDW_Scheme extends TestCase {

    /** Schemes are immutable, so can initialize this once. */
    protected DW_Scheme scheme = null;

    protected double[] data = null;


    public TestDW_Scheme(String name) {
	super(name);
	scheme = getScheme();
    }

    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initializes global resources for tests.
     */
    protected void setUp() {
	data = new double[15];

	for (int i = 0; i < data.length; i++) {
	    data[i] = i+1;
	}
	
    }

    /**
     * Does the opposite of setup. Should be used to free memory, I/O resources, etc.
     */    
    protected void tearDown() {
	data = null;
    }
    
    public static Test suite() {
	return new TestSuite(TestDW_Scheme.class);
    }


    /**
     * Tests normalizeData(), normalizeData(double[], int[]), and normalizeData(double[], double). 
     *
     */
    public void testNormalizeData() {

	double[] tooMuchData = new double[16];

	// should complain if fed too much data
	try {
	    LoggedException.setExceptionLogging(false);
	    scheme.normalizeData(tooMuchData);
	    LoggedException.setExceptionLogging(true);
	    fail("Scheme tried to normalize data array that did not match the number of measurements");
	}
	catch (Exception e) {
	    LoggedException.setExceptionLogging(true);
	}

	double[] normData = scheme.normalizeData(data);

	assertEquals(13, normData.length);

	// normalization constant == 3.0

	assertEquals(2.0 / 3.0, normData[0], 1E-8);
	assertEquals(3.0 / 3.0, normData[1], 1E-8);


	double[] normDataByIndices = scheme.normalizeData(data, new int[] {0, 8});

	for (int i = 0; i < normData.length; i++) {
	    assertEquals(normData[i], normDataByIndices[i], 1E-8);
	}

    }

    public void testGeoMeanZeroMeas() {

	assertEquals(3.0, scheme.geoMeanZeroMeas(data), 1E-8);
	
    }


    /**
     * Only works if getB_Value() also works.
     */
    public void testGetB_Matrix() {

	try {

	    RealMatrix B = scheme.getB_Matrix();

	    java.net.URL url = TestDW_Scheme.class.getResource("/PointSets/Elec007.txt");

	    double[][] points = PointSetToScheme.readPoints(url.getFile());


	    // first measurement is b=0
	    assertEquals(1.0, B.entries[0][0], 1E-8);

	    for (int j = 0; j < 6; j++) {
		assertEquals(0.0, B.entries[0][j+1], 1E-8);
	    }

	    double b = 1000E6;
	
	    for (int i = 0; i < 7; i++) {
		assertEquals(1.0, B.entries[i+1][0], 1E-8);
		
		double ans = -b * points[i][0] * points[i][0];
		
		assertEquals(ans, B.entries[i+1][1], Math.abs(ans) / 1000.0);

		ans = -b * 2.0 * points[i][0] * points[i][1];

		assertEquals(ans, B.entries[i+1][2], Math.abs(ans) / 1000.0);

		ans = -b * 2.0 * points[i][0] * points[i][2];

		assertEquals(ans, B.entries[i+1][3], Math.abs(ans) / 1000.0);

		ans = -b * points[i][1] * points[i][1];

		assertEquals(ans, B.entries[i+1][4], Math.abs(ans) / 1000.0);

		ans = -b * 2.0 * points[i][1] * points[i][2]; 
		
		assertEquals(ans, B.entries[i+1][5], Math.abs(ans) / 1000.0);

		ans = -b * points[i][2] * points[i][2];

		assertEquals(ans, B.entries[i+1][6], Math.abs(ans) / 1000.0);
	    }
	
	    // another b=0
	    assertEquals(1.0, B.entries[8][0], 1E-8);
	
	    for (int j = 0; j < 6; j++) {
		assertEquals(0.0, B.entries[8][j+1], 1E-8);
	    }

	    url = TestDW_Scheme.class.getResource("/PointSets/Elec006.txt");
	    points = PointSetToScheme.readPoints(url.getFile()); 

	    b = 100E6;

	    for (int i = 0; i < 6; i++) {
		assertEquals(1.0, B.entries[i+9][0], 1E-8);
		
		double ans = -b * points[i][0] * points[i][0];
		
		assertEquals(ans, B.entries[i+9][1], Math.abs(ans) / 1000.0);

		ans = -b * 2.0 * points[i][0] * points[i][1];

		assertEquals(ans, B.entries[i+9][2], Math.abs(ans) / 1000.0);

		ans = -b * 2.0 * points[i][0] * points[i][2];

		assertEquals(ans, B.entries[i+9][3], Math.abs(ans) / 1000.0);

		ans = -b * points[i][1] * points[i][1];

		assertEquals(ans, B.entries[i+9][4], Math.abs(ans) / 1000.0);

		ans = -b * 2.0 * points[i][1] * points[i][2]; 
		
		assertEquals(ans, B.entries[i+9][5], Math.abs(ans) / 1000.0);

		ans = -b * points[i][2] * points[i][2];

		assertEquals(ans, B.entries[i+9][6], Math.abs(ans) / 1000.0);
	    }

	    
	}
	catch (java.io.IOException e) {
	    fail(e.toString());
	}
    }


    /**
     * Tests the flipped version for proper flipping of gradient directions, then tries a double flip,
     * ie does scheme.flipX().flipX() produce scheme.
     *
     */
    public void testFlip() {

	try {
	
	    DW_Scheme[] flips = new DW_Scheme[] {scheme.flipX(), scheme.flipY(), scheme.flipZ()};
	
	    for (int j = 0; j < 3; j++) {
		assertEquals(0.0, flips[j].getG_Dir(0)[j], 1E-8);
	    }
	
	    java.net.URL url = TestDW_Scheme.class.getResource("/PointSets/Elec007.txt");
	    
	    double[][] points = PointSetToScheme.readPoints(url.getFile());

	    for (int f = 0; f < 3; f++) {
		for (int i = 0; i < 7; i++) {
		
		    double[] gDir = flips[f].getG_Dir(i+1);
		
		    for (int j = 0; j < 3; j++) {
			if (j == f) {
			    assertEquals(-points[i][j], gDir[j], 1E-6);
			}
			else {
			    assertEquals(points[i][j], gDir[j], 1E-6);
			}
		    }
		}
	    }


	    // now see if everything else gets copied correctly
	    DW_Scheme same = flips[0].flipX();

	    String twiceFlipped = same.toString();

	    assertEquals(scheme.toString(), same.toString());

	}
	catch (java.io.IOException e) {
	    fail(e.toString());
	}

    }


    public void testGradOrder() {
	// just test the first non-zero measurement

	DW_Scheme xzy = scheme.gradOrder(DW_Scheme.gradXZY);

	double[] gDir = scheme.getG_Dir(1);
	double[] xzy1 = xzy.getG_Dir(1);

	assertEquals(gDir[0], xzy1[0], 1E-8);
	assertEquals(gDir[1], xzy1[2], 1E-8);
	assertEquals(gDir[2], xzy1[1], 1E-8);

    }


    public void testGetSubsetScheme() {

	try {

	    java.net.URL url = TestDW_Scheme.class.getResource("/PointSets/Elec007.txt");

	    double[][] points7 = PointSetToScheme.readPoints(url.getFile());

	    url = TestDW_Scheme.class.getResource("/PointSets/Elec006.txt");
	
	    double[][] points6 = PointSetToScheme.readPoints(url.getFile());

	    int[] indices = new int[] {1, 2, 3, 9, 10, 11};

	    DW_Scheme subset = scheme.getSubsetScheme(indices);

	    for (int i = 0; i < indices.length; i++) {
	    
		double[] gDir = scheme.getG_Dir(indices[i]);

		double[] subsetG_Dir = subset.getG_Dir(i);

		for (int j = 0; j < 3; j++) {
		    assertEquals(gDir[j], subsetG_Dir[j], 1E-6);
		}
	    
	    }
	
	    // now get complete scheme out of subset method and ensure that it produces the same
	    // Scheme file

	    indices = new int[15];

	    for (int i = 0; i < indices.length; i++) {
		indices[i] = i;
	    }

	    subset = scheme.getSubsetScheme(indices);

	    assertEquals(scheme.toString(), subset.toString());
	}
	catch (java.io.IOException e) {
	    fail(e.toString());
	}
    }


    /**
     * Doesn't do anything that subclasses can influence, so it doesn't need to run for all subclasses,
     * but I can't figure out an easy way to exclude it.
     */
    public final void testNormalizeGradDirs() {

	double[][] points = new double[3][3];
	
	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {
		points[i][j] = i + j + i*j;
	    }
	}

	double[][] normPoints = DW_Scheme.normalizeGradDirs(points);

	for (int i = 0; i < 3; i++) {
	    Vector3D norm = new Vector3D(points[i]).normalized();
	    assertEquals(norm.x, normPoints[i][0], 1E-8);
	    assertEquals(norm.y, normPoints[i][1], 1E-8);
	    assertEquals(norm.z, normPoints[i][2], 1E-8);
	}

    }


    /**
     * Tests that the toString method gives you what was in the scheme file.
     */
    public void testToString() {

	File schemeFile = new File(getSchemeFile());
	
	
	try {
	    FileInputStream fis = new FileInputStream(schemeFile);
	    
	    byte[] array = new byte[(int)schemeFile.length()];
	    
	    int bytesRead = fis.read(array);
	    
	    while (bytesRead < array.length) {
		fis.read(array, bytesRead, array.length - bytesRead);
	    }
	    
	    String fromFile = new String(array);
	    
	    assertEquals(fromFile, scheme.toString());

	    fis.close();
	    
	}
	catch (IOException e) {
	    fail(e.toString());
	}

    }

    /**
     * Gets the scheme file name for the class.
     *
     */
    protected abstract String getSchemeFile();


    /**
     * Gets the scheme for the class being tested. The test scheme is read from a file
     * in the test/imaging/ directory. The scheme consists of 15 measurements,
     * 1 at b=0 followed by 7 at b=1000, another b=0, and 6 at b=100. The imaging parameters 
     * of all schemes are chosen to produce the same (to within a small tolerance) b-matrix for 
     * all scheme versions.
     *
     */
    protected final DW_Scheme getScheme() {
	return DW_Scheme.readScheme(getSchemeFile());
    }


}

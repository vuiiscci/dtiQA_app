package numerics;

import Jama.Matrix;

import junit.framework.*;
import junit.extensions.*;

/**
 * <dl>
 * <dt>Purpose: General test methods for numerics classes.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This contains some general methods for testing things such as matrix equality with a delta.
 * The tests have tests, so include this class in your test suite.
 * <dd> The coordinate system is right handed, as defined in Foley and Van Dam (red book) page 180.
 * 
 * </dl>
 *
 * @version $Id: GenTestMethods.java,v 1.2 2005/09/26 10:36:45 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class GenTestMethods extends TestCase {


    public GenTestMethods(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }
    
    public static Test suite() {
	return new TestSuite(GenTestMethods.class);
    }


    public void testMatricesEqual() {

	Matrix m = new Matrix(2,2);
	
	m.set(0,0, 193.5107);
	m.set(0,1, 153.10965);
	m.set(1,0, 91.013653);
	m.set(1,1, 18.605);

	Matrix m2 = new Matrix(m.getArrayCopy());

	assertTrue(assertMatricesEqual(m, m2, 0.0000001, true));

	Matrix m3 = new Matrix(m.getArrayCopy());
	
	m3.set(1,1, 18.604);		   

	assertTrue( assertMatricesEqual(m, m3, 0.01, true) );
	assertFalse( assertMatricesEqual(m, m3, 0.0001, true) );
	
	
    }


    public void testRealMatricesEqual() {

	RealMatrix m = new RealMatrix(2,2);
	
	m.entries[0][0] = 193.5107;
	m.entries[0][1] = 153.10965;
	m.entries[1][0] = 91.013653;
	m.entries[1][1] = 18.605;

	RealMatrix m2 = (RealMatrix)m.clone();

	assertTrue( assertMatricesEqual(m, m2, 0.0000001, true) );

	RealMatrix m3 = (RealMatrix)m.clone();
	
	m3.entries[1][1] = 18.604;		   

	assertTrue( assertMatricesEqual(m, m3, 0.01, true) );
	assertFalse( assertMatricesEqual(m, m3, 0.0001, true) );
	
    }



    public void testVectorsEqual() {

	Vector3D v1 = new Vector3D(1.0, 2.0, 3.0);
	Vector3D v2 = new Vector3D(1.1, 2.0, 3.0);	
	Vector3D v3 = new Vector3D(1.0, 2.01, 3.0);
	Vector3D v4 = new Vector3D(1.0, 2.0, 3.001);

	assertTrue(assertVectorsEqual(v1, v1, 0.0000001, true));

	assertFalse(assertVectorsEqual(v1, v2, 0.09, true));
	assertTrue(assertVectorsEqual(v1, v2, 0.11, true));

	assertTrue(assertVectorsEqual(v1, v3, 0.011, true));
	assertFalse(assertVectorsEqual(v1, v3, 0.005, true));


	assertTrue(assertVectorsEqual(v1, v4, 0.0011, true));
	assertFalse(assertVectorsEqual(v1, v4, 0.0005, true));
	
    }


    public void testArraysEqual() {

	double[][] a = new double[3][3];
	double[][] b = new double[3][3];
	
	assertTrue(assertArraysEqual(a,b,0.00000001, true));

	b = new double[3][4];

	assertFalse(assertArraysEqual(a,b,0.00000001, true));

	b = new double[3][3];

	b[1] = new double[5];

	assertFalse(assertArraysEqual(a,b,0.00000001, true));

	b = new double[3][3];

	a[1][1] = 3.14159;
	b[1][1] = 3.14158;

	assertTrue(assertArraysEqual(a,b,0.001, true));

	assertFalse(assertArraysEqual(a,b,0.00000001, true));

	double[] c = new double[4];
	double[] d = new double[5];
	double[] e = new double[5];

	d[4] = 3.14159;
	e[4] = 3.14158;

	assertFalse(assertArraysEqual(c,d, 0.000001, true));
	assertTrue(assertArraysEqual(d,e, 0.001, true));
	assertFalse(assertArraysEqual(d,e, 0.000001, true));

    }




    /**
     * Test if two matrices are equal within the given precision.
     *
     *
     * @param delta maximum difference between corresponding elements.
     *
     */
    public static void assertMatricesEqual(Matrix m1, Matrix m2, double delta) {
	assertMatricesEqual(m1,m2,delta, false);
    }

    private static boolean assertMatricesEqual(Matrix m1, Matrix m2, double delta, boolean test) {
	
	if (m1 == m2) {
	    return true;
	}
	    
	int rows = m1.getRowDimension();
	
	if ( rows != m2.getRowDimension() ) {
	    if (test) {
		return false;
	    }
	    else {
		fail("Expected " + rows + " rows but got " + m2.getRowDimension());
	    }

	}
	
	int cols = m1.getColumnDimension();
	
	if ( cols != m2.getColumnDimension() ) {
	    if (test) {
		return false;
	    }
	    else {
		fail("Expected " + cols + " cols but got " + m2.getColumnDimension());
	    }
	}

	
	for (int i = 0; i < rows; i++) {
	    for (int j = 0; j < cols; j++) {
		if (Math.abs( m1.get(i,j) - m2.get (i,j) ) > delta ) {
		    if (test) {
			return false;
		    }
		    else {
			fail("Expected element [" + i + "][" + j + "] to be " + 
			     m1.get(i,j) + " but was " + m2.get(i,j));
		    }
		}
	    }
	}
	
	return true;
	
    }



    /**
     * Asserts that two matrices are equal within the given precision.
     *
     * @param delta maximum difference between corresponding elements.
     *
     */
    public static void assertMatricesEqual(RealMatrix m1, RealMatrix m2, double delta) {
	assertMatricesEqual(m1, m2, delta, false);
    }


    private static boolean assertMatricesEqual(RealMatrix m1, RealMatrix m2, double delta, boolean test) {
	
	if (m1 == m2) {
	    return true;
	}
	    
	int rows = m1.r;
	
	if ( rows != m2.r ) {
	    if (test) {
		return false;
	    }
	    else {
		fail("Expected " + rows + " rows but got " + m2.r);
		return false;
	    }
	}
	
	int cols = m1.c;
	
	if ( cols != m2.c ) {
	    
	    if (test) {
		return false;
	    }
	    else {
		fail("Expected " + cols + " cols but got " + m2.c);
	    }
	}

	
	for (int i = 0; i < rows; i++) {
	    for (int j = 0; j < cols; j++) {
		if ( Math.abs(m1.entries[i][j] - m2.entries[i][j]) > delta ) {
		    if (test) {
			return false;
		    }
		    else {
			fail("Expected element [" + i + "][" + j + "] to be " + m1.entries[i][j] + 
			     " but was " + m2.entries[i][j]);
		    }
		}
	    }
	}
	
	return true;
	
    }


    /**
     * Asserts that two arrays are equal within the given precision. 
     *
     */
    public static void assertArraysEqual(double[][] a1, double[][] a2, double delta) {
	assertArraysEqual(a1, a2, delta, false);
    }


    private static boolean assertArraysEqual(double[][] a1, double[][] a2, double delta, boolean test) {

	if (a1 == a2) {
	    return true;
	}
	
	if (a1.length != a2.length) {
	    if (test) {
		return false;
	    }
	    else {
		fail("Expected first dimension " + a1.length + " but was " + a2.length);
		return false;
	    }
	}
	
	
	for (int i = 0; i < a1.length; i++) {

	    if (a1[i].length != a2[i].length) {
		if (test) {
		    return false;
		}
		else {
		    fail("Expected second dimension " + a1[i].length + " but was " + a2[i].length);
		    return false;
		}
	    }

	    for (int j = 0; j < a1[i].length; j++) {

		if ( Math.abs(a1[i][j] - a2[i][j]) > delta ) {
		    if (test) {
			return false;
		    }
		    else {
			fail("Expected element [" + i + "][" + j + "] to be " + a1[i][j] + 
			     " but was " + a2[i][j]);
			return false;
		    }
		}
	    }
	}

	return true;
	
    }



    /**
     * Asserts that two arrays are equal within the given precision. 
     *
     */
    public static void assertArraysEqual(double[] a1, double[] a2, double delta) {
	assertArraysEqual(a1, a2, delta, false);
    }
	
    private static boolean assertArraysEqual(double[] a1, double[] a2, double delta, boolean test) {

	if (a1 == a2) {
	    return true;
	}
	
	if (a1.length != a2.length) {
	    if (test) {
		return false;
	    }
	    else {
		fail("expected array length " + a1.length + " but was " + a2.length);
		return false;
	    }
	}
	
	
	for (int i = 0; i < a1.length; i++) {
	    
	    if ( Math.abs(a1[i] - a2[i]) > delta ) {
		if (test) {
		    return false;
		}
		else {
		    fail("Expected element [" + i + "] to be " + a1[i] + " but was " + a2[i]);
		    return false;
		}
	    }
	}
	
	return true;
	
    }


  
    /**
     * Asserts that two vectors are equal within the given precision. 
     *
     */
    public static void assertVectorsEqual(Vector3D v1, Vector3D v2, double delta) {
	assertVectorsEqual(v1, v2, delta, false);
    }

    private static boolean assertVectorsEqual(Vector3D v1, Vector3D v2, double delta, boolean test) {
	
	if ( Math.abs(v1.x - v2.x) > delta) {
	    if (test) {
		return false;
	    }
	    else {
		fail("Expected x to be " + v1.x + " but was " + v2.x);
		return false;
	    }
	}
	if ( Math.abs(v1.y - v2.y) > delta) {
	    if (test) {
		return false;
	    }
	    else {
		fail("Expected x to be " + v1.y + " but was " + v2.y);
		return false;
	    }
	}
	if ( Math.abs(v1.z - v2.z) > delta) {
	    if (test) {
		return false;
	    }
	    else {
		fail("Expected x to be " + v1.z + " but was " + v2.z);
		return false;
	    }
	}
	
	return true;
	
    }
  
}




package numerics;

import junit.framework.*;
import junit.extensions.*;

import Jama.Matrix;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>EigenSystem3D.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>EigenSystem3D</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestEigenSystem3D.java,v 1.1 2005/09/26 10:36:47 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestEigenSystem3D extends TestCase {


    public TestEigenSystem3D(String name) {
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
	return new TestSuite(TestEigenSystem3D.class);
    }


    /**
     * Tests the sort(Jama.EigenvalueDecomposition) method.
     *
     */
    public void testJamaSorting() {
	
	Matrix m = new Matrix(3,3);

	m.set(0,0, 1000.0);
	m.set(1,1, 100.0);
	m.set(2,2, 10.0);

	Vector3D x = new Vector3D(1.0, 0.0, 0.0);
	Vector3D y = new Vector3D(0.0, 1.0, 0.0);	
	Vector3D z = new Vector3D(0.0, 0.0, 1.0);

	EigenSystem3D eig = EigenSystem3D.sort(m.eig());

	assertEquals(1.0, Math.abs(x.dot(eig.eigenvectors[0])), 0.000001);
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[1])), 0.000001);
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(10.0, eig.eigenvalues[2], 0.000001);


	// now change order
	m.set(0,0, 100.0);
	m.set(1,1, 1000.0);
	m.set(2,2, 10.0);

	eig = EigenSystem3D.sort(m.eig());
	
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[0])), 0.000001);
	assertEquals(1.0, Math.abs(x.dot(eig.eigenvectors[1])), 0.000001);
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(10.0, eig.eigenvalues[2], 0.000001);


	
	m.set(0,0, 10.0);
	m.set(1,1, 1000.0);
	m.set(2,2, 100.0);

	eig = EigenSystem3D.sort(m.eig());
	
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[0])), 0.000001);
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[1])), 0.000001);
	assertEquals(1.0, Math.abs(x.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(10.0, eig.eigenvalues[2], 0.000001);

	
	m.set(0,0, 100.0);
	m.set(1,1, 10.0);
	m.set(2,2, 1000.0);

	eig = EigenSystem3D.sort(m.eig());
	
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[0])), 0.000001);
	assertEquals(1.0, Math.abs(x.dot(eig.eigenvectors[1])), 0.000001);
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(10.0, eig.eigenvalues[2], 0.000001);

	// degenerate first second 
	m.set(0,0, 1000.0);
	m.set(1,1, 100.0);
	m.set(2,2, 1000.0);

	eig = EigenSystem3D.sort(m.eig());
	
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(1000.0, eig.eigenvalues[1], 0.000001);
	assertEquals(100.0, eig.eigenvalues[2], 0.000001);



	// degenerate second third
	m.set(0,0, 100.0);
	m.set(1,1, 100.0);
	m.set(2,2, 1000.0);

	eig = EigenSystem3D.sort(m.eig());
	
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[0])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(100.0, eig.eigenvalues[2], 0.000001);

    }


    /**
     * Tests the sort(RealMatrix) method.
     *
     */
    public void testRealMatrixSorting() {
	
	RealMatrix m = new RealMatrix(3,3);

	m.setEntry(0,0, 1000.0);
	m.setEntry(1,1, 100.0);
	m.setEntry(2,2, 10.0);

	Vector3D x = new Vector3D(1.0, 0.0, 0.0);
	Vector3D y = new Vector3D(0.0, 1.0, 0.0);	
	Vector3D z = new Vector3D(0.0, 0.0, 1.0);

	EigenSystem3D eig = EigenSystem3D.sort(m);

	assertEquals(1.0, Math.abs(x.dot(eig.eigenvectors[0])), 0.000001);
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[1])), 0.000001);
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(10.0, eig.eigenvalues[2], 0.000001);


	// now change order
	m.setEntry(0,0, 100.0);
	m.setEntry(1,1, 1000.0);
	m.setEntry(2,2, 10.0);

	eig = EigenSystem3D.sort(m);
	
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[0])), 0.000001);
	assertEquals(1.0, Math.abs(x.dot(eig.eigenvectors[1])), 0.000001);
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(10.0, eig.eigenvalues[2], 0.000001);


	
	m.setEntry(0,0, 10.0);
	m.setEntry(1,1, 1000.0);
	m.setEntry(2,2, 100.0);

	eig = EigenSystem3D.sort(m);
	
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[0])), 0.000001);
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[1])), 0.000001);
	assertEquals(1.0, Math.abs(x.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(10.0, eig.eigenvalues[2], 0.000001);

	
	m.setEntry(0,0, 100.0);
	m.setEntry(1,1, 10.0);
	m.setEntry(2,2, 1000.0);

	eig = EigenSystem3D.sort(m);
	
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[0])), 0.000001);
	assertEquals(1.0, Math.abs(x.dot(eig.eigenvectors[1])), 0.000001);
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(10.0, eig.eigenvalues[2], 0.000001);

	// degenerate first second 
	m.setEntry(0,0, 1000.0);
	m.setEntry(1,1, 100.0);
	m.setEntry(2,2, 1000.0);

	eig = EigenSystem3D.sort(m);
	
	assertEquals(1.0, Math.abs(y.dot(eig.eigenvectors[2])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(1000.0, eig.eigenvalues[1], 0.000001);
	assertEquals(100.0, eig.eigenvalues[2], 0.000001);



	// degenerate second third
	m.setEntry(0,0, 100.0);
	m.setEntry(1,1, 100.0);
	m.setEntry(2,2, 1000.0);

	eig = EigenSystem3D.sort(m);
	
	assertEquals(1.0, Math.abs(z.dot(eig.eigenvectors[0])), 0.000001);

	assertEquals(1000.0, eig.eigenvalues[0], 0.000001);
	assertEquals(100.0, eig.eigenvalues[1], 0.000001);
	assertEquals(100.0, eig.eigenvalues[2], 0.000001);

    }

  
}

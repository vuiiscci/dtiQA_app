package numerics;

import misc.DT;

import Jama.Matrix;

import junit.framework.*;
import junit.extensions.*;
import java.util.Random;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>Rotations.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Rotations</code> with JUnit 3.8.
 * <dd> The coordinate system is right handed, as defined in Foley and Van Dam (red book) page 180.
 * 
 * </dl>
 *
 * @version $Id: TestRotations.java,v 1.1 2005/09/26 10:36:47 ucacpco Exp $
 * @author  Philip Cook
 * @see numerics.Rotations
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestRotations extends TestCase {


    public TestRotations(String name) {
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
	return new TestSuite(TestRotations.class);
    }



    public void testGetRotMat() {

	Random r = new Random(213);

	Vector3D axis = new Vector3D(r.nextGaussian(), r.nextGaussian(), r.nextGaussian()).normalized();

	double theta = Math.PI / 6.0;

	double[] tp = Vector3D.thetaPhi(axis);

	RealMatrix rotAB = Rotations.getRotMat(tp[0], tp[1], theta);
	RealMatrix rotRDouble = Rotations.getRotMat(new double[] {axis.x, axis.y, axis.z}, theta);
	RealMatrix rotRVec = Rotations.getRotMat(axis, theta);

	GenTestMethods.assertMatricesEqual(rotAB, rotRDouble, 0.000001);
	GenTestMethods.assertMatricesEqual(rotAB, rotRVec, 0.000001);

	
    }

    /**
     * If this works, then #getRotMat(Vector3D, double) also works.
     */
    public void testRotateVectorAboutArbAxis() {
	
	Vector3D x = new Vector3D(2.0, 0.0, 0.0); // not mod 1.0, checks for normalization of axis
	Vector3D y = new Vector3D(0.0, 1.0, 0.0);
	Vector3D z = new Vector3D(0.0, 0.0, 1.0);

	double angle = Math.PI / 3.0;

	// should rotate towards -Z axis
	Vector3D xRotY = Rotations.rotateVector(x, y, angle);

	assertEquals(x.mod(), xRotY.mod(), 0.000001);

	assertEquals(x.mod() * Math.cos(angle), xRotY.x, 0.000001);
	assertEquals(x.y, xRotY.y, 0.000001);
	assertEquals(-1.0 * x.mod() * Math.sin(angle), xRotY.z, 0.000001);

	// should rotate towards +Z
	Vector3D xRotMinusY = Rotations.rotateVector(x, y, -1.0 * angle);

	assertEquals(x.mod(), xRotMinusY.mod(), 0.000001);

	assertEquals(x.mod() * Math.cos(angle), xRotMinusY.x, 0.000001);
	assertEquals(x.y, xRotMinusY.y, 0.000001);
	assertEquals(-1.0 * x.mod() * Math.sin(-1.0 * angle), xRotMinusY.z, 0.000001);

	// should rotate X towards Y
	Vector3D xRotZ = Rotations.rotateVector(x, z, angle);

	assertEquals(x.mod(), xRotZ.mod(), 0.000001);

	assertEquals(x.mod() * Math.cos(angle), xRotZ.x, 0.000001);
	assertEquals(x.mod() * Math.sin(angle), xRotZ.y, 0.000001);
	assertEquals(0.0, xRotZ.z, 0.000001);

	
	Vector3D xRotMinusZ = Rotations.rotateVector(x, z, -1.0 * angle);

	assertEquals(x.mod(), xRotZ.mod(), 0.000001);

	assertEquals(x.mod() * Math.cos(angle), xRotMinusZ.x, 0.000001);
	assertEquals(x.mod() * Math.sin(-1.0 * angle), xRotMinusZ.y, 0.000001);
	assertEquals(0.0, xRotMinusZ.z, 0.000001);


	Vector3D yRotX = Rotations.rotateVector(y, x, angle);

	assertEquals(y.mod(), yRotX.mod(), 0.000001);

	assertEquals(0.0, yRotX.x, 0.000001);
	assertEquals(y.mod() * Math.cos(angle), yRotX.y, 0.000001);
	assertEquals(y.mod() * Math.sin(angle), yRotX.z, 0.000001);

	Vector3D yRotMinusX = Rotations.rotateVector(y, Rotations.X_AXIS, -1.0 * angle);

	assertEquals(y.mod(), yRotMinusX.mod(), 0.000001);

	assertEquals(0.0, yRotMinusX.x, 0.000001);
	assertEquals(y.mod() * Math.cos(angle), yRotMinusX.y, 0.000001);
	assertEquals(y.mod() * Math.sin(-1.0 * angle), yRotMinusX.z, 0.000001);

	
	Vector3D yRotZ = Rotations.rotateVector(y, z, angle);

	assertEquals(y.mod(), yRotZ.mod(), 0.000001);


	assertEquals(y.mod() * -1.0 * Math.sin(angle), yRotZ.x, 0.000001);
	assertEquals(y.mod() * Math.cos(angle), yRotZ.y, 0.000001);
	assertEquals(0.0, yRotZ.z, 0.000001);

	Vector3D yRotMinusZ = Rotations.rotateVector(y, z, -1.0 * angle);

	assertEquals(y.mod(), yRotMinusZ.mod(), 0.000001);

	assertEquals(y.mod() * -1.0 * Math.sin(-1.0 * angle), yRotMinusZ.x, 0.000001);
	assertEquals(y.mod() * Math.cos(angle), yRotMinusZ.y, 0.000001);
	assertEquals(0.0, yRotMinusZ.z, 0.000001);


    }

 
    /**
     *
     * If this works, then rotating a vector onto another vector also works, and so does applying a transformation to a tensor.
     */
   public void testRotateTensorOntoVector() {
	
	DT dt = new DT(1000.0, 0.0, 0.0, 100.0, 0.0, 10.0);

	double[][] eigenSys = dt.sortedEigenSystem();

	double[] eigenvalues = new double[3];

	eigenvalues[0] = eigenSys[0][0];
	eigenvalues[1] = eigenSys[0][1];
	eigenvalues[2] = eigenSys[0][2];

	Vector3D rotationVector = Vector3D.vectorFromSPC(1.0, Math.PI / 3.0, Math.PI / 6.0);
	
	DT dtRotated = Rotations.rotateTensor(dt, rotationVector);

	double[] rotatedEVs = new double[3];

	double[][] rotESys = dtRotated.sortedEigenSystem();

	rotatedEVs[0] = rotESys[0][0];
	rotatedEVs[1] = rotESys[0][1];
	rotatedEVs[2] = rotESys[0][2];

	// eigenvalues should be the same
	assertEquals(eigenvalues[0], rotatedEVs[0], 0.00001);
	assertEquals(eigenvalues[1], rotatedEVs[1], 0.00001);
	assertEquals(eigenvalues[2], rotatedEVs[2], 0.00001);

	// e1 of the rotated tensor should lie along the rotationVector
	assertEquals(1.0, Math.abs(rotationVector.dot(new Vector3D(rotESys[1][0], rotESys[2][0], rotESys[3][0]))), 0.0001);


    }


    public void testRotationWithDoubleArg() {

	Vector3D x = new Vector3D(2.0, 0.0, 0.0);
	Vector3D y = new Vector3D(0.0, 1.0, 0.0);
	Vector3D z = new Vector3D(0.0, 0.0, 1.0);

	double angle = Math.PI / 3.0;

	// should rotate towards -Z axis
	RealMatrix rotY = Rotations.getRotMat(y, angle);
	RealMatrix rotYD = Rotations.getRotMat(new double[] {0.0, 1.0, 0.0}, angle);

	GenTestMethods.assertMatricesEqual(rotY, rotYD, 0.000001);

	// should rotate X towards Y
	RealMatrix rotZ = Rotations.getRotMat(z, angle);
	RealMatrix rotZD = Rotations.getRotMat(new double[] {0.0, 0.0, 1.0}, angle);

	GenTestMethods.assertMatricesEqual(rotZ, rotZD, 0.000001);
	
    }

  

}




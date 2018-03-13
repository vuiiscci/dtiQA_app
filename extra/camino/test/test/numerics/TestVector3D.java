package numerics;

import junit.framework.*;
import junit.extensions.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>Vector3D.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Vector3D</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestVector3D.java,v 1.1 2005/09/26 10:36:49 ucacpco Exp $
 * @author  Philip Cook
 * @see numerics.Vector3D
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestVector3D extends TestCase {


    public TestVector3D(String name) {
	super(name);
    }

    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }


      
    public static Test suite() {
	return new TestSuite(TestVector3D.class);
    }


    public void testEquals() {

	Vector3D vector1 = new Vector3D(4.0, 5.0, -6.0);

	Vector3D vector2 = new Vector3D(4.0, 5.0, -6.0);

	Vector3D vector3 = new Vector3D(vector2);

	// these should all be equal
	assertTrue( vector1.equals(vector1) );
	assertTrue( vector1.equals(vector2) );
	assertTrue( vector2.equals(vector3) );

	// Transitivity
	assertTrue( vector1.equals(vector3) );

	// Communitivity
	assertTrue( vector2.equals(vector1) );

	assertFalse( vector3.equals( new Point3D() ) );
	
	// these should be false
	Vector3D vector4 = new Vector3D(4.0, 5.0, 6.0);
	
	assertFalse( vector4.equals(vector3) );

	Vector3D vector5 = new Vector3D(3.0, 5.0, -6.0);
	assertFalse( vector5.equals(vector3) );

	Vector3D vector6 = new Vector3D(4.0, 8.0, -6.0);
	assertFalse( vector6.equals(vector3) );


    }

    
    public void testNegated() {

	Vector3D vector1 = new Vector3D(-4.0, -5.0, 6.0);

	Vector3D vector2 = vector1.negated();

	assertEquals(vector1.x, -1.0 * vector2.x, 0.00001);
	assertEquals(vector1.y, -1.0 * vector2.y, 0.00001);
	assertEquals(vector1.z, -1.0 * vector2.z, 0.00001);
	
	vector2 = vector2.negated();

	assertTrue(vector1.equals(vector2));

    }



    public void testMod() {

	Vector3D vector1 = new Vector3D(-4.0, 4.0, 2.0);

	double mod = vector1.mod();

	// mod should be 6
	assertEquals(6.0, mod, 0.001);

    }



    public void testScaled() {

	Vector3D vector1 = new Vector3D(4.0, -5.0, 6.0);

	Vector3D vector2 = vector1.scaled(2.0);

	assertEquals(2.0 * vector1.mod(), vector2.mod(), 0.0001);

	Vector3D vector3 = vector1.negated();
	Vector3D vector4 = vector1.scaled(-1.0);

	assertTrue( vector3.equals(vector4) );
		     
    }



    public void testDot() {

	Vector3D vector1 = new Vector3D(1.0, 0.0, 0.0);

	Vector3D vector2 = new Vector3D(0.0, 1.0, 0.0);
	
	// dot of orthogonal vectors should be zero
	assertEquals(0.0, vector1.dot(vector2), 0.0001);



	Vector3D vector3 = new Vector3D(4.0, 5.0, -6.0);

	Vector3D vector4 = new Vector3D(4.0, 5.0, -6.0);

	double dot = vector3.dot(vector4);

	// dot product should be |A| * |B|
	assertEquals(dot, vector3.mod() * vector4.mod(), 0.001);

	

    }

    public void testCross() {

	Vector3D vector1 = new Vector3D(2.0, 5.0, 3.0);
	
	Vector3D cross = vector1.cross(vector1);


	// cross product of vector with itself should be 0
	assertEquals(0.0, cross.mod(), 0.0001);
	

	Vector3D vector2 = new Vector3D(2.0, 0.0, 0.0);
	Vector3D vector3 = new Vector3D(0.0, 2.0, 0.0);

	Vector3D cross2 = vector2.cross(vector3);

	// dot product of cross2 should be 0 with both vectors
	assertEquals(0.0, vector2.dot(cross2), 0.0001);
	assertEquals(0.0, vector3.dot(cross2), 0.0001);

	// Actual cross should be 0,0,4
	assertEquals(0.0, cross2.x, 0.0001);
	assertEquals(0.0, cross2.y, 0.0001);
	assertEquals(4.0, cross2.z, 0.0001);
	
	

    }



    public void testConstructFromSPC() {
	

	// z axis unit vector
	Vector3D x = new Vector3D(1.0, 0.0, 0.0);
	Vector3D y = new Vector3D(0.0, 1.0, 0.0);
	Vector3D z = new Vector3D(0.0, 0.0, 1.0);

	Vector3D xy = new Vector3D(1.0, 1.0, 0.0).normalized();
	
	GenTestMethods.assertVectorsEqual(Vector3D.vectorFromSPC(1.0, 0.0, 0.0), z, 0.000001);
	
	GenTestMethods.assertVectorsEqual(Vector3D.vectorFromSPC(10.0, 0.0, 0.0), 
					  z.scaled(10.0), 0.000001);


	GenTestMethods.assertVectorsEqual(Vector3D.vectorFromSPC(10.0, Math.PI, 0.0),
					  z.scaled(-10.0), 0.000001);

	GenTestMethods.assertVectorsEqual(Vector3D.vectorFromSPC(1.0, Math.PI / 2.0, 0.0), 
					  x, 0.000001);


	GenTestMethods.assertVectorsEqual(Vector3D.vectorFromSPC(1.0, Math.PI / 2.0, Math.PI), 
					  x.scaled(-1.0), 0.000001);


	GenTestMethods.assertVectorsEqual(Vector3D.vectorFromSPC(1.0, Math.PI/2.0, Math.PI/2.0),
					  y, 0.000001);

	GenTestMethods.assertVectorsEqual(Vector3D.vectorFromSPC(1.0, Math.PI/2.0, 3.0 * Math.PI/2.0),
					  y.scaled(-1.0), 0.000001);

	
	GenTestMethods.assertVectorsEqual(Vector3D.vectorFromSPC(1.0, 
								 Math.PI / 2.0, Math.PI/4.0), 
					  xy, 0.00001);
	
	GenTestMethods.assertVectorsEqual(Vector3D.vectorFromSPC(1.0, 
								 Math.PI / 2.0, Math.PI *(1.25)), 
					  xy.scaled(-1.0), 0.00001);

	
    }


    public void testThetaPhi() {
	
	testThetaPhi(0.0, 0.0);
	testThetaPhi(Math.PI, 0.0);
	testThetaPhi(Math.PI / 4.0, 0.0);
	testThetaPhi(Math.PI / 4.0, Math.PI / 6.0);
	testThetaPhi(Math.PI * (1.0 / 2.0), Math.PI / 2.0);
	testThetaPhi(Math.PI * (3.0 / 4.0), Math.PI * (3.0 / 2.0));
	testThetaPhi(Math.PI * (1.0 / 2.0), Math.PI);
	testThetaPhi(Math.PI * (3.0 / 5.0), Math.PI * (1.0 + 0.25));
	testThetaPhi(Math.PI * (6.0 / 7.0), Math.PI * (1.0 + 0.5));
    }


    private void testThetaPhi(double theta, double phi) {
	
	double[] thetaPhi;
	
	thetaPhi = Vector3D.thetaPhi(Vector3D.vectorFromSPC(1.0, theta, phi));
	
	assertEquals(theta, thetaPhi[0], 0.000001);
	assertEquals(phi, thetaPhi[1], 0.000001);
	
    }
  
}

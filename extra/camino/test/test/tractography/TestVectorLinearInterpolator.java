package tractography;

import junit.framework.*;
import junit.extensions.*;
import misc.DT;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>VectorLinearInterpolator</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>VectorLinearInterpolator</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id$
 * @author  Philip Cook
 * @see tractography.VectorLinearInterpolator;
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestVectorLinearInterpolator extends TestCase {

    private final double voxelDim = 2.0;


    public TestVectorLinearInterpolator(String name) {
	super(name);
    }

    
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

    
    protected void tearDown() {

    }


    public static Test suite() {
	return new TestSuite(TestVectorLinearInterpolator.class);
    }


    public void testInterpolate() {

	TractographyImage image = null;
	
	VectorLinearInterpolator interpolator = null;

	Vector3D[][][][] vecs = new Vector3D[2][2][2][];

	Vector3D x = new Vector3D(1.0, 0.0, 0.0);
	Vector3D xy = new Vector3D(1.0, 1.0, 0.0).normalized();
	Vector3D yz = new Vector3D(0.0, 1.0, 1.0).normalized();
	Vector3D z = new Vector3D(0.0, 0.0, 1.0);

	vecs[0][0][0] = new Vector3D[] {x, z};
	vecs[0][1][0] = new Vector3D[] {xy.negated(), z};
	vecs[1][0][0] = new Vector3D[] {x.negated(), z.negated()};
	vecs[1][1][0] = new Vector3D[] {xy, z.negated()};

	vecs[0][0][1] = new Vector3D[] {yz, x};
	vecs[0][1][1] = new Vector3D[] {yz, xy};
	vecs[1][0][1] = new Vector3D[] {z, x};
	vecs[1][1][1] = new Vector3D[] {z, xy};
	
	image = new PD_TractographyImage(vecs, new double[] {voxelDim, voxelDim, voxelDim});

	interpolator = new VectorLinearInterpolator(image);
	
	Vector3D interp = interpolator.getTrackingDirection(new Point3D(1.0, 2.0, 1.0), x);

	assertEquals(0.9238795325112867, interp.x, 1E-9);
	assertEquals(0.3826834323650897, interp.y, 1E-9);
	assertEquals(0.0, interp.z, 1E-9);
	

	// test with x negated with respect to previous direction
	interp = interpolator.getTrackingDirection(new Point3D(3.0, 1.5, 1.0), x);

	assertEquals(0.9822902577808736, interp.x, 1E-9);
	assertEquals(0.18736555037889127, interp.y, 1E-9);
	assertEquals(0.0, interp.z, 1E-9);


	// go in the opposite direction
	interp = interpolator.getTrackingDirection(new Point3D(3.0, 1.5, 1.0), x.negated());

	assertEquals(-0.9822902577808736, interp.x, 1E-9);
	assertEquals(-0.18736555037889127, interp.y, 1E-9);
	assertEquals(0.0, interp.z, 1E-9);


	// now check that we choose correct PD in each voxel
	interp = interpolator.getTrackingDirection(new Point3D(1.0, 1.5, 2.0), z);

	assertEquals(0.0, interp.x, 1E-9);
	assertEquals(0.3826834323650897, interp.y, 1E-9);
	assertEquals(0.9238795325112867, interp.z, 1E-9);
	

        // check get first step is correct

        interp = interpolator.getTrackingDirection(new Point3D(1.0, 1.5, 2.1), 1, true);
        
        // should be equivalent to calling the same thing with previous direction == x
	
        Vector3D interpCheck = interpolator.getTrackingDirection(new Point3D(1.0, 1.5, 2.1), x);

        assertEquals(interp, interpCheck);


        interp = interpolator.getTrackingDirection(new Point3D(1.0, 1.5, 2.1), 1, false);
        
        // should be equivalent to calling the same thing with previous direction == x
	
        interpCheck = interpolator.getTrackingDirection(new Point3D(1.0, 1.5, 2.1), x.negated());

        assertEquals(interp, interpCheck);


        interp = interpolator.getTrackingDirection(new Point3D(1.0, 1.5, 2.1), 0, true);
        
        // should be equivalent to calling the same thing with previous direction == x
	
        interpCheck = interpolator.getTrackingDirection(new Point3D(1.0, 1.5, 2.1), yz);
        
        assertEquals(interp, interpCheck);

    }
	
}

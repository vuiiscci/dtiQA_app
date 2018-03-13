package tractography;

import junit.framework.*;
import junit.extensions.*;
import numerics.Point3D;
import Jama.Matrix;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code> CubicVoxelRegion</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code> CubicVoxelRegion</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestCubicVoxelRegion.java,v 1.2 2005/09/26 17:07:04 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestCubicVoxelRegion extends TestCase {


    public TestCubicVoxelRegion(String name) {
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

    public void tearDown() {

    }


    public static Test suite() {
	return new TestSuite(TestCubicVoxelRegion.class);
    }


    public void testContains() {

	
	CubicVoxelRegion region = new CubicVoxelRegion(2, 4, 3, 5, 4, 7, 1.7, 1.7, 2.5);


	for (int k = 0; k < 10; k++) {
	    for (int j = 0; j < 10; j++) {
		for (int i = 0; i < 10; i++) {

		    if (i >= 2 && i <= 4 && j >= 3 && j <= 5 && k >= 4 && k <=7) {
			assertTrue( region.containsVoxel(i,j,k) );
			assertTrue( 
				   region.containsMMPoint( new Point3D( (i + 0.5) * 1.7, 
									(j + 0.5) * 1.7, 
									(k + 0.5) * 2.5 
									) 
							   ) 
				   );
		    }
		    else {
			assertFalse( region.containsVoxel(i,j,k) );
			assertFalse( 
				   region.containsMMPoint( new Point3D( (i + 0.5) * 1.7, 
									(j + 0.5) * 1.7, 
									(k + 0.5) * 2.5 
									) 
							   ) 
				   );

		    }
		}
	    }
	}

	
    }
    


    public void testBorderBehaviour() {
	
	// The test only insists that behaviour at borders be logical and consistent

	// Example: consider the regions A (0,1,0,1,0,1) and B (2,3,0,1,0,1) (voxel dimensions 1.0 mm cubed).
	// The point (2.0, 0.5, 0.5) should be in B and not A
	// The idea is that to be "in" the region, the point must be >= the lower bound of the region, but 
	// < the upper bound. This enforces sensible behaviour, and is consistent with DWITE convenction. 

	// Further, the result should be the same whether we work in voxel or mm coordinates

	Point3D voxelPoint; 
	Point3D mmPoint; 

	CubicVoxelRegion frontLowerLeft = new CubicVoxelRegion(0,1,0,1,0,1, 1.5, 1.5, 1.5);
	CubicVoxelRegion frontLowerRight = new CubicVoxelRegion(2,3,0,1,0,1, 1.5, 1.5, 1.5);
	CubicVoxelRegion rearLowerLeft = new CubicVoxelRegion(0,1,2,3,0,1, 1.5, 1.5, 1.5);
	CubicVoxelRegion frontUpperLeft = new CubicVoxelRegion(0,1,0,1,2,3, 1.5, 1.5, 1.5);


	// X direction

	voxelPoint = new Point3D(3.0, 0.75, 0.75);
	mmPoint = new Point3D(3.0, 0.75, 0.75);

	assertFalse( frontLowerLeft.containsMMPoint(mmPoint) );
	assertTrue( frontLowerRight.containsMMPoint(mmPoint) );

	
	// Y direction

	voxelPoint = new Point3D(0.75, 3.0, 0.75);
	mmPoint = new Point3D(0.75, 3.0, 0.75);

	assertFalse( frontLowerLeft.containsMMPoint(mmPoint) );
	assertTrue( rearLowerLeft.containsMMPoint(mmPoint) );


	// Z direction

	voxelPoint = new Point3D(0.75, 0.75, 3.0);
	mmPoint = new Point3D(0.75, 0.75, 3.0);

	assertFalse( frontLowerLeft.containsMMPoint(mmPoint) );
	assertTrue( frontUpperLeft.containsMMPoint(mmPoint) );


    }

    
    public void testEquals() {

	CubicVoxelRegion region = new CubicVoxelRegion(10, 20, 10, 20, 10, 20, 1.0, 1.0, 1.0);

	CubicVoxelRegion equalRegion = new CubicVoxelRegion(10, 20, 10, 20, 10, 20, 1.0, 1.0, 1.0);

	CubicVoxelRegion diffXMinRegion = new CubicVoxelRegion(9, 20, 10, 20, 10, 20, 1.0, 1.0, 1.0);
	CubicVoxelRegion diffXMaxRegion = new CubicVoxelRegion(10, 21, 10, 20, 10, 20, 1.0, 1.0, 1.0);
	CubicVoxelRegion diffYMinRegion = new CubicVoxelRegion(10, 20, 9, 20, 10, 20, 1.0, 1.0, 1.0);
	CubicVoxelRegion diffYMaxRegion = new CubicVoxelRegion(10, 20, 10, 21, 10, 20, 1.0, 1.0, 1.0);
	CubicVoxelRegion diffZMinRegion = new CubicVoxelRegion(10, 20, 10, 20, 9, 20, 1.0, 1.0, 1.0);
	CubicVoxelRegion diffZMaxRegion = new CubicVoxelRegion(10, 20, 10, 20, 10, 21, 1.0, 1.0, 1.0);

	CubicVoxelRegion diffXVoxDimRegion = new CubicVoxelRegion(10, 20, 10, 20, 10, 20, 1.2, 1.0, 1.0);
	CubicVoxelRegion diffYVoxDimRegion = new CubicVoxelRegion(10, 20, 10, 20, 10, 20, 1.0, 1.2, 1.0);
	CubicVoxelRegion diffZVoxDimRegion = new CubicVoxelRegion(10, 20, 10, 20, 10, 20, 1.0, 1.0, 1.2);

	assertTrue( region.equals(region) );
	assertTrue( region.equals(equalRegion) );

	assertFalse( region.equals(diffXMinRegion) );
	assertFalse( region.equals(diffXMaxRegion) );
	assertFalse( region.equals(diffYMinRegion) );
	assertFalse( region.equals(diffYMaxRegion) );
	assertFalse( region.equals(diffZMinRegion) );
	assertFalse( region.equals(diffZMaxRegion) );

	assertFalse( region.equals(diffXVoxDimRegion) );
	assertFalse( region.equals(diffYVoxDimRegion) );
	assertFalse( region.equals(diffZVoxDimRegion) );



    }

    
    /**
     * Assumes containsVoxel() works.
     */
    public void testSeedPoints() {

	double xVoxelDim = 2.0;
	double yVoxelDim = 3.0;
	double zVoxelDim = 4.0;

	CubicVoxelRegion region = new CubicVoxelRegion(2, 4, 3, 5, 6, 8, xVoxelDim, yVoxelDim, zVoxelDim);
	
	Point3D[] points = region.getSeedPoints();

	int pointCounter = 0;

	for (int k = 0; k < 10; k++) {
	    for (int j = 0; j < 10; j++) {
		for (int i = 0; i < 10; i++) {

		    if (region.containsVoxel(i,j,k)) {
			assertEquals( (i + 0.5) * xVoxelDim, points[pointCounter].x, 1E-8);
			assertEquals( (j + 0.5) * yVoxelDim, points[pointCounter].y, 1E-8);
			assertEquals( (k + 0.5) * zVoxelDim, points[pointCounter].z, 1E-8);
			pointCounter++;
		    }
		    
		}
	    }
	    
	}


    } 



}

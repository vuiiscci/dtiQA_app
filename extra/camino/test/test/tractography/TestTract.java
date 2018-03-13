package tractography;

import junit.framework.*;
import junit.extensions.*;

import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>Tract</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Tract</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestTract.java,v 1.3 2005/11/08 11:47:48 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTract extends TestCase {

    // Array of points used for testing
    Point3D[] points;

    public TestTract(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

	points = new Point3D[10];

	for (int i = 0; i < 10; i++) {
	    points[i] = new Point3D((double)i, 0.0, 0.0);
	}
	
    }

    
    protected void tearDown() {
	points = null;
    }


    public static Test suite() {
	return new TestSuite(TestTract.class);
    }


    /**
     * Make sure that the points (number and value) are preserved when a Tract increases its capacity.
     *
     */
    public void testGrowth() {

	Tract t = new Tract(3, 100.0);

	for (int i = 0; i < points.length; i++) {
	    t.addPoint(points[i]);
	}

	// Adding these points should trigger the growth of the tract, twice
	// But the number of points should be as many as we added
	assertEquals(t.numberOfPoints(), 10);
	
	// check that the points have been copied correctly
	for (int i = 0; i < 10; i++) {
	    assertTrue( points[i].equals( t.getPoint(i) ) );
	}
	
    }

    public void testJoin() {
	
	Tract t1 = new Tract(100, 50.0);
	Tract t2 = new Tract(100, 50.0);

	for (int i = 0; i < points.length; i++) {
	    t1.addPoint(points[i]);
	    t2.addPoint(points[i]);
	}

	t1.joinTract(t2);

	
	// Tracts that are joined are assumed to have a common seed point
	// hence the total number of points in the tract is one less than the sum
	// of the unjoined points
	assertEquals(t1.numberOfPoints(), 2 * points.length - 1);


	for (int i = 0; i < points.length; i++) {
	    assertTrue(  t1.getPoint(i).equals( t2.getPoint(points.length - 1 - i) )  );
	}

	// the contents of points should follow, in its original order
	for (int i = 0; i < points.length - 1; i++) {
	    assertTrue(  t1.getPoint(points.length + i).equals( points[i + 1] )  );
	}


	// test that seed point is in right place
		
	Tract small = new Tract(100, 50.0);
	Tract large = new Tract(100, 50.0);

	for (int i = 0; i < points.length / 2; i++) {
	    small.addPoint(points[i]);
	}

	for (int i = 0; i < points.length; i++) {
	    large.addPoint(points[i]);
	}
	
	small.joinTract(large);

	// omit seed point of one tract so that it is not repeated
	assertEquals(points.length - 1, small.seedPointIndex());
	
	
    }



   public void testToVoxelObjList() {
	
	Tract t = new Tract(100, 100.0);


	t.addPoint(new Point3D(1.8, 2.5, 3.5));
	t.addPoint(new Point3D(1.5, 2.5, 1.5));	
	t.addPoint(new Point3D(1.2, 0.5, 0.5));
	t.addPoint(new Point3D(0.5, 0.5, 0.5));	

	VoxelList list = t.toVoxelList(1.0, 1.0, 2.0);

	assertEquals(0, list.seedPointIndex());

        Voxel[] voxelList = list.getVoxels();
        
	assertEquals(4, voxelList.length);
        
        java.util.Arrays.sort(voxelList);

	assertEquals(0, voxelList[0].x);
	assertEquals(0, voxelList[0].y);
	assertEquals(0, voxelList[0].z);

	assertEquals(1, voxelList[1].x);
	assertEquals(0, voxelList[1].y);
	assertEquals(0, voxelList[1].z);

	assertEquals(1, voxelList[2].x);
	assertEquals(2, voxelList[2].y);
	assertEquals(0, voxelList[2].z);

	assertEquals(1, voxelList[3].x);
	assertEquals(2, voxelList[3].y);
	assertEquals(1, voxelList[3].z);
	    
   }


    public void testPathLengthFromSeed() {

        // tract length smaller than number of points we will have
        // this ensures that displacements are copied correctly on growth.
	Tract t1 = new Tract(6, 100.0);
	Tract t2 = new Tract(6, 100.0);

        Point3D seedPoint = new Point3D(10.0, 10.0, 10.0);
        
        t1.addPoint(seedPoint);
        t2.addPoint(seedPoint);
        
        Point3D currentPos1 = seedPoint;        
        Point3D currentPos2 = seedPoint;

        Vector3D e1 = new Vector3D(0.5, 0.4, 0.2).normalized();

 	for (int i = 0; i < 10; i++) {

            // move twice as far on one side
            // make sure displacement is right
            currentPos1 = currentPos1.displace(e1.scaled(0.5));
            currentPos2 = currentPos2.displace(e1.negated().scaled(0.25));

            t1.addPoint(currentPos1);
            t2.addPoint(currentPos2);
	}

	t1.joinTract(t2);

        for (int i = 0; i < 10; i++) {
            assertEquals(0.25 * (10 - i), t1.pathLengthFromSeed(i), 1E-6);
        }
        
        assertEquals(0.0, t1.pathLengthFromSeed(10), 1E-6);
        
        for (int i = 0; i < 10; i++) {
            assertEquals(0.5 * i, t1.pathLengthFromSeed(10+i), 1E-6);
        }

        
    }


    public void testLength() {

	Tract t1 = new Tract(20, 100.0);

	java.util.Random ran = new java.util.Random(12345);

	double length = 0.0;

	Point3D currentPos = new Point3D(100.0, 100.0, 100.0);

	t1.addPoint(currentPos, 0.0);

        for (int i = 0; i < 10; i++) {

	    double step = ran.nextDouble();

	    length += step;
	    
	    double x = ran.nextGaussian();
	    double y = ran.nextGaussian();
	    double z = ran.nextGaussian();

	    double mod = Math.sqrt(x*x + y*y + z*z);

	    x /= mod;
	    y /= mod;
	    z /= mod;

	    currentPos = currentPos.displace(new Vector3D(x,y,z).scaled(step));

	    t1.addPoint(currentPos, step);
	}

	assertEquals(length, t1.length(), 1E-6);

    }
  

    public void testResample() {

	Tract t = new Tract();
	Tract t2 = new Tract();

	for (int i = 0; i < 6; i++) {
	    
	    t2.addPoint(new Point3D(10 - 2*i, 0.0, 0.0));
	    t.addPoint(new Point3D(10 + 2*i, 0.0, 0.0));
	}

	t.joinTract(t2);

	Tract resampled = t.resample(1.0);

	assertEquals(21, resampled.numberOfPoints());
	
	assertEquals(10, resampled.seedPointIndex());

	for (int i = 0; i < 21; i++) {
	    
	    assertEquals(i, resampled.getPoint(i).x, 1E-6);
	    assertEquals(0.0, resampled.getPoint(i).y, 1E-6);
	    assertEquals(0.0, resampled.getPoint(i).z, 1E-6);
	}


	resampled = t.resample(0.4);

	assertEquals(0.4, resampled.getPoint(1).x, 1E-6);
	assertEquals(0.8, resampled.getPoint(2).x, 1E-6);
	assertEquals(1.2, resampled.getPoint(3).x, 1E-6);
	assertEquals(1.6, resampled.getPoint(4).x, 1E-6);
	assertEquals(2.0, resampled.getPoint(5).x, 1E-6);

	resampled = t.resample(0.9999);
	
	assertEquals(21, resampled.numberOfPoints());

    }


    public void testChop() {
	
	Tract t = new Tract();
	Tract t2 = new Tract();

	for (int i = 0; i < 6; i++) {
	    
	    t.addPoint(new Point3D(5 - i, 0.0, 0.0));
	    t2.addPoint(new Point3D(5 + i, 0.0, 0.0));
	}

	t.joinTract(t2);

	t.chop(4,6);
	
	assertEquals(3, t.numberOfPoints());
	assertEquals(1, t.seedPointIndex());
	
	assertEquals(new Point3D(6.0, 0.0, 0.0), t.getPoint(0));
	assertEquals(new Point3D(5.0, 0.0, 0.0), t.getPoint(1));
	assertEquals(new Point3D(4.0, 0.0, 0.0), t.getPoint(2));
    }


    public void testEquals() {

	Tract t = new Tract();
	Tract t2 = new Tract();

	for (int i = 0; i < 6; i++) {
	    
	    t.addPoint(new Point3D(5 + i, 0.0, 0.0));
	    t2.addPoint(new Point3D(5 + i, 0.0, 0.0));

	}

	assertTrue(t.equals(t2));

	t2.addPoint(new Point3D(11.0, 0.0, 0.0));

	assertFalse(t2.equals(t));
    }


    public void testVoxelPathLengths() {
	
	double[] voxelDims = new double[] {2.0, 2.5, 3.0};

	Tract t = new Tract();
	Tract t2 = new Tract();

	// tract runs from (1.5,0,0) to (8.5,0,0)
	for (int i = 0; i < 5; i++) {
	    
	    t.addPoint(new Point3D(5 + i + 0.5, 1.0, 1.0));
	    t2.addPoint(new Point3D(5 - i + 0.5, 1.0, 1.0));
	}

	t.addPoint(new Point3D(9.5, 3.0, 4.0)); 

	t.joinTract(t2);

	double[] pathLengths = t.getVoxelPathLengths(voxelDims);

	Voxel[] voxels = t.toVoxelList(voxelDims).getVoxels();
	
	assertEquals(voxels.length, pathLengths.length);

	assertEquals(0.5, pathLengths[0], 1E-8);

	for (int i = 1; i < 4; i++) {
	    // straight through these voxels
	    assertEquals(2.0, pathLengths[i], 1E-8);
	}

	// vector from penultimate point in voxel (4,0,0) to last point in voxel (4,1,1)
	Vector3D vec = new Vector3D(t.getPoint(9), t.getPoint(8));

	Vector3D uvec = vec.normalized();

	// path length in voxel (4,0,0) is 1.5 (distance from point 8 to the boundary with (3,0,0)) plus
	// the portion of the last segment that is in voxel (4,0,0)
	assertEquals(1.5 + 2.0 / uvec.z, pathLengths[4], 1E-8);

	// Rest of the segment, including part that is actually in (4,0,1), is assigned to (4,1,1)
	assertEquals(vec.mod() - 2.0 / uvec.z, pathLengths[5], 1E-8);


	Tract tr = t.resample(0.1);

	pathLengths = tr.getVoxelPathLengths(voxelDims);

	voxels = tr.toVoxelList(voxelDims).getVoxels();
	
	assertEquals(voxels.length, pathLengths.length);

	// first part the same as before
	assertEquals(0.5, pathLengths[0], 1E-8);

	for (int i = 1; i < 4; i++) {
	    assertEquals(2.0, pathLengths[i], 1E-8);
	}

	
	// path length in voxel (4,0,0) is 1.5 (distance from point 8 to the boundary with (3,0,0)) plus
	// the portion of the next segment that is in voxel (4,0,0)
	assertEquals(1.5 + 2.0 / uvec.z, pathLengths[4], 1E-8);

	// rest of the above segment is in voxel (4,0,1)
	assertEquals(new Voxel(4,0,1), voxels[5]);

	// point of intersection of tract with voxel boundary between (4,0,0) and (4,0,1)
	Point3D intersectionPointZ = t.getPoint(8).displace(uvec.scaled(2.0 / uvec.z));

	// point of intersection of tract with voxel boundary between (4,0,1) and (4,1,1)
	Point3D intersectionPointY = t.getPoint(8).displace(uvec.scaled(1.5 / uvec.y));

	// segment in voxel (4,0,1)
	assertEquals(new Vector3D(intersectionPointY, intersectionPointZ).mod(), pathLengths[5], 1E-8);

	// segment in (4,1,1) is everything from the intersectionPointY to the end
	assertEquals(new Vector3D(tr.getPoint(tr.numberOfPoints() - 1), intersectionPointY).mod(), 
		     pathLengths[6], 1E-8);

    }



    public void testTransforms() {

        Nifti1Dataset nds = new Nifti1Dataset();

        nds.setPixDims(1.0f, 2.0f, 3.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f);
        
        nds.setQuaternion((short)1, (short)1, new float[] {0.25f, 0.5f, 0.25f}, new float[] {100.0f, 200.0f, 300.0f});


        RealMatrix trans = nds.getVoxelToPhysicalTransform();

        Tract t = new Tract();

        t.addPoint(new Point3D(10.0, 20.0, 30.0));

        t.transformToPhysicalSpace(trans, 2.0, 3.0, 4.0);

        /*

          A =
          
          0.75000    -0.43585     3.66228   100.00000
          1.29057     2.25000    -0.58114   200.00000
          -1.33114     1.93585     1.50000   300.00000
          0.00000     0.00000     0.00000     1.00000
          
          octave-3.4.0:70> A * [10 / 2 - 0.5; 20 / 3 - 0.5; 30 / 4 - 0.5; 1]ans =
          
          126.3232
          215.6146
          316.4476
          1.0000
          
        */
        
        Point3D p = t.getPoint(0);

        assertEquals(126.3232, p.x, 1E-4);
        assertEquals(215.6146, p.y, 1E-4);
        assertEquals(316.4476, p.z, 1E-4);

        t.transformToCaminoSpace(nds.getPhysicalToVoxelTransform(), 2.0, 3.0, 4.0);

        
        p = t.getPoint(0);

        assertEquals(10.0, p.x, 1E-4);
        assertEquals(20.0, p.y, 1E-4);
        assertEquals(30.0, p.z, 1E-4);

        
    }
    
}

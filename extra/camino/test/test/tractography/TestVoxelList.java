package tractography;

import junit.framework.*;
import junit.extensions.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>VoxelList</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>VoxelList</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestVoxelList.java,v 1.2 2005/11/08 11:47:48 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestVoxelList extends TestCase {

  
    public TestVoxelList(String name) {
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

    
    protected void tearDown() {

    }


    public static Test suite() {
	return new TestSuite(TestVoxelList.class);
    }


    public void testToTract() {
	
	Voxel[] voxels = new Voxel[4];

	voxels[0] = new Voxel(0, 0, 0);
	voxels[1] = new Voxel(1, 0, 0);
	voxels[2] = new Voxel(1, 1, 0);
	voxels[3] = new Voxel(1, 1, 1);

	VoxelList list = new VoxelList(voxels, 2, 1.0, 2.0, 3.0, 
				       new Vector3D(0.0, 0.0, 1.0), new Vector3D(0.0, -1.0, 0.0));

	Tract t = list.toTract();

	assertEquals(6, t.numberOfPoints());
	assertEquals(3, t.seedPointIndex());
	
	assertEquals(new Point3D(0.5, 1.0, 1.5), t.getPoint(0));
	assertEquals(new Point3D(1.5, 1.0, 1.5), t.getPoint(1));
	assertEquals(new Point3D(1.5, 2.99, 1.5), t.getPoint(2));
	assertEquals(new Point3D(1.5, 3.0, 1.5), t.getPoint(3));
	assertEquals(new Point3D(1.5, 3.0, 1.51), t.getPoint(4));
	assertEquals(new Point3D(1.5, 3.0, 4.5), t.getPoint(5));

    }

}

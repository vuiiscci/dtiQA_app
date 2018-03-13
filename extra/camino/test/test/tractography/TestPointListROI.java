package tractography;

import junit.framework.*;
import junit.extensions.*;

import imaging.*;
import numerics.*;
import tools.*;


/**
 * Tests for <code>PointListROI</code>.
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestPointListROI extends TestCase {


    private final double xVoxelDim = 2.0;
    private final double yVoxelDim = 3.0;
    private final double zVoxelDim = 4.0;

    private final int xDataDim = 10;
    private final int yDataDim = 10;
    private final int zDataDim = 10;

    

    private final Point3D[] points;


    public TestPointListROI(String name) {
	super(name);

	points = new Point3D[37];

	int pointCounter = 0;

	for (short k = 0; k < zDataDim; k++) {
	    for (short j = 0; j < yDataDim; j++) {
		for (short i = 0; i < xDataDim; i++) {
                    
		    if (i >= 2 && i <= 4 && j >= 3 && j <= 5 && k >= 4 && k <=7) {
                        points[pointCounter++] = new Point3D(xVoxelDim * (i + 0.5), yVoxelDim * (j + 0.5),  
							     zVoxelDim * (k + 0.5));
                    }
                }
            }
        }

	// add duplicate point in 2,3,4
	points[36] = new Point3D(6.0, 10.0, 17.0);
	

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
	return new TestSuite(TestPointListROI.class);
    }


    public void testSeedPoints() {


	PointListROI region = new PointListROI(points, 1);
	
	Point3D[] pointList = region.getSeedPoints();

	int pointCounter = 0;

	Voxel[] voxels = region.getSeedVoxels(xVoxelDim, yVoxelDim, zVoxelDim);

	for (int i = 0; i < points.length; i++) {
	    assertEquals(points[i], pointList[i]);

	    Voxel v = new Voxel( (int)(points[i].x / xVoxelDim), (int)(points[i].y / yVoxelDim), 
				 (int)(points[i].z / zVoxelDim));

	    assertEquals(v, voxels[i]);
	}
	
    } 

    
    public void testReadPoints() {

	FileOutput out = new FileOutput("roiPoints.txt");
	
	out.writeString("1.0 1.5 2.0\n"); // centre of voxel 0,0,0
	out.writeString("3.0 4.5 6.0\n"); // centre of voxel 1,1,1
	
	out.close();

        Nifti1Dataset ih = new Nifti1Dataset();

        ih.setPixDims(1.0f, 2.0f, 3.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f);

        // align Camino space with physical space
        ih.setQuaternion((short)1, (short)1, new float[] {0.0f, 0.0f, 0.0f}, new float[] {1.0f, 1.5f, 2.0f});

        
        // // DEBUG

        // RealMatrix trans = ih.getPhysicalToVoxelTransform();

        // Point3D p = new Point3D(1.0, 1.5, 2.0);
        
        // p = p.transform(trans);

        // System.err.println("\nBORK " + trans + "\n");

        // // /DEBUG


	PointListROI region = PointListROI.readPoints("roiPoints.txt", ih);

	Point3D[] pointList = region.getSeedPoints();

	assertEquals(2, pointList.length);

	assertEquals(1.0, pointList[0].x, 1E-6);
	assertEquals(1.5, pointList[0].y, 1E-6);
	assertEquals(2.0, pointList[0].z, 1E-6);

	assertEquals(3.0, pointList[1].x, 1E-6);
	assertEquals(4.5, pointList[1].y, 1E-6);
	assertEquals(6.0, pointList[1].z, 1E-6);
        
    }

}

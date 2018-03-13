package tractography;

import junit.framework.*;
import junit.extensions.*;
import numerics.Point3D;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>EightNeighbourInterpolator</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>EightNeighbourInterpolator</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestEightNeighbourInterpolator.java,v 1.1 2005/09/26 10:36:53 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestEightNeighbourInterpolator extends TestCase {

  
    public TestEightNeighbourInterpolator(String name) {
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
	return new TestSuite(TestEightNeighbourInterpolator.class);
    }


    public void testInterpolation() {

	int xDataDim = 10;
	int yDataDim = 10;
	int zDataDim = 10;

	double xVoxelDim = 2.0;
	double yVoxelDim = 2.0;
	double zVoxelDim = 2.0;
	
	
	EightNeighbourInterpolator interpolator = new EightNeighbourInterpolator(xDataDim,
										 yDataDim,
										 zDataDim,
										 xVoxelDim,
										 yVoxelDim,
										 zVoxelDim);

	int[] dims = new int[6];
	double[] components = new double[8];
	
	Point3D point = new Point3D(6.0 * xVoxelDim, 4.9 * yVoxelDim, 2.7 * zVoxelDim);

	interpolator.setInterpolationVoxels(point,components,dims); 

	// should be interpolating between
	// (5,4,2), (6,4,2), (5,5,2), (6,5,2)
	// (5,4,3), (6,4,3), (5,5,3), (6,5,3)

	assertEquals(5, dims[0]);
	assertEquals(6, dims[1]);
	assertEquals(4, dims[2]);
	assertEquals(5, dims[3]);
	assertEquals(2, dims[4]);
	assertEquals(3, dims[5]);

	// linear interp components
	assertEquals(0.5 * 0.6 * 0.8, components[0], 0.00000001);
	assertEquals(0.5 * 0.6 * 0.2, components[1], 0.00000001);
	assertEquals(0.5 * 0.4 * 0.8, components[2], 0.00000001);
	assertEquals(0.5 * 0.4 * 0.2, components[3], 0.00000001);
	assertEquals(0.5 * 0.6 * 0.8, components[4], 0.00000001);
	assertEquals(0.5 * 0.6 * 0.2, components[5], 0.00000001);
	assertEquals(0.5 * 0.4 * 0.8, components[6], 0.00000001);
	assertEquals(0.5 * 0.4 * 0.2, components[7], 0.00000001);
	

	// check boundaries
	point = new Point3D((xDataDim - 0.1) * xVoxelDim, 4.9 * yVoxelDim, 2.7 * zVoxelDim);
	
	interpolator.setInterpolationVoxels(point,components,dims); 
	
	assertEquals(xDataDim - 1, dims[0]);
	assertEquals(xDataDim - 1, dims[1]);
	assertEquals(4, dims[2]);
	assertEquals(5, dims[3]);
	assertEquals(2, dims[4]);
	assertEquals(3, dims[5]);

	point = new Point3D(xVoxelDim / 10.0, 4.9 * yVoxelDim, 2.7 * zVoxelDim);
	
	interpolator.setInterpolationVoxels(point,components,dims); 
	
	assertEquals(0, dims[0]);
	assertEquals(0, dims[1]);
	assertEquals(4, dims[2]);
	assertEquals(5, dims[3]);
	assertEquals(2, dims[4]);
	assertEquals(3, dims[5]);

	
	point = new Point3D(6.0 * xVoxelDim, (yDataDim - 0.1) * yVoxelDim, 2.7 * zVoxelDim);

	interpolator.setInterpolationVoxels(point,components,dims); 

	assertEquals(5, dims[0]);
	assertEquals(6, dims[1]);
	assertEquals(yDataDim-1, dims[2]);
	assertEquals(yDataDim-1, dims[3]);
	assertEquals(2, dims[4]);
	assertEquals(3, dims[5]);


	point = new Point3D(5.5 * xVoxelDim, yVoxelDim / 10.0, 2.7 * zVoxelDim);

	interpolator.setInterpolationVoxels(point,components,dims); 

	assertEquals(5, dims[0]);
	assertEquals(6, dims[1]);
	assertEquals(0, dims[2]);
	assertEquals(0, dims[3]);
	assertEquals(2, dims[4]);
	assertEquals(3, dims[5]);


	point = new Point3D(6.0 * xVoxelDim, 4.9 * yVoxelDim, (zDataDim - 0.1) * zVoxelDim);

	interpolator.setInterpolationVoxels(point,components,dims); 

	assertEquals(5, dims[0]);
	assertEquals(6, dims[1]);
	assertEquals(4, dims[2]);
	assertEquals(5, dims[3]);
	assertEquals(zDataDim - 1, dims[4]);
	assertEquals(zDataDim - 1, dims[5]);


	point = new Point3D(5.5 * xVoxelDim, 4.9 * yVoxelDim, zVoxelDim / 10.0);

	interpolator.setInterpolationVoxels(point,components,dims); 

	assertEquals(5, dims[0]);
	assertEquals(6, dims[1]);
	assertEquals(4, dims[2]);
	assertEquals(5, dims[3]);
	assertEquals(0, dims[4]);
	assertEquals(0, dims[5]);


    }

    

}

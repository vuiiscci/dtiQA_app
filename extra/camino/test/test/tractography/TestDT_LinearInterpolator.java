package tractography;

import junit.framework.*;
import junit.extensions.*;
import misc.DT;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>DT_LinearInterpolator</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>DT_LinearInterpolator</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestDT_LinearInterpolator.java,v 1.1 2005/09/26 10:36:52 ucacpco Exp $
 * @author  Philip Cook
 * @see tractography.DT_LinearInterpolator;
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestDT_LinearInterpolator extends TestCase {

    protected final DT_TractographyImage image = Images.getTwoTensorCube();
    protected final DT_LinearInterpolator interpolator = new DT_LinearInterpolator(image);


    public TestDT_LinearInterpolator(String name) {
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
	return new TestSuite(TestDT_LinearInterpolator.class);
    }


    /**
     * Test that interpolation is weighted at voxel centres, rather than edges
     */
    public void testInterpolationCentre() {


	// The point of interpolation is the centre of voxel (1,1,1)
	DT interpolated = interpolator.getDT(Images.getPointAtVoxelCentre(1,1,1), Rotations.X_AXIS );
	
	Images.assertTensorsEqual(this, interpolated, image.getDTs(1,1,1)[0], 1E-12);


    }



    public void testInterpolationInEachDimension() {
	testInterpolationInEachDimension(0, Rotations.X_AXIS);
	testInterpolationInEachDimension(1, Rotations.Y_AXIS);
    }

    
    /**
     * Test linear interpolation in x, y, z directions
     */
    private void testInterpolationInEachDimension(int dtIndex, Vector3D previousDirection) {
	
	
	// interpolate at the centre of the boundary between (0,1,1) and (1,1,1)
       	DT interpByInterpolator = interpolator.getDT( Images.getPointAtVoxel(1, 1.5, 1.5), previousDirection);
	
	// Should be linear interpolation of sphere and ellipsoid
	DT interpByTest = oneDimInterpolate(0,1,1, 1,1,1, dtIndex, 0.5);

	Images.assertTensorsEqual(this, interpByInterpolator, interpByTest, 1E-12);

	// Check the same thing works when we're not halfway between the two
	interpByInterpolator = interpolator.getDT( Images.getPointAtVoxel(1.25, 1.5, 1.5),  previousDirection);

	interpByTest = oneDimInterpolate(0,1,1, 1,1,1, dtIndex, 0.75);

	Images.assertTensorsEqual(this, interpByInterpolator, interpByTest, 1E-12);


	// Now repeat in the Y direction
	// interpolate at the centre of the boundary between (1,0,1) and (1,1,1)
       	interpByInterpolator = interpolator.getDT( Images.getPointAtVoxel(1.5, 1, 1.5), previousDirection);

	interpByTest = oneDimInterpolate(1,0,1, 1,1,1, dtIndex, 0.5);

	Images.assertTensorsEqual(this, interpByInterpolator, interpByTest, 1E-12);


	// Now repeat in the Z direction
	// interpolate at the centre of the boundary between (1,1,0) and (1,1,1)
       	interpByInterpolator = interpolator.getDT( Images.getPointAtVoxel(1.5, 1.5, 1), previousDirection);

	interpByTest = oneDimInterpolate(1,1,0, 1,1,1, dtIndex, 0.5);

	Images.assertTensorsEqual(this, interpByInterpolator, interpByTest, 1E-12);

    }


    /**
     * Test that dataset interpolates properly at boundaries
     */
    public void testInterpolationAtDatasetBoundaries() { 
	
       	DT interpByInterpolator;
	DT interpByTest;
	
	// at boundaries, the data should interpolate with itself


	// interpolate within voxel (0,0,0)
       	interpByInterpolator = interpolator.getDT( Images.getPointAtVoxel(0.25, 0.5, 0.5), Rotations.X_AXIS);
	interpByTest = oneDimInterpolate(0,0,0, 0,0,0, 0, 0.5);

	Images.assertTensorsEqual(this, interpByInterpolator, interpByTest, 1E-12);


       	interpByInterpolator = interpolator.getDT( Images.getPointAtVoxel(0.5, 0.25, 0.5), Rotations.X_AXIS);
	// interpByTest will be the same thing
	Images.assertTensorsEqual(this, interpByInterpolator, interpByTest, 1E-12);


       	interpByInterpolator = interpolator.getDT(Images.getPointAtVoxel(0.5, 0.5, 0.25), Rotations.X_AXIS);
	// interpByTest will be the same thing
	Images.assertTensorsEqual(this, interpByInterpolator, interpByTest, 1E-12);

	
    }       



    /**
     * Interpolate in 1D between two voxels.
     * @param distance the distance along the axis between the two voxel centres, in voxel coordinates 
     * (ie in range 0.0 to 1.0, where 0.0 is the centre of voxel 1).
     *
     */
    private DT oneDimInterpolate(int x1, int y1, int z1, int x2, int y2, int z2, int dtIndex, 
				 double distance) {

	double[] comps = image.getDTs(x1,y1,z1)[dtIndex].getComponents();

	double d1xx = comps[0];
	double d1xy = comps[1];
	double d1xz = comps[2];
	double d1yy = comps[3];
	double d1yz = comps[4];
	double d1zz = comps[5];

	comps = image.getDTs(x2,y2,z2)[dtIndex].getComponents();

	double d2xx = comps[0];
	double d2xy = comps[1];
	double d2xz = comps[2];
	double d2yy = comps[3];
	double d2yz = comps[4];
	double d2zz = comps[5];


	double interpXX = (1.0 - distance) * d1xx + distance * d2xx;

	double interpXY = (1.0 - distance) * d1xy + distance * d2xy;

	double interpXZ = (1.0 - distance) * d1xz + distance * d2xz;
	
	double interpYY = (1.0 - distance) * d1yy + distance * d2yy;
	
	double interpYZ = (1.0 - distance) * d1yz + distance * d2yz;

	double interpZZ = (1.0 - distance) * d1zz + distance * d2zz;

	return new DT(interpXX, interpXY, interpXZ, interpYY, interpYZ, interpZZ);

    }
	
}

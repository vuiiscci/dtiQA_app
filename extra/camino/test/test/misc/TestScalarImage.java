package misc;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;

/**
 * Automated tests for ScalarImage.
 *
 * @version $Id$
 * @author  Philip Cook
 * @see misc.ScalarImage
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestScalarImage extends TestCase {

    
    private ScalarImage image = null;

    private double[][][] data = null;

    private double[] voxelDims = new double[] {1.0, 2.0, 3.0};

    public TestScalarImage(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

	data = new double[2][2][2];

	for (int i = 0; i < 2; i++) {
	    for (int j = 0; j < 2; j++) {	
		for (int k = 0; k < 2; k++) {
		    data[i][j][k] = 1.0 + i * 2.0 + j * 4.0 + k * 16.0;
		}
	    }
	}


	image = new ScalarImage(data, voxelDims);


    }
    
    protected void tearDown() {
	image = null;
	data = null;
    }
    
    public static Test suite() {
	return new TestSuite(TestScalarImage.class);
    }


    public void testNN_Interpolation() {

	image.setInterpolation("nearestneighbour");
	
	double[] xPos = new double[] {0.5, 1.5};
	double[] yPos = new double[] {1.0, 3.0};
	double[] zPos = new double[] {1.5, 4.5};

	
	for (int i = 0; i < 2; i++) {
	    for (int j = 0; j < 2; j++) {	
		for (int k = 0; k < 2; k++) {
		    assertEquals(data[i][j][k], image.valueAt(new Point3D(xPos[i], yPos[j], zPos[k])), 1E-8);
		}
	    }
	}

    }


    public void testNN_DerivInterpolation() {

	image.setInterpolation("nearestneighbour");
	
	double[] xPos = new double[] {0.5, 1.5};
	double[] yPos = new double[] {1.0, 3.0};
	double[] zPos = new double[] {1.5, 4.5};

	
	for (int i = 0; i < 2; i++) {
	    for (int j = 0; j < 2; j++) {	
		for (int k = 0; k < 2; k++) {
                    double testDerivX = 2.0/voxelDims[0];
                    double testDerivY = 4.0/voxelDims[1];
                    double testDerivZ = 16.0/voxelDims[2];

                    double[] deriv = image.derivAt(new Point3D(xPos[i], yPos[j], zPos[k]));

		    assertEquals((i==0)?testDerivX:0.0, deriv[0], 1E-8);
		    assertEquals((j==0)?testDerivY:0.0, deriv[1], 1E-8);
		    assertEquals((k==0)?testDerivZ:0.0, deriv[2], 1E-8);
		}
	    }
	}
    }


    public void testLinearInterpolation() {

	image.setInterpolation("linear");

	Point3D point = new Point3D(voxelDims[0], voxelDims[1] / 2.0, voxelDims[2] / 2.0);

	assertEquals(0.5 * data[1][0][0] + 0.5 * data[0][0][0], image.valueAt(point), 1E-8);

	point = new Point3D(voxelDims[0] / 2, voxelDims[1], voxelDims[2] / 2.0);

	assertEquals(0.5 * data[0][1][0] + 0.5 * data[0][0][0], image.valueAt(point), 1E-8);

	point = new Point3D(voxelDims[0] / 2, voxelDims[1] / 2.0, voxelDims[2]);

	assertEquals(0.5 * data[0][0][1] + 0.5 * data[0][0][0], image.valueAt(point), 1E-8);
	
    }


    public void testLinearDerivInterpolation() {

	image.setInterpolation("linear");

        double testDerivX = 2.0/voxelDims[0];
        double testDerivY = 4.0/voxelDims[1];
        double testDerivZ = 16.0/voxelDims[2];

	Point3D point = new Point3D(voxelDims[0], voxelDims[1] / 2.0, voxelDims[2] / 2.0);
        double[] deriv = image.derivAt(point);
	assertEquals(testDerivX, deriv[0], 1E-8);
	assertEquals(testDerivY, deriv[1], 1E-8);
	assertEquals(testDerivZ, deriv[2], 1E-8);

	point = new Point3D(voxelDims[0] / 2, voxelDims[1], voxelDims[2] / 2.0);
        deriv = image.derivAt(point);
	assertEquals(testDerivX, deriv[0], 1E-8);
	assertEquals(testDerivY, deriv[1], 1E-8);
	assertEquals(testDerivZ, deriv[2], 1E-8);

	point = new Point3D(voxelDims[0] / 2, voxelDims[1] / 2.0, voxelDims[2]);
        deriv = image.derivAt(point);
	assertEquals(testDerivX, deriv[0], 1E-8);
	assertEquals(testDerivY, deriv[1], 1E-8);
	assertEquals(testDerivZ, deriv[2], 1E-8);
	

    }


}

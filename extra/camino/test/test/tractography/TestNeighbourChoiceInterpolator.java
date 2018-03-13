package tractography;

import junit.framework.*;
import junit.extensions.*;

import data.*;
import misc.DT;
import numerics.*;

import java.util.Random;


/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>NeighbourChoiceInterpolator</code>.
 * <BR><BR>
  * </dl>
 *
 * @version $Id: TestPICoNeighbourChoiceInterpolator.java,v 1.3 2006/05/25 17:38:14 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestNeighbourChoiceInterpolator extends TestCase {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestNeighbourChoiceInterpolator(String name) {
	super(name);
	
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestNeighbourChoiceInterpolator.class);
    }


    public void testNeighbourChoice() {

	TractographyImage image = Images.getCube();

	Vector3D[][][] e1s = new Vector3D[2][2][2];

	for (int k = 0; k < 2; k++) {
	    for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 2; i++) {
		    e1s[i][j][k] = image.getPDs(i,j,k)[0];
		}
	    }
	}


	// test in cube
	long seed = 562l;

	Random ran = new Random(seed);

	NeighbourChoiceInterpolator interpolator = 
	    new NeighbourChoiceInterpolator(image, ran);

	// put point such that neighbour selection should be 60-40 between x|y|z and x|y|z+1

	Point3D point = new Point3D(0.9*Images.xVoxelDim, 0.9*Images.yVoxelDim, 0.9*Images.zVoxelDim);

	int[][][] voxelSelections = new int[2][2][2];

	// make 100000 interpolation steps; check that they are distributed approximately correctly
	for (int n = 0; n < 100000; n++) {
	    Vector3D td = interpolator.getTrackingDirection(point, e1s[0][0][0]);
	    for (int k = 0; k < 2; k++) {
		for (int j = 0; j < 2; j++) {
		    for (int i = 0; i < 2; i++) {
			if (Math.abs(e1s[i][j][k].dot(td)) > 0.9999999999) {
			    voxelSelections[i][j][k]++;
			}
		    }
		}
	    }
	}


	// make sure that direction is correct
	Point3D voxelCentre = new Point3D(0.5*Images.xVoxelDim, 0.5*Images.yVoxelDim, 0.5*Images.zVoxelDim);
	assertEquals(-1.0, interpolator.getTrackingDirection(voxelCentre, e1s[0][0][0]).dot(interpolator.getTrackingDirection(voxelCentre, e1s[0][0][0].negated())), 0.00000001);

	int counter = 0;

	for (int k = 0; k < 2; k++) {
	    for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 2; i++) {
		    counter += voxelSelections[i][j][k];
		}
	    }
	}

	assertEquals(100000, counter);

	// assert equality to within 5%

	assertEquals(14400.0, voxelSelections[0][0][1], 720.0);
	assertEquals(9600.0, voxelSelections[0][1][1], 480.0);
	assertEquals(9600.0, voxelSelections[1][0][1], 480.0);
	assertEquals(6400.0, voxelSelections[1][1][1], 320.0);

	assertEquals(21600.0, voxelSelections[0][0][0], 1008.0);
	assertEquals(14400.0, voxelSelections[0][1][0], 720.0);
	assertEquals(14400.0, voxelSelections[1][0][0], 720.0);
	assertEquals(9600.0, voxelSelections[1][1][0], 320.0);

    }



    public void testPDChoice() {

	// test in cube
	long seed = 562l;

	Random ran = new Random(seed);

	TractographyImage image = Images.getTwoTensorCube();
	
	NeighbourChoiceInterpolator interpolator = 
	    new NeighbourChoiceInterpolator(image, ran);

	// place point in centre
	Point3D point = new Point3D(Images.xVoxelDim, Images.yVoxelDim, Images.zVoxelDim);

	Vector3D[] pds = image.getPDs(0,0,0);

	boolean pd1Chosen = false;

        Vector3D td = interpolator.getTrackingDirection(point, pds[0]);
	
        if (Math.abs(td.dot(pds[0])) > 0.99999999) {
            pd1Chosen = true;
        }
	
        assertTrue(pd1Chosen);

    }


   
}

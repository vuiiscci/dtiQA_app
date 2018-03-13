package tractography;

import junit.framework.*;
import junit.extensions.*;

import data.*;
import misc.DT;
import numerics.*;


/**
 * Tests for RK4 tracker
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestRK4FibreTracker extends TestFibreTracker {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestRK4FibreTracker(String name) {
	super(name);
	
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestRK4FibreTracker.class);
    }



 	
    protected FibreTracker getTracker(TractographyImage image) {
	return new RK4FibreTracker(new VectorLinearInterpolator(image), 0.5);
    }



    /**
     * Ensure we get the correct first / next step at various points.
     *
     */
    public void testStep() {

        DT_TractographyImage cube = Images.getCube();

        RK4FibreTracker tracker = new RK4FibreTracker(new DT_LinearInterpolator(cube), 0.5);

        // in the middle
        Point3D point = new Point3D(Images.xVoxelDim, Images.yVoxelDim, Images.zVoxelDim);

        Vector3D step = tracker.getNextStep(point, new Vector3D(1.0, 0.0, 0.0));

        assertEquals(0.5, step.mod(), 1E-6);

        // at the edge

        point = new Point3D(0.001, 0.001, 0.001);

        step = tracker.getNextStep(point, new Vector3D(1.0, 0.0, 0.0));

        assertEquals(0.5, step.mod(), 1E-6);

    }


    public void testBG_Interp() {
        
        DT_TractographyImage cube = new DT_TractographyImage(Images.getCube());
        
        double[][][] scalar = new double[2][2][2];
        
        
        for (int k = 0; k < 2; k++) {
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 2; i++) {
                    scalar[i][j][k] = 1.0;
                }
            }
        }
        
        scalar[1][0][0] = 0.0;
        
        cube.computeIsotropicMask(scalar, 0.5);

        cube.numPDs[1][0][0] = 0;

        RK4FibreTracker tracker = new RK4FibreTracker(new DT_LinearInterpolator(cube), 1.5);

        Point3D point = new Point3D(1.1, 0.5, 0.5);

        Vector3D step = tracker.getNextStep(point, new Vector3D(1.0, 0.0, 0.0));

        assertEquals(1.5, step.mod(), 1E-6);
    }

    
}

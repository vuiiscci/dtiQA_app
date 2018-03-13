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
public class TestEulerFibreTracker extends TestFibreTracker {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestEulerFibreTracker(String name) {
	super(name);
	
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestEulerFibreTracker.class);
    }



 	
    protected FibreTracker getTracker(TractographyImage image) {
	return new EulerFibreTracker(new VectorLinearInterpolator(image), 0.25);
    }



    /**
     * Ensure we get the correct step
     *
     */
    public void testStep() {

        DT_TractographyImage cube = Images.getCube();

        EulerFibreTracker tracker = new EulerFibreTracker(new VectorLinearInterpolator(cube), 0.5);

        // in the middle
        Point3D point = Images.getPointAtVoxel(1.0, 0.5, 0.5);

        Vector3D step = tracker.getNextStep(point, new Vector3D(1.0, 0.0, 0.0));
        
        assertEquals(0.5, step.mod(), 1E-6);
        
        assertEquals(0.48296, step.x, 1E-4);
        assertEquals(0.12941, step.y, 1E-4);
        assertEquals(0.0, step.z, 1E-6);

    }


    
}

package tractography;

import junit.framework.*;
import junit.extensions.*;

import data.*;
import misc.DT;
import numerics.*;

import java.util.Random;


/**
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestNearestNeighbourInterpolator extends TestCase {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestNearestNeighbourInterpolator(String name) {
	super(name);
	
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestNearestNeighbourInterpolator.class);
    }


    public void testInterpolation() {
        
        DT_TractographyImage crossing = Images.getCrossing();
        
        NearestNeighbourInterpolator interpolator = new NearestNeighbourInterpolator(crossing);
        
        Vector3D[] pds = crossing.getPDs(1,1,0);
     
        for (int i = 0; i < 2; i++) { 

            Vector3D v = interpolator.getTrackingDirection(new Point3D(1.9 * Images.xVoxelDim, 1.8 * Images.yVoxelDim, 0.1 * Images.zVoxelDim), pds[i]);
            
            assertEquals(pds[i], v);

            v = interpolator.getTrackingDirection(new Point3D(1.9 * Images.xVoxelDim, 1.8 * Images.yVoxelDim, 0.1 * Images.zVoxelDim), i, true);

            assertEquals(pds[i], v);
        }


    }




}
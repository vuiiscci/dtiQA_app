package tractography;

import junit.framework.*;
import junit.extensions.*;

import data.*;
import misc.DT;
import numerics.*;


/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>FACT_FibreTracker</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>FACT_FibreTracker</code> with JUnit 3.8. 
 * The tests should work for deterministic fibre trackers. 
 * 
 * </dl>
 *
 * @version $Id: TestFACT_FibreTracker.java,v 1.5 2005/11/08 11:47:47 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestFACT_FibreTracker extends TestFibreTracker {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestFACT_FibreTracker(String name) {
	super(name);
	
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestFACT_FibreTracker.class);
    }



 	
    protected FibreTracker getTracker(TractographyImage image) {
	return new FACT_FibreTracker(image);
    }


    /**
     * Tracks along a straight line and checks the resulting tract
     *
     */
    public void testTractComposition() {

	DT_TractographyImage linear = Images.getLinear();

	FibreTracker tracker = getTracker(linear);
	
	TractCollection tc = tracker.trackFromSeed(new Point3D(1.5 * Images.xVoxelDim, 0.5 *  Images.yVoxelDim, 0.5 *  Images.zVoxelDim));
        
	assertEquals(1, tc.numberOfTracts());

	Tract t = tc.getTract(0);

	assertEquals(5, t.numberOfPoints());
        
        boolean[] hitVoxel = new boolean[5];

        for (int i = 0; i < 5; i++) {

            Point3D p = t.getPoint(i);

            int vox = (int)(p.x / linear.xVoxelDim());

            hitVoxel[vox] = true;
        }

        for (int i = 0; i < 5; i++) {
            assertTrue("Tract did not intersect all voxels", hitVoxel[i]);
        }

    }


    /**
     * Ensure we get the correct first / next step at various points.
     *
     */
    public void testStep() {

        DT_TractographyImage crossing = Images.getCrossing();

        FACT_FibreTracker tracker = new FACT_FibreTracker(crossing);

        Vector3D[] pds = crossing.getPDs(1,1,1);

        Vector3D step = tracker.getFirstStep(new Point3D(1.5 * Images.xVoxelDim, 1.2 * Images.yVoxelDim, 1.1 * Images.zVoxelDim), 0, true);

        assertEquals(1.0, pds[0].dot(step.normalized()), 1E-8);

        assertEquals(Images.xVoxelDim / 2.0, step.mod(), 0.05);
       
        step = tracker.getFirstStep(new Point3D(1.5 * Images.xVoxelDim, 1.2 * Images.yVoxelDim, 1.1 * Images.zVoxelDim), 0, false);

        assertEquals(-1.0, pds[0].dot(step.normalized()), 1E-8);

        step = tracker.getFirstStep(new Point3D(1.5 * Images.xVoxelDim, 1.2 * Images.yVoxelDim, 1.1 * Images.zVoxelDim), 1, true);

        assertEquals(1.0, pds[1].dot(step.normalized()), 1E-8);


        for (int i = 0; i < 2; i++) {

            step = tracker.getNextStep(new Point3D(1.5 * Images.xVoxelDim, 1.2 * Images.yVoxelDim, 1.1 * Images.zVoxelDim), new Vector3D(pds[i].x + 0.01, pds[i].y + 0.02, pds[i].z + 0.03).normalized());

            assertEquals(1.0, pds[i].dot(step.normalized()), 1E-8);

        }

        
    }


}

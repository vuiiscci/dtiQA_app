package tractography;

import junit.framework.*;
import junit.extensions.*;

import data.*;
import misc.DT;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>FibreTracker</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>FibreTracker</code> with JUnit 3.8. The tests should
 * work for deterministic fibre trackers. 
 * 
 * </dl>
 *
 * @version $Id: TestFibreTracker.java,v 1.5 2005/11/08 11:47:47 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public abstract class TestFibreTracker extends TestCase {
  

    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestFibreTracker(String name) {
	super(name);
	
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestFibreTracker.class);
    }



    public void testConnectionProbability() {
	
	DT_TractographyImage linear = Images.getLinear();

	// default anisotropy threshold is 0.0

	FibreTracker tracker = getTracker(linear);

  
	double[][][] prob = connectionProbability(tracker, new Point3D(2.5 * Images.xVoxelDim / 2.0, Images.yVoxelDim / 2.0, Images.zVoxelDim / 2.0), 2);

	assertEquals(1.0, prob[0][0][0], 0.00000001);
	assertEquals(1.0, prob[1][0][0], 0.00000001);
	assertEquals(1.0, prob[2][0][0], 0.00000001);
	assertEquals(1.0, prob[3][0][0], 0.00000001);
	assertEquals(1.0, prob[4][0][0], 0.00000001);


	// test snake
	DT_TractographyImage snake = Images.getSnake();

	tracker = getTracker(snake);

	prob = connectionProbability(tracker, new Point3D(1.5 * Images.xVoxelDim / 2.0, Images.yVoxelDim / 2.0, Images.zVoxelDim / 2.0), 2);

	assertEquals(1.0, prob[0][0][0], 0.00000001);
	assertEquals(1.0, prob[1][0][0], 0.00000001);
	assertEquals(1.0, prob[2][0][0], 0.00000001);
	assertEquals(0.0, prob[3][0][0], 0.00000001);
	assertEquals(0.0, prob[4][0][0], 0.00000001);
	assertEquals(0.0, prob[0][1][0], 0.00000001);
	assertEquals(0.0, prob[1][1][0], 0.00000001);
	assertEquals(1.0, prob[2][1][0], 0.00000001);
	assertEquals(1.0, prob[3][1][0], 0.00000001);
	assertEquals(1.0, prob[4][1][0], 0.00000001);

	// test crossing

	tracker = getTracker(Images.getCrossing());

	prob = connectionProbability(tracker, new Point3D(Images.xVoxelDim / 2.0, 1.5 * Images.yVoxelDim, Images.zVoxelDim / 2.0), 2);

	assertEquals(0.0, prob[0][0][0], 0.00000001);
	assertEquals(0.0, prob[1][0][0], 0.00000001);
	assertEquals(0.0, prob[2][0][0], 0.00000001);

	assertEquals(1.0, prob[0][1][0], 0.00000001);
	assertEquals(1.0, prob[1][1][0], 0.00000001);
	assertEquals(1.0, prob[2][1][0], 0.00000001);

	assertEquals(0.0, prob[0][2][0], 0.00000001);
	assertEquals(0.0, prob[1][2][0], 0.00000001);
	assertEquals(0.0, prob[2][2][0], 0.00000001);


	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {
		assertEquals(0.0, prob[i][j][1], 0.00000001);
	    }
	}

    }


    public void testIPThreshold() {
 
	DT_TractographyImage snake = Images.getSnake();
	
	FibreTracker tracker = getTracker(snake);

        tracker.setIP_Threshold(0.99);

        tracker.setCurveCheckInterval(1.0);

	double[][][] prob = connectionProbability(tracker, new Point3D(0.75 * Images.xVoxelDim, Images.yVoxelDim / 2.0, Images.zVoxelDim / 2.0), 2);

	assertEquals(1.0, prob[0][0][0], 0.00000001);
	assertEquals(1.0, prob[1][0][0], 0.00000001);
	assertEquals(0.0, prob[2][0][0], 0.00000001);
	assertEquals(0.0, prob[3][0][0], 0.00000001);
	assertEquals(0.0, prob[4][0][0], 0.00000001);

	assertEquals(0.0, prob[0][1][0], 0.00000001);
	assertEquals(0.0, prob[1][1][0], 0.00000001);
        assertEquals(0.0, prob[2][1][0], 0.00000001);

	assertEquals(0.0, prob[3][1][0], 0.00000001);
	assertEquals(0.0, prob[4][1][0], 0.00000001);


	tracker = getTracker(snake);

        tracker.setIP_Threshold(0.2);

        tracker.setCurveCheckInterval(1.0);

	prob = connectionProbability(tracker, new Point3D(1.5 * Images.xVoxelDim / 2.0, Images.yVoxelDim / 2.0, Images.zVoxelDim / 2.0), 2);

	assertEquals(1.0, prob[0][0][0], 0.00000001);
	assertEquals(1.0, prob[1][0][0], 0.00000001);
	assertEquals(1.0, prob[2][0][0], 0.00000001);
	assertEquals(0.0, prob[3][0][0], 0.00000001);
	assertEquals(0.0, prob[4][0][0], 0.00000001);
	
	assertEquals(0.0, prob[0][1][0], 0.00000001);
	assertEquals(0.0, prob[1][1][0], 0.00000001);
	assertEquals(1.0, prob[2][1][0], 0.00000001);
	assertEquals(1.0, prob[3][1][0], 0.00000001);
	assertEquals(1.0, prob[4][1][0], 0.00000001);
 

    }

    
    public void testBoundsChecking() {

	FibreTracker tracker = getTracker(Images.getLinear());
	
	// point is outside the image
	TractCollection tc = tracker.trackFromSeed(new Point3D(1.5 * Images.xVoxelDim, 
							       1.5 * Images.yVoxelDim, 0.5 * 
							       Images.zVoxelDim));

	assertEquals(1, tc.numberOfTracts());

	Tract t = tc.getTract(0);

	assertEquals(1, t.numberOfPoints());

    }



    public void testIsotropicMasking() {

	DT_TractographyImage linear = Images.getLinear();

        linear.computeIsotropicMask(0.2);

	FibreTracker tracker = getTracker(linear);

	double[][][] prob = connectionProbability(tracker, new Point3D(2.5 * Images.xVoxelDim / 2.0, Images.yVoxelDim / 2.0, Images.zVoxelDim / 2.0), 2);

	assertEquals(0.0, prob[0][0][0], 0.00000001);
	assertEquals(1.0, prob[1][0][0], 0.00000001);
	assertEquals(1.0, prob[2][0][0], 0.00000001);
	assertEquals(1.0, prob[3][0][0], 0.00000001);
	assertEquals(0.0, prob[4][0][0], 0.00000001);

    }


    protected abstract FibreTracker getTracker(TractographyImage image);



    protected double[][][] connectionProbability(FibreTracker tracker, Point3D seed, int iterations) {

	ConnectionProbabilityImage image = new ConnectionProbabilityImage(tracker.image.getDataDims(), tracker.image.getVoxelDims());
        
	for (int i = 0; i < iterations; i++) {
	    TractCollection t = tracker.trackFromSeed(seed);
	    image.processTracts(t);
	}
	
	return image.getConnectionProbabilities();

    }
  
    

}

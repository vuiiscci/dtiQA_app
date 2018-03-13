package misc;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;
import tools.*;

/**
 * Automated tests for DynamicScalarImage.
 *
 * @version $Id$
 * @author  Philip Cook
 * @see misc.DynamicScalarImage
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestDynamicScalarImage extends TestCase {

    /**
     * Test double array.
     */
    private double[] testData = {0.0, 1.01, 1.03, 1.02, 1.01, 1.05, 1.00, 1.09};

    /**
     * Test weights.
     */
    private double[] testWeights = {1.0, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0};


    
    public TestDynamicScalarImage(String name) {
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
	return new TestSuite(TestDynamicScalarImage.class);
    }


    public void testScalarStats() {
	
	int[] dataDims = new int[] {2, 2, 2};

	double[] voxelDims = new double[] {1.0, 2.0, 3.0};
	
	DynamicScalarImage image = new DynamicScalarImage(dataDims, voxelDims);


	for (int i = 0; i < testData.length; i++) {

	    image.addValue(1,1,1, testData[i], testWeights[i]);

	}

	double[][][] result = image.getVoxelStatistic("mean");
	
	assertEquals(ArrayOps.weightedMean(testData, testWeights), result[1][1][1], 1E-8);

	result = image.getVoxelStatistic("max");

	assertEquals(ArrayOps.max(testData), result[1][1][1], 1E-8);

	result = image.getVoxelStatistic("min");

	assertEquals(ArrayOps.min(testData), result[1][1][1], 1E-8);

	result = image.getVoxelStatistic("sum");

	assertEquals(ArrayOps.weightedMean(testData, testWeights) * ArrayOps.sum(testWeights), result[1][1][1], 1E-8);


    }


}

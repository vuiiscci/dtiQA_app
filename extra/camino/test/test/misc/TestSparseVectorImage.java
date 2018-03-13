package misc;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;
import tools.*;

/**
 * Automated tests for SparseVectorImage.
 *
 * @version $Id$
 * @author  Philip Cook
 * @see misc.SparseVectorImage
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestSparseVectorImage extends TestCase {

    /**
     * Test double array.
     */
    private double[] testData = {0.0, 1.01, 1.03, 1.02, 1.01, 1.05, 1.00, 1.09};

    
    public TestSparseVectorImage(String name) {
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
	return new TestSuite(TestSparseVectorImage.class);
    }


    public void testScalarStats() {
	
	int[] dataDims = new int[] {2, 2, 2};

	double[] voxelDims = new double[] {1.0, 2.0, 3.0};
	
	SparseVectorImage image = new SparseVectorImage(dataDims, voxelDims);

	for (int i = 0; i < testData.length; i++) {

	    image.addValue(1,1,1, testData[i]);

	}

	double[][][] result = image.getVoxelStatistic("mean");
	
	assertEquals(ArrayOps.mean(testData), result[1][1][1], 1E-8);

	result = image.getVoxelStatistic("max");

	assertEquals(ArrayOps.max(testData), result[1][1][1], 1E-8);

	result = image.getVoxelStatistic("min");

	assertEquals(ArrayOps.min(testData), result[1][1][1], 1E-8);

        result = image.getVoxelStatistic("sum");

	assertEquals(ArrayOps.sum(testData), result[1][1][1], 1E-8);

	result = image.getVoxelStatistic("median");

	assertEquals(ArrayOps.median(testData), result[1][1][1], 1E-8);

	double[][][] var = image.getVoxelStatistic("var");

	assertEquals(ArrayOps.var(testData, ArrayOps.mean(testData)), 
		     var[1][1][1], 1E-8);

	double[][][] std = image.getVoxelStatistic("std");

        assertEquals(Math.sqrt(var[1][1][1]), std[1][1][1], 1E-8);

    }



    public void testGrowth() {

	// add a large number of points to the image

	
	int[] dataDims = new int[] {2, 2, 2};

	double[] voxelDims = new double[] {1.0, 2.0, 3.0};
	
	SparseVectorImage image = new SparseVectorImage(dataDims, voxelDims);

	
	for (int r = 0; r < image.initialArrayLength; r++) {
	    
	    for (int i = 0; i < testData.length; i++) {
		
		image.addValue(1,1,1, testData[i]);
	    }
	}

	double[][][] result = image.getVoxelStatistic("mean");
	
	assertEquals(ArrayOps.mean(testData), result[1][1][1], 1E-8);

    }

}

package tractography;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;


/**
 * Automated tests for <code>TargetCP_Image</code>.
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 */
public class TestTargetCP_Image extends TestCase {

    // test image is 6x6 block, one slice

    private final int xDataDim = 6;
    private final int yDataDim = 6;
    private final int zDataDim = 1;

    private final double xVoxelDim = 1.0;
    private final double yVoxelDim = 1.0;
    private final double zVoxelDim = 1.0;

    // tracts
    private Tract horizontal = null;
    private Tract vertical = null;

    // seeded in voxel 2 3 0
    private Tract centreSeededHorizontal = null;

    // these should not be modified
    private final int[][][] horizontalEnds = new int[xDataDim][yDataDim][zDataDim];
    private final int[][][] verticalEnds = new int[xDataDim][yDataDim][zDataDim];
    private final int[][][] centre = new int[xDataDim][yDataDim][zDataDim];    

	

    public TestTargetCP_Image(String name) {
	super(name);

	horizontalEnds[0][2][0] = 1;
	horizontalEnds[1][2][0] = 1;
	horizontalEnds[0][3][0] = 1;
	horizontalEnds[1][3][0] = 1;

	horizontalEnds[4][2][0] = 2;
	horizontalEnds[5][2][0] = 2;
	horizontalEnds[4][3][0] = 2;
	horizontalEnds[5][3][0] = 2;


	verticalEnds[2][0][0] = 1;
	verticalEnds[2][1][0] = 1;
	verticalEnds[3][0][0] = 1;
	verticalEnds[3][1][0] = 1;

	verticalEnds[2][4][0] = 2;
	verticalEnds[2][5][0] = 2;
	verticalEnds[3][4][0] = 2;
	verticalEnds[3][5][0] = 2;

	centre[2][2][0] = 1;
	centre[2][3][0] = 1;

	centre[3][2][0] = 1;
	centre[3][3][0] = 1;

    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

	horizontal = new Tract();


	Tract t = new Tract();
	Tract t2 = new Tract();

	for (int i = 0; i < 3; i++) {
	    t2.addPoint(new Point3D(2.5 - i, 3.5, 0.0));
	    t.addPoint(new Point3D(2.5 + i, 3.5, 0.0));
	}

	t.joinTract(t2);

	centreSeededHorizontal = t.resample(0.1);

	//	System.err.println(t);
	//	System.err.println(t.toVoxelList(xVoxelDim, yVoxelDim, zVoxelDim));



	horizontal.addPoint(new Point3D(1.5, 3.5, 0.5));
	horizontal.addPoint(new Point3D(3.5, 3.5, 0.5));
	horizontal.addPoint(new Point3D(5.5, 3.5, 0.5));

	vertical = new Tract();

	vertical.addPoint(new Point3D(3.5, 1.5, 0.5));
	vertical.addPoint(new Point3D(3.5, 3.5, 0.5));
	vertical.addPoint(new Point3D(3.5, 5.5, 0.5));

	horizontal = horizontal.resample(0.1);
	vertical = vertical.resample(0.1);

    }

    
    protected void tearDown() {
	horizontal = null;
	vertical = null;
    }


    public static Test suite() {
	return new TestSuite(TestTargetCP_Image.class);
    }


    
    public void testTargetCP() {

	TargetCP_Image image = new TargetCP_Image(horizontalEnds, xVoxelDim, yVoxelDim, zVoxelDim);

	TractCollection tc = new TractCollection(3, 100.0);

	image.setCountFirstEntry(false);

	tc.addTract(horizontal);
	tc.addTract(vertical);
	
	image.processTracts(tc);
	image.processTracts(tc);

	double[][][] sc = image.getStreamlineCounts();

	for (int i = 0; i < xDataDim; i++) {
	    for (int j = 0; j < yDataDim; j++) {
		if (horizontalEnds[i][j][0] > 0) {
		    assertEquals(2.0, sc[i][j][0]);
		}
		else {
		    assertEquals(0.0, sc[i][j][0]);
		}
	    }
	}

	double[][][] cp = image.getConnectionProbabilities();

	for (int i = 0; i < xDataDim; i++) {
	    for (int j = 0; j < yDataDim; j++) {
		if (horizontalEnds[i][j][0] > 0) {
		    assertEquals(0.5, cp[i][j][0], 1E-6);
		}
		else {
		    assertEquals(0.0, cp[i][j][0], 1E-6);
		}
	    }
	}


	image.reset();

	image.setCountFirstEntry(true);


	tc = new TractCollection(3, 100.0);

	tc.addTract(centreSeededHorizontal);

	// this time count first entry only
	image.processTracts(tc);

	cp = image.getConnectionProbabilities();
	
	for (int i = 0; i < xDataDim; i++) {
	    for (int j = 0; j < yDataDim; j++) {
		if (horizontalEnds[i][j][0] == 1) {
		    // seed closer to target on left side
		    assertEquals(1.0, cp[i][j][0]);
		}
		else {
		    assertEquals(0.0, cp[i][j][0]);
		}
	    }
	}
    }


    public void testGetNonZeroTargetSC() {

	int[][][] testBlock = new int[10][2][1];

        testBlock[0][0][0] = 1;
        testBlock[0][1][0] = 2;

        testBlock[9][0][0] = 3;
        testBlock[9][1][0] = 4;

        TargetCP_Image image = new TargetCP_Image(testBlock, 1.0, 1.0, 1.0);
        
	image.setCountFirstEntry(false);

        TractCollection tc = new TractCollection();

        Tract t = new Tract();

        t.addPoint(new Point3D(0.1, 0.5, 0.5));
        t.addPoint(new Point3D(9.1, 1.5, 0.5));

        tc.addTract(t);

        t = new Tract();

        t.addPoint(new Point3D(0.1, 1.5, 0.5));
        t.addPoint(new Point3D(9.1, 1.5, 0.5));
        
        tc.addTract(t);

        image.processTracts(tc);

        int[][] sc = image.getNonZeroStreamlineCounts();

	assertEquals(3, sc.length);

	assertEquals(1, sc[0][0]);
	assertEquals(1, sc[0][1]);
	
	assertEquals(2, sc[1][0]);
	assertEquals(1, sc[1][1]);
	
	assertEquals(4, sc[2][0]);
	assertEquals(2, sc[2][1]);
	

    }


}




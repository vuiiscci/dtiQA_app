package tractography;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;


/**
 * Automated tests for <code>StreamlineROI_Filter</code>.
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 */
public class TestStreamlineROI_Filter extends TestCase {

    // test image is 3x3 block, one slice

    private final int xDataDim = 6;
    private final int yDataDim = 6;
    private final int zDataDim = 1;

    private final double xVoxelDim = 1.0;
    private final double yVoxelDim = 1.0;
    private final double zVoxelDim = 1.0;

    private Tract horizontal = null;
    private Tract vertical = null;

    // these should not be modified
    private final int[][][] horizontalEnds = new int[xDataDim][yDataDim][zDataDim];
    private final int[][][] verticalEnds = new int[xDataDim][yDataDim][zDataDim];
    private final int[][][] centre = new int[xDataDim][yDataDim][zDataDim];    


    StreamlineROI_Filter filter = null;
	

    public TestStreamlineROI_Filter(String name) {
	super(name);

	horizontalEnds[0][2][0] = 3;
	horizontalEnds[1][2][0] = 3;
	horizontalEnds[0][3][0] = 3;
	horizontalEnds[1][3][0] = 3;

	horizontalEnds[4][2][0] = 5;
	horizontalEnds[5][2][0] = 5;
	horizontalEnds[4][3][0] = 5;
	horizontalEnds[5][3][0] = 5;


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

	horizontal.addPoint(new Point3D(1.5, 3.5, 0.5));
	horizontal.addPoint(new Point3D(2.5, 3.5, 0.5));
	horizontal.addPoint(new Point3D(3.5, 3.5, 0.5));
	horizontal.addPoint(new Point3D(4.5, 3.5, 0.5));
	horizontal.addPoint(new Point3D(5.5, 3.5, 0.5));

	vertical = new Tract();

	vertical.addPoint(new Point3D(3.5, 1.5, 0.5));
	vertical.addPoint(new Point3D(3.5, 2.5, 0.5));
	vertical.addPoint(new Point3D(3.5, 3.5, 0.5));
	vertical.addPoint(new Point3D(3.5, 4.5, 0.5));
	vertical.addPoint(new Point3D(3.5, 5.5, 0.5));

	filter = new StreamlineROI_Filter(new int[] {xDataDim, yDataDim, zDataDim},
                                          new double[] {xVoxelDim, yVoxelDim, zVoxelDim});
    }

    
    protected void tearDown() {
	horizontal = null;
	vertical = null;
        filter = null;
    }


    public static Test suite() {
	return new TestSuite(TestStreamlineROI_Filter.class);
    }



    public void testWaypoints() {
	
	filter.setResampleTracts(false);

	TractCollection tc = new TractCollection(3, 100.0);

	tc.addTract(vertical);
	tc.addTract(horizontal);

	filter.setWaypoints(horizontalEnds);

	tc = filter.processTracts(tc);

	assertEquals(1, tc.numberOfTracts());

	assertTrue(tc.getTract(0).equals(horizontal));

	filter.setWaypoints(verticalEnds);

	tc = new TractCollection(3, 100.0);

	tc.addTract(vertical);
	tc.addTract(horizontal);

	tc = filter.processTracts(tc);
	
	assertEquals(1, tc.numberOfTracts());

	assertTrue(tc.getTract(0).equals(vertical));

    }


    public void testDiscardOnExclusion() {

	filter.setResampleTracts(false);
	filter.setDiscardOnExclusionEntry(true);


	TractCollection tc = new TractCollection(3, 100.0);

	tc.addTract(vertical);
	tc.addTract(horizontal);

	filter.setExclusionROIs(horizontalEnds);

	tc = filter.processTracts(tc);

	assertEquals(1, tc.numberOfTracts());

	assertTrue(tc.getTract(0).equals(vertical));
	
    }


    public void testTerminateOnExclusion() {

	filter.setResampleTracts(false);
	filter.setDiscardOnExclusionEntry(false);

	TractCollection tc = new TractCollection(3, 100.0);

	tc.addTract(vertical);
	
	Tract left = new Tract(10, 100.0);
	Tract right = new Tract(10, 100.0);
	
	left.addPoint(new Point3D(3.5, 3.5, 0.5));
	left.addPoint(new Point3D(2.5, 3.5, 0.5));
	left.addPoint(new Point3D(1.5, 3.5, 0.5));
	left.addPoint(new Point3D(1.25, 3.5, 0.5));

	right.addPoint(new Point3D(3.5, 3.5, 0.5));
	right.addPoint(new Point3D(4.5, 3.5, 0.5));
	right.addPoint(new Point3D(5.5, 3.5, 0.5));
	right.addPoint(new Point3D(6.0, 3.5, 0.5));
	

	// should be horizontal with two extra points, one on each end
	right.joinTract(left);

	tc.addTract(right);
	
	filter.setExclusionROIs(horizontalEnds);

	filter.setDiscardOnExclusionEntry(false);

	tc = filter.processTracts(tc);

	assertEquals(2, tc.numberOfTracts());

	assertTrue(tc.getTract(0).equals(vertical));

	Tract chopped = tc.getTract(1);


	// should just leave horizontal
	assertEquals(4, chopped.numberOfPoints());
	assertEquals(2, chopped.seedPointIndex());

	for (int i = 0; i < 4; i++) {
	    assertTrue(chopped.getPoint(i).equals(horizontal.getPoint(i)));
	}



    }


    public void testLoopCheck() {

	filter.setResampleTracts(false);

	TractCollection tc = new TractCollection(3, 100.0);

	// back in the waypoint
	Tract verticalPlus = new Tract(vertical);

	verticalPlus.addPoint(new Point3D(3.5, 3.5, 0.5));
	verticalPlus.addPoint(new Point3D(3.5, 1.5, 0.5));
	verticalPlus.addPoint(new Point3D(3.5, 0.5, 0.5));

	tc.addTract(verticalPlus);
	
	filter.setWaypoints(centre);

	filter.setTruncateLoops(true);

	tc = filter.processTracts(tc);

	// should take out the last two points
	assertEquals(1, tc.numberOfTracts());

	Tract chopped = tc.getTract(0); //  == verticalPlus

	assertEquals(vertical.numberOfPoints() + 1, chopped.numberOfPoints());

	for (int i = 0; i < vertical.numberOfPoints(); i++) {
	    assertTrue(chopped.getPoint(i).equals(vertical.getPoint(i)));
	}

	assertTrue(chopped.getPoint(vertical.numberOfPoints()).equals(new Point3D(3.5, 3.5, 0.5)));

   }


    public void testLoopDiscard() {

	filter.setResampleTracts(false);

	TractCollection tc = new TractCollection(3, 100.0);

	// back in the waypoint
	vertical.addPoint(new Point3D(3.5, 3.5, 0.5));
	vertical.addPoint(new Point3D(3.5, 1.5, 0.5));

	tc.addTract(vertical);
	tc.addTract(horizontal);
	
	filter.setWaypoints(centre);

	filter.setDiscardLoops(true);

	tc = filter.processTracts(tc);

	// should take out the last point
	assertEquals(1, tc.numberOfTracts());

	Tract h = tc.getTract(0);

	assertEquals(h, horizontal);
   }


    public void testMinTractLength() {
	
	TractCollection tc = new TractCollection(3, 100.0);
	
	filter.setMinTractLength(3.0);
	
	tc.addTract(vertical);

	tc = filter.processTracts(tc);

	assertEquals(1, tc.numberOfTracts());

	filter.setMinTractLength(10.0);
	
	tc = filter.processTracts(tc);

	assertEquals(0, tc.numberOfTracts());
	
    }




    public void testMinTractPoints() {
	

	TractCollection tc = new TractCollection(3, 100.0);
	
	filter.setMinTractPoints(1);
	
	tc.addTract(vertical);

	tc = filter.processTracts(tc);

	assertEquals(1, tc.numberOfTracts());

	filter.setMinTractLength(100);
	
	tc = filter.processTracts(tc);

	assertEquals(0, tc.numberOfTracts());
	
    }


    public void testMaxTractPoints() {
	

	for (int i = 300; i < 302; i++) {

	    filter.setMaxTractPoints(i);

	    TractCollection tc = new TractCollection(3, 100.0);
	    
	    tc.addTract(vertical.resample(1E-2));
	    
	    tc = filter.processTracts(tc);
	    
	    assertEquals(1, tc.numberOfTracts());

	    assertEquals(i, tc.getTract(0).numberOfPoints());
	    
	}


	// test what happens if the seed point is near the end of a tract
	// instead of in the middle
	
	Tract t1 = new Tract(100, 100.0);
	Tract t2 = new Tract(100, 100.0);

	for (int i = 0; i < 10; i++) {
	    t1.addPoint(new Point3D(5.0 - i / 2.0, 0.0, 0.0));
	}

	t2.addPoint(new Point3D(5.0, 0.0, 0.0));
	t2.addPoint(new Point3D(5.1, 0.0, 0.0));
	t2.addPoint(new Point3D(5.2, 0.0, 0.0));

	t2.joinTract(t1);
	

	// t2 now contains 12 points

	filter.setMaxTractPoints(7);
	filter.setResampleTracts(false);

	TractCollection tc = new TractCollection(3, 100.0);
	
	tc.addTract(t2);
	
	tc = filter.processTracts(tc);
	
	assertEquals(1, tc.numberOfTracts());
	
	assertEquals(7, tc.getTract(0).numberOfPoints());
	
	
    }


    public void testMaxTractLength() {
	
	filter.setMaxTractLength(3.0);

	TractCollection tc = new TractCollection(3, 100.0);
	
	tc.addTract(vertical.resample(1E-2));
	
	tc = filter.processTracts(tc);
	
	assertEquals(1, tc.numberOfTracts());
	
	assertEquals(3.0, tc.getTract(0).length(), 1E-2);
	
	


	// test what happens if the seed point is near the end of a tract
	// instead of in the middle
	
	Tract t1 = new Tract(100, 100.0);
	Tract t2 = new Tract(100, 100.0);

	for (int i = 0; i < 10; i++) {
	    t1.addPoint(new Point3D(5.0 - i / 2.0, 0.0, 0.0));
	}

	t2.addPoint(new Point3D(5.0, 0.0, 0.0));
	t2.addPoint(new Point3D(5.1, 0.0, 0.0));
	t2.addPoint(new Point3D(5.2, 0.0, 0.0));

	t2.joinTract(t1);

	// t2 now contains 12 points

	filter.setMaxTractLength(4.0);

	tc = new TractCollection(3, 100.0);
	
	tc.addTract(t2.resample(0.01));
	
	tc = filter.processTracts(tc);
	
	assertEquals(1, tc.numberOfTracts());
	
	assertEquals(4.0, tc.getTract(0).length(), 1E-2);	
    }


    public void testEndPoints() {
	
	// 10 1 1 line, end zones at ends + 1 in the middle
	int[][][] endZones = new int[10][1][1];
	
	Point3D[] points = new Point3D[10];

	for (int i = 0; i < 10; i++) {
	    points[i] = new Point3D(i + 0.5, 0.5, 0.5);
	}

	endZones[0][0][0] = 1;
	endZones[1][0][0] = 2;
	
	endZones[5][0][0] = 3;
	endZones[6][0][0] = 3;

	endZones[8][0][0] = 4;
	endZones[9][0][0] = 5;

	

	filter = new StreamlineROI_Filter(new int[] {10, 1, 1}, new double[] {1.0, 1.0, 1.0});
	filter.setResampleTracts(false);
	filter.setEndZones(endZones);

	// case 1: seed between two end zones

	TractCollection tc = new TractCollection(2, 100.0);

	Tract t1 = new Tract();
	
	t1.addPoint(points[4]);
	t1.addPoint(points[3]);
	t1.addPoint(points[2]);
	t1.addPoint(points[1]);
	t1.addPoint(points[0]);

	Tract t2 = new Tract();

	t2.addPoint(points[4]);
	t2.addPoint(points[5]);
	t2.addPoint(points[6]);
	t2.addPoint(points[7]);
	t2.addPoint(points[8]);


	t2.joinTract(t1);

	tc.addTract(t2);

	TractCollection output = filter.processTracts(tc);

	assertEquals(1, output.numberOfTracts());

	Tract endZoneTrunc = output.getTract(0);

	// should truncate between voxels 1 and 5
	assertEquals(5, endZoneTrunc.numberOfPoints());

	// case 2: seed within one end zone, connects to two others

	tc = new TractCollection(2, 100.0);
	t1 = new Tract();
	t2 = new Tract();
	
	t1.addPoint(points[6]);
	t1.addPoint(points[5]);
	t1.addPoint(points[4]);
	t1.addPoint(points[3]);
	t1.addPoint(points[2]);
	t1.addPoint(points[1]);

	t2.addPoint(points[6]);
	t2.addPoint(points[7]);
	t2.addPoint(points[8]);
	t2.addPoint(points[9]);

	t2.joinTract(t1);

	tc.addTract(t2);

	output = filter.processTracts(tc);

	assertEquals(1, output.numberOfTracts());
	
	endZoneTrunc = output.getTract(0);

	assertEquals(points[6], endZoneTrunc.getPoint(0));
	assertEquals(points[8], endZoneTrunc.getPoint(2));


	// case 3: both ends connect to same end zone, seed in zone

	tc = new TractCollection(2, 100.0);
	t1 = new Tract();
	t2 = new Tract();
	
	t1.addPoint(points[6]);
	t1.addPoint(points[5]);
	t1.addPoint(points[4]);
	t1.addPoint(points[3]);
	t1.addPoint(points[2]);
	t1.addPoint(points[3]);
	t1.addPoint(points[4]);
	t1.addPoint(points[5]);

	tc.addTract(t1);

	output = filter.processTracts(tc);

	assertEquals(0, output.numberOfTracts());
	
	
	// case 4: one end doesn't connect

	tc = new TractCollection(2, 100.0);
	t1 = new Tract();
	t2 = new Tract();
	
	t1.addPoint(points[4]);
	t1.addPoint(points[3]);
	t1.addPoint(points[2]);
	t1.addPoint(points[1]);
	t1.addPoint(points[0]);
	
	tc.addTract(t1);

	output = filter.processTracts(tc);

	assertEquals(0, output.numberOfTracts());
	
	
	// case 5: both ends don't connect, seed in one zone

	tc = new TractCollection(2, 100.0);
	t1 = new Tract();
	t2 = new Tract();
	
	t1.addPoint(points[1]);
	t1.addPoint(points[2]);
	t1.addPoint(points[3]);
	t1.addPoint(points[4]);
	
	tc.addTract(t1);

	output = filter.processTracts(tc);

	assertEquals(0, output.numberOfTracts());


	// case 6: both ends connect to same end zone

	tc = new TractCollection(2, 100.0);
	t1 = new Tract();
	t2 = new Tract();
	

	t1.addPoint(points[4]);
	t1.addPoint(points[3]);
	t1.addPoint(points[2]);
	t1.addPoint(points[1]);


	t2.addPoint(points[4]);
	t2.addPoint(points[3]);
	t2.addPoint(points[2]);
	t2.addPoint(points[1]);

	t2.joinTract(t1);

	tc.addTract(t2);

	output = filter.processTracts(tc);

	assertEquals(0, output.numberOfTracts());

	    
    }

}

package tractography;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;


/**
 * Automated tests for <code>TestStreamlineConnectivityGraph</code>.
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 */
public class TestStreamlineConnectivityGraph extends TestCase {

    // test image is 5x5 block, one slice

    private final int xDataDim = 5;
    private final int yDataDim = 5;
    private final int zDataDim = 1;

    private final double xVoxelDim = 1.0;
    private final double yVoxelDim = 1.0;
    private final double zVoxelDim = 1.0;

    private final Tract t1;
    private final Tract t2;
    private final Tract t3;
    private final Tract t4;
    private final Tract t5;
    private final Tract t6;
    private final Tract t7;

    private int[][][] targets = null;

    StreamlineROI_Filter filter = null;
	

    public TestStreamlineConnectivityGraph(String name) {
	super(name);


        targets = new int[xDataDim][yDataDim][zDataDim];


        targets[0][0][0] = 2;
        targets[0][1][0] = 3;
        targets[0][3][0] = 5;
        targets[0][4][0] = 7;

        targets[4][0][0] = 18;
        targets[4][1][0] = 16;
        targets[4][3][0] = 14;
        targets[4][4][0] = 11;

	t1 = new Tract();

        t1.addPoint(new Point3D(2.5, 2.5, 0.5));
        t1.addPoint(new Point3D(0.5, 0.5, 0.5)); // label 2, node 1
        
        
	t2 = new Tract();
        
        t2.addPoint(new Point3D(2.5, 2.5, 0.5));
        t2.addPoint(new Point3D(4.5, 0.5, 0.5)); // label 18, node 8
                    

        t1.joinTract(t2);

	t3 = new Tract();

        t3.addPoint(new Point3D(2.5, 2.5, 0.5));
        t3.addPoint(new Point3D(0.5, 1.5, 0.5)); // label 3, node 2


	t4 = new Tract();

        t4.addPoint(new Point3D(2.5, 2.5, 0.5));
        t4.addPoint(new Point3D(4.5, 4.5, 0.5)); // label 11, node 5


        t3.joinTract(t4);
        

	t5 = new Tract();
        
        t5.addPoint(new Point3D(0.5, 0.5, 0.5)); // label 2, node 1
        t5.addPoint(new Point3D(2.5, 2.5, 0.5)); 
        t5.addPoint(new Point3D(4.5, 4.5, 0.5)); // label 11, node 5

        t6 = new Tract();

        t6.addPoint(new Point3D(0.5, 0.5, 0.5)); // label 2, node 1
        t6.addPoint(new Point3D(2.5, 2.5, 0.5)); 
        t6.addPoint(new Point3D(4.5, 4.5, 0.5)); // label 11, node 5


        t7 = new Tract(); // self connectivity in a seed region, should get ignored

        t7.addPoint(new Point3D(0.5, 0.5, 0.5)); // label 2, node 1
        t7.addPoint(new Point3D(2.5, 2.5, 0.5)); 
        t7.addPoint(new Point3D(0.5, 0.5, 0.5)); // label 2, node 1

        

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
	return new TestSuite(TestStreamlineConnectivityGraph.class);
    }



    public void testStreamlineCounts() {
        
        StreamlineConnectivityGraph graph = new StreamlineConnectivityGraph(targets, new double[] {xVoxelDim, yVoxelDim, zVoxelDim});

    

        graph.processTract(t1);


        graph.processTract(t1);
        graph.processTract(t3);

        graph.processTract(t2); // only connects one target, should be ignored

        graph.processTract(t5); 
        graph.processTract(t6); // should both update same entry

        graph.processTract(t7); // should be ignored

        RealMatrix matrix = graph.getStreamlineCountMatrix();

        assertEquals(8, matrix.rows());
        assertEquals(8, matrix.columns());

	for (int i = 0; i < 8; i++) {
	    assertEquals(0.0, matrix.entries[i][i], 1E-8);
	}


	for (int i = 1; i < 8; i++) {
	    assertEquals(0.0, matrix.entries[i][7], 1E-8);
	}

	assertEquals(2.0, matrix.entries[0][7]);
	assertEquals(2.0, matrix.entries[7][0]);

	assertEquals(1.0, matrix.entries[1][4]);

	assertEquals(2.0, matrix.entries[0][4]);


    }


    public void testMatrixString() {
 
	StreamlineConnectivityGraph graph = new StreamlineConnectivityGraph(targets, new double[] {xVoxelDim, yVoxelDim, zVoxelDim});

    

        graph.processTract(t1);
        graph.processTract(t1);

        graph.processTract(t3);


	String csv = graph.getStreamlineCountMatrixCSV();
	

	String[] tokens = csv.split("[,\\n]");

	// 8 rows plus header
	assertEquals(72, tokens.length);

	// Labels should be ordered by intensity
	assertEquals("Label_0016", tokens[6]);

	assertEquals(2.0, Double.parseDouble(tokens[15]), 1E-8);

	

    }


    public void testTractStats() {


        double[][][] scalars = new double[xDataDim][yDataDim][zDataDim];

	scalars[4][0][0] = 6.0;
	scalars[0][0][0] = 6.0;

	scalars[1][2][0] = 12.0;


	StreamlineConnectivityGraph graph = new StreamlineConnectivityGraph(targets, new double[] {xVoxelDim, yVoxelDim, zVoxelDim});


	TractStatisticFilter filter = new TractStatisticFilter(scalars, new double[] {xVoxelDim, yVoxelDim, zVoxelDim});

	graph.setTractStatFilter(filter);
	
        graph.processTract(t1);


	// same connectivity as t1 but different scalars
	Tract t = new Tract();

	t.addPoint(new Point3D(0.5, 0.5, 0.5));
	t.addPoint(new Point3D(1.5, 2.5, 0.5));
	t.addPoint(new Point3D(4.5, 0.5, 0.5));

	graph.processTract(t);

	RealMatrix matrix = graph.getTractStatisticMatrix();

	for (int i = 1; i < 7; i++) {
	    for (int j = 1; j < 7; j++) {
		assertEquals(0.0, matrix.entries[0][0], 1E-8);
	    }
	}

	assertEquals(6.0, matrix.entries[0][7], 1E-8);
    }

}

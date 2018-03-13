package tractography;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;

/**
 *  Automated tests for <code>ConnectionProbabilityImage</code>.
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestConnectionProbabilityImage extends TestCase {

    // test image is 3x3 block, one slice

    private final int xDataDim = 3;
    private final int yDataDim = 3;
    private final int zDataDim = 1;

    private final double xVoxelDim = 2.0;
    private final double yVoxelDim = 2.0;
    private final double zVoxelDim = 2.0;

    // tracts sampled every 2 mm (diffusion space)
    private Tract horizontal = null;
    private Tract vertical = null;


    public TestConnectionProbabilityImage(String name) {
	super(name);
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
	horizontal.addPoint(new Point3D(3.5, 3.5, 0.5));
	horizontal.addPoint(new Point3D(5.5, 3.5, 0.5));

	vertical = new Tract();

	vertical.addPoint(new Point3D(3.5, 1.5, 0.5));
	vertical.addPoint(new Point3D(3.5, 3.5, 0.5));
	vertical.addPoint(new Point3D(3.5, 5.5, 0.5));


    }

    
    protected void tearDown() {

	horizontal = null;
	vertical = null;

    }


    public static Test suite() {
	return new TestSuite(TestConnectionProbabilityImage.class);
    }


    public void testCP() {


	ConnectionProbabilityImage image = 
	    new ConnectionProbabilityImage(xDataDim, yDataDim, zDataDim, xVoxelDim, yVoxelDim,
					   zVoxelDim);

	TractCollection tc = new TractCollection(3, 100.0);

	tc.addTract(horizontal);
	tc.addTract(vertical);
	
	image.processTracts(tc);
	image.processTracts(tc);

	double[][][] sc = image.getStreamlineCounts();

	assertEquals(0.0, sc[0][0][0]);
	assertEquals(2.0, sc[0][1][0]);
	assertEquals(0.0, sc[0][2][0]);
	assertEquals(2.0, sc[1][0][0]);
	assertEquals(4.0, sc[1][1][0]);
	assertEquals(2.0, sc[1][2][0]);
	assertEquals(0.0, sc[2][0][0]);
	assertEquals(2.0, sc[2][1][0]);
	assertEquals(0.0, sc[2][2][0]);

	double[][][] cp = image.getConnectionProbabilities();

	assertEquals(0.0, cp[0][0][0], 1E-6);
	assertEquals(0.5, cp[0][1][0], 1E-6);
	assertEquals(0.0, cp[0][2][0], 1E-6);
	assertEquals(0.5, cp[1][0][0], 1E-6);
	assertEquals(1.0, cp[1][1][0], 1E-6);
	assertEquals(0.5, cp[1][2][0], 1E-6);
	assertEquals(0.0, cp[2][0][0], 1E-6);
	assertEquals(0.5, cp[2][1][0], 1E-6);
	assertEquals(0.0, cp[2][2][0], 1E-6);
	

    }


}

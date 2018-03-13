package tractography;

import junit.framework.*;
import junit.extensions.*;

import misc.*;
import numerics.*;
import tools.*;


/**
 * 
 * Automated tests for <code>TractStatisticImage</code>.
 * 
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTractStatisticImage extends TestCase {

    
    private ScalarImage scalars = null;

    private Tract tInterp = null;

    private Tract tFACT = null;

    private double[][][] data = null;

    private double[] voxelDims = null;



    public TestTractStatisticImage(String name) {
	super(name);

	tInterp = new Tract();

	Tract t1 = new Tract();
	Tract t2 = new Tract();

	Point3D point = new Point3D(5.5, 0.5, 0.5);

	for (int i = 0; i < 51; i++) {
	    t1.addPoint(point);
	    point = point.displace(new Vector3D(0.1, 0.0, 0.0));
	}


	point = new Point3D(5.5, 0.5, 0.5);
	
	for (int i = 0; i < 51; i++) {
	    t2.addPoint(point);
	    point = point.displace(new Vector3D(-0.1, 0.0, 0.0));
	}

	t2.joinTract(t1);
	
	tInterp = t2;



	t1 = new Tract();
	t2 = new Tract();
	
	t1.addPoint(new Point3D(5.5, 0.5, 0.5));
	t1.addPoint(new Point3D(6.001, 0.5, 0.5));
	t1.addPoint(new Point3D(7.001, 0.5, 0.5));
	t1.addPoint(new Point3D(8.001, 0.5, 0.5));
	t1.addPoint(new Point3D(9.001, 0.5, 0.5));
	t1.addPoint(new Point3D(10.001, 0.5, 0.5));

	t2.addPoint(new Point3D(5.5, 0.5, 0.5));
	t2.addPoint(new Point3D(4.999, 0.5, 0.5));
	t2.addPoint(new Point3D(3.999, 0.5, 0.5));
	t2.addPoint(new Point3D(2.999, 0.5, 0.5));
	t2.addPoint(new Point3D(1.999, 0.5, 0.5));
	t2.addPoint(new Point3D(0.999, 0.5, 0.5));

	t2.joinTract(t1);

	tFACT = t2;

	
	data = new double[11][2][1];

	

	for (int i = 0; i < 11; i++) {
	    
	    data[i][0][0] = i + 1.0;
	}

	voxelDims = new double[] {1.0, 1.0, 1.0};
	
	scalars = new ScalarImage(data, voxelDims);
	
	
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }



    protected void setUp() {
	
    }

    
    protected void tearDown() {

    }


    public static Test suite() {
	return new TestSuite(TestTractStatisticImage.class);
    }

    
    
    public void testSVO_Stats() {

	// assuming that this works
	TractStatisticFilter filter = new TractStatisticFilter(scalars);
	filter.setInterpolate(false);
	filter.setTractStatistic("mean");
	
	
	TractStatisticImage image = new TractStatisticImage(scalars);
	
	image.setTractStatistic("mean");
	image.setImageStatistic("mean");
	
	image.setInterpolate(false);
	
	image.processTract(tFACT);
	image.processTract(tFACT);
	
	VoxelList voxels = tFACT.toVoxelList(voxelDims);
	
	Voxel seedPointVox = voxels.getVoxel(voxels.seedPointIndex());
	
	double[][][] result = image.getImageStatistic();
	
	double factStat = filter.processTract(tFACT)[0];

	for (int i = 0; i < data.length; i++) {
	    if (i == seedPointVox.x) {
		assertEquals(factStat, result[i][0][0], 1E-8);
	    }
	    else {
		assertEquals(0.0, result[i][0][0], 1E-8);
	    }
	}

	



    }

    

    public void testCountIntersect() {
	
	TractStatisticFilter filter = new TractStatisticFilter(scalars);
	filter.setInterpolate(false);
	filter.setTractStatistic("max");
	
	double factStat = filter.processTract(tFACT)[0];

	TractStatisticImage image = new TractStatisticImage(scalars);
	
	image.setTractStatistic("max");
	image.setImageStatistic("mean");
	image.setCountIntersect(true);
	image.setInterpolate(false);


	Tract t3 = new Tract();

	t3.addPoint(new Point3D(1.5, 1.5, 0.5)); // 1 1 0
	t3.addPoint(new Point3D(1.5, 0.5, 0.5)); // 1 0 0

	image.processTract(tFACT);
	image.processTract(t3);
	
	double[][][] result = image.getImageStatistic();

	for (int i = 0; i < data.length; i++) {
	    if (i == 1) {
		assertEquals((factStat + 2.0 ) / 2.0 , result[i][0][0], 1E-8);
	    }
	    else {
		assertEquals(factStat, result[i][0][0], 1E-8);
	    }
	}
	
	assertEquals(2.0, result[1][1][0], 1E-8);

    }


    public void testNonImageStats() {

	// assuming that this works
	TractStatisticFilter filter = new TractStatisticFilter(scalars.getDataDims(), scalars.getVoxelDims());
	filter.setInterpolate(false);
	filter.setTractStatistic("length");
	
	TractStatisticImage image = new TractStatisticImage(scalars.getDataDims(), scalars.getVoxelDims());
	
	image.setTractStatistic("length");
	image.setImageStatistic("mean");
	
	image.setInterpolate(false);
	
	image.processTract(tFACT);

	double factValue = filter.processTract(tFACT)[0];

	double[][][] result = image.getImageStatistic();
	
	assertEquals(factValue, result[5][0][0], 1E-8);


    }

}

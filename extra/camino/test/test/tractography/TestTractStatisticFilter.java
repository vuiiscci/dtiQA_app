package tractography;

import junit.framework.*;
import junit.extensions.*;

import misc.*;
import numerics.*;
import tools.*;


/**
 * 
 * Automated tests for <code>TractStatisticFilter</code>.
 * 
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTractStatisticFilter extends TestCase {

    
    private ScalarImage scalars = null;

    private Tract tInterp = null;

    private Tract tFACT = null;

    private double[][][] data = null;

    private double[] voxelDims = null;



    public TestTractStatisticFilter(String name) {
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
	return new TestSuite(TestTractStatisticFilter.class);
    }

    
    
    public void testScalarStats() {

	TractStatisticFilter filter = new TractStatisticFilter(scalars);
	
	filter.setInterpolate(false);

	VoxelList voxels = tFACT.toVoxelList(voxelDims);

	double[] values = scalars.valuesAt(voxels.getVoxels());
	
		
	// mean
	filter.setTractStatistic("mean");
	double factStat = filter.processTract(tFACT)[0];

	assertEquals(ArrayOps.mean(values), factStat, 1E-8);
	
	// min
	filter.setTractStatistic("min");
	factStat = filter.processTract(tFACT)[0];

	assertEquals(ArrayOps.min(values), factStat, 1E-8);


	// max
	filter.setTractStatistic("max");

	factStat = filter.processTract(tFACT)[0];
	assertEquals(ArrayOps.max(values), factStat, 1E-8);


	// max
	filter.setTractStatistic("sum");

	factStat = filter.processTract(tFACT)[0];
	assertEquals(ArrayOps.sum(values), factStat, 1E-8);



	// median

	filter.setTractStatistic("median");

	factStat = filter.processTract(tFACT)[0];
	assertEquals(ArrayOps.median(values), factStat, 1E-8);


	// var
	filter.setTractStatistic("var");
	
	factStat = filter.processTract(tFACT)[0];
	assertEquals(ArrayOps.var(values, ArrayOps.mean(values)), factStat, 1E-8);
	

	// meanvar
	filter.setTractStatistic("meanvar");
	
	double[] meanVar = filter.processTract(tFACT);
	
	assertEquals(ArrayOps.mean(values), meanVar[0], 1E-8);
	assertEquals(ArrayOps.var(values, meanVar[0]), meanVar[1], 1E-8);


	filter.setInterpolate(true);

	values = scalars.valuesAt(tInterp.getPoints());


	// mean
	filter.setTractStatistic("mean");

	double intStat = filter.processTract(tInterp)[0];

	assertEquals(ArrayOps.mean(values), intStat, 1E-8);
	
	// min
	filter.setTractStatistic("min");
	intStat = filter.processTract(tInterp)[0];

	assertEquals(ArrayOps.min(values), intStat, 1E-8);


	// max
	filter.setTractStatistic("max");

	intStat = filter.processTract(tInterp)[0];
	assertEquals(ArrayOps.max(values), intStat, 1E-8);



	// median

	filter.setTractStatistic("median");

	intStat = filter.processTract(tInterp)[0];
	assertEquals(ArrayOps.median(values), intStat, 1E-8);


	// var
	filter.setTractStatistic("var");
	
	intStat = filter.processTract(tInterp)[0];
	assertEquals(ArrayOps.var(values, ArrayOps.mean(values)), intStat, 1E-8);
	

	// meanvar
	filter.setTractStatistic("meanvar");
	
	meanVar = filter.processTract(tInterp);
	
	assertEquals(ArrayOps.mean(values), meanVar[0], 1E-8);
	assertEquals(ArrayOps.var(values, meanVar[0]), meanVar[1], 1E-8);

    }


    public void testTractLength() {
	
	// length
	TractStatisticFilter filter = new TractStatisticFilter(scalars.getDataDims(), scalars.getVoxelDims());

	filter.setTractStatistic("length");

	double length = filter.processTract(tFACT)[0];

	assertEquals(9.002, length, 1E-8);

	length = filter.processTract(tInterp)[0];

	assertEquals(10.0, length, 1E-8);
	
    }


    public void testEndPointSep() {
	

	Tract tCurve = new Tract();

	Tract t1 = new Tract();
	Tract t2 = new Tract();

	Point3D point = new Point3D(5.5, 0.5, 0.5);

	for (int i = 0; i < 5; i++) {
	    t1.addPoint(point);
	    point = point.displace(new Vector3D(0.1, 0.2 * i, 0.0));
	}


	point = new Point3D(5.5, 0.5, 0.5);
	
	for (int i = 0; i < 5; i++) {
	    t2.addPoint(point);
	    point = point.displace(new Vector3D(-0.1, 0.2 * i, 0.0));
	}

	t2.joinTract(t1);
	
	tCurve = t2;

        TractStatisticFilter filter = new TractStatisticFilter(scalars.getDataDims(), scalars.getVoxelDims());

	filter.setTractStatistic("endpointsep");

	double ep = filter.processTract(tCurve)[0];

        // endpoint 1 = (5.9, 1.1, 0.5)
        // endpoint 2 = (5.1, 1.1, 0.5)

	assertEquals(0.8, ep, 1E-8);

	
    }

    
    
}

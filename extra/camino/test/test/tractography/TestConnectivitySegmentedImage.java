package tractography;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;


/**
 * Automated tests for <code>ConnectivitySegmentedImage</code>.
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 */
public class TestConnectivitySegmentedImage extends TestCase {

    // test image is 10x2x2 block, with targets at both ends

    private final int xDataDim = 10;
    private final int yDataDim = 2;
    private final int zDataDim = 2;

    private final double xVoxelDim = 1.0;
    private final double yVoxelDim = 1.5;
    private final double zVoxelDim = 2.0;

  
    // seeded in voxel 6, 0, 0, connects to targets 1 and 3
    private Tract tractSeed1_1 = null;

    // seeded in voxel 6, 0, 0, connects to targets 1 and 4
    private Tract tractSeed1_2 = null;


    // seeded in voxel 4, 1, 1, connects to targets 1 and 3
    private Tract tractSeed2_1 = null;

    // seeded in voxel 4, 1, 1, connects to target 3 only
    private Tract tractSeed2_2 = null;
    
    private final Point3D seedPoint1 = new Point3D(6.5, 0.5, 0.5);
    private final Point3D seedPoint2 = new Point3D(4.5, 2.0, 3.0);
    
    private final Voxel seedVoxel1 = new Voxel(6, 0, 0);
    private final Voxel seedVoxel2 = new Voxel(4, 1, 1);


    // these should not be modified
    private final int[][][] targets;
    
    private final PointListROI seeds;

	

    public TestConnectivitySegmentedImage(String name) {
	super(name);

        targets = new int[xDataDim][yDataDim][zDataDim];

        targets[0][0][0] = 1;
        targets[0][0][1] = 2;
        targets[0][1][0] = 1;
        targets[0][1][1] = 2;
        
        targets[9][0][0] = 3;
        targets[9][0][1] = 4;
        targets[9][1][0] = 3;
        targets[9][1][1] = 4;

        seeds = new PointListROI(new Point3D[] {seedPoint1, seedPoint2}, 1);
   
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
       
        tractSeed1_1 = new Tract();

        tractSeed1_1.addPoint(new Point3D(seedPoint1));
        tractSeed1_1.addPoint(new Point3D(9.5, 0.5, 0.5));

        Tract tmp = new Tract();

        tmp.addPoint(seedPoint1);
        tmp.addPoint(new Point3D(0.5, 0.5, 0.5));

        tractSeed1_1.joinTract(tmp);

        tractSeed1_1.resample(0.5);


        tractSeed1_2 = new Tract();

        tractSeed1_2.addPoint(new Point3D(seedPoint1));
        tractSeed1_2.addPoint(new Point3D(9.5, 0.5, 2.5));

        tmp = new Tract();

        tmp.addPoint(seedPoint1);
        tmp.addPoint(new Point3D(0.5, 0.5, 0.5));

        tractSeed1_2.joinTract(tmp);
        tractSeed1_2.resample(0.5);
        
        
        // now seed point 2


        tractSeed2_1 = new Tract();

        tractSeed2_1.addPoint(new Point3D(seedPoint2));
        tractSeed2_1.addPoint(new Point3D(9.5, 0.5, 1.5));

        tmp = new Tract();

        tmp.addPoint(seedPoint2);
        tmp.addPoint(new Point3D(0.5, 0.5, 0.5));

        tractSeed2_1.joinTract(tmp);

        tractSeed2_1.resample(0.5);


        tractSeed2_2 = new Tract();

        tractSeed2_2.addPoint(new Point3D(seedPoint2));
        tractSeed2_2.addPoint(new Point3D(9.5, 2.5, 0.5));

        tmp = new Tract();

        tmp.addPoint(seedPoint1);
        
        // not in any target region
        tmp.addPoint(new Point3D(2.5, 0.5, 0.5));

        tractSeed2_2.joinTract(tmp);
        tractSeed2_2.resample(0.5);
        


    }

    
    protected void tearDown() {

        tractSeed1_1 = null;
        tractSeed1_2 = null;
        tractSeed2_1 = null;
        tractSeed2_2 = null;

    }


    public static Test suite() {
	return new TestSuite(TestConnectivitySegmentedImage.class);
    }


  
    public void testSegmentationCountFirstEntry() {

        ConnectivitySegmentedImage cbs = new ConnectivitySegmentedImage(seeds, targets, xVoxelDim, yVoxelDim, zVoxelDim);

        cbs.setCountFirstEntry(true);

        TractCollection tc = new TractCollection();

        tc.addTract(tractSeed1_1);
        tc.addTract(tractSeed1_1);
        tc.addTract(tractSeed1_1);
        tc.addTract(tractSeed1_2);

        cbs.processTracts(0, tc);

        tc = new TractCollection();

        tc.addTract(tractSeed2_1);
        tc.addTract(tractSeed2_1);
        tc.addTract(tractSeed2_1);
        tc.addTract(tractSeed2_2);
	
        cbs.processTracts(1, tc);


        double[][][] seg = cbs.getSegmentedSeeds();

        
        for (int i = 0; i < xDataDim; i++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int k = 0; k < zDataDim; k++) {
                    
                    Voxel v = new Voxel(i,j,k);
                    
                    if (v.equals(seedVoxel1)) {
                        assertEquals(3.0, seg[i][j][k], 1E-8);
                    }
                    else if (v.equals(seedVoxel2)) {
                        assertEquals(1.0, seg[i][j][k], 1E-8);
                    }
                    else {
                        assertEquals(0.0, seg[i][j][k], 1E-8);
                    }
                }
            }
        }


        // now check probabilities
        double[][][] targetCP = cbs.getMaxTargetCP();

 
        for (int i = 0; i < xDataDim; i++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int k = 0; k < zDataDim; k++) {
                    
                    Voxel v = new Voxel(i,j,k);
                    
                    if (v.equals(seedVoxel1)) {
                        assertEquals(0.75, targetCP[i][j][k], 1E-8);
                    }
                    else if (v.equals(seedVoxel2)) {
                        assertEquals(0.75, targetCP[i][j][k], 1E-8);
                    }
                    else {
                        assertEquals(0.0, targetCP[i][j][k], 1E-8);
                    }
                }
            }
        }


    }





    public void testSegmentationCountAllEntries() {

        ConnectivitySegmentedImage cbs = new ConnectivitySegmentedImage(seeds, targets, xVoxelDim, yVoxelDim, zVoxelDim);

        cbs.setCountFirstEntry(false);

        TractCollection tc = new TractCollection();

        tc.addTract(tractSeed1_1);
        tc.addTract(tractSeed1_1);
        tc.addTract(tractSeed1_1);
        tc.addTract(tractSeed1_2);

        cbs.processTracts(0, tc);

        tc = new TractCollection();

        tc.addTract(tractSeed2_1);
        tc.addTract(tractSeed2_1);
        tc.addTract(tractSeed2_1);
        tc.addTract(tractSeed2_2);
	
        cbs.processTracts(1, tc);


        double[][][] seg = cbs.getSegmentedSeeds();

        
        for (int i = 0; i < xDataDim; i++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int k = 0; k < zDataDim; k++) {
                    
                    Voxel v = new Voxel(i,j,k);
                    
                    if (v.equals(seedVoxel1)) {
                        assertEquals(1.0, seg[i][j][k], 1E-8);
                    }
                    else if (v.equals(seedVoxel2)) {
                        assertEquals(3.0, seg[i][j][k], 1E-8);
                    }
                    else {
                        assertEquals(0.0, seg[i][j][k], 1E-8);
                    }
                }
            }
        }


        // now check probabilities
        double[][][] targetCP = cbs.getMaxTargetCP();

 
        for (int i = 0; i < xDataDim; i++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int k = 0; k < zDataDim; k++) {
                    
                    Voxel v = new Voxel(i,j,k);
                    
                    if (v.equals(seedVoxel1)) {
                        assertEquals(1.0, targetCP[i][j][k], 1E-8);
                    }
                    else if (v.equals(seedVoxel2)) {
                        assertEquals(1.0, targetCP[i][j][k], 1E-8);
                    }
                    else {
                        assertEquals(0.0, targetCP[i][j][k], 1E-8);
                    }
                }
            }
        }


    }



}




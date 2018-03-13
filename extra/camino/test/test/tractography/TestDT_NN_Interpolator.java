package tractography;

import junit.framework.*;
import junit.extensions.*;

import data.*;
import misc.DT;
import numerics.*;

import java.util.Random;


/**
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestDT_NN_Interpolator extends TestCase {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestDT_NN_Interpolator(String name) {
	super(name);
	
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestDT_NN_Interpolator.class);
    }


    public void testInterpolation() {
        
        DT_TractographyImage crossing = Images.getCrossing();
        
        DT_NN_Interpolator interpolator = new DT_NN_Interpolator(crossing);
        
        Vector3D[] pds = crossing.getPDs(1,1,0);

        DT[] dts = crossing.getDTs(1,1,0);
     
        for (int i = 0; i < 2; i++) { 

            Point3D p = new Point3D(1.5 * Images.xVoxelDim, 1.5 * Images.yVoxelDim, 0.5 * Images.zVoxelDim);
            
            DT d = interpolator.getDT(p, pds[i]);
            
            Images.assertTensorsEqual(this, dts[i], d, 1E-16);

        }


    }




}
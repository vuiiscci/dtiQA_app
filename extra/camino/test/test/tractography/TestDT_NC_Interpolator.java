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
public class TestDT_NC_Interpolator extends TestCase {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestDT_NC_Interpolator(String name) {
	super(name);
	
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestDT_NC_Interpolator.class);
    }


    public void testInterpolation() {
        
        DT_TractographyImage crossing = Images.getCrossing();
        
        DT_NC_Interpolator interpolator = new DT_NC_Interpolator(crossing, new Random(12345));
        
        Vector3D[] pds = crossing.getPDs(1,1,0);

        DT[] dts = crossing.getDTs(1,1,0);
     
        for (int i = 0; i < 2; i++) { 

            Point3D p = Images.getPointAtVoxelCentre(1,1,0);
            
            DT d = interpolator.getDT(p, pds[i]);
            
            Images.assertTensorsEqual(this, dts[i], d, 1E-16);

        }


    }




}
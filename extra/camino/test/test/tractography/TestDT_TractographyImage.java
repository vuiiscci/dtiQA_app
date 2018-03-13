package tractography;

import data.*;
import misc.DT;
import numerics.*;

import junit.framework.*;
import junit.extensions.*;



/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>DT_TractographyImage</code>.
 * <BR><BR>
 *
 * </dl>
 *
 * @version $Id: TestDT_TractographyImage.java,v 1.4 2005/11/25 18:50:47 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestDT_TractographyImage extends TestTractographyImage {

    
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {
    }


    public TestDT_TractographyImage(String name) {
	super(name);

    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestDT_TractographyImage.class);
    }


    public void testDTIsotropicMask() {

	DT_TractographyImage crossing = Images.getCrossing();

	// should match voxelClassification

	boolean[][][] isotropic = null;
	
	// should all be isotropic
	crossing.computeIsotropicMask(0.9);
	isotropic = crossing.getIsotropicMask();

	for (int k = 0; k < 2; k++) {
	    assertTrue(isotropic[0][0][k]);
	    assertTrue(isotropic[1][0][k]);
	    assertTrue(isotropic[2][0][k]);
	    
	    assertTrue(isotropic[0][1][k]);
	    assertTrue(isotropic[1][1][k]);
	    assertTrue(isotropic[2][1][k]);
	    
	    assertTrue(isotropic[0][2][k]);
	    assertTrue(isotropic[1][2][k]);
	    assertTrue(isotropic[2][2][k]);
	}


	// background should be isotropic, whatever the FA
	crossing.computeIsotropicMask(0.01);
	isotropic = crossing.getIsotropicMask();

	for (int k = 0; k < 2; k++) {
	    assertTrue(!isotropic[0][0][k]);
	    assertTrue(!isotropic[1][0][k]);
	    assertTrue(!isotropic[2][0][k]);
	    
	    assertTrue(!isotropic[0][1][k]);
	    assertTrue(!isotropic[1][1][k]);
	    assertTrue(!isotropic[2][1][k]);
	    
	    assertTrue(isotropic[0][2][k]);
	    assertTrue(!isotropic[1][2][k]);
	    assertTrue(isotropic[2][2][k]);

	}


    }


    public void testGetMix() {
	DT_TractographyImage crossing = Images.getCrossing();
        
        double[] mix = crossing.getMix(0,0,0);
        assertEquals(1, mix.length);
        assertEquals(1.0, mix[0], 1E-6);
        
        mix = crossing.getMix(0,1,0);
        assertEquals(1, mix.length);
        assertEquals(1.0, mix[0], 1E-6);
        
   
        // check that background voxels don't have nulls
        mix = crossing.getMix(0,2,0);
        assertEquals(0, mix.length);
           
        mix = crossing.getMix(1,0,0);
        assertEquals(1, mix.length);
        assertEquals(1.0, mix[0], 1E-6);
        
        mix = crossing.getMix(1,1,0);
        assertEquals(2, mix.length);
        assertEquals(0.5, mix[0], 1E-6);
        assertEquals(0.5, mix[1], 1E-6);
        
        mix = crossing.getMix(1,2,0);
        assertEquals(1, mix.length);
        assertEquals(1.0, mix[0], 1E-6);
        
        mix = crossing.getMix(2,0,0);
        assertEquals(1, mix.length);
        assertEquals(1.0, mix[0], 1E-6);
        
        mix = crossing.getMix(2,1,0);
        assertEquals(1, mix.length);
        assertEquals(1.0, mix[0], 1E-6);
        
        mix = crossing.getMix(2,2,0);
        assertEquals(0, mix.length);
      
        
        
    }
    
    protected TractographyImage getCrossingImage() {
	return Images.getCrossing();
    }
    
}

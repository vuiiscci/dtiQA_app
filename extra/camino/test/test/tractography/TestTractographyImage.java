package tractography;

import data.*;
import misc.DT;
import numerics.*;

import junit.framework.*;
import junit.extensions.*;



/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>TractographyImage</code>.
 * <BR><BR>
 *
 * </dl>
 *
 * @version $Id: TestTractographyImage.java,v 1.3 2005/10/24 17:46:48 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public abstract class TestTractographyImage extends TestCase {

    
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {
    }


    public TestTractographyImage(String name) {
	super(name);

    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestTractographyImage.class);
    }


    public void testIsotropicMask() {

	TractographyImage crossing = getCrossingImage();

	// background should be isotropic
	
	boolean[][][] isotropic = crossing.getIsotropicMask();

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



	double[][][] anisMap = new double[3][3][2];

	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {
		anisMap[i][j][0] = 0.4;
		anisMap[i][j][1] = 0.8;
	    }
	}

	crossing.computeIsotropicMask(anisMap, 0.5);
	isotropic = crossing.getIsotropicMask();


        for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 2; j++) {
                assertTrue(isotropic[i][j][0]);
            }
        }    
     
        assertTrue(!isotropic[0][0][1]);
        assertTrue(!isotropic[1][0][1]);
        assertTrue(!isotropic[2][0][1]);
        
        assertTrue(!isotropic[0][1][1]);
        assertTrue(!isotropic[1][1][1]);
        assertTrue(!isotropic[2][1][1]);
        
        assertTrue(isotropic[0][2][1]);
        assertTrue(!isotropic[1][2][1]);
        assertTrue(isotropic[2][2][1]);
        
	
    }


    public void testGetPDs() {

	TractographyImage crossing = getCrossingImage();

	Vector3D[] x = crossing.getPDs(0,1,0);

	assertEquals(1, x.length);

	assertEquals(1.0, Math.abs(x[0].dot(new Vector3D(1.0,0.0,0.0))), 0.00000001);

	// try crossing

	Vector3D[] c = crossing.getPDs(1,1,0);
	
	assertEquals(2, c.length);

	assertEquals(1.0, Math.abs(c[0].dot(new Vector3D(1.0,0.0,0.0))), 0.00000001);
	assertEquals(1.0, Math.abs(c[1].dot(new Vector3D(0.0,1.0,0.0))), 0.00000001);
		     

    }

  
    public void testNumberOfPDs() {

	TractographyImage crossing = getCrossingImage();

	int x = crossing.numberOfPDs(0,1,0);

	assertEquals(1, x);

	int c = crossing.numberOfPDs(1,1,0);
	
	assertEquals(2, c);

	int y = crossing.numberOfPDs(0,2,0);

	assertEquals(0, y);
	

    }


    /**
     * @return the crossing image in either PD, PICo, or DT format.
     */
    protected abstract TractographyImage getCrossingImage();



}

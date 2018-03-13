package tractography;

import data.*;
import misc.DT;
import numerics.*;

import junit.framework.*;
import junit.extensions.*;



/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>PICoTractographyImage</code>.
 * <BR><BR>
 *
 * </dl>
 *
 * @version $Id: TestPICoTractographyImage.java,v 1.1 2005/10/24 13:35:47 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestPICoTractographyImage extends TestTractographyImage {

    
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {
    }


    public TestPICoTractographyImage(String name) {
	super(name);

    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestPICoTractographyImage.class);
    }


    public void testPICoParams() {

        PICoTractographyImage crossing = (PICoTractographyImage)getCrossingImage();

        double[] pdfParams = crossing.getPICoPDFParams(0,0,0);

        assertEquals(1, pdfParams.length);

        assertEquals(100.0, pdfParams[0], 1E-6);


        pdfParams = crossing.getPICoPDFParams(1,1,0);

        assertEquals(2, pdfParams.length);

        assertEquals(50.0, pdfParams[0], 1E-6);
        assertEquals(50.0, pdfParams[1], 1E-6);
    }


    public void testGetPDs() {

	TractographyImage crossing = getCrossingImage();

	Vector3D[] x = crossing.getPDs(0,1,0);

	assertEquals(1, x.length);

	assertEquals(1.0, Math.abs(x[0].dot(new Vector3D(-0.9985213194730316, -0.025976007103758057,
							 -0.047753760195206511))), 0.00000001);
	
	// try crossing

	Vector3D[] c = crossing.getPDs(1,1,0);

	assertEquals(2, c.length);

	assertEquals(1.0, Math.abs(c[0].dot(new Vector3D(-0.9927563964927568, 0.09531088185625143,
							 -0.07314761118792366))), 0.00000001);
	
	assertEquals(1.0, Math.abs(c[1].dot(new Vector3D(-0.025908093080063373, 0.9985075537020853, 
							 0.048077395029597354))), 0.00000001);
		     

    }



    /**
     * Make sure that PDs are randomized when necessary and only then
     */
    public void testReset() {
        
        
        PICoTractographyImage image = (PICoTractographyImage)getCrossingImage();
        
	long seed = 562l;

        Vector3D x = image.getPDs(1,1,1)[0];

        Vector3D y = image.getPDs(1,1,1)[0];

        // should be the same
        assertEquals(1.0, x.dot(y), 1E-6);

        
        image.resetRandomization();

        
        Vector3D z = image.getPDs(1,1,1)[0];

        assertTrue(x.dot(y) > z.dot(y));
        
        
    }

    
    protected TractographyImage getCrossingImage() {
	return Images.getCrossingWatson();
    }

    
   
}

package tractography;

import junit.framework.*;
import junit.extensions.*;

import data.*;
import misc.DT;
import numerics.*;
import java.util.Random;

import optimizers.*;


/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>PICoBinghamRandomizer</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>PICoBinghamRandomizer</code> with JUnit 3.8. 
 * </dl>
 *
 * @version $Id: TestPICoBinghamRandomizer.java,v 1.2 2005/10/24 13:34:58 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestPICoBinghamRandomizer extends TestCase {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestPICoBinghamRandomizer(String name) {
	super(name);
   }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestPICoBinghamRandomizer.class);
    }



    public void testOneFibreDistribution() {

	long seed = 973l;

	PICoTractographyImage image = Images.getCrossingBingham();

	PICoBinghamRandomizer randomizer = new PICoBinghamRandomizer(image, new Random(seed));
	
	// sufficient to check the PDFs are initialized properly
	
	BinghamDistribution pdf = (BinghamDistribution)randomizer.getPDFs(0,1,0)[0];

	assertEquals(-200.0, pdf.k1(), 1E-6);
	assertEquals(-100.0, pdf.k2(), 1E-6);
	
	assertEquals(1.0, Math.abs(image.getEigenvectors(0,1,0)[0].dot(pdf.e3())), 1E-6);
    }

    
    public void testTwoFibreDistribution() {
	long seed = 973l;

	PICoTractographyImage image = Images.getCrossingBingham();

	PICoBinghamRandomizer randomizer = new PICoBinghamRandomizer(image, new Random(seed));
	
	// sufficient to check the PDFs are initialized properly
	
	AxialDistribution[] tmp = randomizer.getPDFs(1,1,0);

	BinghamDistribution[] pdfs = new BinghamDistribution[2];

	pdfs[0] = (BinghamDistribution)tmp[0];
	pdfs[1] = (BinghamDistribution)tmp[1];

	assertEquals(-100.0, pdfs[0].k1(), 1E-6);
	assertEquals(-50.0, pdfs[0].k2(), 1E-6);

	assertEquals(-100.0, pdfs[1].k1(), 1E-6);
	assertEquals(-50.0, pdfs[1].k2(), 1E-6);
	
	
	if ( Math.abs(image.getEigenvectors(1,1,0)[0].dot(pdfs[0].e3())) > 0.9 ) {
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[0].dot(pdfs[0].e3())), 1E-6);
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[1].dot(pdfs[1].e3())), 1E-6);
	}
	else {
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[1].dot(pdfs[0].e3())), 1E-6);
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[0].dot(pdfs[1].e3())), 1E-6);
	}

	Vector3D normal = pdfs[0].e3().cross(pdfs[1].e3()).normalized();

	assertEquals(1.0, Math.abs(normal.dot(pdfs[0].e1())), 1E-6);
	assertEquals(1.0, Math.abs(normal.dot(pdfs[1].e1())), 1E-6);

	Vector3D e2_1 = Rotations.rotateVector(normal, pdfs[0].e3(), Math.PI / 2.0);
	Vector3D e2_2 = Rotations.rotateVector(normal, pdfs[1].e3(), Math.PI / 2.0);

	assertEquals(1.0, Math.abs(e2_1.dot(pdfs[0].e2())), 1E-6);
	assertEquals(1.0, Math.abs(e2_2.dot(pdfs[1].e2())), 1E-6);
	
 
    }


}

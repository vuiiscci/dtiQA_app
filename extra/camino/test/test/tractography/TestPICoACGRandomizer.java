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
 * <dt>Purpose: Automated tests for <code>PICoACGRandomizer</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>PICoACGRandomizer</code> with JUnit 3.8. 
 * </dl>
 *
 * @version $Id: TestPICoACGRandomizer.java,v 1.2 2005/10/24 13:34:58 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestPICoACGRandomizer extends TestCase {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestPICoACGRandomizer(String name) {
	super(name);
   }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestPICoACGRandomizer.class);
    }



    public void testOneFibreDistribution() {

	long seed = 973l;

	PICoTractographyImage image = Images.getCrossingACG();

	PICoACGRandomizer randomizer = new PICoACGRandomizer(image, new Random(seed));
	
	// sufficient to check the PDFs are initialized properly
	
	ACG_Distribution pdf = (ACG_Distribution)randomizer.getPDFs(0,1,0)[0];

	double[] sigmas = pdf.sigmas();

	Vector3D[] evecs = pdf.eigenvectors();

	assertEquals(Math.sqrt(2.942083193342585), sigmas[0], 1E-6);
	assertEquals(Math.sqrt(0.08377478758707287), sigmas[1], 1E-6);
	assertEquals(Math.sqrt(0.008377478758707287), sigmas[2], 1E-6);
	
	assertEquals(1.0, Math.abs(image.getEigenvectors(0,1,0)[0].dot(evecs[0])), 1E-6);
    }

    
    public void testTwoFibreDistribution() {
	long seed = 973l;

	PICoTractographyImage image = Images.getCrossingACG();

	PICoACGRandomizer randomizer = new PICoACGRandomizer(image, new Random(seed));
	
	// sufficient to check the PDFs are initialized properly
	
	ACG_Distribution[] pdfs = (ACG_Distribution[])randomizer.getPDFs(1,1,0);

	double[] sigmas = pdfs[0].sigmas();

	assertEquals(Math.sqrt(2.942083193342585), sigmas[0], 1E-6);
	assertEquals(Math.sqrt(0.08377478758707287), sigmas[1], 1E-6);
	assertEquals(Math.sqrt(0.008377478758707287), sigmas[2], 1E-6);

	sigmas = pdfs[1].sigmas();

	assertEquals(Math.sqrt(2.942083193342585), sigmas[0], 1E-6);
	assertEquals(Math.sqrt(0.08377478758707287), sigmas[1], 1E-6);
	assertEquals(Math.sqrt(0.008377478758707287), sigmas[2], 1E-6);

	
	Vector3D[] evecs1 = pdfs[0].eigenvectors();
	Vector3D[] evecs2 = pdfs[1].eigenvectors();
	
	if ( Math.abs(image.getEigenvectors(1,1,0)[0].dot(evecs1[0])) > 0.9 ) {
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[0].dot(evecs1[0])), 1E-6);
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[1].dot(evecs2[0])), 1E-6);
	}
	else {
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[0].dot(evecs2[0])), 1E-6);
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[1].dot(evecs1[0])), 1E-6);
	}

	Vector3D normal = evecs1[0].cross(evecs2[0]).normalized();

	assertEquals(1.0, Math.abs(normal.dot(evecs1[2])), 1E-6);
	assertEquals(1.0, Math.abs(normal.dot(evecs2[2])), 1E-6);

	Vector3D e2_1 = Rotations.rotateVector(normal, evecs1[0], Math.PI / 2.0);
	Vector3D e2_2 = Rotations.rotateVector(normal, evecs2[0], Math.PI / 2.0);

	assertEquals(1.0, Math.abs(e2_1.dot(evecs1[1])), 1E-6);
	assertEquals(1.0, Math.abs(e2_2.dot(evecs2[1])), 1E-6);
	
 
    }


}

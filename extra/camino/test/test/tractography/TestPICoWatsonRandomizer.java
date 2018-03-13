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
 * <dt>Purpose: Automated tests for <code>PICoWatsonRandomizer</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>PICoWatsonRandomizer</code> with JUnit 3.8. 
 * </dl>
 *
 * @version $Id: TestPICoWatsonRandomizer.java,v 1.2 2005/10/24 13:34:58 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestPICoWatsonRandomizer extends TestCase {

    
 
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {

    }


    public TestPICoWatsonRandomizer(String name) {
	super(name);
   }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestPICoWatsonRandomizer.class);
    }



    public void testOneFibreDistribution() {

	long seed = 973l;

	PICoTractographyImage image = Images.getCrossingWatson();

	PICoWatsonRandomizer randomizer = new PICoWatsonRandomizer(image, new Random(seed));

	// sufficient to check the PDFs are initialized properly
	
	WatsonDistribution pdf = (WatsonDistribution)randomizer.getPDFs(0,1,0)[0];

	assertEquals(100.0, pdf.concentration(), 1E-6);
	assertEquals(1.0, Math.abs(image.getEigenvectors(0,1,0)[0].dot(pdf.meanAxis())), 1E-6);
    }

    
    public void testTwoFibreDistribution() {
	
	long seed = 973l;
	
 	PICoTractographyImage image = Images.getCrossingWatson();
	
 	PICoWatsonRandomizer randomizer = new PICoWatsonRandomizer(image, new Random(seed));
	
	

	// sufficient to check the PDFs are initialized properly
	
	WatsonDistribution[] pdfs = (WatsonDistribution[])randomizer.getPDFs(1,1,0);

	assertEquals(50.0, pdfs[0].concentration(), 1E-6);
	assertEquals(50.0, pdfs[1].concentration(), 1E-6);
	
	
	if ( Math.abs(image.getEigenvectors(1,1,0)[0].dot(pdfs[0].meanAxis())) > 0.9 ) {
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[0].dot(pdfs[0].meanAxis())), 1E-6);
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[1].dot(pdfs[1].meanAxis())), 1E-6);
	}
	else {
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[1].dot(pdfs[0].meanAxis())), 1E-6);
	    assertEquals(1.0, Math.abs(image.getEigenvectors(1,1,0)[0].dot(pdfs[1].meanAxis())), 1E-6);
	}
 
    }


}

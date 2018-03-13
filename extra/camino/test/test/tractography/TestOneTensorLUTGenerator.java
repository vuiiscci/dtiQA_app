package tractography;

import data.*;
import imaging.*;
import inverters.ModelIndex;
import misc.DT;
import numerics.*;
import optimizers.*;
import tools.*;

import java.util.Random;

import junit.framework.*;
import junit.extensions.*;



/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>OneTensorLUTGenerator</code>.
 * <BR><BR>
 *
 * </dl>
 *
 * @version $Id: TestOneTensorLUTGenerator.java,v 1.2 2005/09/29 12:44:06 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestOneTensorLUTGenerator extends TestCase {

    private Random ran = null;

    private OneTensorLUTGenerator generator = null;

    
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	// same sequence of random numbers for all tests
	ran = new Random(12345);
	generator = new OneTensorLUTGenerator(DW_Scheme.readScheme("bmx6.scheme"), 20.0, 2100.0E-12, ran);
    }

    
    protected void tearDown() {
	ran = null;
	generator = null;
    }


    public TestOneTensorLUTGenerator(String name) {
	super(name);

    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestOneTensorLUTGenerator.class);
    }

  
    public void testGenerateLUT() {

	double[][][][] lut = generator.generateLUT(1.0, 4.0, 1.0, 100, ModelIndex.LDT, true, false, false);

	assertEquals(4, lut[OneTensorLUTGenerator.WATSON].length);
	assertEquals(4, lut[OneTensorLUTGenerator.WATSON][0].length);
	assertEquals(1, lut[OneTensorLUTGenerator.WATSON][0][0].length);

	double[][] expected = {{-0.856248507098973, 0.0, 0.0, 0.0},
			       {173.0620976580426, -161.42410322540346, 0.0, 0.0},
			       {371.441864627901, 119.85329258808115, -342.1335310733197, 0.0},
			       {406.96067579396345, 218.97405002804373, 62.22365768603128, 
				-367.27010492409426}};

	for (int i = 0; i < 4; i++) {
	    for (int j = 0; j < 4; j++) {
		assertEquals(expected[i][j], lut[OneTensorLUTGenerator.WATSON][i][j][0], 1E-6);
	    }
	}
		
    }

}

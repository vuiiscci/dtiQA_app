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
 * @version $Id: TestTwoTensorLUTGenerator.java,v 1.6 2006/03/08 11:41:40 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTwoTensorLUTGenerator extends TestCase {

    private Random ran = null;

    private TwoTensorLUTGenerator generator = null;

    
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	// same sequence of random numbers for all tests
	ran = new Random(12345);
	generator = new TwoTensorLUTGenerator(DW_Scheme.readScheme("bmx6.scheme"), 
						    40.0, 2100.0E-12, 0.5, Math.PI / 2.0, ran);
    }

    
    protected void tearDown() {
	ran = null;
	generator = null;
    }


    public TestTwoTensorLUTGenerator(String name) {
	super(name);

    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestTwoTensorLUTGenerator.class);
    }

  
    public void testGenerateLUT() {

	double[][][][] lut = generator.generateLUT(0.6, 0.8, 0.1, 200, ModelIndex.POSPOS, true, false, false);


	assertEquals(3, lut[0].length);
	assertEquals(3, lut[0][0].length);
	assertEquals(2, lut[0][0][0].length);

	double[][][] expected = {
	    {{43.312512915518845, 41.39892937150715}, {81.4574928613992, 177.5}, 
	     {116.79034501529844, 339.9847788917596}},
	    {{177.5, 81.4574928613992}, {175.62425427210778, 198.75670236692326}, 
	     {244.21205927636046, 390.77360575766716}},
	    {{339.9847788917596, 116.79034501529844}, {390.77360575766716, 244.21205927636046},
	     {409.17, 482.34}}
	};

	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {

		assertEquals(expected[i][j][0], lut[0][i][j][0], 1.0);
		assertEquals(expected[i][j][1], lut[0][i][j][1], 1.0);
	    }
	   
	}
    }



    public void testFitting() {
        
        Vector3D[] sampleVecs = new Vector3D[500];

        Random ran = new Random(2345);

        Vector3D mu1 = new Vector3D(0.0, 0.0, 1.0);
        Vector3D mu2 = new Vector3D(1.0, 0.0, 0.0);
        

        WatsonDistribution w1 = new WatsonDistribution(mu1, ran);
        WatsonDistribution w2 = new WatsonDistribution(mu2, ran);

        for (int i = 0; i < sampleVecs.length / 2; i++) {
            sampleVecs[2*i] = w1.nextVector(200.0);
            sampleVecs[2*i+1] = w2.nextVector(100.0);
        }
        
        double[] kappas = generator.getWatsonConcentrationParams(sampleVecs, mu1, mu2);

        assertEquals(200.0, kappas[0], 20.0);
        assertEquals(100.0, kappas[1], 10.0);
        

        kappas = generator.getWatsonConcentrationParams(sampleVecs, mu2, mu1);
         
        assertEquals(100.0, kappas[0], 10.0);
        assertEquals(200.0, kappas[1], 20.0);


        kappas = generator.getBinghamConcentrationParams(sampleVecs, mu1, mu2);

        assertEquals(-200.0, kappas[0], 20.0);
        assertEquals(-200.0, kappas[1], 20.0);

        assertEquals(-100.0, kappas[2], 10.0);
        assertEquals(-100.0, kappas[3], 10.0);


        kappas = generator.getBinghamConcentrationParams(sampleVecs, mu2, mu1);

        assertEquals(-100.0, kappas[0], 10.0);
        assertEquals(-100.0, kappas[1], 10.0);

        assertEquals(-200.0, kappas[2], 20.0);
        assertEquals(-200.0, kappas[3], 20.0);


        double[] sigmas = generator.getACGConcentrationParams(sampleVecs, mu1, mu2);

        assertEquals(2.99276, sigmas[0], 1E-4);
        assertEquals(0.00407, sigmas[1], 1E-4);
        assertEquals(0.00317, sigmas[2], 1E-4);
        assertEquals(2.98474, sigmas[3], 1E-4);
        assertEquals(0.00790, sigmas[4], 1E-4);   
        assertEquals(0.00736, sigmas[5], 1E-4);

        sigmas = generator.getACGConcentrationParams(sampleVecs, mu2, mu1);
   
        assertEquals(2.98474, sigmas[0], 1E-4);
        assertEquals(0.00790, sigmas[1], 1E-4);   
        assertEquals(0.00736, sigmas[2], 1E-4);
        assertEquals(2.99276, sigmas[3], 1E-4);
        assertEquals(0.00407, sigmas[4], 1E-4);
        assertEquals(0.00317, sigmas[5], 1E-4);

       
  
    }


}

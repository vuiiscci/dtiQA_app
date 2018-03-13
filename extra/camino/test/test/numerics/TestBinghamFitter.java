package numerics;

import junit.framework.*;
import junit.extensions.*;
import numerics.*;

import java.util.Random;


/**
 * <dl>
 * <dt>Purpose: Automated tests for fitting code. 
 * <BR><BR>
 *
 * 
 * </dl>
 *
 * @version $Id: TestBinghamFitter.java,v 1.1 2005/09/26 10:36:46 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestBinghamFitter extends TestCase {

    private static final long seed = 183842l;
    


    public TestBinghamFitter(String name) {
	super(name);
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    

    public static Test suite() {
	return new TestSuite(TestBinghamFitter.class);
    }



    public void testFitting() {
	
	// example from Mardia
	try {
	    double[] bingPars = BinghamFitter.bngpar(0.365, 0.139);
	    
	    assertEquals(-3.63, bingPars[0], 0.01);
	    assertEquals(-0.68, bingPars[1], 0.01);
	}
	catch (Exception e) {
	    fail("Unexpected exception " + e);
	}
    }




    public void testSampling() {

	
	double[][] trueK1 = {{-50.03613698812956, 0.0, 0.0, 0.0},
			     {-147.97030316143127, -151.41933691407397, 0.0, 0.0},
			     {-443.3758473897254, -430.66091841717787, -478.0454763317299, 0.0},
			     {-994.3105646710299,-1025.6296698434612,-1003.5208469036224,-1022.615630315132}};

	double[][] trueK2 = {{-49.76824018890192, 0.0, 0.0, 0.0},
			     {-50.230150320473584, -143.95845051066314, 0.0, 0.0},
			     {-50.21762865087338, -150.31513370160278, -453.65588502140633, 0.0},
			     {-49.4986382644791,-144.96586784562476,-452.22595048853753,-988.6404600133196}};

	

	// kappas of each distribution
	double[] kappas = {-50.0, -150.0, -450.0, -1000.0};

	int sampleSize = 5000;
	
	Vector3D mu1 = new Vector3D(1.0, 0.0, 0.0);
	Vector3D mu2 = new Vector3D(0.0, 1.0, 0.0);
	Vector3D mu3 = new Vector3D(0.0, 0.0, 1.0);

	java.util.Random ran = new java.util.Random(seed);

	for (int k1 = 0; k1 < kappas.length; k1++) {
	    for (int k2 = 0; k2 < kappas.length; k2++) {
		if (k1 >= k2) {
		    
		    try {

			BinghamDistribution bing = 
			    BinghamDistribution.getBinghamDistribution(new Vector3D[] {mu1, mu2, mu3},
								 kappas[k1], kappas[k2], 
								 ran);
			
			Vector3D[] samples = new Vector3D[sampleSize];
			
			for (int s = 0; s < sampleSize; s++) {
			    samples[s] = bing.nextVector();
			}
		
			EigenSystem3D eig = BinghamFitter.tBarEigenSystem(samples);

			double[] bingPars = BinghamFitter.bngpar(eig.eigenvalues[1], eig.eigenvalues[2]);
			assertEquals("k1, input k1 " + kappas[k1] + " ", trueK1[k1][k2], bingPars[0], 0.000000001);
			assertEquals("k2, input k2 " +  kappas[k2] + " ", trueK2[k1][k2], bingPars[1], 0.000000001);
	
		    }
		    catch (ConvergenceException e) {
			fail("Unexpected exception: " + e);
		    }

		}
		
	    }
	    
	}



    }



    public void testGetBinghamDistribution() {
       
       	Vector3D mu1 = new Vector3D(1.0, 0.0, 0.0);
	Vector3D mu2 = new Vector3D(0.0, 1.0, 0.0);
	Vector3D mu3 = new Vector3D(0.0, 0.0, 1.0);

	EigenSystem3D eig = new EigenSystem3D(new Vector3D[] {mu1, mu2, mu3}, 
					      new double[] {0.97, 0.02, 0.01});

	try {


	    BinghamDistribution b1 = BinghamFitter.getBinghamDistribution(eig, new Random(seed));

	
	    double t2 = eig.eigenvalues[1];
	    double t3 = eig.eigenvalues[2];
	
	    /* L1: */
	    double[] akfc = null;

	    akfc = BinghamFitter.bngpar(t2, t3);
	
	    BinghamDistribution b2 = BinghamDistribution.getBinghamDistribution(new Vector3D[] 
		{mu1, mu2, mu3}, akfc[0], akfc[1], new Random(seed));


	    for (int i = 0; i < 10; i++) {
	    
		Vector3D v1 = b1.nextVector();
		Vector3D v2 = b2.nextVector();

		assertEquals(1.0, v1.dot(v2), 0.000001);
		assertEquals(b1.pdf(v1), b2.pdf(v1), 0.00001);
	    }

	}
	catch (ConvergenceException e) {
	    fail(e.toString());
	}


    }

    
    public void testNormC() {

	try {
	    double[] bngpar = BinghamFitter.bngpar(0.05, 0.01);
	    
	    double bingC = BinghamFitter.bingc(bngpar[0], bngpar[1]);
	    
	    assertEquals(0.0, bngpar[2] - bingC, 1E-6);

	}
	catch (ConvergenceException e) {
	    fail(e.toString());
	}

    }

}

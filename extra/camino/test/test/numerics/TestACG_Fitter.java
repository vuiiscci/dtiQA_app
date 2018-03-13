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
 * @version $Id: TestACG_Fitter.java,v 1.1 2005/09/26 10:36:46 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestACG_Fitter extends TestCase {


    public TestACG_Fitter(String name) {
	super(name);
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    

    public static Test suite() {
	return new TestSuite(TestACG_Fitter.class);
    }


    public void testFitting() {

	Random ran = new Random(3319713);

	Vector3D e1 = new Vector3D(0.0, 0.0, 1.0);
	Vector3D e2 = new Vector3D(1.0, 0.0, 0.0);
	Vector3D e3 = new Vector3D(0.0, 1.0, 0.0);

	double sigma1Sq = 2.9873179232102363;
	double sigma2Sq = 0.009232599188102994;
	double sigma3Sq = 0.003449477601664368;

	ACG_Distribution acg = new ACG_Distribution(new Vector3D[] {e1, e2, e3}, new double[] {sigma1Sq, sigma2Sq, sigma3Sq}, ran);

	Vector3D[] samples = new Vector3D[5000];
	
	for (int i = 0; i < 5000; i++) {
	    samples[i] = acg.nextVector();
	}

	RealMatrix A = ACG_Fitter.findA(samples);

	EigenSystem3D eigA = EigenSystem3D.sort(A);

	assertEquals(eigA.eigenvalues[0], sigma1Sq, sigma1Sq * 0.01);
	assertEquals(eigA.eigenvalues[1], sigma2Sq, sigma2Sq * 0.01);
	assertEquals(eigA.eigenvalues[2], sigma3Sq, sigma2Sq * 0.01);
		
	assertEquals(1.0, Math.abs(eigA.eigenvectors[0].dot(e1)), 0.0005);
	assertEquals(1.0, Math.abs(eigA.eigenvectors[1].dot(e2)), 0.0005);
	assertEquals(1.0, Math.abs(eigA.eigenvectors[2].dot(e3)), 0.0005);

    }



}
		     

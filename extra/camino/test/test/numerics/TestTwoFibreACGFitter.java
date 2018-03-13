package numerics;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;
import optimizers.*;

import java.util.Random;

/**
 * <dl>
 * <dt>Purpose: Automated tests for TwoFibreACGFitter.
 * <BR><BR>
 *
 * 
 * </dl>
 *
 * @version $Id: TestTwoFibreACGFitter.java,v 1.1 2005/09/26 10:36:48 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTwoFibreACGFitter extends TestCase {

    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

    
    protected void tearDown() {

    }


    public TestTwoFibreACGFitter(String name) {
	super(name);
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    

    public static Test suite() {
	return new TestSuite(TestTwoFibreACGFitter.class);
    }

   

    /**
     * Test a simple example. This uses test data with a known solution. It verifies that the optimization
     * works OK and that the parameters are set properly.
     *
     */
    public void testOptimization() {

	Random ran = new Random(19713);

	Vector3D e1_1 = new Vector3D(0.0, 0.0, 1.0);
	Vector3D e2_1 = new Vector3D(1.0, 0.0, 0.0);
	Vector3D e3_1 = new Vector3D(0.0, 1.0, 0.0);

	Vector3D e1_2 = new Vector3D(1.0, 0.0, 0.0);
	Vector3D e2_2 = new Vector3D(0.0, 1.0, 0.0);
	Vector3D e3_2 = new Vector3D(0.0, 0.0, 1.0);

	double sigma1Sq = 2.942083193342585;
	double sigma2Sq = 0.08377478758707287;
	double sigma3Sq = 0.008377478758707287;

	ACG_Distribution acg1 = new ACG_Distribution(new Vector3D[] {e1_1, e2_1, e3_1}, new double[] {sigma1Sq, sigma2Sq, sigma3Sq}, ran);

	ACG_Distribution acg2 = new ACG_Distribution(new Vector3D[] {e1_2, e2_2, e3_2}, new double[] {sigma1Sq, sigma2Sq, sigma3Sq}, ran);

	Vector3D[] samples1 = new Vector3D[1000];
	Vector3D[] samples2 = new Vector3D[1000];

	Vector3D[] sampleVecs = new Vector3D[2000];
	
	for (int i = 0; i < 1000; i++) {
	    
	    samples1[i] = acg1.nextVector();
	    samples2[i] = acg2.nextVector();

	    sampleVecs[2*i] = samples1[i];
	    sampleVecs[2*i+1] = samples2[i];
	}

	

	RealMatrix A = ACG_Fitter.findA(samples1);
	RealMatrix A2 = ACG_Fitter.findA(samples2);

	Vector3D[] samples1Init = new Vector3D[1000];
	Vector3D[] samples2Init = new Vector3D[1000];

	for (int i = 0; i < 1000; i++) {
	    if (Math.abs(e1_1.dot(sampleVecs[2*i])) + Math.abs(e1_2.dot(sampleVecs[2*i+1])) > Math.abs(e1_2.dot(sampleVecs[2*i+1])) + Math.abs(e1_2.dot(sampleVecs[2*i]))) {
		samples1Init[i] = sampleVecs[2*i];
		samples2Init[i] = sampleVecs[2*i+1];
	    }
	    else {
		samples1Init[i] = sampleVecs[2*i+1];
		samples2Init[i] = sampleVecs[2*i];
	    }

	}


	RealMatrix aInit = ACG_Fitter.findA(samples1Init);
	RealMatrix a2Init = ACG_Fitter.findA(samples2Init);

	TwoFibreACGFitter fitter = new TwoFibreACGFitter(sampleVecs);
	double[] params = null;

	
	double[] aParams = new double[6];

	aParams[0] = aInit.entries[0][0];
	aParams[1] = aInit.entries[0][1];
	aParams[2] = aInit.entries[0][2];
	aParams[3] = aInit.entries[1][1];
	aParams[4] = aInit.entries[1][2];		
	aParams[5] = aInit.entries[2][2];

	double[] a2Params = new double[6];
	
	a2Params[0] = a2Init.entries[0][0];
	a2Params[1] = a2Init.entries[0][1];
	a2Params[2] = a2Init.entries[0][2];
	a2Params[3] = a2Init.entries[1][1];
	a2Params[4] = a2Init.entries[1][2];		
	a2Params[5] = a2Init.entries[2][2];

	try {
	    params = fitter.fitEstimatedParams(aParams, a2Params, 1E-8);
	}
	catch (ConjGradMinimizerException e) {
	    fail(e.toString());
	}
	

	RealMatrix aFitted = new RealMatrix(3,3);
	
	aFitted.setEntry(0,0, params[0]);
	aFitted.setEntry(0,1, params[1]);
	aFitted.setEntry(0,2, params[2]);
	aFitted.setEntry(1,0, params[1]);
	aFitted.setEntry(1,1, params[3]);
	aFitted.setEntry(1,2, params[4]);
	aFitted.setEntry(2,0, params[2]);
	aFitted.setEntry(2,1, params[4]);
	aFitted.setEntry(2,2, params[5]);
	
	RealMatrix a2Fitted = new RealMatrix(3,3);
	
	a2Fitted.setEntry(0,0, params[6]);
	a2Fitted.setEntry(0,1, params[7]);
	a2Fitted.setEntry(0,2, params[8]);
	a2Fitted.setEntry(1,0, params[7]);
	a2Fitted.setEntry(1,1, params[9]);
	a2Fitted.setEntry(1,2, params[10]);
	a2Fitted.setEntry(2,0, params[8]);
	a2Fitted.setEntry(2,1, params[10]);
	a2Fitted.setEntry(2,2, params[11]);
	
	EigenSystem3D eigAf = EigenSystem3D.sort(aFitted);
	EigenSystem3D eigA2f = EigenSystem3D.sort(a2Fitted);

	EigenSystem3D eigA = EigenSystem3D.sort(A);
	EigenSystem3D eigA2 = EigenSystem3D.sort(A2);
		

	
	if ( Math.abs(eigAf.eigenvectors[0].dot(eigA.eigenvectors[0])) <  Math.abs(eigAf.eigenvectors[0].dot(eigA2.eigenvectors[0])) ) {
	    // trade
	    EigenSystem3D tmp = eigA2f;
	    eigA2f = eigAf;
	    eigAf = tmp;
	}
	    
	assertEquals(eigA.eigenvalues[0], eigAf.eigenvalues[0], eigA.eigenvalues[0] * 0.05);
	assertEquals(eigA.eigenvalues[1], eigAf.eigenvalues[1], eigA.eigenvalues[1] * 0.05);
	assertEquals(eigA.eigenvalues[2], eigAf.eigenvalues[2], eigA.eigenvalues[2] * 0.05);
	
	assertEquals(eigA2.eigenvalues[0], eigA2f.eigenvalues[0], eigA2.eigenvalues[0] * 0.05);
	assertEquals(eigA2.eigenvalues[1], eigA2f.eigenvalues[1], eigA2.eigenvalues[1] * 0.05);
	assertEquals(eigA2.eigenvalues[2], eigA2f.eigenvalues[2], eigA2.eigenvalues[2] * 0.05);
	
	assertEquals(1.0, Math.abs(eigAf.eigenvectors[0].dot(eigA.eigenvectors[0])), 0.005);
	assertEquals(1.0, Math.abs(eigAf.eigenvectors[1].dot(eigA.eigenvectors[1])), 0.005);
	assertEquals(1.0, Math.abs(eigAf.eigenvectors[2].dot(eigA.eigenvectors[2])), 0.005);
	
	assertEquals(1.0, Math.abs(eigA2f.eigenvectors[0].dot(eigA2.eigenvectors[0])), 0.005);
	assertEquals(1.0, Math.abs(eigA2f.eigenvectors[1].dot(eigA2.eigenvectors[1])), 0.005);
	assertEquals(1.0, Math.abs(eigA2f.eigenvectors[2].dot(eigA2.eigenvectors[2])), 0.005);
	    

    }




}

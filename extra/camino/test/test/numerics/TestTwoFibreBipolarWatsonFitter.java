package numerics;

import junit.framework.*;
import junit.extensions.*;
import numerics.*;
import optimizers.*;
import java.util.Random;

/**
 * <dl>
 * <dt>Purpose: Automated tests for optimizer. 
 * <BR><BR>
 *
 * 
 * </dl>
 *
 * @version $Id: TestTwoFibreBipolarWatsonFitter.java,v 1.3 2006/02/16 21:33:25 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTwoFibreBipolarWatsonFitter extends TestTwoFibreWatsonFitter {


    public TestTwoFibreBipolarWatsonFitter(String name) {
	super(name);
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    

    public static Test suite() {
	return new TestSuite(TestTwoFibreBipolarWatsonFitter.class);
    }


    public void testFObj() {

	double[] initParams = new double[7];

	Random ran = new Random(4082);

	WatsonDistribution w1 = new WatsonDistribution(Rotations.Z_AXIS, 100.0,ran);
	WatsonDistribution w2 = new WatsonDistribution(Rotations.X_AXIS, 50.0,ran);

	Vector3D[] samples = new Vector3D[2000];

	Vector3D[] samples1 = new Vector3D[1000];
	Vector3D[] samples2 = new Vector3D[1000];


	for (int i = 0; i < 1000; i++) {
	    samples1[i] = w1.nextVector();
	    samples2[i] = w2.nextVector();
	}
	
	
	for (int i = 0; i < 1000; i++) {
	    samples[2*i] = samples1[i];
	    samples[2*i+1] = samples2[i];
	}

 	TwoFibreBipolarWatsonFitter fitter = new TwoFibreBipolarWatsonFitter(samples);

	EigenSystem3D eig1 = WatsonFitter.tBarEigenSystem(samples1);

	double[] tp = Vector3D.thetaPhi(eig1.eigenvectors[0]);

	initParams[1] = tp[0];
	initParams[2] = tp[1];
	initParams[3] = Math.sqrt(WatsonFitter.fitKappa(eig1, samples1));

	EigenSystem3D eig2 = WatsonFitter.tBarEigenSystem(samples2);

	tp = Vector3D.thetaPhi(eig2.eigenvectors[0]);

	initParams[4] = tp[0];
	initParams[5] = tp[1];
	initParams[6] = Math.sqrt(WatsonFitter.fitKappa(eig2, samples2)); 


	double f = 0.0;
	for (int i = 0; i < 2000; i++) {
	    f -= Math.log( (WatsonDistribution.pdf(samples[i], Vector3D.vectorFromSPC(1.0, initParams[1], initParams[2]), initParams[3]*initParams[3]) + WatsonDistribution.pdf(samples[i], Vector3D.vectorFromSPC(1.0, initParams[4], initParams[5]), initParams[6]*initParams[6])) / 2.0 );
	}


	double[] dfda = new double[7];
	double[][] d2fda2 = new double[7][7];

	double fObj = fitter.fObj(initParams, dfda , d2fda2);

	for (int i = 0; i < 7; i++) {
	    assertEquals(0.0, dfda[i], 0.00001);
	}

	assertEquals(f, fObj, 0.000000001);

	

    }


    // test other ways of fitting
    public void testFitEstimatedParamsCone() {
	
	try {
	    

	    TwoFibreBipolarWatsonFitter fitter = new TwoFibreBipolarWatsonFitter(samples);
	    
	    fitter.fitEstimatedParams(samples[0], samples[samples.length - 1], 4, 8, 4, 8, 50.0, 50.0, 5);
	    
	    double[] params = fitter.getParameters();
	    
	    double[] retParams = new double[7];
	    
	    for (int i = 0; i < 6; i++) {
		retParams[i] = params[i];
	    }
	    
	    // prop is fixed
	    retParams[6] = alpha;
	    
	    checkOptimizedParams(retParams);
	    
	}
	catch (MarquardtMinimiserException e) {
	    fail(e.toString());
	}
	
    }


    
    public void testFitEstimatedParamsAutoCone() {
    
	try {


	    TwoFibreBipolarWatsonFitter fitter = new TwoFibreBipolarWatsonFitter(samples);

	    fitter.fitEstimatedParams(samples[0], samples[samples.length - 1], 5);
	    
	    double[] params = fitter.getParameters();
	    
	    double[] retParams = new double[7];
	    
	    for (int i = 0; i < 6; i++) {
		retParams[i] = params[i];
	    }
	    
	    // prop is fixed
	    retParams[6] = alpha;
	    
	    checkOptimizedParams(retParams);

	}
	catch (MarquardtMinimiserException e) {
	    fail(e.toString());
	}

    }


   public void testHighConcentration() {
    
	try {

	    double[] initParams = new double[7];
	    
	    Random ran = new Random(4082);
	    
	    WatsonDistribution w1 = new WatsonDistribution(Rotations.Z_AXIS, 800.0,ran);
	    WatsonDistribution w2 = new WatsonDistribution(Rotations.X_AXIS, 600.0,ran);
	    
	    Vector3D[] sampleVecs = new Vector3D[1000];
	    
	    Vector3D[] samples1 = new Vector3D[500];
	    Vector3D[] samples2 = new Vector3D[500];

	    
	    for (int i = 0; i < 500; i++) {
		samples1[i] = w1.nextVector();
		samples2[i] = w2.nextVector();
	    }
	    
	    
	    for (int i = 0; i < 500; i++) {
		sampleVecs[2*i] = samples1[i];
		sampleVecs[2*i+1] = samples2[i];
	    }
	    
	    TwoFibreBipolarWatsonFitter fitter = new TwoFibreBipolarWatsonFitter(sampleVecs);

	    EigenSystem3D eig1 = WatsonFitter.tBarEigenSystem(samples1);
	    
	    double[] tp = Vector3D.thetaPhi(eig1.eigenvectors[0]);
	    
	    initParams[0] = tp[0];
	    initParams[1] = tp[1];
	    initParams[2] = WatsonFitter.fitKappa(eig1, samples1);
	    
	    EigenSystem3D eig2 = WatsonFitter.tBarEigenSystem(samples2);
	    
	    tp = Vector3D.thetaPhi(eig2.eigenvectors[0]);
	    
	    initParams[3] = tp[0];
	    initParams[4] = tp[1];
	    initParams[5] = WatsonFitter.fitKappa(eig2, samples2); 
	    	    
	    fitter.fitEstimatedParams(sampleVecs[0], sampleVecs[1], 5);
	    
	    double[] params = fitter.getParameters();

	    Vector3D mu1 = Vector3D.vectorFromSPC(1.0, params[0], params[1]);
	    Vector3D mu2 = Vector3D.vectorFromSPC(1.0, params[3], params[4]);
	    double kappa1 = params[2];
	    double kappa2 = params[5];

	    if (Math.abs(mu1.dot(eig1.eigenvectors[0])) > Math.abs(mu2.dot(eig1.eigenvectors[0]))) {
		assertEquals(1.0, Math.abs(mu1.dot(eig1.eigenvectors[0])), 0.00001);
		assertEquals(1.0, Math.abs(mu2.dot(eig2.eigenvectors[0])), 0.00001);

		assertEquals(initParams[2], kappa1, 0.0001);
		assertEquals(initParams[5], kappa2, 0.0001);

	    }
	    else {

		assertEquals(1.0, Math.abs(mu1.dot(eig2.eigenvectors[0])), 0.00001);
		assertEquals(1.0, Math.abs(mu2.dot(eig1.eigenvectors[0])), 0.00001);

		assertEquals(initParams[2], kappa2, 0.0001);
		assertEquals(initParams[5], kappa1, 0.0001);

	    }
	    

	}
	catch (MarquardtMinimiserException e) {
	    fail(e.toString());
	}

        // try with both concentrations just under limit

	try {

	    double[] initParams = new double[7];
	    
	    Random ran = new Random(4082);
	    
	    WatsonDistribution w1 = new WatsonDistribution(Rotations.Z_AXIS, 380.0,ran);
	    WatsonDistribution w2 = new WatsonDistribution(Rotations.X_AXIS, 380.0,ran);
	    
	    Vector3D[] sampleVecs = new Vector3D[1000];
	    
	    Vector3D[] samples1 = new Vector3D[500];
	    Vector3D[] samples2 = new Vector3D[500];

	    
	    for (int i = 0; i < 500; i++) {
		samples1[i] = w1.nextVector();
		samples2[i] = w2.nextVector();
	    }
	    
	    
	    for (int i = 0; i < 500; i++) {
		sampleVecs[2*i] = samples1[i];
		sampleVecs[2*i+1] = samples2[i];
	    }
	    
	    TwoFibreBipolarWatsonFitter fitter = new TwoFibreBipolarWatsonFitter(sampleVecs);

	    EigenSystem3D eig1 = WatsonFitter.tBarEigenSystem(samples1);
	    
	    double[] tp = Vector3D.thetaPhi(eig1.eigenvectors[0]);
	    
	    initParams[0] = tp[0];
	    initParams[1] = tp[1];
	    initParams[2] = WatsonFitter.fitKappa(eig1, samples1);
	    
	    EigenSystem3D eig2 = WatsonFitter.tBarEigenSystem(samples2);
	    
	    tp = Vector3D.thetaPhi(eig2.eigenvectors[0]);
	    
	    initParams[3] = tp[0];
	    initParams[4] = tp[1];
	    initParams[5] = WatsonFitter.fitKappa(eig2, samples2); 
	    	    
	    fitter.fitEstimatedParams(sampleVecs[0], sampleVecs[1], 5);
	    
	    double[] params = fitter.getParameters();

	    Vector3D mu1 = Vector3D.vectorFromSPC(1.0, params[0], params[1]);
	    Vector3D mu2 = Vector3D.vectorFromSPC(1.0, params[3], params[4]);
	    double kappa1 = params[2];
	    double kappa2 = params[5];

	    if (Math.abs(mu1.dot(eig1.eigenvectors[0])) > Math.abs(mu2.dot(eig1.eigenvectors[0]))) {
		assertEquals(1.0, Math.abs(mu1.dot(eig1.eigenvectors[0])), 0.00001);
		assertEquals(1.0, Math.abs(mu2.dot(eig2.eigenvectors[0])), 0.00001);

		assertEquals(initParams[2], kappa1, 2.0);
		assertEquals(initParams[5], kappa2, 2.0);

	    }
	    else {

		assertEquals(1.0, Math.abs(mu1.dot(eig2.eigenvectors[0])), 0.00001);
		assertEquals(1.0, Math.abs(mu2.dot(eig1.eigenvectors[0])), 0.00001);

		assertEquals(initParams[2], kappa2, 2.0);
		assertEquals(initParams[5], kappa1, 2.0);

	    }
	    

	}
	catch (MarquardtMinimiserException e) {
	    fail(e.toString());
	}

    }


	 
}


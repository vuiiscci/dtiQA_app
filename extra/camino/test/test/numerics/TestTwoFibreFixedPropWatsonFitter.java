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
 * @version $Id: TestTwoFibreFixedPropWatsonFitter.java,v 1.2 2006/02/16 21:33:25 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTwoFibreFixedPropWatsonFitter extends TestTwoFibreWatsonFitter {


    public TestTwoFibreFixedPropWatsonFitter(String name) {
	super(name);
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    

    public static Test suite() {
	return new TestSuite(TestTwoFibreFixedPropWatsonFitter.class);
    }


   
    public void testFObj() {

	double[] initParams = new double[7];

	Random ran = new Random(4082);

	WatsonDistribution w1 = new WatsonDistribution(Rotations.Z_AXIS, 200.0,ran);
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

 	TwoFibreFixedPropWatsonFitter fitter = new TwoFibreFixedPropWatsonFitter(samples, 0.5);

	EigenSystem3D eig1 = WatsonFitter.tBarEigenSystem(samples1);

	double[] tp = Vector3D.thetaPhi(eig1.eigenvectors[0]);

	initParams[1] = tp[0];
	initParams[2] = tp[1];
	initParams[3] = WatsonFitter.fitKappa(eig1, samples1);

	EigenSystem3D eig2 = WatsonFitter.tBarEigenSystem(samples2);

	tp = Vector3D.thetaPhi(eig2.eigenvectors[0]);

	initParams[4] = tp[0];
	initParams[5] = tp[1];
	initParams[6] = WatsonFitter.fitKappa(eig2, samples2); 


	double f = 0.0;
	for (int i = 0; i < 2000; i++) {
	    f -= Math.log( (WatsonDistribution.pdf(samples[i], Vector3D.vectorFromSPC(1.0, initParams[1], initParams[2]), initParams[3]) + WatsonDistribution.pdf(samples[i], Vector3D.vectorFromSPC(1.0, initParams[4], initParams[5]), initParams[6])) / 2.0 );
	}


	double[] dfda = new double[7];

	double fObj = fitter.fObj(initParams, dfda , new double[7][7]);

	for (int i = 0; i < 7; i++) {
	    assertEquals(0.0, dfda[i], 0.00001);
	}
	
	assertEquals(f, fObj, 0.000000001);

	

    }


    // test other ways of fitting
    public void testFitEstimatedParamsCone() {
	
	try {
	    

	    TwoFibreFixedPropWatsonFitter fitter = new TwoFibreFixedPropWatsonFitter(samples, alpha);
	    
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


	    TwoFibreFixedPropWatsonFitter fitter = new TwoFibreFixedPropWatsonFitter(samples, alpha);

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


	 
}


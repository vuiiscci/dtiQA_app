package numerics;

import junit.framework.*;
import junit.extensions.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for fitting code. 
 * <BR><BR>
 *
 * 
 * </dl>
 *
 * @version $Id: TestWatsonFitter.java,v 1.2 2005/10/29 00:51:14 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestWatsonFitter extends TestCase {


    public TestWatsonFitter(String name) {
	super(name);
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    

    public static Test suite() {
	return new TestSuite(TestWatsonFitter.class);
    }


    public void testFitKappaBP() {

	// dataset B18 from Fisher (1987)
	Vector3D[] samples = new Vector3D[33];

	samples[0] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -5.0 * Math.PI / 180.0,   12.0 * Math.PI / 180.0);
	samples[1] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -1.0 * Math.PI / 180.0,  17.0 * Math.PI / 180.0);
	samples[2] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -10.0 * Math.PI / 180.0,  9.0 * Math.PI / 180.0);
	samples[3] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 8.0 * Math.PI / 180.0,  342.0 * Math.PI / 180.0);
	samples[4] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 0.0 * Math.PI / 180.0,  12.0 * Math.PI / 180.0);
	samples[5] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 8.0 * Math.PI / 180.0,  0.0 * Math.PI / 180.0);
	samples[6] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -7.0 * Math.PI / 180.0,  4.0 * Math.PI / 180.0);
	samples[7] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -13.0 * Math.PI / 180.0,  2.0 * Math.PI / 180.0);
	samples[8] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 5.0 * Math.PI / 180.0,  4.0 * Math.PI / 180.0);
	samples[9] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -7.0 * Math.PI / 180.0,  4.0 * Math.PI / 180.0);
	samples[10] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 5.0 * Math.PI / 180.0,  15.0 * Math.PI / 180.0);

	samples[11] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -1.0 * Math.PI / 180.0,  2.0 * Math.PI / 180.0);
	samples[12] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -12.0 * Math.PI / 180.0,  353.0 * Math.PI / 180.0);
	samples[13] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 11.0 * Math.PI / 180.0,  350.0 * Math.PI / 180.0);
	samples[14] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 3.0 * Math.PI / 180.0,  9.0 * Math.PI / 180.0);
	samples[15] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 15.0 * Math.PI / 180.0,  355.0 * Math.PI / 180.0);
	samples[16] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 3.0 * Math.PI / 180.0,  344.0 * Math.PI / 180.0);
	samples[17] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 0.0 * Math.PI / 180.0,  2.0 * Math.PI / 180.0);
	samples[18] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -9.0 * Math.PI / 180.0,  359.0 * Math.PI / 180.0);
	samples[19] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 0.0 * Math.PI / 180.0,  12.0 * Math.PI / 180.0);
	samples[20] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -11.0 * Math.PI / 180.0,  10.0 * Math.PI / 180.0);
	samples[21] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -13.0 * Math.PI / 180.0,  14.0 * Math.PI / 180.0);

	samples[22] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 0.0 * Math.PI / 180.0,  12.0 * Math.PI / 180.0);
	samples[23] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 1.0 * Math.PI / 180.0,  13.0 * Math.PI / 180.0);
	samples[24] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 3.0 * Math.PI / 180.0,  353.0 * Math.PI / 180.0);
	samples[25] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 10.0 * Math.PI / 180.0,  3.0 * Math.PI / 180.0);
	samples[26] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 12.0 * Math.PI / 180.0,  4.0 * Math.PI / 180.0);
	samples[27] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 15.0 * Math.PI / 180.0,  347.0 * Math.PI / 180.0);
	samples[28] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -6.0 * Math.PI / 180.0,  2.0 * Math.PI / 180.0);
	samples[29] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -13.0 * Math.PI / 180.0,  7.0 * Math.PI / 180.0);
	samples[30] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + -17.0 * Math.PI / 180.0,  354.0 * Math.PI / 180.0);
	samples[31] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 25.0 * Math.PI / 180.0,  4.0 * Math.PI / 180.0);
	samples[32] = Vector3D.vectorFromSPC(1.0, Math.PI / 2.0 + 0.0 * Math.PI / 180.0,  4.0 * Math.PI / 180.0);

	assertEquals(20.0, WatsonFitter.fitKappa(WatsonFitter.tBarEigenSystem(samples), samples), 0.5);

    }


    public void testFitKappaGirdle() {

	
	Vector3D[] samples = new Vector3D[39];
	// dataset B12 from Fisher (1987)
	samples[0] = Vector3D.vectorFromSPC(1.0, 62.0 * Math.PI / 180.0, 128.0 * Math.PI / 180.0);
	samples[1] = Vector3D.vectorFromSPC(1.0, 53.0 * Math.PI / 180.0, 133.0 * Math.PI / 180.0);
	samples[2] = Vector3D.vectorFromSPC(1.0, 52.0 * Math.PI / 180.0, 129.0 * Math.PI / 180.0);
	samples[3] = Vector3D.vectorFromSPC(1.0, 48.0 * Math.PI / 180.0, 124.0 * Math.PI / 180.0);
	samples[4] = Vector3D.vectorFromSPC(1.0, 44.0 * Math.PI / 180.0, 141.0 * Math.PI / 180.0);
	samples[5] = Vector3D.vectorFromSPC(1.0, 37.0 * Math.PI / 180.0, 138.0 * Math.PI / 180.0);
	samples[6] = Vector3D.vectorFromSPC(1.0, 38.0 * Math.PI / 180.0, 135.0 * Math.PI / 180.0);
	samples[7] = Vector3D.vectorFromSPC(1.0, 25.0 * Math.PI / 180.0, 156.0 * Math.PI / 180.0);
	samples[8] = Vector3D.vectorFromSPC(1.0, 15.0 * Math.PI / 180.0, 156.0 * Math.PI / 180.0);
	samples[9] = Vector3D.vectorFromSPC(1.0, 22.0 * Math.PI / 180.0, 130.0 * Math.PI / 180.0);
	samples[10] = Vector3D.vectorFromSPC(1.0, 28.0 * Math.PI / 180.0, 113.0 * Math.PI / 180.0);
	samples[11] = Vector3D.vectorFromSPC(1.0, 22.0 * Math.PI / 180.0, 110.0 * Math.PI / 180.0);
	samples[12] = Vector3D.vectorFromSPC(1.0, 33.0 * Math.PI / 180.0, 106.0 * Math.PI / 180.0);
	samples[13] = Vector3D.vectorFromSPC(1.0, 24.0 * Math.PI / 180.0, 77.0 * Math.PI / 180.0);
	samples[14] = Vector3D.vectorFromSPC(1.0, 8.0 * Math.PI / 180.0, 111.0 * Math.PI / 180.0);
	samples[15] = Vector3D.vectorFromSPC(1.0, 6.0 * Math.PI / 180.0, 122.0 * Math.PI / 180.0);
	samples[16] = Vector3D.vectorFromSPC(1.0, 8.0 * Math.PI / 180.0, 140.0 * Math.PI / 180.0);
	samples[17] = Vector3D.vectorFromSPC(1.0, 11.0 * Math.PI / 180.0, 48.0 * Math.PI / 180.0);
	samples[18] = Vector3D.vectorFromSPC(1.0, 8.0 * Math.PI / 180.0, 28.0 * Math.PI / 180.0);
	samples[19] = Vector3D.vectorFromSPC(1.0, 20.0 * Math.PI / 180.0, 310.0 * Math.PI / 180.0);
	samples[20] = Vector3D.vectorFromSPC(1.0, 25.0 * Math.PI / 180.0, 326.0 * Math.PI / 180.0);
	samples[21] = Vector3D.vectorFromSPC(1.0, 28.0 * Math.PI / 180.0, 332.0 * Math.PI / 180.0);
	samples[22] = Vector3D.vectorFromSPC(1.0, 32.0 * Math.PI / 180.0, 308.0 * Math.PI / 180.0);
	samples[23] = Vector3D.vectorFromSPC(1.0, 34.0 * Math.PI / 180.0, 304.0 * Math.PI / 180.0);
	samples[24] = Vector3D.vectorFromSPC(1.0, 38.0 * Math.PI / 180.0, 304.0 * Math.PI / 180.0);
	samples[25] = Vector3D.vectorFromSPC(1.0, 45.0 * Math.PI / 180.0, 293.0 * Math.PI / 180.0);
	samples[26] = Vector3D.vectorFromSPC(1.0, 47.0 * Math.PI / 180.0, 313.0 * Math.PI / 180.0);
	samples[27] = Vector3D.vectorFromSPC(1.0, 45.0 * Math.PI / 180.0, 319.0 * Math.PI / 180.0);
	samples[28] = Vector3D.vectorFromSPC(1.0, 43.0 * Math.PI / 180.0, 320.0 * Math.PI / 180.0);
	samples[29] = Vector3D.vectorFromSPC(1.0, 59.0 * Math.PI / 180.0, 312.0 * Math.PI / 180.0);
	samples[30] = Vector3D.vectorFromSPC(1.0, 66.0 * Math.PI / 180.0, 317.0 * Math.PI / 180.0);
	samples[31] = Vector3D.vectorFromSPC(1.0, 65.0 * Math.PI / 180.0, 314.0 * Math.PI / 180.0);
	samples[32] = Vector3D.vectorFromSPC(1.0, 70.0 * Math.PI / 180.0, 312.0 * Math.PI / 180.0);
	samples[33] = Vector3D.vectorFromSPC(1.0, 66.0 * Math.PI / 180.0, 311.0 * Math.PI / 180.0);
	samples[34] = Vector3D.vectorFromSPC(1.0, 66.0 * Math.PI / 180.0, 310.0 * Math.PI / 180.0);
	samples[35] = Vector3D.vectorFromSPC(1.0, 72.0 * Math.PI / 180.0, 305.0 * Math.PI / 180.0);
	samples[36] = Vector3D.vectorFromSPC(1.0, 67.0 * Math.PI / 180.0, 301.0 * Math.PI / 180.0);
	samples[37] = Vector3D.vectorFromSPC(1.0, 69.0 * Math.PI / 180.0, 301.0 * Math.PI / 180.0);
	samples[38] = Vector3D.vectorFromSPC(1.0, 82.0 * Math.PI / 180.0, 300.0 * Math.PI / 180.0);

	assertEquals(-38.0, WatsonFitter.fitKappa(WatsonFitter.tBarEigenSystem(samples), samples), 0.5);
	
    }




    /**
     * Test that samples drawn from distribution are correctly distributed. 
     *
     */
    public void testSamplingLargeKappa() {
	
 	double minKappa = -1500.0;

 	double maxKappa = 1500.0;

 	double interval = 500.0;

 	int samples = 5000;

	long seed = 56267934101l;

 	// kappas of each distribution
 	double[] kappas = new double[1 + (int)((maxKappa - minKappa) / interval)];

	for (int i = 0; i < kappas.length; i++) {
	    
	    kappas[i] = minKappa + (maxKappa - minKappa) * (double)(i) / (kappas.length - 1.0);
    
	}


	// run with sphereDistFit.TestWatsonFitting
	double[] trueKappas = {-1474.8427964432933, -993.888911272413, -503.63570534999815, -0.09146032182611312, 498.0639110203141, 984.8447713171712, 1505.5639265498291};
	
	double[] goodness = {178236.53585279398, 271957.71260195354, 95499.3326990485, 77.4453301868365, 0.41600755326637445, 0.6398381056110158, 0.864681061619831};

	double[] kuiper = {0.037152496837928616, 3.230244704573572, 1.186098986179703, 1.1852221178786189, 0.9828461137137623, 0.03785948674957811, 1.9439441499084562};

	testSampling(seed, kappas, samples, trueKappas, goodness, kuiper);
    }


    /**
     * Test that samples drawn from distribution are correctly distributed. 
     *
     */
    public void testSampling() {
	
 	double minKappa = -360.0;

 	double maxKappa = 360.0;

 	double interval = 80.0;

 	int samples = 10000;

	long seed = 56267934101l;

 	// kappas of each distribution
 	double[] kappas = new double[1 + (int)((maxKappa - minKappa) / interval)];

	for (int i = 0; i < kappas.length; i++) {
	    
	    kappas[i] = minKappa + (maxKappa - minKappa) * (double)(i) / (kappas.length - 1.0);
    
	}


	double[] trueKappas = {-359.7886734354719, -279.6606561753043, -192.7829754756786, 
			       -116.03567057619581, -39.919352002631385, 40.15060803487996, 
			       121.01786185027629, 200.94375549162766, 285.84524531451393, 
			       357.26992837164795};
	
	double[] goodness = {57243.016721678956, 579508.7113235587, 3455432.5127393296, 
			     231652.63213093835, 564901.3279569907, 0.7534131046166762, 
			     0.7915366655504337, 0.7665818481103431, 0.6030785646471142, 
			     0.6340134624519244};

	double[] kuiper = {6.512595470185472, 0.498295673526757, 0.6532586877022716, 
			   3.705005976594561, 0.8254997130257996, 2.2953349278915303, 
			   0.31167895829493286, 6.3532719124372266, 1.6364276017829034, 
			   8.90909719986106};

	testSampling(seed, kappas, samples, trueKappas, goodness, kuiper);
    }



    /**
     * Avoids code duplication. 
     *
     */
    private void testSampling(long seed, double[] kappas, int samples, double[] trueKappas, double[] goodness, double[] kuiper) {

	WatsonDistribution[] watsons = new WatsonDistribution[kappas.length];

	java.util.Random ran = new java.util.Random(seed);

	Vector3D mu = new Vector3D(ran.nextGaussian(), ran.nextGaussian(), 
				   ran.nextGaussian()).normalized();

	EigenSystem3D[] tBarEigs = new EigenSystem3D[kappas.length];
	
	Vector3D[][] sampleVecs = new Vector3D[kappas.length][samples];

	
	for (int i = 0; i < kappas.length; i++) {
		watsons[i] = new WatsonDistribution(mu, ran);
	}

	for (int i = 0; i < samples; i++) {
	    for (int n = 0; n < kappas.length; n++) {
		
		Vector3D s = watsons[n].nextVector(kappas[n]);

		sampleVecs[n][i] = s;
	    }
	}


	for (int n = 0; n < kappas.length; n++) {
	    
	    tBarEigs[n] = WatsonFitter.tBarEigenSystem(sampleVecs[n]);

	    double fittedKappa = WatsonFitter.fitKappa(tBarEigs[n], sampleVecs[n]);

	    assertEquals(trueKappas[n], fittedKappa, 0.00001);
	    
	    if (Math.abs(fittedKappa) > 5.0) { // can't expect great accuracy for low K
		if (fittedKappa > 0.0) {
		    assertEquals(1.0, Math.abs(tBarEigs[n].eigenvectors[0].dot(mu)), 0.00001);
		}
		else {
		    assertEquals(1.0, Math.abs(tBarEigs[n].eigenvectors[2].dot(mu)), 0.00001);
		}
	    }
	    	    
	}
	

	// calculate goodness of fit


	for (int n = 0; n < kappas.length; n++) {
   
	    if (trueKappas[n] > 0.0) {
		assertEquals(goodness[n], WatsonFitter.ksTest(sampleVecs[n], 
							      tBarEigs[n].eigenvectors[0], 
							      trueKappas[n]), 1E-6);
	    }
	    else {
		assertEquals(goodness[n], WatsonFitter.ksTest(sampleVecs[n], 
							      tBarEigs[n].eigenvectors[2], 
							      trueKappas[n]), 1E-6);
	    }
	    

	}
	
	
	

	// test rotational symmetry
	
	
	for (int n = 0; n < kappas.length; n++) {
	    
	    if (kappas[n] > 0.0) {

		assertEquals(kuiper[n], WatsonFitter.testBipolarRotSymm(tBarEigs[n], 
									tBarEigs[n].eigenvectors[0], 
									sampleVecs[n]), 0.000000001);
	    }
	    else {
		assertEquals(kuiper[n], WatsonFitter.testGirdleRotSymm(tBarEigs[n], 
								       tBarEigs[n].eigenvectors[2], 
								       sampleVecs[n]), 0.000000001);
	    }
	    
	    
	}
	
	
    }
    
  
 

}
		     

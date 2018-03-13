package numerics;

import junit.framework.*;
import junit.extensions.*;
import java.util.StringTokenizer;
import java.util.Random;

import tools.FileInput;


/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>WatsonDistribution.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>WatsonDistribution</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestWatsonDistribution.java,v 1.1 2005/09/26 10:36:49 ucacpco Exp $
 * @author  Philip Cook
 * @see numerics.WatsonDistribution
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestWatsonDistribution extends TestCase {


    public TestWatsonDistribution(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }
    
    public static Test suite() {
	return new TestSuite(TestWatsonDistribution.class);
    }


    /**
     * Test hyper1F1. If this works, then #GammaFunctions.gammln(double) also works.
     *
     */
    public void testHyper1F1() {

	// from the book

	assertEquals(12.15039235, WatsonDistribution.hyper1F1(2.0, 0.5, 1.0, 1.0e-10), 1.0e-5);
	assertEquals(5.436563657, WatsonDistribution.hyper1F1(2.0, 1.0, 1.0, 1.0e-10), 1.0e-5);
	assertEquals(3.545117704, WatsonDistribution.hyper1F1(2.0, 1.5, 1.0, 1.0e-10), 1.0e-5);
	
	assertEquals(277.5103024, WatsonDistribution.hyper1F1(2.0, 0.5, 3.0, 1.0e-10), 1.0e-5);
	assertEquals(80.34214769, WatsonDistribution.hyper1F1(2.0, 1.0, 3.0, 1.0e-10), 1.0e-5);
	assertEquals(35.9550392, WatsonDistribution.hyper1F1(2.0, 1.5, 3.0, 1.0e-10), 1.0e-5);

	assertEquals(3823.379408, WatsonDistribution.hyper1F1(2.0, 0.5, 5.0, 1.0e-10), 1.0e-5);
	assertEquals(890.4789546, WatsonDistribution.hyper1F1(2.0, 1.0, 5.0, 1.0e-10), 1.0e-5);
	assertEquals(323.5090268, WatsonDistribution.hyper1F1(2.0, 1.5, 5.0, 1.0e-10), 1.0e-5);

	assertEquals(3.346503339e24, WatsonDistribution.hyper1F1(2.0, 0.5, 50.0, 1.0e-10), 1.0e15);
	assertEquals(2.64419982e23, WatsonDistribution.hyper1F1(2.0, 1.0, 50.0, 1.0e-10), 1.0e15);
	assertEquals(3.281522692e22, WatsonDistribution.hyper1F1(2.0, 1.5, 50.0, 1.0e-10), 1.0e15);


	// from Mathematica
	assertEquals(0.053, WatsonDistribution.hyper1F1(1.0, 1.5, -10.0, 1.0e-10), 0.001);

	assertEquals(0.017, WatsonDistribution.hyper1F1(1.0, 1.5, -30.0, 1.0e-10), 0.001);

	assertEquals(0.010, WatsonDistribution.hyper1F1(1.0, 1.5, -50.0, 1.0e-10), 0.001);

	

	assertEquals(7.15873e37, WatsonDistribution.hyper1F1(2.0, 0.5, 80.0, 1.0e-10), 1.0e32);
	assertEquals(4.4879e36, WatsonDistribution.hyper1F1(2.0, 1.0, 80.0, 1.0e-10), 1.0e32);
	assertEquals(4.41931e35, WatsonDistribution.hyper1F1(2.0, 1.5, 80.0, 1.0e-10), 1.0e30);


    }


   /**
     * Test logHyper1F1
     *
     */
    public void testLogHyper1F1() {

	assertEquals(WatsonDistribution.hyper1F1(0.5, 1.5, 10.0, 1.0e-10), Math.exp(WatsonDistribution.logHyper1F1(0.5, 1.5, 10.0, 1.0e-10)), WatsonDistribution.hyper1F1(0.5, 1.5, 10.0, 1.0e-10) / 1000000.0 );

	assertEquals(WatsonDistribution.hyper1F1(0.5, 1.5, -10.0, 1.0e-10), Math.exp(WatsonDistribution.logHyper1F1(0.5, 1.5, -10.0, 1.0e-10)), WatsonDistribution.hyper1F1(0.5, 1.5, -10.0, 1.0e-10) / 1000000.0 );


	assertEquals(WatsonDistribution.hyper1F1(0.5, 1.5, 200.0, 1.0e-10), Math.exp(WatsonDistribution.logHyper1F1(0.5, 1.5, 200.0, 1.0e-10)), WatsonDistribution.hyper1F1(0.5, 1.5, 200.0, 1.0e-10) / 1000000.0 );

	assertEquals(WatsonDistribution.hyper1F1(0.5, 1.5, 500.0, 1.0e-10), Math.exp(WatsonDistribution.logHyper1F1(0.5, 1.5, 500.0, 1.0e-10)), WatsonDistribution.hyper1F1(0.5, 1.5, 500.0, 1.0e-10) / 1000000.0 );


	

    }




    public void testLogPDF() {

	Vector3D mu = new Vector3D(0.0, 0.0, 1.0);

	java.util.Random r = new java.util.Random();

	for (int i = -700; i < 700; i += 10) {
	    
	    Vector3D x = new Vector3D(r.nextGaussian(), r.nextGaussian(), r.nextGaussian()).normalized();

	    double kappa = (double)i;

	    double expLogPDF = Math.exp( kappa * Math.pow(mu.dot(x), 2.0) - WatsonDistribution.logHyper1F1(0.5, 1.5, kappa, 1.0e-10) );

	    double logPDF = WatsonDistribution.logPDF(mu, x, kappa);

	    assertEquals(WatsonDistribution.pdf(mu, x, kappa), expLogPDF, WatsonDistribution.pdf(mu, x, kappa) / 100);
	    assertEquals(Math.exp(WatsonDistribution.logPDF(mu, x, kappa)), WatsonDistribution.pdf(mu, x, kappa), WatsonDistribution.pdf(mu, x, kappa) / 100);

	}

	// crank it up, check it doesn't blow
	for (int i = 700; i < 2000; i += 10) {
	    
	    Vector3D x = new Vector3D(r.nextDouble(), r.nextDouble(), r.nextDouble()).normalized();

	    double kappa = (double)i;

	    double expLogPDF = Math.exp( kappa * Math.pow(mu.dot(x), 2.0) - WatsonDistribution.logHyper1F1(0.5, 1.5, kappa, 1.0e-10) );

	    assertFalse(Double.isNaN(expLogPDF));
	    assertFalse(Double.isInfinite(expLogPDF));
	    

	}


    }


    /**
     * Test that the PDF integrates correctly.
     *
     */
    public void testPDF() {
	
	FileInput in = new FileInput("./test/numerics/SHREWD_ZCW700.txt");
	
	double sumBP = 0.0;
	double sumG = 0.0;
	double sumStatic = 0.0;
	double sumUniform = 0.0;
	
	double[] weight = new double[700];
	double[] phi = new double[700];
	double[] theta = new double[700];


	for (int i = 0; i < 700; i++) {
	    StringTokenizer tokens = new StringTokenizer(in.readString(), " ");
	    phi[i] = Double.parseDouble(tokens.nextToken());
	    theta[i] = Double.parseDouble(tokens.nextToken());
	    weight[i] = Double.parseDouble(tokens.nextToken());

	}

	in.close();


	// arbitrary axis
	Vector3D mu = new Vector3D(0.4, 0.3, 0.2);

	mu = mu.normalized();

	assertEquals(1.0, mu.mod(), 0.000001);

	WatsonDistribution wd = new WatsonDistribution(mu, new java.util.Random());
	
	// making this large causes problems as the elements will not be well spaced enough to 
	// sample the distribution accurately
	double kappa = 80.0;

	for (int i = 0; i < 700; i++) {
	    
	    Vector3D v = Vector3D.vectorFromSPC(1.0, theta[i], phi[i]);

	    double bpPDF = wd.pdf(v, kappa);
	    double staticPDF = WatsonDistribution.pdf(mu, v, kappa);

	    sumStatic += staticPDF * weight[i];
	    sumBP += bpPDF * weight[i];
	    sumG += wd.pdf(v, -1.0 * kappa) * weight[i];
	    sumUniform += 1.0 * weight[i];

	    assertEquals(bpPDF, staticPDF, 0.00000001);
	  
	}


	assertEquals(1.0, sumUniform, 0.001);
	assertEquals(1.0, sumBP, 0.001);
	assertEquals(1.0, sumG, 0.001);
	assertEquals(1.0, sumStatic, 0.001);

	
    }


    public void testSumLogPDF() {

	java.util.Random ran = new java.util.Random(3197278l);
	
	Vector3D mu = new Vector3D(ran.nextGaussian(), ran.nextGaussian(), ran.nextGaussian()).normalized();

	WatsonDistribution w = new WatsonDistribution(mu, ran);

	// get some samples
	Vector3D[] samples = new Vector3D[1000];

	double l = 0.0;
	
	double kappa = 200.0;

	for (int i = 0; i < 1000; i++) {
	    samples[i] = w.nextVector(kappa);
	    l += WatsonDistribution.logPDF(mu, samples[i], kappa);
	}
	
	assertEquals(l, WatsonDistribution.sumLogPDF(mu, samples, kappa), 0.000000001);
	
    }



    public void testSampling() {
	
	Vector3D[] expected = new Vector3D[10];

	Random ran = new Random(123456);
	
	Vector3D mu = new Vector3D(ran.nextGaussian(), ran.nextGaussian(), ran.nextGaussian()).normalized();
	
	WatsonDistribution d = new WatsonDistribution(mu, ran);
	
	Vector3D[] samples = new Vector3D[10];

	for (int i = 0; i < samples.length; i++) {
	    samples[i] = d.nextVector(10.0);
	}

	expected[0] = new Vector3D(0.09558127483699945, 0.0073072147250661055, 0.9953948083617352);
	expected[1] = new Vector3D(-0.18615524926290153, 0.311794885874182, 0.9317350333193298);
	expected[2] = new Vector3D(0.1700413539289699, -0.0174797913653887, 0.985281885983816);
	expected[3] = new Vector3D(0.3123825633645738, -0.7879135167218196, 0.530668846148738);
	expected[4] = new Vector3D(-0.08495998198990422, 0.26130117803246333, -0.9615110482043991);
	expected[5] = new Vector3D(0.3360625900639824, 3.9416738016787023E-4, -0.9418395724281096);
	expected[6] = new Vector3D(-0.02111916942370834, 0.07314043785515555, -0.9970980177661616);
	expected[7] = new Vector3D(-0.4415143919340837, 0.26022252665628975, 0.8586904438362483);
	expected[8] = new Vector3D(0.043537330441317613, -0.5958430725000211, 0.8019199048606902);
	expected[9] = new Vector3D(0.022177390702294814, 0.05929407838105004, 0.9979941761406124);


	for (int i = 0; i < samples.length; i++) {
	    GenTestMethods.assertVectorsEqual(samples[i], expected[i], 0.00000001);
	}

	for (int i = 0; i < samples.length; i++) {
	    samples[i] = d.nextVector(-10.0);
	}

	expected[0] = new Vector3D(0.01771743233387237, 0.9822335256046605, 0.18682449991832883);
	expected[1] = new Vector3D(-0.6582801864813183, -0.567979060862669, -0.4940313578181872);
	expected[2] = new Vector3D(0.8956772300685346, 0.44416899290551937, -0.02182212817426349);
	expected[3] = new Vector3D(-0.9517289625411491, -0.2237970139924375, -0.21006398641465615);
	expected[4] = new Vector3D(-0.5912759231402993, -0.8052259161835031, -0.04476613252253153);
	expected[5] = new Vector3D(-0.8510845878945809, -0.47228579466398274, 0.2293494111763137);
	expected[6] = new Vector3D(-0.7580840318548429, -0.6437642723910174, 0.10428884043640296);
	expected[7] = new Vector3D(-0.5882375033448793, -0.8086779497001403, -0.004075944965374223);
	expected[8] = new Vector3D(-0.9210599249571205, 0.3067171123749345, 0.23994421771396202);
	expected[9] = new Vector3D(0.9337234020912559, 0.2972331336928285, -0.19953213430992045);

	for (int i = 0; i < samples.length; i++) {
	    GenTestMethods.assertVectorsEqual(samples[i], expected[i], 0.00000001);
	}


	// test large k
	for (int i = 0; i < samples.length; i++) {
	    samples[i] = d.nextVector(900.0);
	}

	expected[0] = new Vector3D(0.0029221120977905926, -0.08878122763022357, 0.9960468637977624);
	expected[1] = new Vector3D(0.06222346376882666, -0.12153035133652171, 0.9906354598239615);
	expected[2] = new Vector3D(-0.03748594158956315, 0.15018849547560187, -0.9879464661660162);
	expected[3] = new Vector3D(0.0204966023497901, -0.09347722677626405, 0.9954104165450217);
	expected[4] = new Vector3D(-0.0010367536324464022, 0.13633591789364802, -0.9906621233468057);
	expected[5] = new Vector3D(0.039452745907439385, -0.1099933841307397, 0.9931490000437149);
	expected[6] = new Vector3D(-0.0034356269691688396, -0.15163599597561678, 0.9884304331574436);
	expected[7] = new Vector3D(0.041758027192315805, -0.11539896970483754, 0.9924411040238449);
	expected[8] = new Vector3D(-1.414567317282872E-4, -0.10868454327410469, 0.9940762797910895);
	expected[9] = new Vector3D(-0.025120398291198432, 0.09203626665698635, -0.995438743072388);


    }


    /**
     * Test that duplicated code is consistent.
     */
    public void testSamplingVersusStatic() {

	long seed = 24921l;
	
	java.util.Random ran = new java.util.Random(seed);

	double kappa = 10.0 + 20.0 * ran.nextDouble();
	
	Vector3D mu = new Vector3D(ran.nextGaussian(), ran.nextGaussian(), ran.nextGaussian()).normalized();
    	WatsonDistribution w = new WatsonDistribution(mu, new java.util.Random(seed));
		
	java.util.Random staticRan = new java.util.Random(seed);
	
	for (int i = 0; i < 20; i++) {
	    Vector3D sample = w.nextVector(kappa);
	    double[] tp = Vector3D.thetaPhi(mu);
	    Vector3D staticSample = WatsonDistribution.nextVector(tp[0], tp[1], kappa, staticRan);

	    assertEquals(1.0, sample.dot(staticSample), 0.0000001);
	}    


	// try now for negative kappa 
	kappa = -1.0 * kappa;
	
	w = new WatsonDistribution(mu, new java.util.Random(seed));
	staticRan = new java.util.Random(seed);
	
	for (int i = 0; i < 20; i++) {
	    Vector3D sample = w.nextVector(kappa);

	    double[] tp = Vector3D.thetaPhi(mu);
	    Vector3D staticSample = WatsonDistribution.nextVector(tp[0], tp[1], kappa, staticRan);
	    
	    assertEquals(1.0, sample.dot(staticSample), 0.0000001);
	}    


	// now for large kappa
	kappa = 800.0;
	
	w = new WatsonDistribution(mu, new java.util.Random(seed));
	staticRan = new java.util.Random(seed);
	
	for (int i = 0; i < 20; i++) {
	    Vector3D sample = w.nextVector(kappa);

	    double[] tp = Vector3D.thetaPhi(mu);
	    Vector3D staticSample = WatsonDistribution.nextVector(tp[0], tp[1], kappa, staticRan);

	    assertEquals(1.0, sample.dot(staticSample), 0.0000001);
	}    


    }


    
}

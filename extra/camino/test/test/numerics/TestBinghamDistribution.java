package numerics;

import junit.framework.*;
import junit.extensions.*;
import java.util.StringTokenizer;
import java.util.Random;

import tools.FileInput;


/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>BinghamDistribution.java</code>.
 * <BR><BR>
 *
 * </dl>
 *
 * @version $Id: TestBinghamDistribution.java,v 1.2 2006/06/28 16:20:08 ucacpco Exp $
 * @author  Philip Cook
 * @see numerics.BinghamDistribution
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestBinghamDistribution extends TestCase {


    public TestBinghamDistribution(String name) {
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
	return new TestSuite(TestBinghamDistribution.class);
    }





    public void testPDF() {

	Vector3D[] axes = new Vector3D[3];
	
	axes[0] = new Vector3D(1.0, 0.0, 0.0);
	axes[1] = new Vector3D(0.0, 1.0, 0.0);
	axes[2] = new Vector3D(0.0, 0.0, 1.0);

	double[] bingPars = null;

	try {
	    bingPars = BinghamFitter.bngpar(0.02, 0.01);
	}
	catch (ConvergenceException e) {
	    fail(e.toString());
	}
	// test watson BP like values
	BinghamDistribution b = new BinghamDistribution(axes, bingPars[0], bingPars[1], 
							bingPars[2], new java.util.Random());

 	
	FileInput in = new FileInput("./test/numerics/SHREWD_ZCW700.txt");
	
	double sumB = 0.0;
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

	for (int i = 0; i < 700; i++) {
	    
	    Vector3D v = Vector3D.vectorFromSPC(1.0, theta[i], phi[i]);

	    sumStatic += BinghamDistribution.pdf(axes, bingPars, v) * weight[i];
	    sumB += b.pdf(v) * weight[i];
	    sumUniform += 1.0 * weight[i];
	  
	}


	assertEquals(1.0, sumUniform, 0.0001);
	assertEquals(1.0, sumB, 0.0001);
	assertEquals(1.0, sumStatic, 0.0001);


    }	
	

    public void testSampling() {

	Random ran = new Random(1827l);

	Vector3D e1 = new Vector3D(0.0, 0.0, 1.0);
	Vector3D e2 = new Vector3D(0.0, 1.0, 0.0);
	Vector3D e3 = new Vector3D(1.0, 0.0, 0.0);


	Vector3D[] expected = new Vector3D[] { 
	    new Vector3D(-0.03625305424623798, 0.10318035476753543, -0.9940017758776213),
	    new Vector3D(-0.08963764783922124, 0.08568006912937227, 0.9922822269111937),
	    new Vector3D(0.029136127313843182, -0.09453888894341106, -0.9950947113528927),
	    new Vector3D(-0.019140442467471895, 0.03831804442684211, 0.9990822643473619),
	    new Vector3D(-0.009918169581551834, -0.02922804573862027, -0.9995235621307047),
	    new Vector3D(0.22139306332226294, -0.2189011817581943, -0.9502985763104407),
	    new Vector3D(-0.01755652291845943, -0.08480775797127829, 0.9962426474965317),
	    new Vector3D(-0.08697589023846451, -0.02067775270315838, -0.9959957957041654),
	    new Vector3D(0.05512535419138733, 0.08390898998213352, -0.9949474743550295),
	    new Vector3D(0.12785367562788752, 0.09134608307524013, 0.9875775061914129)};

	
	BinghamDistribution bingDist = null;
	
	try {
	    bingDist = BinghamDistribution.getBinghamDistribution(new Vector3D[] {e1, e2, e3}, -100.0, -50.0, ran);
	}
	catch (ConvergenceException e) {
	    fail("Couldn't fit Bingham distribution");
	}


	for(int i = 0; i < expected.length; i++) {
	    Vector3D sample = bingDist.nextVector();
	    assertEquals(1.0, expected[i].dot(sample), 0.000001);
	}


    }


    public void testStaticSampling() {

	Random ran = new Random(1827l);

	Vector3D e1 = new Vector3D(0.0, 0.0, 1.0);
	Vector3D e2 = new Vector3D(0.0, 1.0, 0.0);
	Vector3D e3 = new Vector3D(1.0, 0.0, 0.0);


	Vector3D[] expected = new Vector3D[] { 
	    new Vector3D(-0.03625305424623798, 0.10318035476753543, -0.9940017758776213),
	    new Vector3D(-0.08963764783922124, 0.08568006912937227, 0.9922822269111937),
	    new Vector3D(0.029136127313843182, -0.09453888894341106, -0.9950947113528927),
	    new Vector3D(-0.019140442467471895, 0.03831804442684211, 0.9990822643473619),
	    new Vector3D(-0.009918169581551834, -0.02922804573862027, -0.9995235621307047),
	    new Vector3D(0.22139306332226294, -0.2189011817581943, -0.9502985763104407),
	    new Vector3D(-0.01755652291845943, -0.08480775797127829, 0.9962426474965317),
	    new Vector3D(-0.08697589023846451, -0.02067775270315838, -0.9959957957041654),
	    new Vector3D(0.05512535419138733, 0.08390898998213352, -0.9949474743550295),
	    new Vector3D(0.12785367562788752, 0.09134608307524013, 0.9875775061914129)};

	
	BinghamDistribution bingDist = null;
	
	try {
	    bingDist = BinghamDistribution.getBinghamDistribution(new Vector3D[] {e1, e2, e3}, -100.0, -50.0, ran);
	}
	catch (ConvergenceException e) {
	    fail("Couldn't fit Bingham distribution");
	}

	
	for(int i = 0; i < expected.length; i++) {
	    Vector3D sample = BinghamDistribution.nextVector(new Vector3D[] {e1, e2, e3}, -100.0, -50.0, bingDist.logNormC(), ran);
	    assertEquals(1.0, expected[i].dot(sample), 0.000001);
	}


    }


    public void testUniform() {

	Vector3D[] axes = new Vector3D[3];
	
	axes[0] = new Vector3D(1.0, 0.0, 0.0);
	axes[1] = new Vector3D(0.0, 1.0, 0.0);
	axes[2] = new Vector3D(0.0, 0.0, 1.0);


	Random ran = new Random(1827l);

        BinghamDistribution bingDist = null;
	
	try {
	    bingDist = BinghamDistribution.getBinghamDistribution(axes, 0.0, 0.0, ran);
	}
	catch (ConvergenceException e) {
	    fail("Couldn't fit Bingham distribution");
	}

        Vector3D[] samples = new Vector3D[100];
       

        for (int i = 0; i < 100; i++) {

            samples[i] = bingDist.nextVector();
        }

        EigenSystem3D eig = BinghamFitter.tBarEigenSystem(samples);
        
        try {
            double[] bingPars = BinghamFitter.bngpar(eig.eigenvalues[1], eig.eigenvalues[2]);
            assertEquals(0.0, bingPars[0], 1.0);
            assertEquals(0.0, bingPars[1], 1.0);
        }
        catch (ConvergenceException e) {
            fail("Unexpected exception: " + e);
        }


    }
    
}

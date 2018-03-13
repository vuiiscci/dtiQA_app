package numerics;

import junit.framework.*;
import junit.extensions.*;
import java.util.StringTokenizer;
import java.util.Random;

import tools.FileInput;


/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>ACG_Distribution.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>ACG_Distribution</code> 
 * with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestACG_Distribution.java,v 1.1 2005/09/26 10:36:45 ucacpco Exp $
 * @author  Philip Cook
 * @see numerics.ACGDistribution
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestACG_Distribution extends TestCase {


    public TestACG_Distribution(String name) {
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
	return new TestSuite(TestACG_Distribution.class);
    }

    
    
    
    /**
     * Test that the PDF integrates correctly.
     *
     */
    public void testPDF() {
	
	FileInput in = new FileInput("./test/numerics/SHREWD_ZCW700.txt");
	
	double sum = 0.0;
	double sumLog = 0.0;

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


	RealMatrix A = new RealMatrix(3,3);

	A.entries[0][0] = 0.795142802508015;
	A.entries[0][1] = -0.1095595767890877;
	A.entries[0][2] = -1.2628682075769877;
	A.entries[1][0] = -0.1095595767890877;
	A.entries[1][1] = 0.05602781014156483;
	A.entries[1][2] = 0.1824041650043007;
	A.entries[2][0] = -1.2628682075769877;
	A.entries[2][1] = 0.1824041650043007;
	A.entries[2][2] = 2.1488293873504145;

	ACG_Distribution gd = new ACG_Distribution(A, new Random(0));

	for (int i = 0; i < 700; i++) {
	    
	    Vector3D v = Vector3D.vectorFromSPC(1.0, theta[i], phi[i]);

	    sum += gd.pdf(v) * weight[i];
	    sumLog += Math.exp(gd.logPDF(v)) * weight[i];
	  
	}


	assertEquals(1.0, sum, 0.005);
	assertEquals(1.0, sumLog, 0.005);

	
    }


    /**
     * Test that both constructors provide the same result
     */
    public void testConstruction() {

	long seed = 4280l;
	
	RealMatrix A = new RealMatrix(3,3);

	A.entries[0][0] = 0.7879112249346113;
	A.entries[0][1] = -0.1119998446803741;
	A.entries[0][2] = -1.3044290212559813;
	A.entries[1][0] = -0.1119998446803741;
	A.entries[1][1] = 0.02383457776978681;
	A.entries[1][2] = 0.1869852745385591;
	A.entries[2][0] = -1.3044290212559813;
	A.entries[2][1] = 0.1869852745385591;
	A.entries[2][2] = 2.188254197295593;

	Jama.Matrix aj = new Jama.Matrix(3,3);
	
	aj.set(0,0, A.entries[0][0]);
	aj.set(0,1, A.entries[0][1]);
	aj.set(0,2, A.entries[0][2]);
	aj.set(1,0, A.entries[1][0]);
	aj.set(1,1, A.entries[1][1]);
	aj.set(1,2, A.entries[1][2]);
	aj.set(2,0, A.entries[2][0]);
	aj.set(2,1, A.entries[2][1]);
	aj.set(2,2, A.entries[2][2]);
	
	EigenSystem3D eigA = EigenSystem3D.sort(aj.eig());

	ACG_Distribution gd1 = new ACG_Distribution(A, new Random(seed));

	ACG_Distribution gd2 = new ACG_Distribution(eigA.eigenvectors, eigA.eigenvalues, new Random(seed));

	Vector3D[] samples1 = new Vector3D[100];
	Vector3D[] samples2 = new Vector3D[100];


	for (int i = 0; i < 100; i++) {
	    samples1[i] = gd1.nextVector();
	    samples2[i] = gd2.nextVector();
	}
	
	assertEquals(WatsonFitter.fitKappa(samples1), WatsonFitter.fitKappa(samples2), 0.00001);
	
    } 


    public void testSampling() {

	long seed = 4280l;
	
	RealMatrix A = new RealMatrix(3,3);

	A.entries[0][0] = 0.7879112249346113;
	A.entries[0][1] = -0.1119998446803741;
	A.entries[0][2] = -1.3044290212559813;
	A.entries[1][0] = -0.1119998446803741;
	A.entries[1][1] = 0.02383457776978681;
	A.entries[1][2] = 0.1869852745385591;
	A.entries[2][0] = -1.3044290212559813;
	A.entries[2][1] = 0.1869852745385591;
	A.entries[2][2] = 2.188254197295593;

	
	ACG_Distribution gd = new ACG_Distribution(A, new Random(seed));

	Vector3D[] expectedSamples = new Vector3D[10];

	expectedSamples[0] = new Vector3D (0.6426094580525922,  0.25585999687348465,  -0.7222110123926763);
	expectedSamples[1] = new Vector3D (0.5482274182919723,  -0.03112799424967646, -0.8357498105335778);
	expectedSamples[2] = new Vector3D (0.5450384246591186,  -0.06502575210480813, -0.8358856184959222);
	expectedSamples[3] = new Vector3D (0.4597321387161732,  0.11696797821641747,  -0.8803208805335796);
	expectedSamples[4] = new Vector3D (0.19401366587957258, 0.013232954129011421, -0.9809095709477959);
	expectedSamples[5] = new Vector3D (0.5270749753288874,  -0.07480593858276531, -0.8465199595607916);
	expectedSamples[6] = new Vector3D (-0.4113839263709455, 0.17784824694846468,  0.8939425407603009);
	expectedSamples[7] = new Vector3D (-0.4959829245282889, 0.07881363556869021,  0.8647481421921711);
	expectedSamples[8] = new Vector3D (-0.5102124672535415, 0.04417872998070172, 0.8589129630389485);
	expectedSamples[9] = new Vector3D (-0.5371838083146823, -0.06259726483949038, 0.8411391909304603);

	for (int i = 0; i < 10; i++) {
	    GenTestMethods.assertVectorsEqual(expectedSamples[i], gd.nextVector(), 0.0000001);
	}

    }
}

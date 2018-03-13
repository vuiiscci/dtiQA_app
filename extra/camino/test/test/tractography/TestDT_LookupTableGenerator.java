package tractography;

import data.*;
import imaging.*;
import inverters.*;
import misc.DT;
import numerics.*;
import optimizers.*;
import tools.*;

import java.util.Random;

import junit.framework.*;
import junit.extensions.*;



/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>DT_LookupTableGenerator</code>.
 * <BR><BR>
 *
 * </dl>
 *
 * @version $Id: TestDT_LookupTableGenerator.java,v 1.1 2005/09/26 10:36:52 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestDT_LookupTableGenerator extends TestCase {

    private Random ran = null;

    private DT_LookupTableGenerator generator = null;

    
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	// same sequence of random numbers for all tests
	ran = new Random(12345);
	generator = new OneTensorLUTGenerator(null, 0.0, 0.0, ran);
    }

   
    protected void tearDown() {
	ran = null;
    }


    public TestDT_LookupTableGenerator(String name) {
	super(name);

    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestDT_LookupTableGenerator.class);
    }

    /**
     * Tests <code>getTensor(double fa, double trace, Random ran)</code>. If this works
     * then {@link #getFullySpecifiedTensor(double, double, double, java.util.Random) 
     * getFullySpecifiedTensor } also works.
     *
     *
     */
    public void testGetTensorFromFA() {

	double trace = 2.1E-9;
	double fa = 0.75;

	DT tensor = generator.getTensor(fa, trace, ran);
	
	assertEquals(0.75, tensor.fa(), 1E-6);
	assertEquals(2.1E-9, tensor.trace(), 1E-10);

	double[][] seig = tensor.sortedEigenSystem();

	Vector3D e1 = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);
	Vector3D e2 = new Vector3D(seig[1][1], seig[2][1], seig[3][1]);
	Vector3D e3 = new Vector3D(seig[1][2], seig[2][2], seig[3][2]);

	assertEquals(1.0, Math.abs(e1.dot(generator.e1Vec)), 1E-6);

	// cylindrical symmetry 
	assertEquals(seig[0][2], seig[0][1], seig[0][2] / 1000.0);
    }


    /**
     * tests <code>protected DT getTensor(double x, double y, double trace, Random ran)</code>.
     *
     */
    public void testGetTensorFromEV() {

	double trace = 2.1E-9;
	double x = 5.0;
	double y = 2.0;

	DT tensor = generator.getTensor(x, y, trace, ran);

	assertEquals(2.1E-9, tensor.trace(), 1E-10);
	
	double[][] seig = tensor.sortedEigenSystem();

	Vector3D e1 = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);
	Vector3D e2 = new Vector3D(seig[1][1], seig[2][1], seig[3][1]);
	Vector3D e3 = new Vector3D(seig[1][2], seig[2][2], seig[3][2]);

	assertEquals(1.0, Math.abs(e1.dot(generator.e1Vec)), 1E-6);
	assertEquals(1.0, Math.abs(e2.dot(generator.e2Vec)), 1E-6);
	assertEquals(1.0, Math.abs(e3.dot(generator.e3Vec)), 1E-6);

	assertEquals(5.0, seig[0][0] / seig[0][2], 1E-6);
	assertEquals(2.0, seig[0][1] / seig[0][2], 1E-6);
    }


    /**
     * Tests <code>getTensor(Vector3D e1, Vector3D e2, Vector3D e3, double fa, 
     * double trace, Random ran)</code>. If this works, then 
     * {@link getTensor(numerics.Vector3D, numerics.Vector3D, numerics.Vector3D, double, double, 
     double, java.util.Random) getTensor(numerics.Vector3D, numerics.Vector3D, 
     * numerics.Vector3D, double, double, double, java.util.Random) } also works.
     *
     */
    public void testGetTensorWithCustomEvecs() {
	
	RealMatrix rot = Rotations.randomRotMat(ran);

	Vector3D customE1 = Rotations.rotateVector(generator.e1Vec, rot);
	Vector3D customE2 = Rotations.rotateVector(generator.e2Vec, rot);
	Vector3D customE3 = Rotations.rotateVector(generator.e3Vec, rot);

	double trace = 2.1E-9;
	double fa = 0.75;
	
	DT tensor = generator.getTensor(customE1, customE2, customE3, fa, trace, ran);
		
	assertEquals(2.1E-9, tensor.trace(), 1E-10);
	assertEquals(0.75, tensor.fa(), 1E-6);
	
	double[][] seig = tensor.sortedEigenSystem();

	Vector3D e1 = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);
	Vector3D e2 = new Vector3D(seig[1][1], seig[2][1], seig[3][1]);
	Vector3D e3 = new Vector3D(seig[1][2], seig[2][2], seig[3][2]);

	assertEquals(1.0, Math.abs(e1.dot(customE1)), 1E-6);
	
	// cylindrical symmetry 
	assertEquals(seig[0][2], seig[0][1], seig[0][2] / 1000.0);
		     
    }

  
    /**
     * Relies on WatsonFitter to work properly.
     * @see numerics.WatsonFitter
     */
    public void testGetNoisyTensors() {
	double trace = 2.1E-9;
	double fa = 0.8;

	int samples = 100;
	
	DW_Scheme imPars = DW_Scheme.readScheme("bmx6.scheme");

	DT noiseFree = generator.getTensor(fa, trace, ran);

	DT[][] noisy = generator.getNoisyTensors(new DT[] {noiseFree}, new double[] {1.0}, ModelIndex.LDT, 
						imPars, 10, samples, ran);

	// sanity checks
	assertEquals(1, noisy.length);
	assertEquals(samples, noisy[0].length);
	
	Vector3D[] sampleVecs = new Vector3D[samples];
	
	double[] faArr = new double[samples];

	for (int i = 0; i < samples; i++) {
	    double[][] eig = noisy[0][i].sortedEigenSystem();
	    sampleVecs[i] = new Vector3D(eig[1][0], eig[2][0], eig[3][0]);

	    faArr[i] = noisy[0][i].fa();
	}
	
	EigenSystem3D tbarEig = WatsonFitter.tBarEigenSystem(sampleVecs);
	
	double kappa = WatsonFitter.fitKappa(tbarEig, sampleVecs);
	
	assertEquals(1.0, Math.abs(tbarEig.eigenvectors[0].dot(generator.e1Vec)), 1E-4);
	assertEquals(154.45, kappa, 0.01);
	assertEquals(0.771, ArrayOps.mean(faArr), 1E-3);
	assertEquals(0.001905, ArrayOps.var(faArr, ArrayOps.mean(faArr)), 1E-6);
		     
	
	// again with snr = 20
	noisy = generator.getNoisyTensors(new DT[] {noiseFree}, new double[] {1.0}, ModelIndex.LDT, imPars, 
					  20, samples, ran);
	
	// sanity checks
	assertEquals(1, noisy.length);
	assertEquals(samples, noisy[0].length);
	
	sampleVecs = new Vector3D[samples];
	
	faArr = new double[samples];
	
	for (int i = 0; i < samples; i++) {
	    double[][] eig = noisy[0][i].sortedEigenSystem();
	    sampleVecs[i] = new Vector3D(eig[1][0], eig[2][0], eig[3][0]);
	    
	    faArr[i] = noisy[0][i].fa();
	}
		     
	kappa = WatsonFitter.fitKappa(sampleVecs);

	assertEquals(558.51, kappa, 0.01);
	assertEquals(0.7967, ArrayOps.mean(faArr), 1E-3);
	assertEquals(6.862E-4, ArrayOps.var(faArr, ArrayOps.mean(faArr)), 1E-6);


	RealMatrix rot = Rotations.getRotMat(generator.e3Vec, Math.PI / 2.0);

	// again with two tensors
	DT noiseFreeRot = noiseFree.transform(rot);

	noisy = generator.getNoisyTensors(new DT[] {noiseFree, noiseFreeRot}, new double[] {0.5, 0.5}, 
					  ModelIndex.POSPOS, imPars, 30.0, samples, ran);
	
	// sanity checks
	assertEquals(2, noisy.length);
	assertEquals(samples, noisy[0].length);
	
	sampleVecs = new Vector3D[2*samples];
	
	for (int i = 0; i < samples; i++) {
	    double[][] eig1 = noisy[0][i].sortedEigenSystem();
	    sampleVecs[2*i] = new Vector3D(eig1[1][0], eig1[2][0], eig1[3][0]);

	    double[][] eig2 = noisy[1][i].sortedEigenSystem();
	    sampleVecs[2*i+1] = new Vector3D(eig2[1][0], eig2[2][0], eig2[3][0]);

	}
	
	TwoFibreBipolarWatsonFitter fitter = new TwoFibreBipolarWatsonFitter(sampleVecs);
	
	try {
	    fitter.fitEstimatedParams(sampleVecs[0], sampleVecs[1], 5);
	    
	}
	catch (MarquardtMinimiserException e) {
	    fail(e.toString());
	}
	
	double[] kappas = fitter.getKappas();

	Vector3D noiseFreeRotE1 = Rotations.rotateVector(generator.e1Vec, 
							 generator.e3Vec, Math.PI / 2.0);

	assertEquals(262.33, kappas[0], 0.01);
	assertEquals(212.49, kappas[1], 0.01);

		     
    }


}

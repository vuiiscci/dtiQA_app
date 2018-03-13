package numerics;

import junit.framework.*;
import junit.extensions.*;
import numerics.*;

import Jama.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for fitting code. 
 * <BR><BR>
 *
 * 
 * </dl>
 *
 * @version $Id: TestSphericalDistributionFitter.java,v 1.1 2005/09/26 10:36:47 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestSphericalDistributionFitter extends TestCase {


    public TestSphericalDistributionFitter(String name) {
	super(name);
    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    

    public static Test suite() {
	return new TestSuite(TestSphericalDistributionFitter.class);
    }


    public void testTBarEigenSystem() {

	Vector3D[] samples = new Vector3D[10];

	/* from matlab

	samples =
	
	0.5110    0.7124    0.4810
	0.8118    0.4032    0.4224
	0.6410    0.3635    0.6760
	0.6204    0.7822    0.0579
	0.4500    0.5340    0.7158
	0.7988    0.1372    0.5857
	0.9144    0.0112    0.4047
	0.4518    0.7728    0.4457
	0.9147    0.2070    0.3471
	0.3124    0.5395    0.7819

	*/
	
	samples[0] = new Vector3D(0.5110, 0.7124, 0.4810);
	samples[1] = new Vector3D(0.8118, 0.4032, 0.4224);
	samples[2] = new Vector3D(0.6410, 0.3635, 0.6760);
	samples[3] = new Vector3D(0.6204, 0.7822, 0.0579);
	samples[4] = new Vector3D(0.4500, 0.5340, 0.7158);
	samples[5] = new Vector3D(0.7988, 0.1372, 0.5857);
	samples[6] = new Vector3D(0.9144, 0.0112, 0.4047);
	samples[7] = new Vector3D(0.4518, 0.7728, 0.4457);
	samples[8] = new Vector3D(0.9147, 0.2070, 0.3471);
	samples[9] = new Vector3D(0.3124, 0.5395, 0.7819);
	
	EigenSystem3D eig = SphericalDistributionFitter.tBarEigenSystem(samples);

	
	/*
	  [eVecs, eVals] = eig(transpose(samples) * samples)

	  eVecs =
	  
	  0.39357625506323   0.60421229409026   0.69283853466777
	  0.37096839135509  -0.79396364102022   0.48166813196777
	  -0.84111841263036  -0.06744805711134   0.53662889926474

	  
	  eVals / 10

	  ans =
	  
	  0.0495         0         0
	  0    0.0944         0
	  0         0    0.8562
  


	*/
	
	// long format

	assertEquals(1.0, eig.eigenvectors[0].dot(new Vector3D(0.69283853466777, 0.48166813196777, 0.53662889926474)), 0.01);
	assertEquals(1.0, eig.eigenvectors[1].dot(new Vector3D(0.60421229409026, -0.79396364102022, -0.06744805711134)), 0.01);
	assertEquals(1.0, eig.eigenvectors[2].dot(new Vector3D(0.39357625506323, 0.37096839135509, -0.84111841263036)), 0.01);

	assertEquals(0.8562, eig.eigenvalues[0], 0.0001);
	assertEquals(0.0944, eig.eigenvalues[1], 0.0001);
	assertEquals(0.0495, eig.eigenvalues[2], 0.0001);


    }




 
    


}
		     

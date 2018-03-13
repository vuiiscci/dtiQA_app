package misc;

import junit.framework.*;
import junit.extensions.*;

import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>DT.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestDT.java,v 1.6 2006/06/30 14:16:33 ucacpco Exp $
 * @author  Philip Cook
 * @see misc.DT
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestDT extends TestCase {


    public TestDT(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	// any class variables you declare should be initialized here. This is called before each test
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestDT.class);
    }


    public void testTransform() {

	DT dt = new DT(1000.0, 0.0, 0.0, 200.0, 0.0, 100.0);
	
	// quick sanity check
	DT dtRot = dt.transform(RealMatrix.identity(3));

	assertEquals(dt.dxx, dtRot.dxx, 0.000001);
	assertEquals(dt.dxy, dtRot.dxy, 0.000001);
	assertEquals(dt.dxz, dtRot.dxz, 0.000001);
	assertEquals(dt.dyy, dtRot.dyy, 0.000001);
	assertEquals(dt.dyz, dtRot.dyz, 0.000001);
	assertEquals(dt.dzz, dtRot.dzz, 0.000001);

	RealMatrix rot = Rotations.getRotMat(0.0, 0.0, Math.PI / 3.0);

	dtRot = dt.transform(rot);

	// shouldn't change tensor trace or FA 
	assertEquals(dt.trace(), dtRot.trace(), 0.000001);
	assertEquals(dt.fa(), dtRot.fa(), 0.000001);

    }



    public void testSortedEigenSystem() {
	
	for (int order = 0; order < 6; order++) {
	    
	    DT d = null;

	    switch(order) {
	    case 0: d = new DT(1600.0, 0.0, 0.0, 100.0, 0.0, 90.0); break;
	    case 1: d = new DT(1600.0, 0.0, 0.0, 90.0, 0.0, 100.0); break;
	    case 2: d = new DT(100.0, 0.0, 0.0, 1600.0, 0.0, 90.0); break;
	    case 3: d = new DT(100.0, 0.0, 0.0, 90.0, 0.0, 1600.0); break;
	    case 4: d = new DT(90.0, 0.0, 0.0, 100.0, 0.0, 1600.0); break;
	    case 5: d = new DT(90.0, 0.0, 0.0, 1600.0, 0.0, 100.0); break;
	    }


	    double[][] eig = d.eigenSystem();

	    double[][] seig = d.sortedEigenSystem();

	    assertEquals(1600.0, seig[0][0], 0.000001);
	    assertEquals(100.0, seig[0][1], 0.000001);
	    assertEquals(90.0, seig[0][2], 0.000001);


	}

    }
 


    public void testPPD() {
	DT dt1 = new DT(1.7, 0.0, 0.0, 0.3, 0.0, 0.1);
	DT dt2 = new DT(0.3, 0.0, 0.0, 1.7, 0.0, 0.1);

        RealMatrix jac1 = new RealMatrix(3,3);
        jac1.entries[0][0] = 1.0;
        jac1.entries[1][0] = Math.cos(Math.PI/6);
        jac1.entries[1][1] = 1.0;
        jac1.entries[2][2] = 1.0;
	
        RealMatrix jac2 = new RealMatrix(3,3);
        jac2.entries[0][0] = 1.0;
        jac2.entries[2][0] = Math.cos(Math.PI/6);
        jac2.entries[1][1] = 1.0;
        jac2.entries[2][2] = 1.0;
	
        DT dtReor11 = dt1.ppd(jac1);
        DT dtReor12 = dt1.ppd(jac2);
        DT dtReor21 = dt2.ppd(jac1);
        DT dtReor22 = dt2.ppd(jac2);

        double[] comps = dtReor11.getComponents();
        assertEquals(comps[0], 1.1, 0.00001);
        assertEquals(comps[1], 0.69282, 0.00001);
        assertEquals(comps[2], 0.0, 0.00001);
        assertEquals(comps[3], 0.9, 0.00001);
        assertEquals(comps[4], 0.0, 0.00001);
        assertEquals(comps[5], 0.1, 0.00001);

        comps = dtReor12.getComponents();
        assertEquals(comps[0], 1.01429, 0.00001);
        assertEquals(comps[1], 0.0, 0.00001);
        assertEquals(comps[2], 0.79179, 0.00001);
        assertEquals(comps[3], 0.3, 0.00001);
        assertEquals(comps[4], 0.0, 0.00001);
        assertEquals(comps[5], 0.78571, 0.00001);

        comps = dtReor21.getComponents();
        assertEquals(comps[0], 0.3, 0.00001);
        assertEquals(comps[1], 0.0, 0.00001);
        assertEquals(comps[2], 0.0, 0.00001);
        assertEquals(comps[3], 1.7, 0.00001);
        assertEquals(comps[4], 0.0, 0.00001);
        assertEquals(comps[5], 0.1, 0.00001);

        comps = dtReor22.getComponents();
        assertEquals(comps[0], 0.21429, 0.00001);
        assertEquals(comps[1], 0.0, 0.00001);
        assertEquals(comps[2], 0.09897, 0.00001);
        assertEquals(comps[3], 1.7, 0.00001);
        assertEquals(comps[4], 0.0, 0.00001);
        assertEquals(comps[5], 0.18571, 0.00001);

    }


   public void testDegenerateSecondThirdEigenSystemSorting() {

	for (int order = 0; order < 3; order++) {
	    
	    DT d = null;

	    switch(order) {
	    case 0: d = new DT(1600.0, 0.0, 0.0, 100.0, 0.0, 100.0); break;
	    case 1: d = new DT(100.0, 0.0, 0.0, 1600.0, 0.0, 100.0); break;
	    case 2: d = new DT(100.0, 0.0, 0.0, 100.0, 0.0, 1600.0); break;
	    }
	    
	    double[][] eig = d.eigenSystem();


	    double[][] seig = d.sortedEigenSystem();

	    assertEquals(1600.0, seig[0][0], 0.000001);
	    assertEquals(100.0, seig[0][1], 0.000001);
	    assertEquals(100.0, seig[0][2], 0.000001);

	}

   }


   public void testDegenerateFirstSecondEigenSystemSorting() {

	for (int order = 0; order < 3; order++) {
	    
	    DT d = null;

	    switch(order) {
	    case 0: d = new DT(1600.0, 0.0, 0.0, 1600.0, 0.0, 100.0); break;
	    case 1: d = new DT(100.0, 0.0, 0.0, 1600.0, 0.0, 1600.0); break;
	    case 2: d = new DT(100.0, 0.0, 0.0, 1600.0, 0.0, 1600.0); break;
	    }
	    double[][] eig = d.eigenSystem();

	    double[] lams = new double[3];

	    lams[0] = eig[0][0];
	    lams[1] = eig[0][1];
	    lams[2] = eig[0][2];

	    double[][] seig = d.sortedEigenSystem();

	    assertEquals(1600.0, seig[0][0], 0.000001);
	    assertEquals(1600.0, seig[0][1], 0.000001);
	    assertEquals(100.0, seig[0][2], 0.000001);

	}

   }



   

    public void testFA() {
	DT d01 = new DT(781.0, 0.0, 0.0, 660.0, 0.0, 660.0);
	DT d09 = new DT(1773.0, 0.0, 0.0, 164.0, 0.0, 164.0);

	d01 = d01.scale(1E-12);
	d09 = d09.scale(1E-12);

	/* Matlab code
		 
	>> d = [1143, 0.0, 0.0; 0.0, 479, 0.0; 0.0, 0.0, 479]

	d =

	1143           0           0
	0         479           0
	0           0         479

	>> r = [cos(pi/6), -sin(pi/6), 0.0; sin(pi/6), cos(pi/6), 0.0; 0.0, 0.0, 1.0]

	r =

	0.86602540378444  -0.50000000000000                  0
	0.50000000000000   0.86602540378444                  0
	0                  0   1.00000000000000

	>> transpose(r)*d*r

	ans =

	1.0e+02 *

	9.77000000000000  -2.87520434056434                  0
	-2.87520434056434   6.45000000000000                  0
	0                  0   4.79000000000000


		   
	*/		


	DT d05 = new DT(1143.0E-12, 0.0, 0.0, 479.0E-12, 0.0, 479.0E-12);
		
	DT d05Rotated = new DT(977.0E-12, -287.5E-12, 0.0, 645.0E-12, 0.0, 479.0E-12);


	assertEquals(0.1, d01.fa(), 1E-3);
	assertEquals(0.5, d05.fa(), 1E-3);
	assertEquals(0.5, d05Rotated.fa(), 1E-3);
	assertEquals(0.9, d09.fa(), 1E-3);


        // test isotropic

	DT d00 = new DT(700.0E-12, 0.0, 0.0, 700.0E-12, 0.0, 700.0E-12);
        
	assertEquals(0.0, d00.fa(), 1E-3);
    }


    /**
     * If this passes, then determinant is also OK.
     */
    public void testInverse() {
	DT d = new DT(1000.0, 100.0, 0.0, 200.0, 0.0, 50.0);

	/*
	  d =
	  
	  1000        100        0
	  100         200        0
	  0           0          50
	  
	  >> inv(d)
	  
	  ans =
	  
	  0.00105263157895   -0.00052631578947   0
	  -0.00052631578947  0.00526315789474    0
	  0                  0                   0.02000000000000
	  
	*/
	DT inv = d.inverse();

	assertEquals(0.00105263157895, inv.dxx, 0.0000001);
	assertEquals(-0.00052631578947, inv.dxy, 0.0000001);
	assertEquals(0.0, inv.dxz, 0.0000001);
	assertEquals(0.00526315789474, inv.dyy, 0.0000001);
	assertEquals(0.0, inv.dyz, 0.0000001);
	assertEquals(0.02, inv.dzz, 0.0000001);

	

    }

    public void testContract() {
	
	/* Matlab
	   
	D =

	1.0e+02 *
	
	9.77000000000000  -2.87520400000000   0.10000000000000
	-2.87520400000000   6.45000000000000   0.05000000000000
	0.10000000000000   0.05000000000000   4.79000000000000
	

	>> [0.5, 0.2, 0.2]*D*transpose([0.5, 0.2, 0.2])
	
	ans =
	
	2.341059200000000e+02
	
	>> 
	   
	*/

	DT d = new DT(977.0, -287.5204, 10.0000, 645.0000 , 5.0000, 479.0000);
	
	assertEquals(2.3410592E02, d.contractBy(new double[] {0.5, 0.2, 0.2}), 1E-6);

    }

    public void testITransform() {

	DT dt = new DT(1000.0, 0.0, 0.0, 200.0, 0.0, 100.0);
	
	RealMatrix rot = Rotations.getRotMat(0.0, 0.0, Math.PI / 3.0);

	DT dtRot = dt.transform(rot);

	RealMatrix iRot = Rotations.getRotMat(0.0, 0.0, -Math.PI / 3.0);

	// iTransform by iRot should give same result
	GenTestMethods.assertArraysEqual(dt.iTransform(iRot).getComponents(), dtRot.getComponents(), 1E-6);
	

    }

    public void testDT_FromEig() {
      DT d05Rotated = new DT(977.0E-3, -287.5E-3, 0.0, 645.0E-3, 0.0, 479.0E-3);

      double[][] eig = d05Rotated.eigenSystem();
     
      DT recon = DT.dtFromEig(eig);

      double[] comps = d05Rotated.getComponents();
      double[] recomps = recon.getComponents();
 
      for (int i = 0; i < 6; i++) {
	assertEquals(comps[i], recomps[i], 1E-8);
      }

    }

    
}

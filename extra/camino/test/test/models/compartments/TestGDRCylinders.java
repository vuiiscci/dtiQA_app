package models.compartments;

import junit.framework.*;
import junit.extensions.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import imaging.*;
import numerics.*;
import models.compartments.*;


/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>GDRCylinders.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>GDRCylinders</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/GDRCylinders.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestGDRCylinders extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestGDRCylinders(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
        ip = DW_Scheme.readScheme("ActiveAxG140_PM.scheme1");
    }
    
    public static Test suite() {
	return new TestSuite(TestGDRCylinders.class);
    }


    public void testGetSignalScheme() {

	GDRCylinders gdr_cyl = new GDRCylinders();

	// k, b, diff, theta, phi
	double[] params = {1.8,2E-6,6E-10,1.5,-1.5};

	RealMatrix signals = gdr_cyl.getSignals(params, ip);
	RealMatrix jac = gdr_cyl.getJacobian(params, ip);
	
//	int[] indices = {1,3,97,191,285};
//	for(int i=0; i<indices.length; i++)
//	{
//		System.out.println("assertEquals("+signals.entries[indices[i]][0]+", signals.entries["+indices[i]+"][0], 0.0000001);");
//	}
//	System.out.println("");
//	System.out.println("");
//	for(int j=0; j<params.length; j++)
//	{
//		for(int i=0; i<indices.length; i++)
//		{
//			
//				System.out.println("assertEquals("+jac.entries[indices[i]][j]+", jac.entries["+indices[i]+"]["+j+"], 0.0000001);");
//		}
//		System.out.println("");
//	}
	
	assertEquals(1.0, signals.entries[1][0], 0.0000001);
	assertEquals(0.7684692993774136, signals.entries[3][0], 0.0000001);
	assertEquals(0.8172519900088466, signals.entries[97][0], 0.0000001);
	assertEquals(0.8445129800664916, signals.entries[191][0], 0.0000001);
	assertEquals(5.114815173991515E-4, signals.entries[285][0], 0.0000001);


	assertEquals(0.0, jac.entries[1][0], 0.0000001);
	assertEquals(-0.16196724772453308, jac.entries[3][0], 0.0000001);
	assertEquals(-0.18403860926628113, jac.entries[97][0], 0.0000001);
	assertEquals(-0.15170319378376007, jac.entries[191][0], 0.0000001);
	assertEquals(-2.0811643480556086E-5, jac.entries[285][0], 0.0000001);

	assertEquals(0.0, jac.entries[1][1], 0.0000001);
	assertEquals(-128455.734375, jac.entries[3][1], 0.0000001);
	assertEquals(-145775.65625, jac.entries[97][1], 0.0000001);
	assertEquals(-125048.6640625, jac.entries[191][1], 0.0000001);
	assertEquals(-17.560148239135742, jac.entries[285][1], 0.0000001);

	assertEquals(0.0, jac.entries[1][2], 0.0000001);
	assertEquals(-1.1231728E8, jac.entries[3][2], 0.0000001);
	assertEquals(-1.865934E7, jac.entries[97][2], 0.0000001);
	assertEquals(-2.1362158E7, jac.entries[191][2], 0.0000001);
	assertEquals(-6430687.5, jac.entries[285][2], 0.0000001);

	assertEquals(0.0, jac.entries[1][3], 0.0000001);
	assertEquals(0.18440648913383484, jac.entries[3][3], 0.0000001);
	assertEquals(0.22742344439029694, jac.entries[97][3], 0.0000001);
	assertEquals(-0.45611709356307983, jac.entries[191][3], 0.0000001);
	assertEquals(-9.111716644838452E-4, jac.entries[285][3], 0.0000001);

	assertEquals(0.0, jac.entries[1][4], 0.0000001);
	assertEquals(-0.3775908648967743, jac.entries[3][4], 0.0000001);
	assertEquals(-0.05939415842294693, jac.entries[97][4], 0.0000001);
	assertEquals(-0.059447240084409714, jac.entries[191][4], 0.0000001);
	assertEquals(0.0012815413065254688, jac.entries[285][4], 0.0000001);

    }	

}



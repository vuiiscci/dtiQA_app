package models.compartments;

import junit.framework.*;
import junit.extensions.*;
import java.util.Random;
import imaging.*;
import numerics.*;
import models.compartments.*;


/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>Tensor.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Tensor</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/Tensor.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTensor extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestTensor(String name) {
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
	return new TestSuite(TestTensor.class);
    }


    public void testGetSignalScheme() {

	Tensor tensor = new Tensor();

	// diffusivity, theta, phi, diffperp1, diffperp2, alpha
	double[] params = {1.7E-9, 1,1, 0.6E-9, 0.3E-9 , 1};

	RealMatrix signals = tensor.getSignals(params, ip);
	RealMatrix jac = tensor.getJacobian(params, ip);
//
//	System.err.println("signals");
//	for(int i=0; i<signals.entries.length; i++)
//	    System.err.print(i + "," + signals.entries[i][0] + " ");
//	System.err.println();
//	System.err.println("jac");
//	for(int i=0; i<jac.entries.length; i++) {
//	    System.err.print(i + ":");
//	    for(int j=0; j<jac.entries[0].length; j++)
//		System.err.print(" " + jac.entries[i][j]);
//	    System.err.println();
//	}
	

        assertEquals(1.0, signals.entries[1][0], 0.0000001);
        assertEquals(0.23576346649633492, signals.entries[3][0], 0.0000001);
        assertEquals(0.21035501939737827, signals.entries[97][0], 0.0000001);
        assertEquals(0.05573023545626891, signals.entries[191][0], 0.0000001);
        assertEquals(2.32943017215995E-7, signals.entries[285][0], 0.0000001);
         
              
             
        assertEquals(0.0, jac.entries[1][0], 0.0000001);
        assertEquals(-1.109166E8, jac.entries[3][0], 0.0000001);
        assertEquals(-1.43499264E8, jac.entries[97][0], 0.0000001);
        assertEquals( -6.428806E7, jac.entries[191][0], 0.0000001);
        assertEquals(-1869.43017578125, jac.entries[285][0], 0.0000001);

        assertEquals(0.0, jac.entries[1][1], 0.0000001);
        assertEquals(0.17173750698566437, jac.entries[3][1], 0.0000001);
        assertEquals(0.501343846321106, jac.entries[97][1], 0.0000001);
        assertEquals(0.20181480050086975, jac.entries[191][1], 0.0000001);
        assertEquals(-3.0380965654330794E-6, jac.entries[285][1], 0.0000001);

        assertEquals(0.0, jac.entries[1][2], 0.0000001);
        assertEquals(0.4450336694717407, jac.entries[3][2], 0.0000001);
        assertEquals(0.14114096760749817, jac.entries[97][2], 0.0000001);
        assertEquals(-0.04922059550881386, jac.entries[191][2], 0.0000001);
        assertEquals(-2.303100473000086E-6, jac.entries[285][2], 0.0000001);
        
        assertEquals(0.0, jac.entries[1][3], 0.0000001);
        assertEquals(-1.63975264E8, jac.entries[3][3], 0.0000001);
        assertEquals(-1.7030586E7, jac.entries[97][3], 0.0000001);
        assertEquals(-6.3919192E7, jac.entries[191][3], 0.0000001);
        assertEquals(-61.78946304321289, jac.entries[285][3], 0.0000001);
        
        
        assertEquals(0.0, jac.entries[1][4], 0.0000001);
        assertEquals(-1.79059056E8 , jac.entries[3][4], 0.0000001);
        assertEquals(-2.45890704E8, jac.entries[97][4], 0.0000001);
        assertEquals(-4.42152E7, jac.entries[191][4], 0.0000001);
        assertEquals(-1141.66943359375, jac.entries[285][4], 0.0000001);
        
        assertEquals(0.0, jac.entries[1][5], 0.0000001);
        assertEquals(-0.10281112045049667, jac.entries[3][5], 0.0000001);
        assertEquals(0.0388275571167469, jac.entries[97][5], 0.0000001);
        assertEquals(0.03189735859632492, jac.entries[191][5], 0.0000001);
        assertEquals(-1.5935923158849619E-7, jac.entries[285][5], 0.0000001);

    }	

}



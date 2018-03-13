package models.compartments;

import junit.framework.*;
import junit.extensions.*;
import java.util.Random;
import imaging.*;
import numerics.*;
import models.compartments.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>Stick.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Stick</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/Stick.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestStick extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestStick(String name) {
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
	return new TestSuite(TestStick.class);
    }


    public void testGetSignalScheme() {

	Stick stick = new Stick();

	// diffusivity, theta, phi.
	double[] params = {0.6E-9, 1, 1};

	RealMatrix signals = stick.getSignals(params, ip);
	RealMatrix jac = stick.getJacobian(params, ip);

	/*
	System.err.println("signals");
	for(int i=0; i<signals.entries.length; i++)
	    System.err.print(i + "," + signals.entries[i][0] + " ");
	System.err.println();
	System.err.println("jac");
	for(int i=0; i<jac.entries.length; i++) {
	    System.err.print(i + ":");
	    for(int j=0; j<jac.entries[0].length; j++)
		System.err.print(" " + jac.entries[i][j]);
	    System.err.println();
	}
	*/

        assertEquals(1.0, signals.entries[1][0], 0.0000001);
        assertEquals(0.7540671766779545, signals.entries[3][0], 0.0000001);
        assertEquals(0.6641104457977804, signals.entries[97][0], 0.0000001);
        assertEquals(0.5005054179624501, signals.entries[191][0], 0.0000001);
        assertEquals(0.008105800699469212, signals.entries[285][0], 0.0000001);

        assertEquals(0.0, jac.entries[1][0], 0.0000001);
        assertEquals(-3.5475570567748845E8, jac.entries[3][0], 0.0000001);
        assertEquals(-4.5304154599254364E8, jac.entries[97][0], 0.0000001);
        assertEquals(-5.773637433198513E8, jac.entries[191][0], 0.0000001);
        assertEquals(-6.505141936961775E7, jac.entries[285][0], 0.0000001);

        assertEquals(0.0, jac.entries[1][1], 0.0000001);
        assertEquals(0.1754802645525041, jac.entries[3][1], 0.0000001);
        assertEquals(0.700022375020228, jac.entries[97][1], 0.0000001);
        assertEquals(0.8567589291988403, jac.entries[191][1], 0.0000001);
        assertEquals(-0.043664577383441906, jac.entries[285][1], 0.0000001);

        assertEquals(0.0, jac.entries[1][2], 0.0000001);
        assertEquals(0.6124213117528668, jac.entries[3][2], 0.0000001);
        assertEquals(0.19093722473182717, jac.entries[97][2], 0.0000001);
        assertEquals(-0.22793591102463992, jac.entries[191][2], 0.0000001);
        assertEquals(-0.03778395669529067, jac.entries[285][2], 0.0000001);

    }	

}




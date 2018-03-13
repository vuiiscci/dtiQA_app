package models.compartments;

import junit.framework.*;
import junit.extensions.*;
import java.util.Random;
import imaging.*;
import numerics.*;
import models.compartments.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>Zeppelin.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Zeppelin</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/Zeppelin.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestZeppelin extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestZeppelin(String name) {
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
	return new TestSuite(TestZeppelin.class);
    }


    public void testGetSignalScheme() {

	Zeppelin zep = new Zeppelin();

	// diffusivity parallel, theta, phi, diffusivity perpendicular.
	double[] params = {0.6E-9, 1, 1, 0.2E-9};

	RealMatrix signals = zep.getSignals(params, ip);
	RealMatrix jac = zep.getJacobian(params, ip);

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
        assertEquals(0.5636794235774016, signals.entries[3][0], 0.0000001);
        assertEquals(0.5172216462024697, signals.entries[97][0], 0.0000001);
        assertEquals(0.33952712515376043, signals.entries[191][0], 0.0000001);
        assertEquals(0.002884442732501836, signals.entries[285][0], 0.0000001);

        assertEquals(0.0, jac.entries[1][0], 0.0000001);
        assertEquals(-2.6518657471346638E8, jac.entries[3][0], 0.0000001);
        assertEquals(-3.5283723618424445E8, jac.entries[97][0], 0.0000001);
        assertEquals(-3.91665394423582E8, jac.entries[191][0], 0.0000001);
        assertEquals(-2.3148495848400276E7, jac.entries[285][0], 0.0000001);

        assertEquals(0.0, jac.entries[1][1], 0.0000001);
        assertEquals(0.08744986957080238, jac.entries[3][1], 0.0000001);
        assertEquals(0.3634603525939952, jac.entries[97][1], 0.0000001);
        assertEquals(0.3874655308293724, jac.entries[191][1], 0.0000001);
        assertEquals(-0.010358670090805499, jac.entries[285][1], 0.0000001);

        assertEquals(0.0, jac.entries[1][2], 0.0000001);
        assertEquals(0.3051976469931967, jac.entries[3][2], 0.0000001);
        assertEquals(0.09913698976028225, jac.entries[97][2], 0.0000001);
        assertEquals(-0.10308303275324444, jac.entries[191][2], 0.0000001);
        assertEquals(-0.008963593960724277, jac.entries[285][2], 0.0000001);

        assertEquals(0.0, jac.entries[1][3], 0.0000001);
        assertEquals(-8.201416276208899E8, jac.entries[3][3], 0.0000001);
        assertEquals(-6.464675022210473E8, jac.entries[97][3], 0.0000001);
        assertEquals(-6.587922605704966E8, jac.entries[191][3], 0.0000001);
        assertEquals(-1.4901726597001161E7, jac.entries[285][3], 0.0000001);

    }	

}




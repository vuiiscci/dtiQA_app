package models.compartments;

import junit.framework.*;
import junit.extensions.*;
import java.util.Random;
import imaging.*;
import numerics.*;
import models.compartments.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>Ball.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Ball</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/Ball.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestBall extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestBall(String name) {
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
	return new TestSuite(TestBall.class);
    }


    public void testGetSignalScheme() {

	Ball ball = new Ball();

	// diffusivity.
	double[] params = {0.6E-9};

	RealMatrix signals = ball.getSignals(params, ip);
	RealMatrix jac = ball.getJacobian(params, ip);

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
        assertEquals(0.314975280267681, signals.entries[3][0], 0.0000001);
        assertEquals(0.31372523552404236, signals.entries[97][0], 0.0000001);
        assertEquals(0.15624490434750724, signals.entries[191][0], 0.0000001);
        assertEquals(3.652531079616311E-4, signals.entries[285][0], 0.0000001);

        assertEquals(0.0, jac.entries[1][0], 0.0000001);
        assertEquals(-6.064644910100057E8, jac.entries[3][0], 0.0000001);
        assertEquals(-6.06136879069766E8, jac.entries[97][0], 0.0000001);
        assertEquals(-4.834036625242481E8, jac.entries[191][0], 0.0000001);
        assertEquals(-4818248.54770468, jac.entries[285][0], 0.0000001);
    }	

}




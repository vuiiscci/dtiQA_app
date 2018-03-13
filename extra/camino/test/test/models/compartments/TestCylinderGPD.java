package models.compartments;

import junit.framework.*;
import junit.extensions.*;
import java.util.Random;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>CylinderGPD.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>CylinderGPD</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/CylinderGPD.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestCylinderGPD extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestCylinderGPD(String name) {
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
	return new TestSuite(TestCylinderGPD.class);
    }


    public void testGetSignalScheme() {

	CylinderGPD cyl = new CylinderGPD();

	// diffusivity, theta, phi, radius.
	double[] params = {0.6E-9, 0, 0, 1E-6};

	RealMatrix signals = cyl.getSignals(params, ip);

	double[] G = {0.122725454400000, -0.027900169200000, 0.061411741200000};
	double sig = cyl.getSignal(params, 0.01015, 0.0167, G);
	//System.err.println("sig[3] direct: " + sig);
        assertEquals(0.798653708946276, sig, 0.0000001);

	//System.err.println("signals");
	//for(int i=0; i<signals.entries.length; i++)
	//    System.err.print(i + "," + signals.entries[i][0] + " ");
	//System.err.println();

        assertEquals(1.0, signals.entries[1][0], 0.0000001);
        assertEquals(0.798653708946276, signals.entries[3][0], 0.0000001);
        assertEquals(0.3397753813046995, signals.entries[97][0], 0.0000001);
        assertEquals(0.17730539148848956, signals.entries[191][0], 0.0000001);
        assertEquals(0.9712186281584984, signals.entries[285][0], 0.0000001);
    }	

    public void testGetSignalDirect() {

	CylinderGPD cyl = new CylinderGPD();

	// diffusivity, theta, phi, radius.
	double[] params = {0.6E-9, 0, 0, 1E-6};

	double[] G1 = {0, 0.05, 0};
	double[] G2 = {0, 0.5, 0};
	double[] G3 = {0, 0, 0.05};

	double t11 = cyl.getSignal(params, 0.005, 0.04, G1);
	double t21 = cyl.getSignal(params, 0.005, 0.4, G1);
	double t31 = cyl.getSignal(params, 0.015, 0.04, G1);
	double t41 = cyl.getSignal(params, 0.015, 0.4, G1);

	//System.err.println("t11: " + t11 + " t21: " + t21 + " t31: " + t31 + " t41: " + t41);

	assertEquals(0.9998039457158182, t11, 0.0000001);
	assertEquals(0.9998039457158182, t21, 0.0000001);
	assertEquals(0.9993692719918621, t31, 0.0000001);
	assertEquals(0.9993692719918621, t41, 0.0000001);

	double t12 = cyl.getSignal(params, 0.005, 0.04, G2);
	double t22 = cyl.getSignal(params, 0.005, 0.4, G2);
	double t32 = cyl.getSignal(params, 0.015, 0.04, G2);
	double t42 = cyl.getSignal(params, 0.015, 0.4, G2);

	//System.err.println("t12: " + t12 + " t22: " + t22 + " t32: " + t32 + " t42: " + t42);

	assertEquals(0.9805836233631441, t12, 0.0000001);
	assertEquals(0.9805836233631441, t22, 0.0000001);
	assertEquals(0.938856437595154, t32, 0.0000001);
	assertEquals(0.938856437595154, t42, 0.0000001);

	double t13 = cyl.getSignal(params, 0.005, 0.04, G3);
	double t23 = cyl.getSignal(params, 0.005, 0.4, G3);
	double t33 = cyl.getSignal(params, 0.015, 0.04, G3);
	double t43 = cyl.getSignal(params, 0.015, 0.4, G3);

	//System.err.println("t13: " + t13 + " t23: " + t23 + " t33: " + t33 + " t43: " + t43);

	assertEquals(0.9022407948835627, t13, 0.0000001);
	assertEquals(0.34335397764297787, t23, 0.0000001);
	assertEquals(0.4294050722878022, t33, 0.0000001);
	assertEquals(7.188591047255512E-5, t43, 0.0000001);


	// Radius 0.1um
	params[3] = 1E-7;
        t11 = cyl.getSignal(params, 0.005, 0.04, G1);
        t21 = cyl.getSignal(params, 0.005, 0.4, G1);
        t31 = cyl.getSignal(params, 0.015, 0.04, G1);
        t41 = cyl.getSignal(params, 0.015, 0.4, G1);

        //System.err.println("t11: " + t11 + " t21: " + t21 + " t31: " + t31 + " t41: " + t41);

        assertEquals(0.9999999782786391, t11, 0.0000001);
        assertEquals(0.9999999782786391, t21, 0.0000001);
        assertEquals(0.9999999347932096, t31, 0.0000001);
        assertEquals(0.9999999347932096, t41, 0.0000001);

	t12 = cyl.getSignal(params, 0.005, 0.04, G2);
	t22 = cyl.getSignal(params, 0.005, 0.4, G2);
        t32 = cyl.getSignal(params, 0.015, 0.04, G2);
        t42 = cyl.getSignal(params, 0.015, 0.4, G2);

        //System.err.println("t12: " + t12 + " t22: " + t22 + " t32: " + t32 + " t42: " + t42);

        assertEquals(0.9999978278662395, t12, 0.0000001);
        assertEquals(0.9999978278662395, t22, 0.0000001);
        assertEquals(0.999993479342007, t32, 0.0000001);
        assertEquals(0.999993479342007, t42, 0.0000001);

        t13 = cyl.getSignal(params, 0.005, 0.04, G3);
	t23 = cyl.getSignal(params, 0.005, 0.4, G3);
        t33 = cyl.getSignal(params, 0.015, 0.04, G3);
        t43 = cyl.getSignal(params, 0.015, 0.4, G3);

        //System.err.println("t13: " + t13 + " t23: " + t23 + " t33: " + t33 + " t43: " + t43);

        assertEquals(0.9022407948835627, t13, 0.0000001);
        assertEquals(0.34335397764297787, t23, 0.0000001);
        assertEquals(0.4294050722878022, t33, 0.0000001);
        assertEquals(7.188591047255512E-5, t43, 0.0000001);


	// Radius 10um
	params[3] = 1E-5;
        t11 = cyl.getSignal(params, 0.005, 0.04, G1);
        t21 = cyl.getSignal(params, 0.005, 0.4, G1);
        t31 = cyl.getSignal(params, 0.015, 0.04, G1);
        t41 = cyl.getSignal(params, 0.015, 0.4, G1);

        //System.err.println("t11: " + t11 + " t21: " + t21 + " t31: " + t31 + " t41: " + t41);

        assertEquals(0.9429172882166851, t11, 0.0000001);
        assertEquals(0.897869656998364, t21, 0.0000001);
        assertEquals(0.6302459688524225, t31, 0.0000001);
        assertEquals(0.40442838173427287, t41, 0.0000001);

	t12 = cyl.getSignal(params, 0.005, 0.04, G2);
	t22 = cyl.getSignal(params, 0.005, 0.4, G2);
        t32 = cyl.getSignal(params, 0.015, 0.04, G2);
        t42 = cyl.getSignal(params, 0.015, 0.4, G2);

        //System.err.println("t12: " + t12 + " t22: " + t22 + " t32: " + t32 + " t42: " + t42);

        assertEquals(0.0028013014885340136, t12, 0.0000001);
        assertEquals(2.0957018990393086E-5, t22, 0.0000001);
        assertEquals(8.933212284708102E-21, t32, 0.0000001);
        assertEquals(4.832396541236077E-40, t42, 0.0000001);

        t13 = cyl.getSignal(params, 0.005, 0.04, G3);
	t23 = cyl.getSignal(params, 0.005, 0.4, G3);
        t33 = cyl.getSignal(params, 0.015, 0.04, G3);
        t43 = cyl.getSignal(params, 0.015, 0.4, G3);

        //System.err.println("t13: " + t13 + " t23: " + t23 + " t33: " + t33 + " t43: " + t43);

        assertEquals(0.9022407948835627, t13, 0.0000001);
        assertEquals(0.34335397764297787, t23, 0.0000001);
        assertEquals(0.4294050722878022, t33, 0.0000001);
        assertEquals(7.188591047255512E-5, t43, 0.0000001);


	// Radius 20um
	params[3] = 2E-5;
        t11 = cyl.getSignal(params, 0.005, 0.04, G1);
        t21 = cyl.getSignal(params, 0.005, 0.4, G1);
        t31 = cyl.getSignal(params, 0.015, 0.04, G1);
        t41 = cyl.getSignal(params, 0.015, 0.4, G1);

        //System.err.println("t11: " + t11 + " t21: " + t21 + " t31: " + t31 + " t41: " + t41);

        assertEquals(0.9214804699235097, t11, 0.0000001);
        assertEquals(0.6802931891711229, t21, 0.0000001);
        assertEquals(0.5160007413414626, t31, 0.0000001);
        assertEquals(0.03356777238152799, t41, 0.0000001);

	t12 = cyl.getSignal(params, 0.005, 0.04, G2);
	t22 = cyl.getSignal(params, 0.005, 0.4, G2);
        t32 = cyl.getSignal(params, 0.015, 0.04, G2);
        t42 = cyl.getSignal(params, 0.015, 0.4, G2);

        //System.err.println("t12: " + t12 + " t22: " + t22 + " t32: " + t32 + " t42: " + t42);

        assertEquals(2.8093995673086095E-4, t12, 0.0000001);
        assertEquals(1.8604255687867848E-17, t22, 0.0000001);
        assertEquals(1.8409100053599087E-29, t32, 0.0000001);
        assertEquals(3.9106789455025054E-148, t42, 0.0000001);

        t13 = cyl.getSignal(params, 0.005, 0.04, G3);
	t23 = cyl.getSignal(params, 0.005, 0.4, G3);
        t33 = cyl.getSignal(params, 0.015, 0.04, G3);
        t43 = cyl.getSignal(params, 0.015, 0.4, G3);

        //System.err.println("t13: " + t13 + " t23: " + t23 + " t33: " + t33 + " t43: " + t43);

        assertEquals(0.9022407948835627, t13, 0.0000001);
        assertEquals(0.34335397764297787, t23, 0.0000001);
        assertEquals(0.4294050722878022, t33, 0.0000001);
        assertEquals(7.188591047255512E-5, t43, 0.0000001);


    }
  

}




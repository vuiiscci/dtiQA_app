package sphfunc;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for ImagSH
 * <BR>
 * </dl>
 *
 * @version $Id: TestGreatCircleIntegrals.java,v 1.1 2006/01/05 12:28:58 ucacmgh Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestGreatCircleIntegrals extends TestCase {

    private static final double Rt2 = 1.414213562;


    /**
     * a few normal directions to test things in
     */
    double[] u=new double[] {1.0/Rt2, 1.0/Rt2, 0.0};



    public TestGreatCircleIntegrals(String name) {
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
    
    protected void tearDown() {
	// does the opposite of setup. 
	// Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestGreatCircleIntegrals.class);
    }
    

    public void testBasisFunctions(){

	RealSH rsh = new RealSH(0, 0);
	ImagSH ish = new ImagSH(0, 0);
	TuchRBF trbf = new TuchRBF(u, 1.0);


	assertEquals(0.5328, trbf.greatCircleIntegral(u), 1E-4);

	assertEquals(1.7724, rsh.greatCircleIntegral(u), 1E-4);
	assertEquals(0.0, ish.greatCircleIntegral(u), 1E-4);

	rsh= new RealSH(2, 0);
	ish= new ImagSH(2, 0);

	assertEquals(0.9908, rsh.greatCircleIntegral(u), 1E-4);
	assertEquals(0.0, ish.greatCircleIntegral(u), 1E-4);

	rsh= new RealSH(2, 2);
	ish= new ImagSH(2, 2);

	assertEquals(0.0, rsh.greatCircleIntegral(u), 1E-4);
	assertEquals(-1.2135, ish.greatCircleIntegral(u), 1E-4);

	rsh= new RealSH(4, 0);
	ish= new ImagSH(4, 0);

	assertEquals(0.7477, rsh.greatCircleIntegral(u), 1E-4);
	assertEquals(0.0, ish.greatCircleIntegral(u), 1E-4);

	rsh= new RealSH(4, 3);
	ish= new ImagSH(4, 3);

	assertEquals(0.0, rsh.greatCircleIntegral(u), 1E-4);
	assertEquals(0.0, ish.greatCircleIntegral(u), 1E-4);
    }



}

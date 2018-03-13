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
 * <dt>Purpose: Automated tests for <code>SphereGPD.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>SphereGPD</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/SphereGPD.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestSphereGPD extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestSphereGPD(String name) {
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
	return new TestSuite(TestSphereGPD.class);
    }


    public void testGetSignalScheme() {

	SphereGPD sphere = new SphereGPD();

	// diff, R
	double[] params = {1e-9, 1e-6};

	RealMatrix signals = sphere.getSignals(params, ip);
	RealMatrix jac = sphere.getJacobian(params, ip);
	
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
	assertEquals(0.9987279691746842, signals.entries[3][0], 0.0000001);
	assertEquals(0.9987254109212748, signals.entries[97][0], 0.0000001);
	assertEquals(0.9991701146908735, signals.entries[191][0], 0.0000001);
	assertEquals(0.9977557531825973, signals.entries[285][0], 0.0000001);


	assertEquals(0.0, jac.entries[1][0], 0.0000001);
	assertEquals(1241061.5, jac.entries[3][0], 0.0000001);
	assertEquals(1243330.0, jac.entries[97][0], 0.0000001);
	assertEquals(803539.5, jac.entries[191][0], 0.0000001);
	assertEquals(2212876.0, jac.entries[285][0], 0.0000001);

	assertEquals(0.0, jac.entries[1][1], 0.0000001);
	assertEquals(-5025.634765625, jac.entries[3][1], 0.0000001);
	assertEquals(-5036.296875, jac.entries[97][1], 0.0000001);
	assertEquals(-3266.304931640625, jac.entries[191][1], 0.0000001);
	assertEquals(-8907.755859375, jac.entries[285][1], 0.0000001);

    }	

}



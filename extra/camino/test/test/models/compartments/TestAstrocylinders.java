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
 * <dt>Purpose: Automated tests for <code>Astrocylinders.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Astrocylinders</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/Astrocylinders.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestAstrocylinders extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestAstrocylinders(String name) {
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
	return new TestSuite(TestAstrocylinders.class);
    }


    public void testGetSignalScheme() {

	Astrocylinders astrocyl = new Astrocylinders();

	// diff, R
	double[] params = {1e-9, 1e-6};

	RealMatrix signals = astrocyl.getSignals(params, ip);
	RealMatrix jac = astrocyl.getJacobian(params, ip);
	
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
	assertEquals(0.6059389250061481, signals.entries[3][0], 0.0000001);
	assertEquals(0.605145403176657, signals.entries[97][0], 0.0000001);
	assertEquals(0.4968023945468477, signals.entries[191][0], 0.0000001);
	assertEquals(0.24316735924445906, signals.entries[285][0], 0.0000001);


	assertEquals(0.0, jac.entries[1][0], 0.0000001);
	assertEquals(-2.29351808E8, jac.entries[3][0], 0.0000001);
	assertEquals(-2.2943608E8, jac.entries[97][0], 0.0000001);
	assertEquals(-2.25299936E8, jac.entries[191][0], 0.0000001);
	assertEquals(-1.20793976E8, jac.entries[285][0], 0.0000001);

	assertEquals(0.0, jac.entries[1][1], 0.0000001);
	assertEquals(-3864.534423828125, jac.entries[3][1], 0.0000001);
	assertEquals(-3869.310791015625, jac.entries[97][1], 0.0000001);
	assertEquals(-2180.40087890625, jac.entries[191][1], 0.0000001);
	assertEquals(-3312.765380859375, jac.entries[285][1], 0.0000001);

    }	

}



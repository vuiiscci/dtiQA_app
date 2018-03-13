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
 * <dt>Purpose: Automated tests for <code>Astrosticks.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Astrosticks</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/Astrosticks.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestAstrosticks extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestAstrosticks(String name) {
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
	return new TestSuite(TestAstrosticks.class);
    }


    public void testGetSignalScheme() {

	Astrosticks astrosticks = new Astrosticks();

	// diff
	double[] params = {1e-9};

	RealMatrix signals = astrosticks.getSignals(params, ip);
	RealMatrix jac = astrosticks.getJacobian(params, ip);
	
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
	assertEquals(0.6069207641976696, signals.entries[3][0], 0.0000001);
	assertEquals(0.6061283729910426, signals.entries[97][0], 0.0000001);
	assertEquals(0.49735908819272595, signals.entries[191][0], 0.0000001);
	assertEquals(0.24400413128496698, signals.entries[285][0], 0.0000001);


	assertEquals(0.0, jac.entries[1][0], 0.0000001);
	assertEquals(-2.30553184E8, jac.entries[3][0], 0.0000001);
	assertEquals(-2.3063928E8, jac.entries[97][0], 0.0000001);
	assertEquals(-2.26016128E8, jac.entries[191][0], 0.0000001);
	assertEquals(-1.22001136E8, jac.entries[285][0], 0.0000001);



    }	

}



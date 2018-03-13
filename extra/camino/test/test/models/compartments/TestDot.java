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
 * <dt>Purpose: Automated tests for <code>Dot.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Dot</code> with JUnit.
 * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Danny
 * @see models/compartments/Dot.java
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestDot extends TestCase {

    // Set of imaging parameters.
    private DW_Scheme ip;


    public TestDot(String name) {
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
	return new TestSuite(TestDot.class);
    }


    public void testGetSignalScheme() {

	Dot dot = new Dot();

	double[] params = {};

	RealMatrix signals = dot.getSignals(params, ip);
	RealMatrix jac = dot.getJacobian(params, ip);
	
//	int[] indices = {1,3,97,191,285};
//	for(int i=0; i<indices.length; i++)
//	{
//		System.out.println("assertEquals("+signals.entries[indices[i]][0]+", signals.entries["+indices[i]+"][0], 0.0000001);");
//	}

	

	

		assertEquals(1.0, signals.entries[1][0], 0.0000001);
		assertEquals(1, signals.entries[3][0], 0.0000001);
		assertEquals(1, signals.entries[97][0], 0.0000001);
		assertEquals(1, signals.entries[191][0], 0.0000001);
		assertEquals(1, signals.entries[285][0], 0.0000001);
         
              
             
		

    }	

}



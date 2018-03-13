package optimizers;

import junit.framework.*;
import junit.extensions.*;
import java.util.Random;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>ConjGradMinimizer.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>ConjGradMinimizer</code> with JUnit 3.8.
  * 
 * </dl>
 *
 * @version $Id: TestConjGradMinimizer.java,v 1.1 2005/09/26 10:36:50 ucacpco Exp $
 * @author  Philip Cook
 * @see optimizers.ConjGradMinimizer
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestConjGradMinimizer extends TestCase {


    public TestConjGradMinimizer(String name) {
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
    
    public static Test suite() {
	return new TestSuite(TestConjGradMinimizer.class);
    }


    public void testSimpleFunction() {

	ConjGradMinimizer minimizer = new ConjGradMinimizer() {
		
		public double fObj(double[] atry) {
		    return atry[0] * atry[0] + atry[1] * atry[1] + atry[2] * atry[2]; 
		}

		public double[] dfObj(double[] atry) {
		    double[] der = new double[4];
		    der[1] = 2.0 * atry[1];
		    der[2] = 2.0 * atry[2];
		    der[3] = 2.0 * atry[3]; 

		    return der;
		}

		
	    };
	
	minimizer.init(3);
	
	double[] atry = new double[] {0.0, 100.0, 10000.0, 1000.0};

	try {
	    minimizer.minimise(atry, 1E-9);
	    
	}
	catch (ConjGradMinimizerException e) {
	    fail(e.toString());
	}

	
	assertEquals(0.0, atry[1], 0.0000001);
	assertEquals(0.0, atry[2], 0.0000001);
	assertEquals(0.0, atry[3], 0.0000001);


    }
  

}




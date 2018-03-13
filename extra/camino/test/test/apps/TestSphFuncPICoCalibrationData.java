package apps;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import imaging.*;

import java.util.Random;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>SphFuncPICoCalibrationData</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestSphFuncPICoCalibrationData.java,v 1.1 2006/06/30 14:16:33 ucacpco Exp $
 * @author  Philip Cook
 * @see apps.SphFuncPICoCalibrationData
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestSphFuncPICoCalibrationData extends TestCase {



    public TestSphFuncPICoCalibrationData(String name) {
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
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestSphFuncPICoCalibrationData.class);
    }


  
    public void testOneTensorBlock() {

        // make sure that we get the expected block of single tensor data

        double[] fa = new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
        
        DT[] block = SphFuncPICoCalibrationData.getOneTensorBlock(fa, new Random(123));

        // want cylindrically symmetric tensors with FA as described in array

        assertEquals(fa.length, block.length);

        for (int i = 0; i < fa.length; i++) {
            assertEquals(fa[i], block[i].fa(), 1E-6);

            double[][] seig = block[i].sortedEigenSystem();

            assertEquals(seig[0][1], seig[0][2], 1E-6);
        }
        
    }


    public void testTwoTensorBlock() {

        double[] fa = new double[] {0.4, 0.5, 0.6, 0.7, 0.8};

        double[] mix = new double[] {0.4, 0.5, 0.6};

        double dt2RotAngle = Math.PI / 6.0;
              
        DT[][] block = SphFuncPICoCalibrationData.getTwoTensorBlock(fa, mix, dt2RotAngle, new Random(123));

        double expectedE1Dot = Math.cos(Math.PI / 2.0 - dt2RotAngle);

        assertEquals(mix.length * fa.length * (fa.length + 1) / 2, block.length);

        int counter = 0;

        for (int i = 0; i < fa.length; i++) {
            for (int j = 0; j <= i; j++) {
                for (int k = 0; k < mix.length; k++) { 

                    assertEquals(counter, mix.length * (i * (i+1) / 2) + mix.length * j + k);

                    assertEquals(fa[i], block[counter][0].fa(), 1E-6);
                    assertEquals(fa[j], block[counter][1].fa(), 1E-6);
                    
                    double[][] seig = block[counter][0].sortedEigenSystem();
                    
                    assertEquals(seig[0][1], seig[0][2], 1E-6);

                    seig = block[counter][1].sortedEigenSystem();

                    assertEquals(seig[0][1], seig[0][2], 1E-6);

                    counter++;
                }
                
            }

        }


    
    }
    
}

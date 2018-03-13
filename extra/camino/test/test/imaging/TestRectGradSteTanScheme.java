package imaging;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 *
 * Automated tests for .
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see imaging.RectGradSteTanScheme
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestRectGradSteTanScheme extends TestStejskalTannerScheme {

    public TestRectGradSteTanScheme(String name) {
	super(name);
    }

    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    
    public static Test suite() {
	return new TestSuite(TestRectGradSteTanScheme.class);
    }


    public void testGetGradImpulse() {

	// total echo time TE = DELTA + delta + (time before and after the DW gradients)
	// (time before and after each DW gradient) = (TE - DELTA - delta) / 2.0 = P

	// time (gradient)  
	//
	// 0---(0)---P---(|g|)---P+delta---(0)---P+DELTA---(|g|)---P+DELTA+delta---(0)---TE


	//  if after TE, grad strength is zero
        double[] impulse = stScheme.getGradImpulse(1, 10000.0, 10001.0);

        for (int j = 0; j < 3; j++) {
            assertEquals(0.0, impulse[j], 1E-8);
        }

	for (int i = 0; i < 15; i++) {
	    
	    double delta = stScheme.getDelta(i);
	    double DELTA = stScheme.getDELTA(i);
	    double TE = stScheme.getTE(i);

	    double pad = (TE - DELTA - delta) / 2.0; 

            double[] gDir = stScheme.getG_Dir(i);

	    // close to beginning of measurement
            impulse = stScheme.getGradImpulse(i, 0.75 * pad, pad / 2.0);

            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, impulse[j], 1E-8);
            }
            

	    // first DW
            impulse = stScheme.getGradImpulse(i, pad + 0.75 * delta, pad + delta / 2.0);
            

            for (int j = 0; j < 3; j++) {
                assertEquals(stScheme.getModG(i) * delta * 0.25 * gDir[j], impulse[j], 1E-8);
            }
                

	    // between DWs
            impulse = stScheme.getGradImpulse(i, pad + DELTA, pad + delta);
       
            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, impulse[j], 1E-8);
            }

	    // between DWs with part of the second one
            impulse = stScheme.getGradImpulse(i, pad + DELTA + delta / 2.0, pad + DELTA - (DELTA - delta) / 2.0);

            for (int j = 0; j < 3; j++) {
                assertEquals(stScheme.getModG(i) * delta * -0.5 * gDir[j], impulse[j], 1E-8);
            }
            
            
	    impulse = stScheme.getGradImpulse(i, pad + DELTA + 0.75 * delta, pad + DELTA + delta / 2.0);
                
            // second DW
            for (int j = 0; j < 3; j++) {
                assertEquals(stScheme.getModG(i) * delta * -0.25 * gDir[j], impulse[j], 1E-8);
            }

            
            impulse = stScheme.getGradImpulse(i, TE - pad / 2.0, TE - pad);
            
	    // end
            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, impulse[j], 1E-8);
            }
	    
	}
        

    }
    


    public String getSchemeFile() {
	return "test/imaging/v1_test.scheme";
    }

   
  


}

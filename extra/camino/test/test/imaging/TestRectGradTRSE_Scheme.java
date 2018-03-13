package imaging;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 *
 * Automated tests for RectGradTRSE_Scheme.
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see imaging.RectGradTRSE_Scheme
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestRectGradTRSE_Scheme extends TestTRSE_Scheme {

    public TestRectGradTRSE_Scheme(String name) {
	super(name);
    }

    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    
    public static Test suite() {
	return new TestSuite(TestRectGradTRSE_Scheme.class);
    }


    public void testGetGradImpulse() {

	double[] impulse = trseScheme.getGradImpulse(1, 100001.0, 100000.0);

        for (int j = 0; j < 3; j++) {
            assertEquals(0.0, impulse[j], 1E-8);
        }

	for (int i = 0; i < 15; i++) {

	    double t_del1 = trseScheme.getT_Del1(i);
	    double del1 = trseScheme.getDel1(i);
	    double t_del2 = trseScheme.getT_Del2(i);
	    double del2 = trseScheme.getDel2(i);
	    double t_del3 = trseScheme.getT_Del3(i);
	    double del3 = trseScheme.getDel3(i);
	    double t_del4 = trseScheme.getT_Del4(i);
	    double del4 = trseScheme.getDel4(i);

            double[] gDir = trseScheme.getG_Dir(i);
	    
	    // before first pulse
            impulse = trseScheme.getGradImpulse(i, t_del1 - t_del1 / 4.0, t_del1 - t_del1 / 2.0);

            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, impulse[j], 1E-8);
            }


            impulse = trseScheme.getGradImpulse(i, t_del1 + del1, t_del1 + del1 / 2.0);

            for (int j = 0; j < 3; j++) {
                assertEquals(trseScheme.getModG(i) * (del1 / 2.0) * gDir[j], impulse[j], 1E-8);
            }
            

	    // include 1ms after del1 ends, one before
            impulse = trseScheme.getGradImpulse(i, t_del1 + del1 + 1E-3, t_del1 + del1 - 1E-3);


	    for (int j = 0; j < 3; j++) {
                assertEquals(trseScheme.getModG(i) * 1E-3 * gDir[j], impulse[j], 1E-8);
            }

            
            impulse = trseScheme.getGradImpulse(i, t_del2, t_del1 + del1);

            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, impulse[j], 1E-8);
            }


            impulse = trseScheme.getGradImpulse(i, t_del2 + del2, t_del2 + del2 / 2.0);
            
            for (int j = 0; j < 3; j++) {
                assertEquals(trseScheme.getModG(i) * (del2 / 2.0) * gDir[j], impulse[j], 1E-8);
            }

            impulse = trseScheme.getGradImpulse(i, t_del3 + del3, t_del3 + del3 / 2.0);

            for (int j = 0; j < 3; j++) {
                assertEquals(-1.0 * trseScheme.getModG(i) * (del3 / 2.0) * gDir[j], impulse[j] , 1E-8);
            }
            
            impulse = trseScheme.getGradImpulse(i, t_del4, t_del3 + del3);
            
            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, impulse[j], 1E-8);
            }
            
            impulse = trseScheme.getGradImpulse(i, t_del4 + del4, t_del4 + del4 / 2.0);

            for (int j = 0; j < 3; j++) {
                assertEquals(-1.0 * trseScheme.getModG(i) * (del4 / 2.0) * gDir[j], impulse[j], 1E-8);
            }

            impulse = trseScheme.getGradImpulse(i, t_del4 + del4 + 1.0, t_del4 + del4);

            for (int j = 0; j < 3; j++) {
                assertEquals(0.0, impulse[j], 1E-8);
            }
	    
	}


    }

	
    protected String getSchemeFile() {
	return "test/imaging/v3_test.scheme";
    }
	
 


}

package imaging;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 *
 * Automated tests for TRSE_Scheme.
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see imaging.StejskalTannerScheme
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public abstract class TestTRSE_Scheme extends TestDW_Scheme {
    
    protected final TRSE_Scheme trseScheme;


    public TestTRSE_Scheme(String name) {
	super(name);
	trseScheme = (TRSE_Scheme)scheme;
    }

    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    
    public static Test suite() {
	return new TestSuite(TestTRSE_Scheme.class);
    }



    public void testGetG() {

	for (int i = 0; i < 15; i++) {

	    double[] g = trseScheme.getG(i);

	    double[] gDir = trseScheme.getG_Dir(i);

	    double modG = trseScheme.getModG(i);

	    assertEquals(modG * gDir[0], g[0], 1E-8);
	    assertEquals(modG * gDir[1], g[1], 1E-8);
	    assertEquals(modG * gDir[2], g[2], 1E-8);

	}
    }


    public void testGetNonZeroGs() {

	double[][] nzg = trseScheme.getNonZeroGs();

	for (int i = 0; i < 7; i++) {
	    double[] g = trseScheme.getG(i+1);

	    for (int j = 0; j < 3; j++) {
		assertEquals(g[j], nzg[i][j], 1E-8);
	    }
	}

	for (int i = 0; i < 6; i++) {
	    double[] g = trseScheme.getG(i+9);

	    for (int j = 0; j < 3; j++) {
		assertEquals(g[j], nzg[i+7][j], 1E-8);
	    }
	}
 
    }


    /**
     * Tests the underlying getNonZeroParam(double[]) method. 
     *
     */
    public void testNonZeroModGs() {
	double[] nzModG = trseScheme.getNonZeroModGs();

	assertEquals(13, nzModG.length);

	for (int i = 0; i < 7; i++) {
	    assertEquals(0.039887, nzModG[i], 1E-8);
	}

	for (int i = 0; i < 6; i++) {
	    assertEquals(0.012614, nzModG[i+7], 1E-8);
	}

	
    }


}

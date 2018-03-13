package imaging;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 *
 * Automated tests for StejskalTannerScheme.
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see imaging.StejskalTannerScheme
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public abstract class TestStejskalTannerScheme extends TestDW_Scheme {

    protected final StejskalTannerScheme stScheme;

    public TestStejskalTannerScheme(String name) {
	super(name);
	stScheme = (StejskalTannerScheme)scheme;
    }

    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    
    public static Test suite() {
	return new TestSuite(TestStejskalTannerScheme.class);
    }


    public void testGetQs() {

	for (int i = 0; i < 15; i++) {

	    double[] q = stScheme.getQ(i);
	    
	    double[] gDir = stScheme.getG_Dir(i);

	    double modQ = DW_Scheme.GAMMA * stScheme.getDelta(i) * stScheme.getModG(i);

	    assertEquals(modQ * gDir[0], q[0], 1E-8);
	    assertEquals(modQ * gDir[1], q[1], 1E-8);
	    assertEquals(modQ * gDir[2], q[2], 1E-8);

	}
    }


    public void testGetG() {

	for (int i = 0; i < 15; i++) {

	    double[] g = stScheme.getG(i);

	    double[] gDir = stScheme.getG_Dir(i);

	    double modG = stScheme.getModG(i);

	    assertEquals(modG * gDir[0], g[0], 1E-8);
	    assertEquals(modG * gDir[1], g[1], 1E-8);
	    assertEquals(modG * gDir[2], g[2], 1E-8);

	}
    }


    public void testGetNormNonZeroQs() {
	double norm = 0.0;

	for (int i = 0; i < 7; i++) {
	    norm += stScheme.getModQ(i+1);
	}

	for (int i = 0; i < 6; i++) {
	    norm += stScheme.getModQ(i+9);
	}
	
	norm /= 13.0;

	double[][] nzq = stScheme.getNormNonZeroQs();

	for (int i = 0; i < 7; i++) {
	    double[] q = stScheme.getQ(i+1);

	    for (int j = 0; j < 3; j++) {
		assertEquals(q[j] / norm, nzq[i][j], 1E-8);
	    }
	}

	for (int i = 0; i < 6; i++) {
	    double[] q = stScheme.getQ(i+9);

	    for (int j = 0; j < 3; j++) {
		assertEquals(q[j] / norm, nzq[i+7][j], 1E-8);
	    }
	}
 

    }


    public void testGetNonZeroQs() {

	double[][] nzq = stScheme.getNonZeroQs();

	for (int i = 0; i < 7; i++) {
	    double[] q = stScheme.getQ(i+1);

	    for (int j = 0; j < 3; j++) {
		assertEquals(q[j], nzq[i][j], 1E-8);
	    }
	}

	for (int i = 0; i < 6; i++) {
	    double[] q = stScheme.getQ(i+9);

	    for (int j = 0; j < 3; j++) {
		assertEquals(q[j], nzq[i+7][j], 1E-8);
	    }
	}
 

    }


    public void testGetNonZeroGs() {

	double[][] nzg = stScheme.getNonZeroGs();

	for (int i = 0; i < 7; i++) {
	    double[] g = stScheme.getG(i+1);

	    for (int j = 0; j < 3; j++) {
		assertEquals(g[j], nzg[i][j], 1E-8);
	    }
	}

	for (int i = 0; i < 6; i++) {
	    double[] g = stScheme.getG(i+9);

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
	double[] nzModG = stScheme.getNonZeroModGs();

	assertEquals(13, nzModG.length);

	for (int i = 0; i < 7; i++) {
	    assertEquals(0.022138, nzModG[i], 1E-8);
	}

	for (int i = 0; i < 6; i++) {
	    assertEquals(0.007001, nzModG[i+7], 1E-8);
	}

	
    }


}

package apps;

import junit.framework.*;
import junit.extensions.*;

import imaging.*;
import numerics.*;


/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>RGB_ScalarImage</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestSphFuncPICoCalibrationData.java,v 1.1 2006/06/30 14:16:33 ucacpco Exp $
 * @author  Philip Cook
 * @see apps.RGB_ScalarImage
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestRGB_ScalarImage extends TestCase {


    private Vector3D[][][][] vectors;

    private ImageHeader ih;


    public TestRGB_ScalarImage(String name) {
	super(name);
    }



    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    public static Test suite() {
	return new TestSuite(TestRGB_ScalarImage.class);
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

	vectors = new Vector3D[2][2][2][];

	// 8 vectors

	// gray
	vectors[0][0][0] = new Vector3D[] {new Vector3D(1.0, 1.0, -1.0).normalized()};

	// blue
	vectors[0][0][1] = new Vector3D[] {new Vector3D(0.0, 0.0, 1.0)};

	// green
	vectors[0][1][0] = new Vector3D[] {new Vector3D(0.0, -1.0, 0.0)};

	// cyan (but darker)
	vectors[0][1][1] = new Vector3D[] {new Vector3D(0.0, 1.0, 1.0).normalized()};

	// red
	vectors[1][0][0] = new Vector3D[] {new Vector3D(-1.0, 0.0, 0.0)};

	// magenta (but darker)
	vectors[1][0][1] = new Vector3D[] {new Vector3D(1.0, 0.0, 1.0).normalized()};

	// yellow (but darker)
	vectors[1][1][0] = new Vector3D[] {new Vector3D(1.0, 1.0, 0.0).normalized()};

	// and lastly, a mixture of red and green, which should also be yellow
	vectors[1][1][1] = new Vector3D[] {new Vector3D(1.0, 0.0, 0.0), 
					   new Vector3D(0.0, 1.0, 0.0)};
	
        Nifti1Dataset nds = new Nifti1Dataset();

        nds.setDims(3, 2, 2, 2, 1, 0, 0, 0);
        nds.setPixDims(1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f);

        nds.setQuaternion((short)1, (short)1, new float[] {1.0f, 0.0f, 0.0f}, new float[] {0.5f, 0.5f, 0.5f});

        ih = nds;


    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    

    public void testRGB_Mapping() {

	double[][][] scalars = new double[2][2][2];

	for (int k = 0; k < 2; k++) {
	    for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 2; i++) {
		    scalars[i][j][k] = 1.0;
		}
	    }
	}
	
	RGB_ScalarImage image = new RGB_ScalarImage(vectors, ih, scalars, 0.0, 1.0);
	
	assertEquals(256 * (147 * 256 + 147) + 147, image.rgbIndex(0,0,0));
	assertEquals(256 * (000 * 256 + 000) + 255, image.rgbIndex(0,0,1));
	assertEquals(256 * (000 * 256 + 255) + 000, image.rgbIndex(0,1,0));
	assertEquals(256 * (000 * 256 + 180) + 180, image.rgbIndex(0,1,1));
	assertEquals(256 * (255 * 256 + 000) + 000, image.rgbIndex(1,0,0));
	assertEquals(256 * (180 * 256 + 000) + 180, image.rgbIndex(1,0,1));
	assertEquals(256 * (180 * 256 + 180) + 000, image.rgbIndex(1,1,0));
	assertEquals(256 * (180 * 256 + 180) + 000, image.rgbIndex(1,1,1));


	// test gamma

	image.setRGB_Gamma(1.5);

	assertEquals(256 * (151 * 256 + 000) + 151, image.rgbIndex(1,0,1));

    }

    public void testScalarMapping() {

	// test normalized scalars
	double[][][] scalars = new double[2][2][2];
	
	for (int k = 0; k < 2; k++) {
	    for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 2; i++) {
		    scalars[i][j][k] = 1.0 * i + 2.0 * j + 4.0 * k;
		}
	    }
	}

	// test calculate scalar range
	RGB_ScalarImage image = new RGB_ScalarImage(vectors, ih, scalars, 
						    0.0, 0.0);

	for (int k = 0; k < 2; k++) {
	    for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 2; i++) {
		    assertEquals(scalars[i][j][k] / 7.0, image.normScalarVol[i][j][k], 1E-6);
		}
	    }
	}

	image.setScalarGamma(0.5);
	image.setRGB_Gamma(0.0);

	assertEquals((int)(255 * Math.sqrt(1.0 / 7.0)), image.rgbIndex(1,0,0) & 0xff);


	// should chop one value off each end
	image.calculateScalarRange(0.13);

	assertEquals(1.0, image.minScalarValue, 1E-6);
	assertEquals(6.0, image.maxScalarValue, 1E-6);

	
    }

    
    
    
}

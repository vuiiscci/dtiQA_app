package imaging;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 *
 * Automated tests for B_VectorScheme. This class tests the reading of legacy
 * V0 scheme files by the B_Vector class.
 *
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * @see imaging.B_VectorScheme
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestQ_VectorScheme extends TestDW_Scheme {

    public TestQ_VectorScheme(String name) {
	super(name);
    }

    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    
    public static Test suite() {
	return new TestSuite(TestQ_VectorScheme.class);
    }


    public void testToString() {
	// this scheme format is deprecated, toString returns scheme in B_Vector format
    }


    protected String getSchemeFile() {
	return "test/imaging/v0_test.scheme";
    }
}

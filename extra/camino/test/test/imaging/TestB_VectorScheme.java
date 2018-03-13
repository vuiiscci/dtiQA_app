package imaging;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 *
 * Automated tests for B_VectorScheme. Note: the old scheme V0 is tested in
 * a separate class. This class tests the reading of BVECTOR scheme files.
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
public class TestB_VectorScheme extends TestDW_Scheme {

    public TestB_VectorScheme(String name) {
	super(name);
    }

    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

    
    public static Test suite() {
	return new TestSuite(TestB_VectorScheme.class);
    }

    
    protected String getSchemeFile() {
	return "test/imaging/v2_test.scheme";
    }

 

}

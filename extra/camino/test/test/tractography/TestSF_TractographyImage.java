package tractography;

import data.*;
import misc.DT;
import numerics.*;

import junit.framework.*;
import junit.extensions.*;



/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>SF_TractographyImage</code>.
 * <BR><BR>
 *
 * </dl>
 *
 * @version $Id: TestSF_TractographyImage.java,v 1.1 2005/09/26 10:36:54 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestSF_TractographyImage extends TestTractographyImage {

    
    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

    }

   
    protected void tearDown() {
    }


    public TestSF_TractographyImage(String name) {
	super(name);

    }
    
    public static void main(String[] args) {
	junit.textui.TestRunner.run(suite());
    }

   
    public static Test suite() {
	return new TestSuite(TestSF_TractographyImage.class);
    }

    
    protected TractographyImage getCrossingImage() {
	return Images.getCrossingSF();
    }
    
}

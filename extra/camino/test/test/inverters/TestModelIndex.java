package inverters;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * Automated tests for ModelIndex.
 *
 * @version $Id$
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestModelIndex extends TestCase {




    public TestModelIndex(String name) {
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
	return new TestSuite(TestModelIndex.class);
    }


    /**
     * Non-exhaustive test of numeric indices.
     *
     */
    public void testNumericIndices() {

	ModelIndex[] indices = ModelIndex.getIndices(new String[] {"22"});

	assertEquals(2, indices.length);

	assertEquals(ModelIndex.CYLCYL_EQ, indices[0]);

	assertEquals(ModelIndex.NLDT_POS, indices[1]);


	indices = ModelIndex.getIndices(new String[] {"3"});

	assertEquals(1, indices.length);

	assertEquals(ModelIndex.NLDT_POS_HESS, indices[0]);

    }


    /**
     * Non-exhaustive test of text indices.
     *
     */
    public void testTextIndices() {
	
	ModelIndex[] indices = ModelIndex.getIndices(new String[] {"pospos", "ldt"});
	
	assertEquals(2, indices.length);
	
	assertEquals(ModelIndex.POSPOS, indices[0]);
	
	assertEquals(ModelIndex.LDT, indices[1]);
	
	
	indices = ModelIndex.getIndices(new String[] {"nldt_pos"});
	
	assertEquals(1, indices.length);
	
	assertEquals(ModelIndex.NLDT_POS, indices[0]);
	
    }



    /**
     * Non-exhaustive test of classified model indices.
     *
     */
    public void testClassifiedModelIndices() {
	
	ModelIndex[][] indices = ModelIndex.getClassifiedModelIndices(new String[] {"1", "2", "22"});
	
	assertEquals(3, indices.length);
	
	assertEquals(ModelIndex.LDT, indices[0][0]);
	
	assertEquals(ModelIndex.NLDT_POS, indices[1][0]);
		
	assertEquals(ModelIndex.CYLCYL_EQ, indices[2][0]);
	
	assertEquals(ModelIndex.NLDT_POS, indices[2][1]);
	
	indices = ModelIndex.getClassifiedModelIndices(new String[] {"dt", "dt", "cylcyl"});

	assertEquals(3, indices.length);
	
	assertEquals(ModelIndex.LDT_ALIAS, indices[0][0]);
	
	assertEquals(ModelIndex.LDT_ALIAS, indices[1][0]);
		
	assertEquals(ModelIndex.CYLCYL, indices[2][0]);
	
	assertEquals(ModelIndex.LDT_ALIAS, indices[2][1]);

    }


}

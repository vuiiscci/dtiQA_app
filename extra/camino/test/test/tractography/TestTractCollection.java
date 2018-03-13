package tractography;

import junit.framework.*;
import junit.extensions.*;
import numerics.Point3D;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>TractCollection</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>TractCollection</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestTractCollection.java,v 1.2 2005/11/08 11:47:48 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTractCollection extends TestCase {

    // Array of points used for testing
    Tract[] tracts;

    public TestTractCollection(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

	tracts = new Tract[10];
	
	for (int i = 0; i < 10; i++) {

	    tracts[i] = new Tract(20, 100.0);
	    
	    for (int j = 0; j < 10; j++) {
		tracts[i].addPoint( new Point3D((double)j, 0.0, 0.0) );
	    }

	}

    }

    
    protected void tearDown() {

	tracts = null;
    }


    public static Test suite() {
	return new TestSuite(TestTractCollection.class);
    }



    /**
     * Make sure that the tracts (number and value) are preserved when a TractCollection increases its capacity.
     *
     */
    public void testGrowth() {

	TractCollection tc = new TractCollection(3, 100.0);

	for (int i = 0; i < tracts.length; i++) {
	    tc.addTract(tracts[i]);
	}

	// Adding these tracts should trigger the growth of the TractCollection, twice
	// But the number of tracts should be as many as we added
	assertEquals(10, tc.numberOfTracts());
	
	// check that the tracts have been copied correctly
	for (int i = 0; i < 10; i++) {
	    assertTrue( tracts[i].equals( tc.getTract(i) ) );
	}
	

	
    }

    
    public void testAddTractCollection() {


	// sufficient capacity
	TractCollection tc1 = new TractCollection(21, 100.0);

	// insufficient capacity
	TractCollection tc2 = new TractCollection(11, 100.0);
	
	TractCollection addThis = new TractCollection(11, 100.0);

	for (int i = 0; i < tracts.length; i++) {
	    tc1.addTract(tracts[i]);
	    tc2.addTract(tracts[i]);
	    addThis.addTract(tracts[i]);	    
	}

	tc1.addTractCollection(addThis);

	tc2.addTractCollection(addThis);

	for (int i = 0; i < 2 * tracts.length; i++) {
	    assertEquals(tc1.getTract(i), tracts[i % tracts.length]);
	    assertEquals(tc1.getTract(i), tc2.getTract(i));
	}

	

    }

    
}

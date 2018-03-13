package numerics;

import junit.framework.*;
import junit.extensions.*;

/**
 * Automated tests for <code>SymmetricMatrix</code>.
 *
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 */
public class TestSymmetricMatrix extends TestCase {


    public TestSymmetricMatrix(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	// any class variables you declare should be initialized here. This is called before each test
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestSymmetricMatrix.class);
    }



    public void testToRealMatrix() {

        SymmetricMatrix mat = new SymmetricMatrix(3);

        mat.set(0,0, 1.0);
        mat.set(0,1, 2.0);
        mat.set(2,0, 3.0);

        mat.set(1,1, 4.0);
        mat.set(2,1, 5.0);
        mat.set(2,2, 6.0);

        RealMatrix r = mat.toRealMatrix();

        assertEquals(1.0, r.entries[0][0], 1E-8);
        assertEquals(2.0, r.entries[0][1], 1E-8);
        assertEquals(3.0, r.entries[0][2], 1E-8);
        assertEquals(r.entries[0][1], r.entries[1][0], 1E-8);
        assertEquals(4.0, r.entries[1][1], 1E-8);
        assertEquals(5.0, r.entries[1][2], 1E-8);
        assertEquals(r.entries[0][2], r.entries[2][0], 1E-8);
        assertEquals(r.entries[1][2], r.entries[2][1], 1E-8);
        assertEquals(6.0, r.entries[2][2], 1E-8);
        
    }


    public void testAdd() {

        SymmetricMatrix mat = new SymmetricMatrix(2);
        
        mat.set(0,0, 1.0);
        mat.set(0,1, 2.0);

        mat.add(0, 1, 2.0);
        mat.add(1, 0, -1.0);


        assertEquals(3.0, mat.get(0,1), 1E-8);
        assertEquals(3.0, mat.get(1,0), 1E-8);

    }



    public void testScale() {

        SymmetricMatrix mat = new SymmetricMatrix(2);
        
        mat.set(0,0, 1.0);
        mat.set(0,1, 2.0);
        mat.add(1, 1, 3.0);

        mat.scale(2.0);

        assertEquals(2.0, mat.get(0,0), 1E-8);
        assertEquals(4.0, mat.get(0,1), 1E-8);
        assertEquals(6.0, mat.get(1,1), 1E-8);

        

        mat.scale(0,0,4.0);

        assertEquals(8.0, mat.get(0,0), 1E-8);
    }
  
}

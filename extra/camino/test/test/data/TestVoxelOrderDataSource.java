package data;

import junit.framework.*;
import junit.extensions.*;
import misc.*;
import data.*;
import imaging.*;
import numerics.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>VoxelOrderDataSource.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestVoxelOrderDataSource.java,v 1.2 2005/10/29 00:51:14 ucacpco Exp $
 * @author  Daniel Alexander
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestVoxelOrderDataSource extends TestCase {


    public TestVoxelOrderDataSource(String name) {
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
	return new TestSuite(TestVoxelOrderDataSource.class);
    }


    public void testVoxelOrderDataSource() {

	String[] files = new String[2];

	files[0] = "FourVoxelsFromG.Bfloat";
	files[1] = "FourVoxelsFromG.Bfloat.gz";

	for (int i = 0; i < 2; i++) {

	    
	    VoxelOrderDataSource v = new VoxelOrderDataSource(files[i], 66, "float");
	    
	    try {
		
		// Test the third voxel of four.
		double[] data = v.nextVoxel();
		data = v.nextVoxel();
		data = v.nextVoxel();
		
		assertEquals(data.length, 66);
		
		assertEquals("Failed while reading file " + i, data[0], 450.0, 0.0);
		assertEquals(data[10], 178.0, 0.0);
		assertEquals(data[20], 280.0, 0.0);
		assertEquals(data[30], 224.0, 0.0);
		assertEquals(data[40], 205.0, 0.0);
		assertEquals(data[50], 137.0, 0.0);
		assertEquals(data[60], 372.0, 0.0);
		
	    } catch(Exception e) {
		fail("Failed to read voxel and got exception " + e);
	    }
	}
	
    }


    public void testVoxelOrderDW_DataSource6() {


	String[] files = new String[2];

	files[0] = "FourVoxelsFromG.Bfloat";
	files[1] = "FourVoxelsFromG.Bfloat.gz";

	for (int i = 0; i < 2; i++) {

	    
	    VoxelOrderDW_DataSource v = 
		new VoxelOrderDW_DataSource(files[i], "float", DW_Scheme.readScheme("bmx6.scheme"));

	    try {
		
		// Test the third voxel of four.
		double[] data = v.nextVoxel();
		data = v.nextVoxel();
		data = v.nextVoxel();
		
		assertEquals(data.length, 66);
		
		assertEquals(data[0], 450.0, 0.0);
		assertEquals(data[10], 178.0, 0.0);
		assertEquals(data[20], 280.0, 0.0);
		assertEquals(data[30], 224.0, 0.0);
		assertEquals(data[40], 205.0, 0.0);
		assertEquals(data[50], 137.0, 0.0);
		assertEquals(data[60], 372.0, 0.0);
		
	    } catch(Exception e) {
		fail("Failed to read voxel.");
	    }
	    
	}

    }

 

 
}

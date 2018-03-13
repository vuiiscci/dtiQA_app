package tractography;

import junit.framework.*;
import junit.extensions.*;

import imaging.*;
import numerics.*;

import java.io.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>TractSource</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>TractSource</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestTractSource.java,v 1.2 2005/11/08 11:47:48 ucacpco Exp $
 * @author  Philip Cook
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestTractSource extends TestCase {

    private Voxel[] voxels = null;
    private VoxelList list = null;
    private Tract tract = null;

    public TestTractSource(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {

	voxels = new Voxel[4];
	
	voxels[0] = new Voxel(0, 0, 0);
	voxels[1] = new Voxel(1, 0, 0);
	voxels[2] = new Voxel(1, 1, 0);
	voxels[3] = new Voxel(1, 1, 1);

	list = new VoxelList(voxels, 2, 1.0, 2.0, 3.0, new Vector3D(0.0, 0.0, 1.0), 
				       new Vector3D(0.0, -1.0, 0.0));

	tract = list.toTract();

    }

    
    protected void tearDown() {
	voxels = null;
	list = null;
	tract = null;
    }


    public static Test suite() {
	return new TestSuite(TestTractSource.class);
    }


    /**
     * Test reading and writing of raw tracts.
     *
     */
    public void testRawTractIO() {

	String filename = "tractSourceTracts.Bfloat";
       
        Nifti1Dataset nds = new Nifti1Dataset();
        
        nds.setPixDims(1.0f, 2.0f, 3.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f);
        
        nds.setQuaternion((short)1, (short)1, new float[] {0.25f, 0.5f, 0.25f}, new float[] {100.0f, 200.0f, 300.0f});

        Tract phys = new Tract(tract);

        phys.transformToPhysicalSpace(nds.getVoxelToPhysicalTransform(), 2.0, 3.0, 4.0);

	try {

	    // try with 1 to 3 tracts in the file
	    for (int i = 1; i < 3; i++) {

		FileOutputStream fout = new FileOutputStream(filename);

		DataOutputStream dout = new DataOutputStream(fout);

		for (int t = 0; t < i; t++) {
		    phys.writeRaw(dout);
		}
		
		dout.close();
		
		TractSource source = new TractSource(filename, nds);
		
		for (int t = 0; t < i; t++) {
		    Tract read = source.nextTract();
	    
		    assertEquals(tract.numberOfPoints(), read.numberOfPoints());
		    
		    for (int p = 0; p < tract.numberOfPoints(); p++) {
			assertEquals(tract.getPoint(p).x, read.getPoint(p).x, 1E-4);
			assertEquals(tract.getPoint(p).y, read.getPoint(p).y, 1E-4);
			assertEquals(tract.getPoint(p).z, read.getPoint(p).z, 1E-4);
		    }
		    
		}

		assertFalse(source.more());
	    }

	}
	catch(IOException e) {
	    fail(e.toString());
	}
    }



    
}

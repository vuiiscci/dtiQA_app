package imaging;

import data.*;
import numerics.*;
import tools.*;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 * Tests for Nifti1Dataset
 * 
 * @author  Philip Cook
 *
 *
 *
 */
public class TestNifti1Dataset extends TestImageHeader {

    public TestNifti1Dataset(String name) {
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
    }


    public static Test suite() {
	return new TestSuite(TestNifti1Dataset.class);
    }


    /**
     * Writes a header to disk and then reads it back in.
     */ 
    public void testHeaderIO() {
	Nifti1Dataset nds = new Nifti1Dataset();
        
        // set some stuff

        String fileRoot = "./test/imaging/testNiftiHeader";
        
        nds.setFilename(fileRoot, true, false);

        nds.setDims(5, 1, 2, 3, 1, 6, 0, 0);

        nds.setPixDims(1.0f, 2.0f, 3.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f);
        
        nds.setDataType("short");

        // check basic I/O functions
        nds.setQuaternion((short)1, (short)1, new float[] {0.25f, 0.5f, 0.25f}, new float[] {100.0f, 200.0f, 300.0f});
        
	try {
    
            nds.writeHeader();
	    
	    Nifti1Dataset read = (Nifti1Dataset)Nifti1Dataset.readHeader(nds.getHeaderFilename());
            
            assertEquals(6, read.components());
	    assertEquals(1, read.xDataDim());
            assertEquals(2, read.yDataDim());
            assertEquals(3, read.zDataDim());

            assertEquals(2.0, read.xVoxelDim(), 1E-6);
            assertEquals(3.0, read.yVoxelDim(), 1E-6);
            assertEquals(4.0, read.zVoxelDim(), 1E-6);

            RealMatrix ndsMat = nds.getVoxelToPhysicalTransform();

            RealMatrix readMat = read.getVoxelToPhysicalTransform();

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    assertEquals(ndsMat.entries[i][j], readMat.entries[i][j], 1E-6);
                }
            } 

            // test ImageHeader methods

            assertTrue(nds.sameSpace(read));

            nds.setQuaternion((short)1, (short)1, new float[] {0.2f, 0.6f, 0.25f}, new float[] {100.0f, 200.0f, 300.0f});

            assertFalse(nds.sameSpace(read));

            
	}
	catch (IOException e) {
	    fail(e.toString());
	}
	

    }



    public void testVoxelToPhysicalTransform() {

	Nifti1Dataset nds = new Nifti1Dataset();
        
        nds.setDims(5, 1, 2, 3, 1, 6, 0, 0);

        nds.setPixDims(1.0f, 2.0f, 3.0f, 4.0f, 0.0f, 0.0f, 0.0f, 0.0f);
        
        nds.setDataType("short");

        nds.setQuaternion((short)1, (short)-1, new float[] {0.25f, 0.5f, 0.25f}, new float[] {100.0f, 200.0f, 300.0f});

        // Using the Nifti website formula and Octave, the quaternion * diag(pixdims) should translate to
        double[][] transform = new double[][] {
            {0.75000,  -0.43585,   -3.66228,  100.0},
            {1.29057,   2.25000,  0.58114,  200.0},
            {-1.33114,   1.93585,  -1.50000,  300.0},
            {0.0,       0.0,      0.0,      1.0}
        }; 


        RealMatrix ndsTrans = nds.getVoxelToPhysicalTransform();

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(transform[i][j], ndsTrans.entries[i][j], 1E-5);
            }
        } 


        // test ImageHeader physicalToVoxelTransform
        // compare to Octave solution
            
      
        double[][] invTrans = new double[][] {
            {0.18750,    0.32264,    -0.33278,    16.55717},
            {-0.04843,     0.25000,     0.21510,  -109.68578},
            {-0.22889,     0.03632,    -0.09375,    43.74991},
            {0.00000,     0.00000,     0.00000,     1.00000}
        };
            
        
        RealMatrix ndsInv = nds.getPhysicalToVoxelTransform();

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(invTrans[i][j], ndsInv.entries[i][j], 1E-3);
            }
         }        

  
        // now set sform
        RealMatrix sformMat = new RealMatrix(4,4);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) { 
                sformMat.entries[i][j] = i + j + i * j;
            }
        }
        
        sformMat.entries[3][0] = 0.0;
        sformMat.entries[3][1] = 0.0;
        sformMat.entries[3][2] = 0.0;
        sformMat.entries[3][3] = 1.0;

        nds.setSform((short)1, sformMat);

        // this should now be the transform
        ndsTrans = nds.getVoxelToPhysicalTransform();

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(sformMat.entries[i][j], ndsTrans.entries[i][j], 1E-6);
            }
        } 

    }
    

    
    public void testScaledDataSource() {

        int[] dataDims = new int[] {2, 3, 2};

        double[] voxelDims = new double[] {1.0, 2.0, 3.0};
        
        double[][][] data = new double[dataDims[0]][dataDims[1]][dataDims[2]];

        for (int k = 0; k < dataDims[2]; k++) {
            for (int j = 0; j < dataDims[1]; j++) {
                for (int i = 0; i < dataDims[0]; i++) {
                    data[i][j][k] = 4.0 * i + 2.0 * j + 8.0 * k + 4.0;
                }
            }
        }

        String fileRoot = "./test/imaging/testScaledDataSource";

        ImageHeader ih = getImageHeader(fileRoot, dataDims, voxelDims);
        
        
        // ImageHeader methods write without scale terms
        Nifti1Dataset written = (Nifti1Dataset)ih.writeScalarImage(data, fileRoot);
        
        written.setScale(2.0, 3.0);
                   
        DataSource scaled = written.getImageDataSource();

        for (int k = 0; k < dataDims[2]; k++) {
            for (int j = 0; j < dataDims[1]; j++) {
                for (int i = 0; i < dataDims[0]; i++) {
                    assertEquals(2.0 * data[i][j][k] + 3.0, scaled.nextVoxel()[0], 1E-6);
                }
            }
        }  
        
    }




    /**
     * Gets a header template for testing.
     *
     * @return a header with gzip compression, and the required file root, dimensions and scale parameters.
     *
     */
    protected ImageHeader getImageHeader(String fileRoot, int[] dataDims, double[] voxelDims) {

        Nifti1Dataset hdr = new Nifti1Dataset();

        hdr.setFilename(fileRoot, true, true);

        hdr.setDataType("float");

        hdr.setDims(3, dataDims[0],  dataDims[1],  dataDims[2], 1, 0, 0, 0);

        hdr.setPixDims(1.0f, (float)voxelDims[0], (float)voxelDims[1], (float)voxelDims[2], 0.0f, 0.0f, 0.0f, 0.0f);
        
        return hdr;
        
    }

    
    /**
     * 
     * @return true if tensors are written in upper triangular format
     */
    protected boolean tensorIsUpperTriangular() {
        return false;
    }
    


}

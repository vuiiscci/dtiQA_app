package imaging;

import data.*;
import tools.*;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 *
 * Base class for testing ImageHeader I/O
 *
 * @author  Philip Cook
 *
 *
 */
public abstract class TestImageHeader extends TestCase {
    
    public TestImageHeader(String name) {
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
	return new TestSuite(TestImageHeader.class);
    }


  
    /**
     * Test writing a scalar image and reading it back in with readVolumeData
     *
     */
    public void testWriteScalarImage() {

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

        String fileRoot = "./test/imaging/testImageHeader";

        ImageHeader ih = getImageHeader(fileRoot, dataDims, voxelDims);
        
        
        // ImageHeader methods write without scale terms
        ImageHeader written = ih.writeScalarImage(data, fileRoot);
        
        double[][][][] read = written.readVolumeData();
     
        for (int k = 0; k < data[0][0].length; k++) {
            for (int j = 0; j < data[0].length; j++) {
                for (int i = 0; i < data.length; i++) {
                    assertEquals(data[i][j][k], read[i][j][k][0], 1E-6);
                }
            }
        }  
        
    }


    /**
     * Test vector image write and read via getImageDataSource
     *
     */
    public void testWriteVectorImage() {
        
        int[] dataDims = new int[] {2, 3, 2};

        double[] voxelDims = new double[] {1.0, 2.0, 3.0};
        
        int numComponents = 3;

        double[][][][] data = new double[dataDims[0]][dataDims[1]][dataDims[2]][numComponents];

        for (int n = 0; n < numComponents; n++) {
            for (int k = 0; k < dataDims[2]; k++) {
                for (int j = 0; j < dataDims[1]; j++) {
                    for (int i = 0; i < dataDims[0]; i++) {
                        data[i][j][k][n] = 4.0 * i + 2.0 * j + 8.0 * k + 3.0 * n;
                    }
                }
            }
        }

        String fileRoot = "./test/imaging/testImageHeader";

        ImageHeader ih = getImageHeader(fileRoot, dataDims, voxelDims);
        
        // ImageHeader methods write without scale terms
        ImageHeader written = ih.writeVectorImage(data, fileRoot);
        
        DataSource ds = written.getImageDataSource();
        
        for (int k = 0; k < dataDims[2]; k++) {
            for (int j = 0; j < dataDims[1]; j++) {
                for (int i = 0; i < dataDims[0]; i++) {
                    
                    double[] voxel = ds.nextVoxel();
                    
                    for (int n = 0; n < numComponents; n++) {
                        assertEquals(data[i][j][k][n], voxel[n], 1E-6);
                    }
                }
            }
        }  
        
    }
 

    /**
     * Test writing of a tensor image
     *
     */
    public void testWriteTensorImage() {
        
        int[] dataDims = new int[] {2, 3, 2};

        double[] voxelDims = new double[] {1.0, 2.0, 3.0};

        int numComponents = 6;
        
        double[][][][] data = new double[dataDims[0]][dataDims[1]][dataDims[2]][numComponents];

        for (int n = 0; n < numComponents; n++) {
            for (int k = 0; k < dataDims[2]; k++) {
                for (int j = 0; j < dataDims[1]; j++) {
                    for (int i = 0; i < dataDims[0]; i++) {
                        data[i][j][k][n] = 4.0 * i + 2.0 * j + 8.0 * k + 3.0 * n;
                    }
                }
            }
        }

        String fileRoot = "./test/imaging/testImageHeader";

        ImageHeader ih = getImageHeader(fileRoot, dataDims, voxelDims);
       
        ImageHeader written = ih.writeTensorImage(data, fileRoot);

        double[][][][] tensors = written.readVolumeData();
        
        for (int k = 0; k < dataDims[2]; k++) {
            for (int j = 0; j < dataDims[1]; j++) {
                for (int i = 0; i < dataDims[0]; i++) {
                    
                    assertEquals(data[i][j][k][0], tensors[i][j][k][0], 1E-6);
                    assertEquals(data[i][j][k][1], tensors[i][j][k][1], 1E-6);

                    if (tensorIsUpperTriangular()) {
                        assertEquals(data[i][j][k][2], tensors[i][j][k][2], 1E-6);
                        assertEquals(data[i][j][k][3], tensors[i][j][k][3], 1E-6);
                    }
                    else {
                        assertEquals(data[i][j][k][2], tensors[i][j][k][3], 1E-6);
                        assertEquals(data[i][j][k][3], tensors[i][j][k][2], 1E-6);
                    }

                    assertEquals(data[i][j][k][4], tensors[i][j][k][4], 1E-6);
                    assertEquals(data[i][j][k][5], tensors[i][j][k][5], 1E-6);

                }
            }
        }  
        
    }
 

    


    /**
     * Gets a header template for testing. Header should have a floating point data type
     *
     * @param fileRoot the root of the image file. The extension varies by image format. Calling the write methods
     * on the header with the same root should overwrite the original file - this is needed to test data scaling
     * in the header.
     *
     * @param dataDims dimensions of the image
     *
     * @param voxelDims dimensions of the voxels
     *
     *
     *
     */
    protected abstract ImageHeader getImageHeader(String fileRoot, int[] dataDims, double[] voxelDims);

    
    /**
     * 
     * @return true if tensors are written in upper triangular format, false otherwise
     */
    protected abstract boolean tensorIsUpperTriangular();
    


    // RGB Image I/O tested via apps.TestRGB_ScalarImage
   
}

    
    

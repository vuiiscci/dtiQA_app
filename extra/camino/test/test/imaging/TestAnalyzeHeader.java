package imaging;

import data.*;
import tools.*;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>AnalyzeHeader</code>.
 * <BR>
 *
 * </dl>
 *
 * @version $Id: TestAnalyzeHeader.java,v 1.4 2005/11/24 12:07:35 ucacpco Exp $
 * @author  Philip Cook
 * @see imaging.AnalyzeHeader;
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestAnalyzeHeader extends TestImageHeader {

    public TestAnalyzeHeader(String name) {
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
	return new TestSuite(TestAnalyzeHeader.class);
    }


    public void testHeaderIO() {
	AnalyzeHeader ah = new AnalyzeHeader();

	ah.datatype = AnalyzeHeader.DT_FLOAT;

	ah.width = 2;
	ah.height = 3;
	ah.depth = 4;
	ah.vox_offset = 0;
	ah.nImages = 1;

	ah.intelByteOrder = false; 
	ah.glmin = 0;
	ah.glmax = 100;
    
	// depth == zVoxelDim
	ah.pixelWidth = 1.0f;
	ah.pixelHeight = 2.0f;
	ah.pixelDepth = 3.0f;
    
	// multiply pixel values by this
	ah.scaleSlope = 2.0f;
	ah.scaleInter = 1.0f;
    
	// bits per pixel
	ah.bitpix = 32;
    
	ah.description = "test header";

	ah.centre[0] = 50;
	ah.centre[1] = 60;
	ah.centre[2] = 70;

	AnalyzeHeader ahRead = null;

	try {

	    ah.setHeaderFile("./test/imaging/testHeader.hdr");
	    ah.writeHeader();
	    
	    ahRead = AnalyzeHeader.readHeader("./test/imaging/testHeader.hdr");
	
	}
	catch (java.io.IOException e) {

	    fail(e.toString());
	}
    
	assertEquals(ah.width, ahRead.width);
	assertEquals(ah.height, ahRead.height);
	assertEquals(ah.depth, ahRead.depth);

	assertEquals(ah.vox_offset, ahRead.vox_offset, 1E-6);
	assertEquals(ah.nImages, ahRead.nImages);

	assertEquals(ah.intelByteOrder, ahRead.intelByteOrder);
	assertEquals(ah.glmin, ahRead.glmin);
	assertEquals(ah.glmax, ahRead.glmax);

	assertEquals(ah.pixelWidth, ahRead.pixelWidth, 1E-6);
	assertEquals(ah.pixelHeight, ahRead.pixelHeight, 1E-6);
	assertEquals(ah.pixelDepth, ahRead.pixelDepth, 1E-6);
    
	assertEquals(ah.scaleSlope, ahRead.scaleSlope, 1E-6);
	assertEquals(ah.scaleInter, ahRead.scaleInter, 1E-6);
	assertEquals(ah.bitpix, ahRead.bitpix);

	assertTrue(ah.description.equals(ahRead.description));

	assertEquals(ah.centre[0], ahRead.centre[0]);
	assertEquals(ah.centre[1], ahRead.centre[1]);
	assertEquals(ah.centre[2], ahRead.centre[2]);
    


    }

    public void testCopyConstruction() {

	AnalyzeHeader ah = new AnalyzeHeader();

	ah.datatype = AnalyzeHeader.DT_FLOAT;

	ah.width = 2;
	ah.height = 3;
	ah.depth = 4;
	ah.vox_offset = 0;
	ah.nImages = 1;

	ah.intelByteOrder = false; 
	ah.glmin = 0;
	ah.glmax = 100;
    
	// depth == zVoxelDim
	ah.pixelWidth = 1.0f;
	ah.pixelHeight = 2.0f;
	ah.pixelDepth = 3.0f;
    
	// multiply pixel values by this
	ah.scaleSlope = 2.0f;
	ah.scaleInter = 1.0f;
    
	// bits per pixel
	ah.bitpix = 32;
    
	ah.description = "test header";

	ah.centre[0] = 50;
	ah.centre[1] = 60;
	ah.centre[2] = 70;

	AnalyzeHeader copyRead = null;

	try {

	    AnalyzeHeader copy = new AnalyzeHeader(ah, "./test/imaging/testCopy");
	    
	    copy.writeHeader();

	    copyRead = AnalyzeHeader.readHeader("./test/imaging/testCopy.hdr");
	    
	}
	catch (java.io.IOException e) {
	    
	    fail(e.toString());
	}
    
	assertEquals(ah.width, copyRead.width);
	assertEquals(ah.height, copyRead.height);
	assertEquals(ah.depth, copyRead.depth);

	assertEquals(ah.vox_offset, copyRead.vox_offset, 1E-6);
	assertEquals(ah.nImages, copyRead.nImages);

	assertEquals(ah.intelByteOrder, copyRead.intelByteOrder);
	assertEquals(ah.glmin, copyRead.glmin);
	assertEquals(ah.glmax, copyRead.glmax);

	assertEquals(ah.pixelWidth, copyRead.pixelWidth, 1E-6);
	assertEquals(ah.pixelHeight, copyRead.pixelHeight, 1E-6);
	assertEquals(ah.pixelDepth, copyRead.pixelDepth, 1E-6);
    
	assertEquals(ah.scaleSlope, copyRead.scaleSlope, 1E-6);
	assertEquals(ah.scaleInter, copyRead.scaleInter, 1E-6);
	assertEquals(ah.bitpix, copyRead.bitpix);

	assertTrue(ah.description.equals(copyRead.description));

	assertEquals(ah.centre[0], copyRead.centre[0]);
	assertEquals(ah.centre[1], copyRead.centre[1]);
	assertEquals(ah.centre[2], copyRead.centre[2]);

	
    }

 
    public void testReadWithOffset() {
	
	// Test whether the voxel offset in the header is correctly applied

	// first test, positive offset (applied once) 

	// data is a 4D volume, each volume is a 3x3x3 cube
	byte[] data4D = new byte[81];

	for (int i = 0; i < data4D.length; i++) {
	    data4D[i] = (byte)i;
	}

	AnalyzeHeader ah = new AnalyzeHeader();
	
	ah.datatype = AnalyzeHeader.DT_UNSIGNED_CHAR;
	ah.bitpix = 8;

	ah.width = 3;	
	ah.height = 3;	
	ah.depth = 3;	

	ah.vox_offset = 2.0f;

	ah.nImages = 3;

	byte[] offset = new byte[2];

	offset[0] = (byte)127;
	offset[1] = (byte)127;

	try {
	    ah.setHeaderFile("./test/imaging/testOffset.hdr");
	    
	    FileOutputStream fid = new FileOutputStream("./test/imaging/testOffset.img");
	    
	    fid.write(offset);

	    fid.write(data4D);

	    fid.close();

	    double[][][][] vol = ah.readVolumeData();

	    assertEquals(3, vol.length);
	    assertEquals(3, vol[0].length);
	    assertEquals(3, vol[0][0].length);
	    assertEquals(3, vol[0][0][0].length);
	    
	    int counter = 0;

	    for (int n = 0; n < 3; n++) {
		for (int k = 0; k < 3; k++) {
		    for (int j = 0; j < 3; j++) {
			for (int i = 0; i < 3; i++) {
			    assertEquals(counter++, vol[i][j][k][n], 1E-6);
			}
		    }
		}  
	    }		
	    
	    
	}
	catch (IOException e) {
	    fail(e.toString());
	}

	try {

	    // second test, negative offset, applied to each volume
	    ah.vox_offset = -2.0f;
	    
	    ah.setHeaderFile("./test/imaging/testNegOffset.hdr");
	    	    
	    FileOutputStream fid = new FileOutputStream("./test/imaging/testNegOffset.img");
	    
	    int counter = 0;
	    
	    for (int v = 0; v < 3; v++) {

		fid.write(offset);

		for (int i = 0; i < 27; i++) {
		    fid.write(data4D[counter++]);
		}
	    }
	    
	    fid.close();

	    double[][][][] vol = ah.readVolumeData();

	    assertEquals(3, vol.length);
	    assertEquals(3, vol[0].length);
	    assertEquals(3, vol[0][0].length);
	    assertEquals(3, vol[0][0][0].length);
	    
	    counter = 0;

	    for (int n = 0; n < 3; n++) {
		for (int k = 0; k < 3; k++) {
		    for (int j = 0; j < 3; j++) {
			for (int i = 0; i < 3; i++) {
			    assertEquals(counter++, vol[i][j][k][n], 1E-6);
			}
		    }
		}  
	    }		
	    

	    // test from data source with 4D neg offset
	    DataSource ds = ah.getImageDataSource();

	    counter = 0;

	    for (int k = 0; k < 3; k++) {
		for (int j = 0; j < 3; j++) {
		    for (int i = 0; i < 3; i++) {
			
			double[] voxel = ds.nextVoxel();

			assertEquals(3, voxel.length);
			
			assertEquals(counter, voxel[0], 1E-6);
			assertEquals(counter + 27, voxel[1], 1E-6);
			assertEquals(counter + 54, voxel[2], 1E-6);
			
			counter++;

		    }
		}
	    }

	}
	catch (IOException e) {
	    fail(e.toString());
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
     * @param scaleSlope slope term of the header scaling, if the image supports scaling (override #supportsScale() to 
     * indicate this).
     *
     * @param scaleInter intercept term of the header scaling, if the image supports scaling (override #supportsScale() to 
     * indicate this).
     *
     *
     */
    protected ImageHeader getImageHeader(String fileRoot, int[] dataDims, double[] voxelDims) {
        
        AnalyzeHeader ah = new AnalyzeHeader(fileRoot, false);

        ah.setDataType("double");

        ah.width = (short)dataDims[0];
        ah.height = (short)dataDims[1];
        ah.depth = (short)dataDims[2];

        ah.pixelWidth = (float)voxelDims[0];
        ah.pixelHeight = (float)voxelDims[1];
        ah.pixelDepth = (float)voxelDims[2];

        return ah;

    }

    
    /**
     * 
     * @return true if tensors are written in upper triangular format, false otherwise
     */
    protected boolean tensorIsUpperTriangular() {
        return false;
    }
    




}

    
    

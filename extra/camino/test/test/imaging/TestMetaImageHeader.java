package imaging;

import data.*;
import numerics.*;
import tools.*;

import junit.framework.*;
import junit.extensions.*;

import java.io.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>MetaImageHeader</code>.
 * <BR>
 *
 * </dl>
 *
 * @version $Id $
 * @author  Philip Cook
 * @see imaging.MetaImageHeader;
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestMetaImageHeader extends TestImageHeader {

    public TestMetaImageHeader(String name) {
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
	return new TestSuite(TestMetaImageHeader.class);
    }


    /**
     * Writes a header to disk and then reads it back in.
     */ 
    public void testHeaderIO() {
	MetaImageHeader mh = new MetaImageHeader();

	mh.transformation = Rotations.randomRotMat(1234);

	mh.offset = new double[] {1.0, 2.0, 3.0};

	mh.rotationCentre = new double[] {4.0, 5.0, 6.0};

	mh.anatomicalOrientation = MetaImageHeader.AnatomicalOrientation.RPI;

	mh.spacing = new double[] {2.0, 2.0, 3.0};

	mh.dimSize = new int[] {128, 128, 60};

	mh.channels = 1;

	mh.dataType = MetaImageHeader.DataType.FLOAT;

	mh.dataFile = "test.mhd";

	MetaImageHeader read = null;

	try {
    
	    mh.writeHeader("test.mha");
	    
	    read = MetaImageHeader.readHeader("test.mha");
	    
	}
	catch (IOException e) {
	    fail(e.toString());
	}
	
	assertHeadersEqual(mh, read);

    }

    

    /**
     * Writes a header using MetaImage, then write an image to the same IO stream, and read the volume
     * back in, first the header, then the image.
     */ 
    public void testWriteImage() {
	
	MetaImageHeader mh = new MetaImageHeader();

	mh.spacing = new double[] {2.0, 2.0, 2.0};

	mh.dimSize = new int[] {2, 2, 2};

	mh.channels = 1;

	mh.dataType = MetaImageHeader.DataType.FLOAT;

	mh.dataFile = "LOCAL";

	try {
	    DataOutputStream dout = 
		new DataOutputStream(new BufferedOutputStream(new FileOutputStream("test2.mhd"), 1024));
	    	    
	    mh.writeHeader(dout);

	    for (int i = 0; i < 8; i++) {
		dout.writeFloat((float)i);
	    }

	    dout.close();

	    DataInputStream din = 
		new DataInputStream(new BufferedInputStream(new FileInputStream("test2.mhd"), 1024));

	    
	    MetaImageHeader read = MetaImageHeader.readHeader(din);

	    assertHeadersEqual(mh, read);

	    for (int i = 0; i < 8; i++) {
		assertEquals((float)i, din.readFloat(), 1E-8);
	    }

	    din.close();
	    

	}
	catch (IOException e) {
	    fail(e.toString());
	}

    }



    /** 
     * Do a one-pass reading of a volume, both from LOCAL data and a separate .mha / .mhd combination.
     *
     */
    public void testReadVolumeData() {

	MetaImageHeader mh = new MetaImageHeader();

	mh.spacing = new double[] {2.0, 2.0, 2.0};

	mh.dimSize = new int[] {2, 2, 2};

	mh.channels = 1;

	mh.dataType = MetaImageHeader.DataType.FLOAT;

	mh.dataFile = "LOCAL";

	try {
	    DataOutputStream dout = 
		new DataOutputStream(new BufferedOutputStream(new FileOutputStream("test3.mhd"), 1024));
	    	    
	    mh.writeHeader(dout);

	    for (int i = 0; i < 8; i++) {
		dout.writeFloat((float)i);
	    }

	    dout.close();

	    MetaImageHeader min = MetaImageHeader.readHeader("test3.mhd");

	    double[][][][] vol = min.readVolumeData();

	    for (int i = 0; i < 8; i++) {
		assertEquals((float)i, vol[i % 2][(i / 2) % 2][i / 4][0], 1E-8);
	    }

	    // try separate hdr / img pair
   
	    mh.dataFile = "test4.raw";

	    mh.writeHeader("test4.mha");

	    dout = 
		new DataOutputStream(new BufferedOutputStream(new FileOutputStream("test4.raw"), 1024));

	    for (int i = 0; i < 8; i++) {
		dout.writeFloat((float)i);
	    }

	    dout.close();

	    min = MetaImageHeader.readHeader("test4.mha");
	    vol = min.readVolumeData();
	    
	    for (int i = 0; i < 8; i++) {
		assertEquals((float)i, vol[i % 2][(i / 2) % 2][i / 4][0], 1E-8);
	    }


	    // try reading some little-endian data
	    mh.dataFile = "LOCAL";
	    
	    mh.intelByteOrder = true;

	    LEFilterOutputStream leOut = 
		new LEFilterOutputStream(new BufferedOutputStream(new FileOutputStream("test5.mhd"), 1024));

	    mh.writeHeader(leOut);

	    
	    for (int i = 0; i < 8; i++) {
		leOut.writeFloat((float)i);
	    }

	    leOut.close();
	    
	    min = MetaImageHeader.readHeader("test5.mhd");
	    vol = min.readVolumeData();
	    
	    for (int i = 0; i < 8; i++) {
		assertEquals((float)i, vol[i % 2][(i / 2) % 2][i / 4][0], 1E-8);
	    }

	}
	catch (IOException e) {
	    fail(e.toString());
	}

    }
    

    public void testGetImageDataSource() {

	MetaImageHeader mh = new MetaImageHeader();

	mh.spacing = new double[] {2.0, 2.0, 2.0};

	mh.dimSize = new int[] {2, 2, 2};

	mh.channels = 1;

	mh.dataType = MetaImageHeader.DataType.FLOAT;

	mh.dataFile = "LOCAL";

	mh.headerFile = "test6.mhd";
	mh.localDataFile = "test6.mhd";

	try {
	    DataOutputStream dout = 
		new DataOutputStream(new BufferedOutputStream(new FileOutputStream("test6.mhd"), 1024));
	    	    
	    mh.writeHeader(dout);

	    for (int i = 0; i < 8; i++) {
		dout.writeFloat((float)i);
	    }

	    dout.close();

	    DataSource ds = mh.getImageDataSource();

	    for (int i = 0; i < 8; i++) {
		assertEquals((float)i, (float)ds.nextVoxel()[0], 1E-6);
	    }

	}
	catch (IOException e) {
	    fail(e.toString());
	}

    }


    public void testReadVolume() {
	MetaImageHeader mh = new MetaImageHeader();
	
	mh.spacing = new double[] {2.0, 2.0, 2.0};
	
	mh.dimSize = new int[] {2, 2, 2};
	
	mh.channels = 2;
	
	mh.dataType = MetaImageHeader.DataType.FLOAT;
	
	mh.dataFile = "LOCAL";

	mh.localDataFile = "test7.mhd";

	try {
	    DataOutputStream dout = 
		new DataOutputStream(new BufferedOutputStream(new FileOutputStream("test7.mhd"), 1024));
	    	    
	    mh.writeHeader(dout);

	    for (int i = 0; i < 8; i++) {
		dout.writeFloat((float)i);
		dout.writeFloat((float)(2 * i));
	    }

	    dout.close();

	    mh = MetaImageHeader.readHeader("test7.mhd");

	    double[][][][] c = mh.readVolumeData();

	    for (int i = 0; i < 8; i++) {
		assertEquals((float)i, c[i % 2][(i / 2) % 2][i / 4][0], 1E-8);
		assertEquals((float)(2*i), c[i % 2][(i / 2) % 2][i / 4][1], 1E-8);
	    }
	    
	}
	catch (IOException e) {
	    fail(e.toString());
	}

	
    }




    private void assertHeadersEqual(MetaImageHeader mh, MetaImageHeader read) { 

	assertTrue(mh.intelByteOrder == read.intelByteOrder);
	
	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {
		assertEquals(mh.transformation.entries[i][j], read.transformation.entries[i][j], 1E-8); 
	    }
	}
	
	for (int i = 0; i < 3; i++) {
	    assertEquals(mh.offset[i], read.offset[i], 1E-8);
	    assertEquals(mh.rotationCentre[i], read.rotationCentre[i], 1E-8);
	    assertEquals(mh.spacing[i], read.spacing[i], 1E-8);
	    assertEquals(mh.dimSize[i], read.dimSize[i]);
	}
	    
	assertTrue(mh.anatomicalOrientation == read.anatomicalOrientation);
	
	assertTrue(mh.channels == read.channels);
	
	assertTrue(mh.dataType == read.dataType);
	    
	assertTrue(mh.dataFile.equals(read.dataFile));
    }



    


    
    /**
     * Gets a header template for testing.
     *
     * @return a header with gzip compression, and the required file root, dimensions and scale parameters.
     *
     */
    protected ImageHeader getImageHeader(String fileRoot, int[] dataDims, double[] voxelDims) {

      	MetaImageHeader mh = new MetaImageHeader();

	mh.spacing = voxelDims;
        
	mh.dimSize = dataDims;

	mh.channels = 1;

	mh.dataType = MetaImageHeader.DataType.FLOAT;

	mh.dataFile = "LOCAL";

        mh.headerFile = fileRoot + ".mhd";

        mh.localDataFile = mh.headerFile;

        return mh;
        
    }

    
    /**
     * 
     * @return true if tensors are written in upper triangular format
     */
    protected boolean tensorIsUpperTriangular() {
        return true;
    }


}

package data;

import imaging.*;
import misc.LoggedException;
import tools.*;

import java.io.*;
import java.util.logging.Logger;
import java.util.zip.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Base class for data sources of data from files and input streams.
 * 
 * <dt>Description:
 * 
 * <dd>Contains some basic methods for controlling data types.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 */
public abstract class ExternalDataSource implements DataSource {

    private static final Logger logger = Logger.getLogger("data.ExternalDataSource");
	
    /**
     * Size of read in file buffer.
     */
    public static int FILEBUFFERSIZE = 32000000;

    /**
     * The number of values in each voxel.
     */
    protected int numComponents;

    /**
     * The data type.
     */
    protected int datatype;

    /**
     * Codes for different possible data types.
     */
    public static final int BYTE = 0;

    public static final int SHORT = 1;

    public static final int INT = 2;

    public static final int LONG = 3;

    public static final int FLOAT = 4;

    public static final int DOUBLE = 5;

    // unsigned 8 bit char
    public static final int UBYTE = 6;

    // SPM2 additional types
    public static final int USHORT = 7;
    public static final int UINT = 8;

    /**
     * The stream to read data from.
     */
    protected EndianNeutralDataInputStream dataIn;

    /**
     * Default constructor required by inherited classes.
     */
    protected ExternalDataSource() {
    }

    /**
     * Initializes the input stream from a filename. If filename is null, the
     * input stream in the standard input.
     * 
     * @param filename
     *            The name of the data file.
     *
     * @param intelByteOrder true if the input stream is little-endian.
     *
     * @param offset number of bytes into the data file to skip before reading data.
     */
    protected void initFileInput(String filename, boolean intelByteOrder, int offset) {

        if (filename == null) {
        	logger.info("reading data from standard input");
            dataIn = 
		new EndianNeutralDataInputStream(new BufferedInputStream(System.in, FILEBUFFERSIZE), intelByteOrder);
        }
        else {
	    
	    if (filename.endsWith(".gz") || filename.endsWith(".zip")) {
		initCompressedFileInput(filename, intelByteOrder, offset);
		return;
	    }

	    
            try {
		
                // Open the file for reading. This can fail if the file
                // is inaccessible.
                FileInputStream fileIn = new FileInputStream(filename);
                dataIn = 
		    new EndianNeutralDataInputStream(new BufferedInputStream(fileIn, FILEBUFFERSIZE), intelByteOrder);		

		int bytesSkipped = 0;
		
		while (bytesSkipped < offset) {
		    bytesSkipped += dataIn.skipBytes(offset - bytesSkipped);
		}
		
            }
            catch (Exception e) {
                throw new LoggedException(e);
            }


	}

    }

    
    /**
     * Initializes the input stream from a gzipped file.
     * 
     * @param filename
     *            The name of the data file.
     *
     * @param intelByteOrder true if the input stream is little-endian.
     *
     * @param offset number of bytes into the data file to skip before reading data.
     */
    protected void initCompressedFileInput(String filename, boolean intelByteOrder, int offset) {

        if (filename == null) {
            throw new LoggedException("Compressed input on stdin is not supported. " + 
				      "Try piping output from gunzip -c <file> or specify the input file.");
        }
        else {
            try {
                
                // Open the file for reading. This can fail if the file
                // is inaccessible.
                FileInputStream fileIn = new FileInputStream(filename);
		
                
		if (filename.endsWith(".gz")) {
                    
                    InflaterInputStream in = 
                        new GZIPInputStream(fileIn, FILEBUFFERSIZE/2);
                        
                    // Buffer some more on the decompressed side for speed
                    dataIn = new EndianNeutralDataInputStream
                        (new BufferedInflaterInputStream(in, FILEBUFFERSIZE/2), intelByteOrder);
                    
		}
		else if (filename.endsWith(".zip")) {
                    
		    ZipInputStream in = new ZipInputStream(new BufferedInputStream(fileIn, FILEBUFFERSIZE/2));
		    in.getNextEntry();
		    dataIn = new EndianNeutralDataInputStream
			(new BufferedInflaterInputStream(in, FILEBUFFERSIZE/2), intelByteOrder);
                    
		}
		
		
		int bytesSkipped = 0;
		
		while (bytesSkipped < offset) {
		    bytesSkipped += dataIn.skipBytes(offset - bytesSkipped);
		}
		
	    }
            catch (Exception e) {
                throw new RuntimeException(e);
            }
	
	}
    }


    /**
     * Converts a String indicating a primitive data type to an integer code for
     * that data type..
     * 
     * @param String
     *            The data type as a String: "byte", "float", etc.
     * 
     * @return The int code for the data type.
     */
    public static int getDataTypeCode(String s) throws DataSourceException {
        if (s.equals("byte")) {
            return BYTE;
        }
        if (s.equals("ubyte")) {
            return UBYTE;
        }
        if (s.equals("short")) {
            return SHORT;
        }
        if (s.equals("int")) {
            return INT;
        }
        if (s.equals("ushort")) {
            return USHORT;
        }
        if (s.equals("uint")) {
            return UINT;
        }
        if (s.equals("long")) {
            return LONG;
        }
        if (s.equals("float")) {
            return FLOAT;
        }
        if (s.equals("double")) {
            return DOUBLE;
        }

        throw new DataSourceException("Unrecognized data type String: " + s);
    }


    /**
     * Converts the data type integer code to a string
     * 
     * @param typeCode integer representing the data type.
     * 
     * @return String indicating the data type, eg "double", "byte", etc..
     */
    public static String typeString(int typeCode) {
        if (typeCode == BYTE) {
            return "byte";
        }
        if (typeCode == UBYTE) {
            return "ubyte";
        }
        if (typeCode == SHORT) {
            return "short";
        }
        if (typeCode == INT) {
            return "int";
        }
        if (typeCode == USHORT) {
            return "ushort";
        }
        if (typeCode == UINT) {
            return "uint";
        }
        if (typeCode == LONG) {
            return "long";
        }
        if (typeCode == FLOAT) {
            return "float";
        }
        if (typeCode == DOUBLE) {
            return "double";
        }

        return "unknown type";
    }


  
    /**
     * Gets a data source from some file. If the file is a type known to ImageHeader, then ImageHeader methods are called
     * to initialize the data source. Otherwise the file is read as raw binary. Raw data must be in voxel order. Image files
     * can be in whatever order is required by the image format.
     * 
     *
     * @param filename the name of the file to read. If this is a supported image (eg .nii.gz), then the remaining arguments
     * are used for checking the correctness of the input file, see below.
     *
     * @param components number of components per voxel. Necessary for raw data. For image data, set to 0 for an unspecified 
     * number of components, or the required number of components.
     *
     * @param type the type of the raw data, such as "double". Ignored for image data.
     *
     * @see imaging.ImageHeader
     *
     */
    public static DataSource getDataSource(String filename, int components, String type) {
        
        if (ImageHeader.imageExists(filename)) {
            
            try {
                ImageHeader header = ImageHeader.readHeader(filename); 

                //
                // This is kind of a hack because it assumes DT data. Tricky to avoid though without getting into
                // format specific issues like checking intent codes, so we bravely rely on users to supply correct data
                // and not pass a random 6-component data set instead of an 8-component dataset that is not a DT
                //
                if (header.components() == 6 && components == 8) {
                    logger.info("Importing diffusion tensors from file: " + filename);
                    return new DT_DataSource(header.getImageDataSource(), header.lowerTriangularSymmMatrix());
                }
                if (components != header.components()) {
                    // With the exception of allowing diffusion tensors without the extra bits, enforce correct number
                    // of components
                    throw new LoggedException("Expected " + components + " components but file has " + header.components());
                }

                return header.getImageDataSource();
            }
            catch (IOException e) {
                throw new LoggedException(e);
            }
            
        }
        else {
            return new VoxelOrderDataSource(filename, components, type);
        }

    }



}

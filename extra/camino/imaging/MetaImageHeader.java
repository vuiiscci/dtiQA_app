package imaging;

import data.*;
import misc.*;
import numerics.*;
import tools.*;

import java.io.*;
import java.util.logging.Logger;
import java.util.zip.*;

/**
 * Meta image header, provides support for a common ITK file format. 
 * Currently only supports 3D binary images (images may have multiple components).
 * <p>
 * Compressed data is read as a ZIP archive, assumed to contain one entry. This feature
 * is untested, so use with caution.
 * <p>
 *
 * 
 * @author Philip Cook
 * @version $Id$
 */
public final class MetaImageHeader extends ImageHeader {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.imaging.MetaImageHeader");


    /** ObjectType, which for us is always Image. */
    protected final ObjectType objectType = ObjectType.IMAGE;

    /** NDims, number of dimensions. Always 3. */
    protected final int nDims = 3;
    
    /** BinaryData, which is always true. */
    protected final boolean binaryData = true;
    
    /** CompressedData, if true, data is compressed with zlib. */
    protected boolean compressedData = false;



    /** 
     * BinaryDataByteOrderMSB = !<code>intelByteOrder</code>. 
     * Field named for consistency with AnalyzeHeader. The Meta specification also mentions
     * the seemingly redundant ElementByteOrderMSB. ITK examples contain either field, but 
     * never both. Currently we read either but only write BinaryDataByteOrderMSB.
     * 
     */
    protected boolean intelByteOrder = false;

    /** TransformMatrix, a 3x3 rotation matrix. */
    protected RealMatrix transformation = RealMatrix.identity(3);
    
    /** Offset */
    protected double[] offset = new double[] {0.0, 0.0, 0.0};

    /** CenterOfRotation */
    protected double[] rotationCentre = new double[] {0.0, 0.0, 0.0};


    
    /**
     * AnatomicalOrientation field. From the Meta IO 
     * <a href="http://www.itk.org/Wiki/images/2/27/MetaIO-Introduction.pdf">spec</a>:
     * <blockquote>
     * "Specify anatomic ordering of the axis.  Use only [R|L] | 
     * [A|P] | [S|I] per axis.   For example, if the three letter code for (column index, 
     * row index, slice index is) ILP, then the origin is at the superior, right, anterior 
     * corner of the volume, and therefore the axes run from superior to inferior, from 
     * right to left, from anterior to posterior."
     * </blockquote>
     * <p>
     * Camino reads files in the order of left-right, posterior-anterior, inferior-superior.
     * The origin of the image (voxel 0,0,0) is left, posterior, inferior. The orientation 
     * could therefore be described as LPI, but according to the above it's RAS. Because SNAP
     * expects the codes to match the position of the origin, we'll default to LPI.
     *
     */
    protected AnatomicalOrientation anatomicalOrientation = AnatomicalOrientation.LPI;


    /** Voxel dimensions. */
    protected double[] spacing = new double[] {1.0, 1.0, 1.0};

    /** Number of voxels in each dimension. */
    protected int[] dimSize = new int[] {1, 1, 1};


    /**
     * ElementNumberOfChannels, number of components.
     */
    protected int channels = 1;


    /**
     * ElementType, ie datatype. 
     */ 
    protected DataType dataType = DataType.DOUBLE;


    /**
     * Header file name, includes all path information.
     */
    protected String headerFile = "";


    /**
     * Data file. Should be "LOCAL" if the data is in the same file as the header.
     * Note that the variable localDataFile, which should be an absolute path, will be used
     * to read the image within Camino. This file, if not "LOCAL", is specified relative to the 
     * header (usually in the same directory).
     *
     */
    protected String dataFile = "LOCAL";


    /**
     * This variable contains the absolute path to the file containing the data (which may be the
     * same as the file containing the header). It is set when readHeader is called. 
     *
     * This is the actual file that gets used to read data.
     */
    protected String localDataFile = "";


    /** 
     * Offset of data in elementDataFile. Set when readHeader is called. 
     */
    protected int dataByteOffset = 0;


    
    /**
     * Data type enum. We use a subset of allowed Meta types that correspond to Java data types.
     * UCHAR is an 8 bit unsigned byte, CHAR is an 8 bit signed byte, ie a Java <code>byte</code>. 
     * The SHORT, INT, FLOAT, LONG and DOUBLE types correspond to the Java types of the respective name.
     */
    public enum DataType {
	
	// numeric values defined for consistency with Analyze
	UCHAR("MET_UCHAR", "ubyte", 2),
            CHAR("MET_CHAR", "byte", 132),
            SHORT("MET_SHORT", "short", 4),
            USHORT("MET_USHORT", "ushort", 130),
            INT("MET_INT", "int", 8),
            UINT("MET_UINT", "uint", 136),
            FLOAT("MET_FLOAT", "float", 16),
            DOUBLE("MET_DOUBLE", "double", 64),
            LONG("MET_LONG", "long", 199);


        DataType(String name, String caminoName, int i) {
	    typeName = name;
            caminoTypeName = caminoName;
	    index = i;
	}


	public String toString() {
	    return typeName;
	}


        /**
         * Gets the data type given a string representation of its name, or its Camino equivalent
         *
         */
	public static DataType getDataType(String s) {

	    for (DataType type : DataType.values()) {
		if (s.equals(type.typeName) || s.equals(type.caminoTypeName)) {
		    return type;
		}
	    }
	    
	    // use of enum should make this impossible
	    throw new LoggedException("Unsupported Meta data type " + s);
	}


	public String getCaminoDataType() {
            return caminoTypeName;
        }
	    
	private final String typeName;
        private final String caminoTypeName;
	private final int index;

    }



    /**
     * Object type enum. Only "Image" is allowed.
     */
    public enum ObjectType {

	// only support image for now, others exist
	IMAGE("Image");
	
	// note that enum.name is set automatically to the instance name, eg IMAGE, 
	// and cannot be modified.
        ObjectType(String s) {
	    typeName = s;
	}

	public String toString() {
	    return typeName;
	}


	public static ObjectType getObjectType(String s) {
	    
	    for (ObjectType type : ObjectType.values()) {
		if (s.equals(type.typeName)) {
		    return type;
		}
	    }

	    // use of enum should make this impossible
	    throw new LoggedException("Unsupported Meta object type: '" + s + "'");
	    
	}	

	private final String typeName;

    }


    /**
     * Supported anatomical orientations [L|R] [A|P] [I|S].
     *
     */
    public enum AnatomicalOrientation {
	// default string value is the variable name
	LAS,
	LAI,
	LPS,
	LPI,
	RAS,
	RAI,
	RPS,
	RPI;

	public static AnatomicalOrientation getAnatomicalOrientation(String s) {
	   
	    for (AnatomicalOrientation orient : AnatomicalOrientation.values()) {
		if (s.equals(orient.name())) {
		    return orient;
		}
	    }
	   
	    // use of enum should make this impossible
	    throw new LoggedException("Unsupported Anatomical Orientation " + s);
	   
	}
       
    }

    


    /**
     * Sets all values to their default.
     *
     */
    public MetaImageHeader() {

	
    }


    /**
     * Construct a header that is a copy of an existing header, with a new file name.
     *
     * @param mh an existing header from which to copy information and settings.
     * @param fileRoot a file root including path information if any. The extension of the 
     * header / data file(s) will be determined by the settings in <code>mh</code>.
     *
     */
    protected MetaImageHeader(MetaImageHeader mh, String fileRoot) {
        
    
        compressedData = mh.compressedData;
        intelByteOrder = mh.intelByteOrder;
        
        transformation = (RealMatrix)mh.transformation.clone();
        
        offset = new double[3];
        rotationCentre = new double[3];
        spacing = new double[3];
        dimSize = new int[3];

        for (int i = 0; i < 3; i++) {
            offset[i] = mh.offset[i];
            rotationCentre[i] = mh.rotationCentre[i];
            spacing[i] = mh.spacing[i];
            dimSize[i] = mh.dimSize[i];
        }

        anatomicalOrientation = mh.anatomicalOrientation;
        channels = mh.channels;
        dataType = mh.dataType;

        headerFile = fileRoot + ".mhd";

        if (mh.dataFile.equals("LOCAL")) {
            dataFile = "LOCAL";
            localDataFile = headerFile;
        }
        else {

	    // no path information
	    String fileRootNoPath = "";

            String prefix = "";

            // Cygwin uses / but Windows Java will report \ as separator. Thus we always use /
	    String slashie = "/";
	    
	    int index = fileRoot.lastIndexOf(slashie);
	    
	    if (index > -1) {
		fileRootNoPath = fileRoot.substring(index+1);
	    
                mh.localDataFile = fileRoot + ".raw";

                dataFile = fileRootNoPath + ".raw";
            }
            else {
                dataFile = fileRoot + ".raw";  
                localDataFile = fileRoot + ".raw";  
            }
        }
        
        // dataByteOffset set when reading / writing a header, it will change 
        // whenever anything is changed in the header, so don't set it here

    }


    public void setDataDims(int[] dims) {
        dimSize[0] = dims[0];
        dimSize[1] = dims[1];
        dimSize[2] = dims[2];
    }


    public void setVoxelDims(double[] dims) {
        spacing[0] = dims[0];
        spacing[1] = dims[1];
        spacing[2] = dims[2];
    }
    

    /**
     * Read the text header from a file. Adds either ".mha" or ".mhd" if omitted.
     *
     */
    public static MetaImageHeader readHeader(String file) throws IOException {

	if (!( file.endsWith(".mha") || file.endsWith(".mhd") )) {

	    // look for file, header first
	    File f = new File(file + ".mha");

	    if (f.exists()) {
		file += ".mha";
	    }
	    else {
		f = new File(file + ".mhd");
		
		if (f.exists()) {
		    file += ".mhd";
		}	    
		else {
		    throw new FileNotFoundException("Can't find file " + file + 
						    ".mha or " + file + ".mhd");
		}
	    }
	}

	DataInputStream din = 
	    new DataInputStream(new BufferedInputStream(new FileInputStream(file), 1024));
	
	MetaImageHeader mh = readHeader(din);

	din.close();

        mh.headerFile = file;

	if (mh.dataFile.equals("LOCAL")) {
	    mh.localDataFile = file;
	}
	else {
	    // attach path to file
	    String prefix = "";

            // Cygwin uses / but Windows Java will report \ as separator. Thus we always use /
	    String slashie = "/";
	    
	    int index = file.lastIndexOf(slashie);
	    
	    if (index > -1) {
		prefix = file.substring(0, index+1);
	    }

	    mh.localDataFile = prefix + mh.dataFile;
	}

	return mh;
    }
    
    
    /**
     * Read the text header from an input stream. 
     *
     */
    public static MetaImageHeader readHeader(DataInput din) throws IOException {

	MetaImageHeader mh = new MetaImageHeader();

	boolean gotDataFile = false;

	//  Unix newline, ie one character
	byte newLine = new String("\n").getBytes()[0];


	// Windows uses \r\n for a new line
	byte cr = new String("\r").getBytes()[0];


	// count how many bytes are in the header
	int hdrBytes = 0;
	

	while (!gotDataFile) {

	    String headerKey = "";
	    String headerValue = "";

	    byte[] lineBytes = new byte[2048];
	    
	    // bytes read from this line
	    int bytesRead = 0;

	    readLine: 
	    while (true) { // EOFException thrown at end of file, which we should not hit
		
		byte b = din.readByte();
		
		lineBytes[bytesRead++] = b; 
		
		if (b == newLine) {
		    hdrBytes += bytesRead;
		    break readLine;
		}
	    }

	    int lineLength = bytesRead - 1; 

 	    // remove \r from line written under Windows
 	    if (lineBytes[bytesRead - 1] == cr) {
 		lineLength -= 1;
 	    }

	    String line = new String(lineBytes, 0, lineLength, "US-ASCII");

	    // split on [spaces]=[spaces] where there might be zero spaces
	    // this actually splits on = and gets rid of spaces
	    // leading and trailing spaces removed first by the call to trim
	    String[] headerKeyValue = line.trim().split("\\s*=\\s*");

	    if (headerKeyValue.length != 2) {
		throw new LoggedException("Header line\n" + line + "\nis not a valid Meta header " +
					  " field. Meta fields are of the form \"Key = Value\"");
	    }

	    mh.setHeaderValue(headerKeyValue[0].trim(), headerKeyValue[1].trim());

	    // ElementDataFile is always last
	    if (headerKeyValue[0].equals("ElementDataFile")) {
		gotDataFile = true;
	    }
	}

	// if data file is local, set dataByteOffset
	if (mh.dataFile.equals("LOCAL")) {
	    mh.dataByteOffset = hdrBytes;
	}

	return mh;
	    
    }

   
    /**
     * Parses header values from their String representation. Warns, but does not throw an Exception,
     * if there is an unrecognized key.
     *
     */
    public void setHeaderValue(String key, String value) {

	if (key.equalsIgnoreCase("ObjectType")) {
	    ObjectType type = ObjectType.getObjectType(value);
	    
	    if (type != ObjectType.IMAGE) {
		throw new LoggedException("Only image type is supported");
	    }

	}
	else if (key.equalsIgnoreCase("NDims")) {
	    int dims = Integer.parseInt(value);
	    
	    if (dims != 3) {
		throw new LoggedException("only 3D images are supported.");
	    }
	}
	else if (key.equalsIgnoreCase("BinaryData")) {
	    boolean binary = Boolean.parseBoolean(value);

	    if (!binary) {
		throw new LoggedException("ASCII data is not supported.");
	    }
	}
	else if (key.equalsIgnoreCase("CompressedData")) {
	    compressedData = Boolean.parseBoolean(value);
	}
	else if (key.equalsIgnoreCase("BinaryDataByteOrderMSB")) {
	    intelByteOrder = !(Boolean.parseBoolean(value));
	}
	else if (key.equalsIgnoreCase("ElementByteOrderMSB")) {
	    // it is unclear what the difference is between ElementByteOrderMSB and BinaryDataByteOrderMSB
	    // if both are specified we go with the later one
	    intelByteOrder = !(Boolean.parseBoolean(value));
	}
	else if (key.equalsIgnoreCase("TransformMatrix")) {
	    String[] values = value.split("\\s+");

	    transformation = new RealMatrix(3,3);

	    for (int i = 0; i < 9; i++) {
		transformation.entries[i / 3][i % 3] = Double.parseDouble(values[i]);
	    }
	}
	else if (key.equalsIgnoreCase("Offset")) {

	    String[] values = value.split("\\s+");

	    offset = new double[3];

	    for (int i = 0; i < 3; i++) {
		offset[i] = Double.parseDouble(values[i]);
	    }
	}
	else if (key.equalsIgnoreCase("Position")) {

	    String[] values = value.split("\\s+");

	    offset = new double[3];

	    for (int i = 0; i < 3; i++) {
		offset[i] = Double.parseDouble(values[i]);
	    }
	}
	else if (key.equalsIgnoreCase("Origin")) {

	    String[] values = value.split("\\s+");

	    offset = new double[3];

	    for (int i = 0; i < 3; i++) {
		offset[i] = Double.parseDouble(values[i]);
	    }
	}
	else if (key.equalsIgnoreCase("CenterOfRotation")) {

	    String[] values = value.split("\\s+");

	    rotationCentre = new double[3];

	    for (int i = 0; i < 3; i++) {
		rotationCentre[i] = Double.parseDouble(values[i]);
	    }
	}
	else if (key.equalsIgnoreCase("AnatomicalOrientation")) {
	    anatomicalOrientation = AnatomicalOrientation.getAnatomicalOrientation(value);
	}
	else if (key.equalsIgnoreCase("ElementSpacing")) {
	    
	    String[] values = value.split("\\s+");

	    spacing = new double[3];

	    for (int i = 0; i < 3; i++) {
		spacing[i] = Double.parseDouble(values[i]);
	    }
	    
	}
	else if (key.equalsIgnoreCase("DimSize")) {
	    
	    String[] values = value.split("\\s+");

	    dimSize = new int[3];

	    for (int i = 0; i < 3; i++) {
		dimSize[i] = Integer.parseInt(values[i]);
	    }
	    
	}
	else if (key.equalsIgnoreCase("ElementNumberOfChannels")) {
	    channels = Integer.parseInt(value);
	}
	else if (key.equalsIgnoreCase("ElementType")) {
	    dataType = DataType.getDataType(value);
	}
	else if (key.equalsIgnoreCase("ElementDataFile")) {
	    dataFile = value;
	}
	else {
	    logger.warning("Unrecognized header field " + key + " ignored");
	}

    }
	


    /**
     * Writes header to a file. Expects either a .mha or .mhd extension, which is
     * added if omitted. The Meta convention is to store everything in a .mhd file,
     * or use two files: header.mha and data.mhd.
     *
     */
    public void writeHeader(String file) throws IOException {

	if (!file.endsWith(".mha") && !file.endsWith(".mhd")) {
	    
	    if (dataFile.equalsIgnoreCase("LOCAL")) {
		file = file + ".mhd";
	    }
	    else {
		file = file + ".mha";
	    }
	    
	}

	DataOutputStream dout = 
	    new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file), 1024));

	writeHeader(dout);

	dout.close();
    }


    public void writeHeader(OutputStream out) throws IOException {
        
        int byteCounter = 0;

        byte[][] hdrBytes = new byte[14][];

        hdrBytes[0] = new String("ObjectType = " + objectType + "\n").getBytes("US-ASCII");

	hdrBytes[1] = new String("NDims = " + nDims + "\n").getBytes("US-ASCII");

	hdrBytes[2] = new String("BinaryData = " + titleCaseBool(binaryData) +"\n").getBytes("US-ASCII");
	
	hdrBytes[3] = new String("CompressedData = " + titleCaseBool(compressedData) + "\n").getBytes("US-ASCII");
	
	hdrBytes[4] = new String("BinaryDataByteOrderMSB = " + titleCaseBool(!intelByteOrder) + "\n").getBytes("US-ASCII");

	hdrBytes[5] = new String("TransformMatrix = " + flatMatrixString(transformation) + "\n").getBytes("US-ASCII");

	hdrBytes[6] = new String("Offset = " + offset[0] + " " + offset[1] + " " + offset[2] + "\n").getBytes("US-ASCII");
	
	hdrBytes[7] = new String("CenterOfRotation = " + rotationCentre[0] + " " + rotationCentre[1] + " " +
			      rotationCentre[2] + "\n").getBytes("US-ASCII");
	
	hdrBytes[8] = new String("AnatomicalOrientation = " + anatomicalOrientation + "\n").getBytes("US-ASCII");
	
	hdrBytes[9] = 
            new String("ElementSpacing = " + spacing[0] + " " + spacing[1] + " " + spacing[2] + "\n").getBytes("US-ASCII");
	
	hdrBytes[10] = 
            new String("DimSize = " + dimSize[0] + " " + dimSize[1] + " " + dimSize[2] + "\n").getBytes("US-ASCII");
	
	hdrBytes[11] = new String("ElementNumberOfChannels = " + channels + "\n").getBytes("US-ASCII");

	hdrBytes[12] = new String("ElementType = " + dataType + "\n").getBytes("US-ASCII");

	hdrBytes[13] = new String("ElementDataFile = " + dataFile + "\n").getBytes("US-ASCII");
	
        for (int i = 0; i < 14; i++) {
            out.write(hdrBytes[i]);

            byteCounter += hdrBytes[i].length;
        }

        if (dataFile.equals("LOCAL")) {
            dataByteOffset = byteCounter;
        }

    }

    

    /**
     * Writes the header and data to disk. 
     *
     * @param data is assumed to have the right dimensions to match the header. The output data type
     * is determined by the header.
     */
    private void writeImage(double[][][][] data) throws IOException {

        FileOutputStream hdrFos = new FileOutputStream(headerFile);

        // may be the same object, or not
        OutputStream hdrOut, dsOut;

        
        if (dataFile.equals("LOCAL")) {
            hdrOut = new BufferedOutputStream(hdrFos, ExternalDataSource.FILEBUFFERSIZE);
            
            if (compressedData) {
                dsOut = new ZipOutputStream(hdrOut);
            }
            else {
                dsOut = hdrOut;
            }

        }
        else {
            // don't need much buffer to write header
            hdrOut = new BufferedOutputStream(hdrFos, 1024*8);


            BufferedOutputStream dataOut = 
                new BufferedOutputStream(new FileOutputStream(localDataFile), ExternalDataSource.FILEBUFFERSIZE);

            if (compressedData) {
                dsOut = new ZipOutputStream(dataOut);
            }
            else {
                dsOut = dataOut;
            }
            
        }
        
        
        writeHeader(hdrOut);

        writeData(data, dsOut);

        dsOut.close();
        
        // don't close an already closed stream, it will throw IOException
        if (!dataFile.equals("LOCAL")) {
            hdrOut.close();
        }


    }


    /**
     *
     * Write data to an output stream using current header settings for data type, scaling etc.
     *
     * @param data the 4D data array.
     * @param out an OutputStream to write the data. It should be buffered.
     * @exception IOException
     * 
     */
    private void writeData(double data[][][][], OutputStream out) throws IOException {
	
	int i,j,k,n;
	EndianCorrectOutputStream ecs;

        int XDIM = dimSize[0];
        int YDIM = dimSize[1];
        int ZDIM = dimSize[2];
        
	ecs = new EndianCorrectOutputStream(out, !intelByteOrder);


	switch (dataType) {


            // need unsigned byte for RGB images
	case UCHAR:
            for (k=0; k<ZDIM; k++) {
                for (j=0; j<YDIM; j++) {
                    for (i=0; i<XDIM; i++) {
                        for (n=0; n<channels; n++) { 
                            short s = (short)Math.round(data[i][j][k][n]);
                            ecs.write( s & 0x00ff );
                        }
                    }
                }
            }
        
            break;


	case CHAR:
            for (k=0; k<ZDIM; k++) {
                for (j=0; j<YDIM; j++) {
                    for (i=0; i<XDIM; i++) {
                        for (n=0; n<channels; n++) { 
                            ecs.write((byte)Math.round(data[i][j][k][n]));
                        }
                    }
                }
            }
        
            break;

            
	case SHORT:
            for (k=0; k<ZDIM; k++) {
                for (j=0; j<YDIM; j++) {
                    for (i=0; i<XDIM; i++) {
                        for (n=0; n<channels; n++) { 
                            ecs.writeShortCorrect((short)Math.round(data[i][j][k][n]));
                        }
                    }
                }
            }
        
            break;

            
        case INT:
            for (k=0; k<ZDIM; k++) {
                for (j=0; j<YDIM; j++) {
                    for (i=0; i<XDIM; i++) {
                        for (n=0; n<channels; n++) { 
                            ecs.writeIntCorrect((int)Math.round(data[i][j][k][n]));
                        }
                    }
                }
            }
        
            break;
            
        case LONG:
            for (k=0; k<ZDIM; k++) {
                for (j=0; j<YDIM; j++) {
                    for (i=0; i<XDIM; i++) {
                        for (n=0; n<channels; n++) { 
                            ecs.writeLongCorrect(Math.round(data[i][j][k][n]));
                        }
                    }
                }
            }
        
            break;

        case FLOAT:
            for (k=0; k<ZDIM; k++) {
                for (j=0; j<YDIM; j++) {
                    for (i=0; i<XDIM; i++) {
                        for (n=0; n<channels; n++) { 
                            ecs.writeFloatCorrect((float)(data[i][j][k][n]));
                        }
                    }
                }
            }
        
            break;            

       
        case DOUBLE:
            for (k=0; k<ZDIM; k++) {
                for (j=0; j<YDIM; j++) {
                    for (i=0; i<XDIM; i++) {
                        for (n=0; n<channels; n++) { 
                            ecs.writeDoubleCorrect(data[i][j][k][n]);
                        }
                    }
                }
            }
        
            break;            


	default:
	    throw new IOException("Sorry, cannot yet write MetaIO datatype " + dataType);

	}

	return;
    }



    /**
     * Gets the Camino data type string
     */
    public String caminoDataTypeString() {

        return dataType.getCaminoDataType();
    }


    
    /**
     * @return "True" if b, "False" if !b.
     */
    private static String titleCaseBool(boolean b) {
	if (b) {
	    return "True";
	}
	else {
	    return "False";
	}
    }
    

    /**
     * @return a String representation of the matrix <code>m</code>, with all
     * rows concatenated into a single line.
     */
    private static String flatMatrixString(RealMatrix m) {
	String line = "";

	double[][] matrix = m.entries;

	for (int r = 0; r < m.rows(); r++) {
	    for (int c = 0; c < m.columns(); c++) {
		line += matrix[r][c] + " ";
	    }
	}

	return line.trim();
    }


    // ImageHeader implementation methods

    public int xDataDim() {
	return dimSize[0];
    }

    public int yDataDim() {
	return dimSize[1]; 
    }

    public int zDataDim() {
	return dimSize[2];
    }

    public int[] getDataDims() {
	return new int[] {dimSize[0], dimSize[1], dimSize[2]};
    }

    public double xVoxelDim() {
	return spacing[0];
    }

    public double yVoxelDim() {
	return spacing[1];
    }

    public double zVoxelDim() {
	return spacing[2];
    }

    public double[] getVoxelDims() {
	return new double[] {spacing[0], spacing[1], spacing[2]};
    }

    public int components() {
	return channels > 0 ? channels : 1;
    }



    /**
     * @return a VoxelOrderDataSource for the image data.
     *
     */
    public DataSource getImageDataSource() {


		
        // Need to skip over appropriate number of header bytes
        // Since there is no standard size, we calculate it on the fly
        if (dataFile.equals("LOCAL")) {
            try {
                MetaImageHeader mh = readHeader(headerFile);
                dataByteOffset = mh.dataByteOffset;
            }
            catch (IOException e) {
                throw new LoggedException(e);
            }
        }
        
	if (compressedData) {
	    
	    logger.warning("Compressed Meta IO is untested, use with caution");

	    try {
		// set up an input stream, then skip header bytes
		FileInputStream fin = new FileInputStream(localDataFile);
		
		int bytesSkipped = 0;
		
		while (bytesSkipped < dataByteOffset) {
		    fin.skip(dataByteOffset - bytesSkipped);
		}
		
		// now buffer on the I/O side
		ZipInputStream zin = new ZipInputStream(new BufferedInputStream(fin, ExternalDataSource.FILEBUFFERSIZE/2));
		
		// assume there is one entry (ie one compressed file in the ZIP archive)
		zin.getNextEntry();
		
		// and buffer on the decompressed side for speed
		EndianNeutralDataInputStream dataIn = 
		    new EndianNeutralDataInputStream(new BufferedInflaterInputStream(zin, 
										     ExternalDataSource.FILEBUFFERSIZE/2), 
						     intelByteOrder);
		
		
		return new VoxelOrderDataSource(dataIn, channels, caminoDataTypeString());

	    }
	    catch (IOException e) {
		throw new LoggedException(e);
	    }
	}
	else {
	    VoxelOrderDataSource dataSource = 
		new VoxelOrderDataSource(localDataFile, channels, 
					 caminoDataTypeString(), intelByteOrder, dataByteOffset);
	    
	    return dataSource;
	}
    }
  


    /**
     * Reads a data array, assuming voxel-ordered data for multi-channel data.
     *
     * @return a 4D array with extent [xDataDim][yDataDim][zDataDim][components()].
     */
    public double[][][][] readVolumeData() {
	
	int xDataDim = dimSize[0];
	int yDataDim = dimSize[1];
	int zDataDim = dimSize[2];

	double[][][][] data = new double[xDataDim][yDataDim][zDataDim][channels];
	
	try {
	    
	    DataSource din = getImageDataSource();
	    
	    for (int k = 0; k < zDataDim; k++) { 
		for (int j = 0; j < yDataDim; j++) {
		    for (int i = 0; i < xDataDim; i++) {
			double[] nextVoxel = din.nextVoxel();

			for (int c = 0; c < channels; c++) {
			    data[i][j][k][c] = nextVoxel[c];
			}
		    }
		}
	    }
	    
	}
	catch(DataSourceException e) {
	    throw new LoggedException(e);
	} 

	return data;
	
    }



    /**
     * Will use TransformMatrix + Offset, because it's not clear what coordinates CenterOfRotation uses
     * 
     *
     *
     */
    public RealMatrix getVoxelToPhysicalTransform() {
        
        logger.warning("Basing transformation on TransformMatrix and Offset only (ignoring CenterOfRotation)");

        RealMatrix trans = new RealMatrix(4,4);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                trans.entries[i][j] = transformation.entries[i][j];
            }
        }

        trans.entries[0][3] = offset[0];
        trans.entries[1][3] = offset[1];
        trans.entries[2][3] = offset[2];
        trans.entries[3][3] = 1.0;

        return trans;

    }


    /**
     * Get the header file name associated with the header object.
     *
     */
    public String getHeaderFilename() {
        return headerFile;
    }

    /**
     * Get the data file name associated with the header object. May be the same as the header
     * file name.
     *
     */
    public String getDataFilename() {
	
        if (dataFile.equals("LOCAL")) {
	    return localDataFile;
	}

        return dataFile;

    }





    /**
     * Use the current header to write a scalar image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, with  fields altered where necessary 
     * for writing the data.
     *
     * @param data will be checked against header data dimensions. The data array should be in the 
     * same voxel space as this image, to ensure a correct definition of physical space.
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings. 
     *
     * @return the image header associated with the output data set
     */
    public ImageHeader writeScalarImage(double[][][] data, String fileRoot) {

        if (data.length != dimSize[0] || data[0].length != dimSize[1] || data[0][0].length != dimSize[2]) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + dimSize[0] +                                      " " + dimSize[1] + " " + dimSize[2]);
        }

        MetaImageHeader mh = new MetaImageHeader(this, fileRoot);
        
        mh.setDataTypeToSigned();

        mh.channels = 1;

        double[][][][] fourD = new double[dimSize[0]][dimSize[1]][dimSize[2]][1];

        for (int i = 0; i < dimSize[0]; i++) {
            for (int j = 0; j < dimSize[1]; j++) {
                for (int k = 0; k < dimSize[2]; k++) { 
                    fourD[i][j][k][0] = data[i][j][k];
                }
            }
        }



        try {
            mh.writeImage(fourD);
        }
        catch(IOException e) {
            throw new LoggedException(e.getMessage());
        }

        return mh;
       
    }


    /**
     * Use the current header to write a vector image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, The new header is a copy of this one, 
     * with fields altered where necessary for writing Vector data. 
     *
     * @param data will be checked against header data dimensions. The data array should be in the 
     * same voxel space as this image, to ensure a correct definition of physical space.
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings. 
     *
     * @return the image header associated with the output data set
     */
    public ImageHeader writeVectorImage(double[][][][] data, String fileRoot) {

        if (data.length != dimSize[0] || data[0].length != dimSize[1] || data[0][0].length != dimSize[2]) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + dimSize[0] +                                      " " + dimSize[1] + " " + dimSize[2]);
        }

        MetaImageHeader mh = new MetaImageHeader(this, fileRoot);
        
        mh.channels = data[0][0][0].length;

        try {
            mh.writeImage(data);
        }
        catch(IOException e) {
            throw new LoggedException(e.getMessage());
        }


        return mh;


    }


    /**
     * Use the current header to write an RGB image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, with fields altered where necessary 
     * for writing RGB data. 
     *
     * @param red the red channel intensity, should be normalized between 0 and 255.
     * @param green the green channel intensity, should be normalized between 0 and 255.
     * @param blue the blue channel intensity, should be normalized between 0 and 255.
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings. 
     *
     * @return the image header associated with the output data set
     */
    public ImageHeader writeRGB_Image(int[][][] red, int[][][] green, int[][][] blue, String fileRoot) {

        if (red.length != dimSize[0] || red[0].length != dimSize[1] || red[0][0].length != dimSize[2]) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + dimSize[0] +                                      " " + dimSize[1] + " " + dimSize[2]);
        }

        MetaImageHeader mh = new MetaImageHeader(this, fileRoot);
        
        mh.channels = 3;

        // RGB is the one case where we write unsigned data
        mh.dataType = DataType.UCHAR;

        int xDataDim = dimSize[0];
	int yDataDim = dimSize[1];
	int zDataDim = dimSize[2];

	double[][][][] data = new double[xDataDim][yDataDim][zDataDim][mh.channels];
	
        for (int k = 0; k < zDataDim; k++) { 
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
                    data[i][j][k][0] = red[i][j][k];
                    data[i][j][k][1] = green[i][j][k];
                    data[i][j][k][2] = blue[i][j][k];
                }
            }
        }
        
        try { 
            mh.writeImage(data);
        }
        catch(IOException e) {
            throw new LoggedException(e.getMessage());
        }


        return mh;


    }


    /**
     * 
     * Use the current header to write a tensor image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, with fields altered where necessary 
     * for writing tensor data. 
     *
     * @param data will be checked against header data dimensions. The data array should be in the 
     * same voxel space as this image, to ensure a correct definition of physical space. data[i][j][k] 
     * contains six components in upper-triangular order. The data will be written to disk in upper triangular order.
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings. 
     *
     * @return the image header associated with the output data set
     */
    public ImageHeader writeTensorImage(double[][][][] data, String fileRoot) {

        if (data.length != dimSize[0] || data[0].length != dimSize[1] || data[0][0].length != dimSize[2]) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + dimSize[0] +                                      " " + dimSize[1] + " " + dimSize[2]);
        }

        MetaImageHeader mh = new MetaImageHeader(this, fileRoot);

        mh.channels = 6;

        mh.setDataTypeToSigned();

        try {
            mh.writeImage(data);
        }
        catch(IOException e) {
            throw new LoggedException(e.getMessage());
        }

        return mh;

    }


    /**
     * Set the data type of the image.
     *
     * @param dataType a Camino data type string.
     *
     */
    public void setDataType(String type) {
        
        dataType = DataType.getDataType(type);

    }



    /**
     * Sets gzip compression. This may alter the file name(s) of the image.
     *
     */
    public void setGzip(boolean gz) {

        compressedData = gz;

        if (gz) {
            logger.warning("MetaIO uses ZLIB compression instead of gzip");
        }

    }

    
    /**
     * Set the data type of the image to a signed type (we don't currently write unsigned)
     *
     *
     */
    private void setDataTypeToSigned() {

	if (dataType == DataType.UCHAR) {
            dataType = DataType.CHAR;
	}
	else if (dataType == DataType.USHORT) {
            dataType = DataType.SHORT;
	}
	else if (dataType == DataType.UINT) {
            dataType = DataType.INT;
	}

    }


    public boolean lowerTriangularSymmMatrix() {
        return false;
    }
    
    
}

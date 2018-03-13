package imaging;

import misc.*;
import numerics.RealMatrix;
import tools.*;
import data.*;

import java.io.*;
import java.util.logging.Logger;
import java.util.zip.*;

/**
 * Analyze image header, provides support for the most common header fields. Some fields,
 * such as patient_id, are not currently supported. The class also does some nonstandard things.
 * The field image_dimension.funused1 is the image scale. The intensity of each pixel in the 
 * associated .img file is (image value from file) * scale. 
 * Also, the origin of the Talairach coordinates
 * (midline of the anterior commisure) are encoded in the field
 * data_history.originator. These changes are included for compatibility with SPM / MRICro.
 * <p>
 * The class also stores the PICo seed point in the fields unused8, unused9 and unused10. 
 * 
 * 
 *
 * 
 * @author Philip Cook
 * @version $Id$
 */
public final class AnalyzeHeader extends ImageHeader {

    // buffer size for I/O. Larger buffers create an overhead when creating many files
    private static int BUFFERSIZE = 1024 * 1024 * 16;

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.imaging.AnalyzeHeader");


    public static final short DT_NONE = 0;

    public static final short DT_UNKNOWN = 0;

    public static final short DT_BINARY = 1;

    /** 8-bit unsigned char (0-255) datatype. */
    public static final short DT_UNSIGNED_CHAR = 2;

    /**
     * 16-bit signed integer (-32768-32767).
     */
    public static final short DT_SIGNED_SHORT = 4;

    // 32 bit signed int
    public static final short DT_SIGNED_INT = 8;

    // 32 bit
    public static final short DT_FLOAT = 16;

    // 64 bit (real float, imag float)
    public static final short DT_COMPLEX = 32;

    // 64 bit
    public static final short DT_DOUBLE = 64;

    // 24 bit
    public static final short DT_RGB = 128;

    // SPM2 data types
    // if you write ushort the value will be cast to short
    // if you write uint the value will be cast to int
    
    /** 8-bit signed char, equivalent to Java byte. */
    public static final short DT_SIGNED_CHAR = 130;

    /** 16-bit unsigned short datatype. */
    public static final short DT_UNSIGNED_SHORT = 132;

    /** 32-bit unsigned int datatype. */
    public static final short DT_UNSIGNED_INT = 136;

    // end SPM2

    public static final short DT_ALL = 255;


    protected short datatype = DT_UNKNOWN;

    protected short width = 0; // image dimensions

    protected short height = 0; // y

    protected short depth = 0; // z

    protected float vox_offset = 0.0f; // byte offset in img file. Voxels begin after

    // $vox_offset bytes

    // ITK will write this as zero, may need to check for this
    protected short nImages = 1; // for 4D images

    protected boolean intelByteOrder = false;

    protected int glmin = 0; // min max grey levels

    protected int glmax = 0;

    protected float pixelWidth = 1.0f; // voxel dimensions

    protected float pixelHeight = 1.0f;

    protected float pixelDepth = 1.0f;

    // multiply pixel values by this
    protected float scaleSlope = 1.0f;
    
    // and add this
    protected float scaleInter = 0.0f;

    // bits per pixel
    protected short bitpix = 1;

    protected String description = "";

    /** Centre of the Talairach coordinate system */
    protected short[] centre = new short[3];

    /** Name of the header file that this object was read from, if any. */
    protected String hdrFile = null;
    

    /** Name of the image file that goes with this header. */
    protected String dataFile = null;

    /** If true, read / write compressed output to/from a .img.gz file */
    protected boolean gzip = false;
    

    /**
     * Creates an AnalyzeHeader object with all of its fields set to their
     * default value.
     *
     * @param fileRoot the output file root. We will write fileRoot.hdr and fileRoot.img[.gz].
     * @param gzipImg if true, the output image will be written to fileRoot.img.gz 
     *
     */
    public AnalyzeHeader(String fileRoot, boolean gzipImg) {
        hdrFile = hdrFile + ".hdr";

        gzip = gzipImg;
        
        if (gzip) {
            dataFile = fileRoot + ".img.gz";
        }
        else {
            dataFile = fileRoot + ".img";
        }


    }


    protected AnalyzeHeader() {
        
    }


    /** Copy constructor. */
    protected AnalyzeHeader(AnalyzeHeader ah, String fileRoot) {

	intelByteOrder = ah.intelByteOrder; 
	
	width = ah.width;
	height = ah.height;
	depth = ah.depth;
	nImages = ah.nImages;
	
	datatype = ah.datatype;
	bitpix = ah.bitpix;

	pixelWidth = ah.pixelWidth;
	pixelHeight = ah.pixelHeight;
	pixelDepth = ah.pixelDepth;

	vox_offset = ah.vox_offset;
	scaleSlope = ah.scaleSlope;
	scaleInter = ah.scaleInter;

	glmax = ah.glmax;
	glmin = ah.glmin;

	// data_history 

	description = ah.description;

	centre[0] = ah.centre[0];
	centre[1] = ah.centre[1];
	centre[2] = ah.centre[2];

	hdrFile = fileRoot + ".hdr";

        gzip = ah.gzip;

        if (gzip) {
            dataFile = fileRoot + ".img.gz";
        }
        else {
            dataFile = fileRoot + ".img";
        }


    }


   

    /**
     * Reads a header from the given file. The file should be either name.hdr (the actual header) or
     * name (where name.hdr exists). 
     *
     * @throws IOException if the header cannot be read.
     */
    public static AnalyzeHeader readHeader(String hdrFile) throws IOException {

        if ( !(hdrFile.endsWith(".hdr")) ) {
            hdrFile = hdrFile + ".hdr";
        }
        
        FileInputStream filein = new FileInputStream(hdrFile);
        DataInputStream input = new DataInputStream(filein);
	
        AnalyzeHeader ah = new AnalyzeHeader();

        int i;

        int byteCounter = 0;

        //  header_key
        // we ignore this
        ah.readInt(input); // sizeof_hdr
        byteCounter += 4;

        for (i = 0; i < 10; i++)
            input.read(); // data_type
        for (i = 0; i < 18; i++)
            input.read(); // db_name

        byteCounter += 28;

        ah.readInt(input); // extents
        ah.readShort(input); // session_error
        input.readByte(); // regular
        input.readByte(); // hkey_un0

        byteCounter += 8;

        // image_dimension

        short endian = ah.readShort(input); // dim[0] should always be 4, may be up to 7 for NIFTI files
        if ((endian < 0) || (endian > 15)) {
            ah.intelByteOrder = true;
        }

        byteCounter += 2;

        ah.width = ah.readShort(input);
        byteCounter += 2; // dim[1]
        ah.height = ah.readShort(input);
        byteCounter += 2; // dim[2]
        ah.depth = ah.readShort(input);
        byteCounter += 2; // dim[3]
        ah.nImages = ah.readShort(input);
        byteCounter += 2; // dim[4]
        for (i = 0; i < 3; i++)
            ah.readShort(input);
        byteCounter += 6; // dim[5-7]
	
	ah.readShort(input);
        byteCounter += 2; // unused8  
	ah.readShort(input);
	byteCounter += 2; // unused9  
	ah.readShort(input);
        byteCounter += 2; // unused10 
        ah.readShort(input);
        byteCounter += 2; // unused11
        ah.readShort(input);
        byteCounter += 2; // unused12
        ah.readShort(input);
        byteCounter += 2; // unused13
        ah.readShort(input);
        byteCounter += 2; // unused14

        ah.datatype = ah.readShort(input);
        byteCounter += 2; // datatype
        ah.bitpix = ah.readShort(input);
        byteCounter += 2; // bitpix
        ah.readShort(input);
        byteCounter += 2; // dim_un0
        ah.readFloat(input);
        byteCounter += 4; // pixdim[0]
        ah.pixelWidth = ah.readFloat(input);
        byteCounter += 4; // pixdim[1]
        ah.pixelHeight = ah.readFloat(input);
        byteCounter += 4; // pixdim[2]
        ah.pixelDepth = ah.readFloat(input);
        byteCounter += 4; // pixdim[3]
        for (i = 0; i < 4; i++)
            ah.readFloat(input);
        byteCounter += 16; // pixdim[4-7]
        ah.vox_offset = ah.readFloat(input);
        byteCounter += 4; // vox_offset
        ah.scaleSlope = ah.readFloat(input);
        byteCounter += 4; // funused1 (SPM scale slope)
        ah.scaleInter = ah.readFloat(input);
        byteCounter += 4; // funused2 (SPM scale intercept)
        ah.readFloat(input);
        byteCounter += 4; // funused3
        ah.readFloat(input);
        byteCounter += 4; // cal_max
        ah.readFloat(input);
        byteCounter += 4; // cal_min
        ah.readInt(input);
        byteCounter += 4; // compressed
        ah.readInt(input);
        byteCounter += 4; // verified
        ah.glmax = ah.readInt(input);
        byteCounter += 4; //(int) s.max // glmax
        ah.glmin = ah.readInt(input);
        byteCounter += 4; //(int) s.min // glmin


        // data_history

        byte[] desc = new byte[80];

        input.read(desc, 0, 80);
        byteCounter += 80; // descrip

        // desc may be shorter than 80, and is null terminated

        int charCounter = 0;

        while (charCounter < 80 && desc[charCounter] != '\0') {
            charCounter++;
        }

        byte[] actualDesc = new byte[charCounter];

        for (i = 0; i < charCounter; i++) {
            actualDesc[i] = desc[i];
        }

        ah.description = new String(actualDesc);

        for (i = 0; i < 24; i++)
            input.read();
        byteCounter += 24; // aux_file
        input.read();
        byteCounter += 1; // orient
        ah.centre[0] = ah.readShort(input);
        byteCounter += 2; // originator
        ah.centre[1] = ah.readShort(input);
        byteCounter += 2; // originator
        ah.centre[2] = ah.readShort(input);
        byteCounter += 2; // originator

        for (i = 0; i < 4; i++)
            input.read();
        byteCounter += 4; // originator, unused remainder
        for (i = 0; i < 10; i++)
            input.read();
        byteCounter += 10; // generated
        for (i = 0; i < 10; i++)
            input.read();
        byteCounter += 10; // scannum
        for (i = 0; i < 10; i++)
            input.read();
        byteCounter += 10; // patient_id
        for (i = 0; i < 10; i++)
            input.read();
        byteCounter += 10; // exp_date
        for (i = 0; i < 10; i++)
            input.read();
        byteCounter += 10; // exp_time
        for (i = 0; i < 3; i++)
            input.read();
        byteCounter += 3; // hist_un0
        ah.readInt(input);
        byteCounter += 4; // views
        ah.readInt(input);
        byteCounter += 4; // vols_added
        ah.readInt(input);
        byteCounter += 4; // start_field
        ah.readInt(input);
        byteCounter += 4; // field_skip
        ah.readInt(input);
        byteCounter += 4; // omax
        ah.readInt(input);
        byteCounter += 4; // omin
        ah.readInt(input);
        byteCounter += 4; // smax
        ah.readInt(input);
        byteCounter += 4; // smin

        if (byteCounter != 348) {
            logger.warning("Header is " + byteCounter + " bytes.\nAnalyze expects 348 bytes");
        }
        
        if (ah.pixelWidth < 0.0f || ah.pixelHeight < 0.0f || ah.pixelDepth < 0.0f) {
            logger.warning("Negative pixel dimensions detected - using absolute values. Physical space " + 
                           "may be misaligned, convert to nii if possible");

            ah.pixelWidth = (float)Math.abs(ah.pixelWidth);
            ah.pixelHeight = (float)Math.abs(ah.pixelHeight);
            ah.pixelDepth = (float)Math.abs(ah.pixelDepth);
        }

	ah.checkBitpix();
	
        input.close();
        filein.close();

	ah.hdrFile = hdrFile;
	ah.dataFile = getIMG_File(hdrFile);

        if (ah.dataFile != null && ah.dataFile.endsWith(".gz")) {
            ah.gzip = true;
        }

        return ah;
    }


    /**
     * Check that the bits per pixel agree with the Java data type.
     */
    private void checkBitpix() {

	int existingBitpix = bitpix;

	switch (datatype) {

        case DT_BINARY:
	    bitpix = 1;
            break;
        case DT_UNSIGNED_CHAR:
        case DT_SIGNED_CHAR:
	    bitpix = 8;
            break;
        case DT_SIGNED_SHORT:
        case DT_UNSIGNED_SHORT:
	    bitpix = 16;
	    break;
        case DT_SIGNED_INT:
        case DT_FLOAT:
        case DT_UNSIGNED_INT:
            bitpix = 32;
	    break;
        case DT_RGB:
            bitpix = 24;
	    break;
        case DT_COMPLEX:
        case DT_DOUBLE:
	    bitpix = 64;
            break;
	}	

    }


    /**
     * Writes header to the file. Will modify bitpix if necessary so that the definition
     * of each datatype agrees with the Java specification, eg bitpix == 32 for float. Call
     * <code>setHeaderFile</code> to set the file name of the header.
     *
     * 
     */
    public void writeHeader() throws IOException {
	
        FileOutputStream fileout = new FileOutputStream(hdrFile);
        DataOutputStream output = new DataOutputStream(fileout);

	writeHeader(output);

	output.close();

    }


    /**
     * Writes the header to the output stream.
     *
     */
    private void writeHeader(DataOutputStream out) throws IOException {

	// enforce correct bitpix
	checkBitpix();

	DataOutput output = null;

	if (intelByteOrder) {
	    output = new tools.LEFilterOutputStream(out);
	}
	else {
	    output = out;
	}

        //     header_key

        int byteCounter = 0;

        output.writeInt(348);
        byteCounter += 4; // sizeof_hdr

        int i;
        for (i = 0; i < 10; i++)
            output.write(0);
        byteCounter += 10; // data_type
        for (i = 0; i < 18; i++)
            output.write(0);
        byteCounter += 18; // db_name
        output.writeInt(16384);
        byteCounter += 4; // extents
        output.writeShort(0);
        byteCounter += 2; // session_error
        output.writeByte((int) 'r');
        byteCounter += 1; // regular
        output.writeByte(0);
        byteCounter += 1; // hkey_un0


        // image_dimension

        output.writeShort((short) 4);
        byteCounter += 2; // dim[0]
        output.writeShort(width);
        byteCounter += 2; // dim[1]
        output.writeShort(height);
        byteCounter += 2; // dim[2]
        output.writeShort(depth);
        byteCounter += 2; // dim[3]
        output.writeShort(nImages);
        byteCounter += 2; // dim[4]
        for (i = 0; i < 3; i++)
            output.writeShort(0);
        byteCounter += 6; // dim[5-7]

        output.writeShort(0);
        byteCounter += 2; // unused8
        output.writeShort(0);
        byteCounter += 2; // unused9
        output.writeShort(0);
        byteCounter += 2; // unused10
        output.writeShort(0);
        byteCounter += 2; // unused11
        output.writeShort(0);
        byteCounter += 2; // unused12
        output.writeShort(0);
        byteCounter += 2; // unused13
        output.writeShort(0);
        byteCounter += 2; // unused14

        //	    output.writeBytes ( "mm\0\0" ); // vox_units
        //	    for (i = 0; i < 8; i++) output.write( 0 ); // cal_units[8]
        //	    output.writeShort( 0 ); // unused1
        output.writeShort(datatype);
        byteCounter += 2; // datatype
        output.writeShort(bitpix);
        byteCounter += 2; // bitpix
        output.writeShort(0);
        byteCounter += 2; // dim_un0

        output.writeFloat(0);
        byteCounter += 4; // pixdim[0]
        output.writeFloat(pixelWidth);
        byteCounter += 4; // pixdim[1]
        output.writeFloat(pixelHeight);
        byteCounter += 4; // pixdim[2]
        output.writeFloat(pixelDepth);
        byteCounter += 4; // pixdim[3]
        for (i = 0; i < 4; i++)
            output.writeFloat(0);
        byteCounter += 16; // pixdim[4-7]

        output.writeFloat(vox_offset);
        byteCounter += 4; // vox_offset
        output.writeFloat(scaleSlope);
        byteCounter += 4; // funused1 (SPM scale)
        output.writeFloat(scaleInter);
        byteCounter += 4; // funused2 (SPM scale int)
        output.writeFloat(0);
        byteCounter += 4; // funused3
        output.writeFloat(0);
        byteCounter += 4; // cal_max
        output.writeFloat(0);
        byteCounter += 4; // cal_min
        output.writeInt(0);
        byteCounter += 4; // compressed
        output.writeInt(0);
        byteCounter += 4; // verified

        output.writeInt(glmax);
        byteCounter += 4; // glmax
        output.writeInt(glmin);
        byteCounter += 4; // glmin

        // image_dimension");

        // data_history
        byte[] desc = description.getBytes();

        int descLength = 0;

        if (desc.length < 80) {
            for (i = 0; i < desc.length; i++) {
                byteCounter++;
                output.write(desc[i]); // descrip
            }
            for (i = 0; i < 80 - desc.length; i++) {
                output.write('\0');
                byteCounter++;
            }

        }
        else {
            for (i = 0; i < 79; i++) {
                output.write(desc[i]);
            }
            output.write('\0');
            byteCounter += 80; // descrip
        }

        for (i = 0; i < 24; i++)
            output.write(0);
        byteCounter += 24; // aux_file
        output.write(0);
        byteCounter += 1; // orient

        // SPM behaviour
        for (i = 0; i < 3; i++)
            output.writeShort(centre[i]);
        byteCounter += 6; // originator
        for (i = 0; i < 4; i++)
            output.write(0);
        byteCounter += 4; // originator
        for (i = 0; i < 10; i++)
            output.write(0);
        byteCounter += 10; // generated
        for (i = 0; i < 10; i++)
            output.write(0);
        byteCounter += 10; // scannum
        for (i = 0; i < 10; i++)
            output.write(0);
        byteCounter += 10; // patient_id
        for (i = 0; i < 10; i++)
            output.write(0);
        byteCounter += 10; // exp_date
        for (i = 0; i < 10; i++)
            output.write(0);
        byteCounter += 10; // exp_time
        for (i = 0; i < 3; i++)
            output.write(0);
        byteCounter += 3; // hist_un0

        output.writeInt(0);
        byteCounter += 4; // views
        output.writeInt(0);
        byteCounter += 4; // vols_added
        output.writeInt(0);
        byteCounter += 4; // start_field
        output.writeInt(0);
        byteCounter += 4; // field_skip
        output.writeInt(0);
        byteCounter += 4; // omax
        output.writeInt(0);
        byteCounter += 4; // omin
        output.writeInt(0);
        byteCounter += 4; // smax
        output.writeInt(0);
        byteCounter += 4; // smin

        if (byteCounter != 348) {
            System.err.println("WARNING: Header is " + byteCounter + " bytes");
            System.err.println("Analyze expects 348 bytes");
        }

    }


    public int readInt(DataInputStream input) throws IOException {
        if (!intelByteOrder) return input.readInt();
        byte b1 = input.readByte();
        byte b2 = input.readByte();
        byte b3 = input.readByte();
        byte b4 = input.readByte();
        return ((((b4 & 0xff) << 24) | ((b3 & 0xff) << 16) | ((b2 & 0xff) << 8) | (b1 & 0xff)));
    }


    public short readShort(DataInputStream input) throws IOException {
        if (!intelByteOrder) return input.readShort();
        byte b1 = input.readByte();
        byte b2 = input.readByte();
        return ((short) (((b2 & 0xff) << 8) | (b1 & 0xff)));
    }


    public float readFloat(DataInputStream input) throws IOException {
        if (!intelByteOrder) return input.readFloat();
        int orig = readInt(input);
        return (Float.intBitsToFloat(orig));
    }


    // This needs moving to an actual application

    public static void main(String[] args) {

        AnalyzeHeader ah = new AnalyzeHeader();

        if (args.length == 0) {
            System.err
                    .println("\nTo write a header: AnalyzeHeader\n"
                            + "\t-voxeldims [x y z] voxel dimensions\n"
                            + "\t-datadims [x y z] data dimensions\n"
                            + "\t-nimages [number] number of images in the img file. Default 1.\n"
                            + "\t-datatype [char | short | int | float | complex | double]\n data type of image\n"
                            + "\t-offset [value] a float value.\n"
                            + "\t-gl [min max] greylevel. Stored as int in the header\n"
                            + "\t-scaleslope [value] SPM scale factor. Default 1.0\n"
                            + "\t-scaleint [value] SPM scale intercept Default 0.0\n"
                            + "\t-description [String] no spaces, max length 79 bytes. Default blank.\n"
                            + "\t-centre [x y z] voxel specifying origin of Talairach coordinate system for SPM, default [0 0 0]\n"
			     + "\t-picoseed [x y z] voxel specifying the seed (for PICo maps), default [0 0 0].\n"
                            +

                            "\t-outputfile [filename] output file name. Extension .hdr will be added if missing.\n"
			    + "\t-intelbyteorder output little-endian header. Default is big-endian\n\n"
                            + "To read a header: AnalyzeHeader -readheader [filename]");
            System.exit(0);
        }

        String filename = "";

        CL_Initializer.CL_init(args);


	// got to do this first
        for (int i = 0; i < args.length; i++) {
	    if (args[i].equals("-initfromheader")) {
		
                try {
		    ah = readHeader(args[i + 1]);
                }
                catch (IOException e) {
                    throw new LoggedException(e);
                }

		CL_Initializer.markAsParsed(i,2);
            }

	}

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-readheader")) {
                try {
		    AnalyzeHeader header = readHeader(args[i + 1]);

                    System.out.println(header);
		    System.out.println(getArgs(header, "all"));
                }
                catch (IOException e) {
                    throw new LoggedException(e);
                }

		System.exit(0);

            }
	    else if (args[i].equals("-printimagedims")) {

		try {
		    ah = readHeader(args[i + 1]);
		    
		    String dims = "-datadims " + ah.width + " " + ah.height + " " + ah.depth + 
			" -voxeldims " + ah.pixelWidth + " " + ah.pixelHeight + " " + ah.pixelDepth;
		    
		    System.out.print(dims);
		    System.exit(0);
		    
		}
                catch (IOException e) {
                    throw new LoggedException(e);
                }

	    }
	    else if (args[i].equals("-printintelbyteorder")) {
		try {
		    ah = readHeader(args[i + 1]);
		    System.out.print(ah.intelByteOrder ? 1 : 0);
		    System.exit(0);
		}
		catch (IOException e) {
                    throw new LoggedException(e);
                }    
	    }
	    else if (args[i].equals("-printbigendian")) {
		try {
		    ah = readHeader(args[i + 1]);
		    System.out.print(ah.intelByteOrder ? 0 : 1);
		    System.exit(0);
		}
		catch (IOException e) {
                    throw new LoggedException(e);
                }    
	    }
	    else if (args[i].equals("-printprogargs")) {

		try {
		    ah = readHeader(args[i + 1]);
		    
		    System.out.print(getArgs(ah, args[i+2]));
		    System.exit(0);
		    
		}
                catch (IOException e) {
                    throw new LoggedException(e);
                }

	    }
	    else if (args[i].equals("-initfromheader")) {
		// already caught
	    }
            else if (args[i].equals("-centre")) {
                ah.centre[0] = Short.parseShort(args[i + 1]);
                ah.centre[1] = Short.parseShort(args[i + 2]);
                ah.centre[2] = Short.parseShort(args[i + 3]);

		CL_Initializer.markAsParsed(i,4);
            }
	    else if (args[i].equals("-intelbyteorder")) {
		ah.intelByteOrder = true;
		CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-networkbyteorder")) {
		ah.intelByteOrder = false;
		CL_Initializer.markAsParsed(i);
	    }
            else if (args[i].equals("-datadims")) {
                ah.width = Short.parseShort(args[i + 1]);
                ah.height = Short.parseShort(args[i + 2]);
                ah.depth = Short.parseShort(args[i + 3]);

		CL_Initializer.markAsParsed(i,4);
            }
            else if (args[i].equals("-datatype")) {
                if (args[i + 1].equals("binary")) {
                    ah.datatype = DT_BINARY;
                }
                if (args[i + 1].equals("char")) {
                    ah.datatype = DT_UNSIGNED_CHAR;
                }
                if (args[i + 1].equals("byte")) {
                    ah.datatype = DT_SIGNED_CHAR;
                }
                if (args[i + 1].equals("short")) {
                    ah.datatype = DT_SIGNED_SHORT;
                }
                if (args[i + 1].equals("int")) {
                    ah.datatype = DT_SIGNED_INT;
                }
                if (args[i + 1].equals("ushort")) {
                    ah.datatype = DT_UNSIGNED_SHORT;
                }
                if (args[i + 1].equals("uint")) {
                    ah.datatype = DT_UNSIGNED_INT;
                }
                if (args[i + 1].equals("float")) {
                    ah.datatype = DT_FLOAT;
                }
                if (args[i + 1].equals("complex")) {
                    ah.datatype = DT_COMPLEX;
                }
                if (args[i + 1].equals("double")) {
                    ah.datatype = DT_DOUBLE;
                }

		CL_Initializer.markAsParsed(i,2);

            }
            else if (args[i].equals("-description")) {

		int j = 0;

		ah.description = args[i + 1];		    

		while (i + j + 2 < args.length && !args[i + j + 2].startsWith("-")) {
		    ah.description += " " + args[i + j + 2];		    
		    j++;
		}

		CL_Initializer.markAsParsed(i, j + 2);

            }
            else if (args[i].equals("-gl")) {
                ah.glmin = Short.parseShort(args[i + 1]);
                ah.glmax = Short.parseShort(args[i + 2]);

		CL_Initializer.markAsParsed(i, 3);
            }
            else if (args[i].equals("-nimages")) {
                ah.nImages = Short.parseShort(args[i + 1]);
		CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-offset")) {
                ah.vox_offset = Float.parseFloat(args[i + 1]);
		CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-scale")) {
                ah.scaleSlope = Float.parseFloat(args[i + 1]);
		CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-scaleslope")) {
                ah.scaleSlope = Float.parseFloat(args[i + 1]);
		CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-scaleinter")) {
                ah.scaleInter = Float.parseFloat(args[i + 1]);
		CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-voxeldims")) {
                ah.pixelWidth = Float.parseFloat(args[i + 1]);
                ah.pixelHeight = Float.parseFloat(args[i + 2]);
                ah.pixelDepth = Float.parseFloat(args[i + 3]);
		CL_Initializer.markAsParsed(i, 4);
            }

        }

	CL_Initializer.checkParsing(args);

        switch (ah.datatype) {

        case DT_BINARY:
            logger.fine("Writing header data type as binary bit");
            ah.bitpix = 1;
            break;
        case DT_UNSIGNED_CHAR:
            logger.fine("Writing header data type as 8 bit char");
            ah.bitpix = 8;
            break;
        case DT_SIGNED_SHORT:
            logger.fine("Writing header data type as short");
            ah.bitpix = 16;
            break;
        case DT_SIGNED_INT:
            logger.fine("Writing header data type as int");
            ah.bitpix = 32;
            break;
        case DT_FLOAT:
            logger.fine("Writing header data type as float");
            ah.bitpix = 32;
            break;
        case DT_RGB:
            logger.fine("Writing header data type as RGB");
            ah.bitpix = 24;
            break;
        case DT_COMPLEX:
            logger.fine("Writing header data type as COMPLEX");
            ah.bitpix = 64;
            break;
        case DT_DOUBLE:
            logger.fine("Writing header data type as double");
            ah.bitpix = 64;
            break;
        case DT_SIGNED_CHAR:
            logger.fine("Writing header data type as byte");
            ah.bitpix = 8;
            break;
        case DT_UNSIGNED_SHORT:
            logger.fine("Writing header data type as ushort");
            ah.bitpix = 16;
            break;
        case DT_UNSIGNED_INT:
            logger.fine("Writing header data type as uint");
            ah.bitpix = 32;
            break;

        // don't let people create headers with silly datatypes

        default:
            throw new LoggedException("Unrecognized datatype " + ah.datatype);

        }

        try {
            if (OutputManager.outputFile != null) {

                if (OutputManager.outputFile.endsWith(".hdr")) {
                    ah.setHeaderFile(OutputManager.outputFile);
		    ah.writeHeader();
                }
                else {
                    ah.setHeaderFile(OutputManager.outputFile + ".hdr");
		    ah.writeHeader();
                }
            }
            else {
                OutputManager om = new OutputManager();
                ah.writeHeader(om.getOutputStream());
		om.close();
            }


        }
        catch (java.io.IOException e) {
            throw new LoggedException("Can't write header");
        }



    }

    public String toString() {

        StringBuffer buffer = new StringBuffer();

	String datatypeString = "Unknown";

        switch (datatype) {

        case DT_BINARY:
            datatypeString = "binary bit";
            break;
        case DT_UNSIGNED_CHAR:
            datatypeString = "char (8-bit)";
            break;
        case DT_SIGNED_SHORT:
	    datatypeString = "short";
	    break;
        case DT_SIGNED_INT:
            datatypeString = "int";
	    break;
        case DT_FLOAT:
	    datatypeString = "float";
            break;
        case DT_RGB:
	    datatypeString = "RGB";
	    break;
        case DT_COMPLEX:
	    datatypeString = "Complex";
            break;
        case DT_DOUBLE:
	    datatypeString = "double";
            break;
        case DT_SIGNED_CHAR:
            datatypeString = "byte (8-bit signed)";
            break;
        case DT_UNSIGNED_SHORT:
            datatypeString = "ushort (16-bit unsigned)";
            break;
        case DT_UNSIGNED_INT:
            datatypeString = "uint (32-bit unsigned)";
            break;

	}

        buffer.append("datatype\t: " + datatypeString + " (" + datatype + ")\n");
        buffer.append("width\t\t: " + width + "\n");
        buffer.append("height\t\t: " + height + "\n");
        buffer.append("depth\t\t: " + depth + "\n");

        buffer.append("pixelWidth\t: " + pixelWidth + "\n");
        buffer.append("pixelHeight\t: " + pixelHeight + "\n");
        buffer.append("pixelDepth\t: " + pixelDepth + "\n");

        buffer.append("voxel_offset\t: " + vox_offset + "\n");
        buffer.append("nImages\t\t: " + nImages + "\n");

        buffer.append("intelByteOrder\t: " + intelByteOrder + "\n");

        buffer.append("glmin\t\t: " + glmin + "\n");
        buffer.append("glmax\t\t: " + glmax + "\n");

        buffer.append("scale slope\t: " + scaleSlope + "\n");
        buffer.append("scale inter\t: " + scaleInter + "\n");

	// need to put scale int here

        buffer.append("bitpix\t\t: " + bitpix + "\n");

        buffer.append("description\t: " + description + "\n");

        buffer.append("origin\t\t: " + centre[0] + " " + centre[1] + " " + centre[2]
		      + "\n");

        return buffer.toString();

    }


    /**
     * Get Camino program arguments for this header. If <code>prog</code> is the name of a 
     * supported Camino program, the returned String is args for that program only, 
     * without any surrounding context. Otherwise, a list of programs and args is returned.
     *
     * @param prog Use "all" to get general information. Other supported program names are
     * "shredder", "vcthreshselect", "pdview", "track". Anything else returns an empty String.
     *
     */
    public static String getArgs(AnalyzeHeader ah, String prog) {


	String datatypeString = "Unknown";
	
	switch (ah.datatype) {
	    
	case DT_UNSIGNED_CHAR:
	    datatypeString = "char";
	    break;
	case DT_SIGNED_CHAR:
	    datatypeString = "byte";
	    break;
	case DT_SIGNED_SHORT:
	    datatypeString = "short";
	    break;
	case DT_SIGNED_INT:
	    datatypeString = "int";
	    break;
	case DT_UNSIGNED_SHORT:
	    datatypeString = "ushort";
	    break;
	case DT_UNSIGNED_INT:
	    datatypeString = "uint";
	    break;
	case DT_FLOAT:
	    datatypeString = "float";
	    break;
	case DT_DOUBLE:
	    datatypeString = "double";
	    break;

	default: 
	    return "Camino cannot read files of this data type";
	}

	String shredderArgs = "0 " + ah.bitpix / 8 + " 0";

	if (ah.intelByteOrder) {
	    shredderArgs = "0 -" + ah.bitpix / 8 + " 0";
	}

	String scanner2voxelArgs = "-voxels " + ah.width * ah.height * ah.depth + 
	    " -inputdatatype " + datatypeString;

	if (ah.nImages > 1) {
	    scanner2voxelArgs +=  " -components " + ah.nImages;
	}
	
	String pdviewArgs = "-datadims " + ah.width + " " + ah.height + " " + ah.depth + 
	    " -voxeldims " + ah.pixelWidth + " " + ah.pixelHeight + " " + ah.pixelDepth;

	String vcthreshselectArgs = "-datadims " + ah.width + " " + ah.height + " " + ah.depth;

	String trackArgs = "-datadims " + ah.width + " " + ah.height + " " + ah.depth + 
	    " -voxeldims " + ah.pixelWidth + " " + ah.pixelHeight + " " + ah.pixelDepth;

	if (prog.equals("all")) {

	    StringBuffer buffer = new StringBuffer();
	    
	    if (ah.intelByteOrder) {
		buffer.append("This file is in little-endian byte order.\n" + 
			      "Convert to network byte order with the following command:\n\tshredder " + shredderArgs + " < [file] > [output].\n\n");
	    }
	    
	    
	    buffer.append("\nCamino args computed from this header:\n");
	    buffer.append("Program\t\targs\n");
	    buffer.append("scanner2voxel\t");
	    buffer.append(scanner2voxelArgs);
	    buffer.append("\n\n");
	    
	    buffer.append("pdview\t\t");
	    buffer.append(pdviewArgs);
	    buffer.append("\n\n");
	    
	    buffer.append("vcthreshselect\t");
	    buffer.append(vcthreshselectArgs);
	    buffer.append("\n\n");
	    
	    buffer.append("track\t\t");
	    buffer.append(trackArgs);
	    buffer.append("\n\n");
	    
	    
	    return buffer.toString();
	    
	}
	else if (prog.equals("shredder")) {
	    return shredderArgs;
	}
	else if (prog.equals("scanner2voxel")) {
	    return scanner2voxelArgs;
	}
	else if (prog.equals("pdview")) {
	    return pdviewArgs;
	}
	else if (prog.equals("vcthreshselect")) {
	    return vcthreshselectArgs;
	}
	else if (prog.equals("track")) {
	    return trackArgs;
	}
	else {
	    return "";
	}
	
    

    }


    /** 
     * Creates an AnalyzeHeader object, from specified dimensions. Bitpix is set automatically, nImages is set to 1. 
     * All other parameters take their default value.
     * 
     */
    public static AnalyzeHeader getHeader(String fileRoot, int[] dataDims, double[] voxelDims, short datatype) {

	AnalyzeHeader ah = new AnalyzeHeader(fileRoot, false);

	ah.datatype = datatype;

	// enforce correct bitpix
	ah.checkBitpix();

	ah.width = (short)dataDims[0];
	ah.height = (short)dataDims[1];
	ah.depth = (short)dataDims[2];

	ah.nImages = 1;
	
	ah.pixelWidth = (float)voxelDims[0]; // voxel dimensions
	ah.pixelHeight = (float)voxelDims[1];
	ah.pixelDepth = (float)voxelDims[2];

	return ah;
  
    }



    /**
     * Creates an AnalyzeHeader object, with voxel and data dims taken from data. GL is also set from
     * the data. Number of images is assumed to be 1.
     * 
     * @param fileRoot the root of the image name (without .hdr / .img)
     * @param data the data dimensions are calculated from this array.
     * @param voxelDims <code>{x y z}</code> the voxel dimensions.
     * @param datatype one of the Analyze data types.
     *
     */
    public static AnalyzeHeader getHeader(String fileRoot, double[][][] data, double[] voxelDims, short datatype) {

	AnalyzeHeader ah = getHeader(fileRoot, new int[] {data.length, data[0].length, data[0][0].length}, 
				     voxelDims, datatype);
	
	double min = Double.MAX_VALUE;
	double max = -min;

	for (int k = 0; k < data[0][0].length; k++) {
	    for (int j = 0; j < data[0].length; j++) {
		for (int i = 0; i < data.length; i++) {
		    if (data[i][j][k] < min) {
			min = data[i][j][k];
		    }
		    if (data[i][j][k] > max) {
			max = data[i][j][k];
		    }
		}
	    }
	}

	ah.glmin = (int)(min);
	ah.glmax = (int)(max);

	return ah;
  
    }

  
     
    /**
     * Writes file.hdr and file.img. Values will be scaled down if the scaleSlope is nonzero. 
     * Writes in network byte order.
     * @param gzip if true, compress .img part of image.
     *
     */
    private void writeImage(double[][][][] data) throws IOException {

        writeHeader();
        
        FileOutputStream fout = null;
        EndianCorrectOutputStream dout = null;
        
        fout = new FileOutputStream(dataFile);

        if (gzip) {
            dout = new EndianCorrectOutputStream(new BufferedOutputStream(new GZIPOutputStream(fout, BUFFERSIZE), BUFFERSIZE), 
                                                 !intelByteOrder);
        }
        else {
            dout = new EndianCorrectOutputStream(new BufferedOutputStream(fout, BUFFERSIZE), !intelByteOrder);
        }
        
        	
        for (int n = 0; n < nImages; n++) {
            for (int k = 0; k < depth; k++) {
                for (int j = 0; j < height; j++) {
                    for (int i = 0; i < width; i++) {
                        writeValue(data[i][j][k][n], dout, this);
                    }
                }
            }
            
        }
        dout.close();
            
	
    }

    
    private static final void writeValue(double value, EndianCorrectOutputStream out, AnalyzeHeader ah) 
	throws IOException {
	
	double scaledData = value;

	if (ah.scaleSlope != 0.0) {
	    scaledData = (value - ah.scaleInter) / ah.scaleSlope;
	}
	
	if (ah.datatype == DT_DOUBLE) {
	    out.writeDoubleCorrect(scaledData);
	}
	else if (ah.datatype == DT_FLOAT) {
	    out.writeFloatCorrect((float)(scaledData));
	}
	else if (ah.datatype == DT_SIGNED_INT) {
	    out.writeIntCorrect((int)Math.round(scaledData));
	}
	else if (ah.datatype == DT_SIGNED_SHORT) {
	    out.writeShortCorrect((short)Math.round(scaledData));
	}
	else if (ah.datatype == DT_SIGNED_CHAR) {
	    out.writeByte((byte)Math.round(scaledData));
	}		    
        else {
            throw new LoggedException("Cannot write Analyze data type " + ah.datatype);
        }
	
    }
    



    /**
     * Sets the header file for this object.
     *
     */
    public void setHeaderFile(String file) {
	String root = file;
	
	if (root.endsWith(".hdr")) {	    
	    root = root.substring(0, root.length() - 4);
	}
	
	hdrFile = root + ".hdr";

	dataFile = getIMG_File(hdrFile);
    }


    /**
     * @return Camino data type.
     *
     * @see data.ExternalDataSource
     */
    public String caminoDataTypeString() {
	
	if (datatype == DT_UNSIGNED_CHAR) {
	    return "ubyte";
	}
	if (datatype == DT_SIGNED_CHAR) {
	    return "byte";
	}
	if (datatype == DT_SIGNED_SHORT) {
	    return "short";
	}
	if (datatype == DT_SIGNED_INT) {
	    return "int";
	}
	if (datatype == DT_UNSIGNED_SHORT) {
	    return "ushort";
	}
	if (datatype == DT_UNSIGNED_INT) {
	    return "uint";
	}
	if (datatype == DT_FLOAT) {
	    return "float";
	}
	if (datatype == DT_DOUBLE) {
	    return "double";
	}

	throw new LoggedException("Header does not have a supported Camino data type");
    }


    private void setDataTypeToSigned() {
        
	if (datatype == DT_UNSIGNED_CHAR) {
	    datatype = DT_SIGNED_CHAR;
	}
	else if (datatype == DT_UNSIGNED_SHORT) {
            datatype = DT_SIGNED_SHORT;
	}
	else if (datatype == DT_UNSIGNED_INT) {
            datatype = DT_SIGNED_INT;
	}
        
    }



    /**
     *
     * @param hdrFile the path to the header file, with or without the ".hdr". 
     *
     * @return The image file corresponding to the header. Looks for .img and root.img.gz, 
     * where root is <code>hdrFile</code> without the ".hdr" extension.
     * If neither root.img nor root.img.gz can be found, the method returns <code>null</code>. 
     *
     */
    public static String getIMG_File(String hdrFile) {

	String root = hdrFile;

	if (root.endsWith(".hdr")) {	    
	    root = root.substring(0, root.length() - 4);
	}

	String imgFile = root + ".img";

	if ( !(new File(imgFile).exists()) ) {
	    imgFile = root + ".img.gz";
	    
	    if ( !(new File(imgFile).exists()) ) {
		return null;
	    }
	    
	}

	return imgFile;

    }



    /**
     * Gets the image root, such that root.hdr and root.[img | img.gz] exist. If this is not possible, then
     * the method returns null.
     * <p>
     * If fileName exists and ends in [.hdr | .img | .img.gz], the method returns the file root, where
     * both root.hdr and root.[img | img.gz] exist and fileName is equal to root.[hdr | img | img.gz].
     * <BR>
     * If fileName does not end in [.hdr | .img | .img.gz], the method returns fileName if fileName.hdr 
     * and fileName.[img | img.gz] exist.
     *
     * @param fileName path to a file, corresponding to an Analyze image.
     *
     * @return root where root.hdr is the header to the image and getIMG_File(root) returns the 
     * name of the image data file. If the header / image pair cannot be found, the method 
     * returns <code>null</code>.
     */
    public static String getImageRoot(String fileName) {

	String root = null;

        if ( (fileName.endsWith(".hdr")) ) {
            root = fileName.substring(0, fileName.length() - 4);
        }
        else if ( (fileName.endsWith(".img")) ) {
            root = fileName.substring(0, fileName.length() - 4);
        }
	else if ( (fileName.endsWith(".img.gz")) ) {
            root = fileName.substring(0, fileName.length() - 7);
        }
	else {
	    root = fileName;
	}

	if (new File(root + ".hdr").exists()) {
	    if ( new File(root + ".img").exists() || new File(root + ".img.gz").exists() ) {
		return root;
	    }
	    else { 
		return null;
	    }
	}
	else {
	    return null;
	}
    }
    

   

    // ImageHeader implementation methods

    public int xDataDim() {
	return width;
    }

    public int yDataDim() {
	return height; 
    }

    public int zDataDim() {
	return depth;
    }

    public int[] getDataDims() {
	return new int[] {width, height, depth};
    }

    public double xVoxelDim() {
	return pixelWidth;
    }

    public double yVoxelDim() {
	return pixelHeight;
    }

    public double zVoxelDim() {
	return pixelDepth;
    }

    public double[] getVoxelDims() {
	return new double[] {pixelWidth, pixelHeight, pixelDepth};
    }

    public double[] getOrigin() {
	return new double[] {centre[0], centre[1], centre[2]};
    }


    /** 
     * @return nImages, or 1 if nImages == 0.
     *
     */
    public int components() {
	return nImages > 0 ? nImages : 1;
    }



    /**
     * Gets the transformation matrix R that provides a transformation from a = (i, j, k, 1) in 
     * voxel space to b = (x, y, z, 1) in physical space.
     * <p>
     * The matrix is simply a scaling by pixdim and a translation, defined by centre.
     *
     */
    public RealMatrix getVoxelToPhysicalTransform() {

        logger.warning("Using Analyze header, physical space transformation may be unreliable");


        // centre is indexed from (1,1,1) in Matlab fashion
        double[] trans = new double[3];

        trans[0] = -(centre[0] - 1) * pixelWidth;
        trans[1] = -(centre[1] - 1) * pixelHeight;
        trans[2] = -(centre[2] - 1) * pixelDepth;

        RealMatrix R = new RealMatrix(4,4);

        R.entries[0][0] = pixelWidth;
        R.entries[1][1] = pixelHeight;
        R.entries[2][2] = pixelDepth;

        R.entries[0][3] = trans[0];
        R.entries[1][3] = trans[1];
        R.entries[2][3] = trans[2];
        R.entries[3][3] = 1.0;


        return R;
     
   
    }



    /**
     * Gets a data source from a Analyze header. Use this method to get a data source
     * in order to read data directly. Searches for an image corresponding to the header file of
     * this object in the order .img, .img.gz. 
     * 
     *
     * @return a data source that returns the raw image data, re-ordered into voxel order if necessary. 
     * If the image is scalar, a VoxelOrderDataSource is returned, so you can limit the memory consumption
     * by reducing buffer size in ExternalDataSource.
     */
    public DataSource getImageDataSource() {

	if (dataFile == null) {
	    dataFile = getIMG_File(hdrFile);
	    
	    if (dataFile == null) {
		throw new DataSourceException("Can't find image data for " + hdrFile);
	    }
	}

	// some programs set scale slope to 0
	double localScaleSlope = scaleSlope == 0.0 ? 1.0 : scaleSlope;

	double localScaleInt = scaleSlope == 0.0 ? 0.0 : scaleInter;

	if (components() == 1) {
	    // return voxel order data source if we have a 3D image 
	    // can save some memory this way if we are reading from a large number
	    // of input images, as is the case with reading DW data
	    if (localScaleInt != 0.0 || localScaleSlope != 1.0) {
		return new VoxelOrderScaledDataSource(dataFile, 1, caminoDataTypeString(), intelByteOrder, 
						      (int)Math.abs(vox_offset), localScaleSlope, localScaleInt);
	    }
	    else {
		return new VoxelOrderDataSource(dataFile, 1, caminoDataTypeString(), intelByteOrder, 
						(int)Math.abs(vox_offset));
	    }
	}
	
	// read 4D data, scale, then create a data source

	int voxels = width * height * depth;

	double[][] data = new double[voxels][components()];

	VoxelOrderDataSource vo = 
	    new VoxelOrderDataSource(dataFile, 1, caminoDataTypeString(), intelByteOrder, (int)Math.abs(vox_offset));


	for (int n = 0; n < components(); n++) {

	    if (vox_offset < 0.0f) {
		// offset applied to each volume

		int byteOffset = (bitpix / 8) * n * voxels + (n + 1) * (int)Math.abs(vox_offset);

		vo = new VoxelOrderDataSource(dataFile, 1, caminoDataTypeString(), intelByteOrder, byteOffset);
	    }

	    
	    for (int v = 0; v < voxels; v++) {
		data[v][n] = vo.nextVoxel()[0] * localScaleSlope + localScaleInt;
	    }
	}

	ScannerOrderDataSource source = new ScannerOrderDataSource(data);

	return source;
    }


    /**
     * Reads a 4D Analyze volume and returns the result as double.
     *
     * @return a 4D array with dimensions [width][height][depth][components()].
     */
    public double[][][][] readVolumeData() {
	
	// check that image file exists
	if (dataFile == null) {
	    dataFile = getIMG_File(hdrFile);
	    
	    if (dataFile == null) {
		throw new DataSourceException("Can't find image data for " + hdrFile);
	    }
	}
	
	// some programs set scale slope to 0
	double localScaleSlope = scaleSlope == 0.0 ? 1.0 : scaleSlope;
	
	double localScaleInt = scaleSlope == 0.0 ? 0.0 : scaleInter;


	if (dataFile == null) {
	    throw new LoggedException("Can't find image data for " + hdrFile);
	}
	
	// do it the long way rather than calling getImageDataSource, in order to save memory

        try {

            double[][][][] vol = null;
            
            int xDataDim = width;
            int yDataDim = height;
            int zDataDim = depth;
            
            vol = new double[xDataDim][yDataDim][zDataDim][components()];


	    VoxelOrderDataSource ds = new VoxelOrderDataSource(dataFile, 1, caminoDataTypeString(), 
							       intelByteOrder, (int)Math.abs(vox_offset));
    
	    for (int n = 0; n < components(); n++) {
		
		if (vox_offset < 0.0f) {
		    // offset applied to each volume

		    int byteOffset = (bitpix / 8) * n * width * height * depth + (n + 1) * (int)Math.abs(vox_offset);

		    ds = new VoxelOrderDataSource(dataFile, 1, caminoDataTypeString(), 
						  intelByteOrder, byteOffset);
		}
		

		for (int k = 0; k < zDataDim; k++) {
		    for (int j = 0; j < yDataDim; j++) {
			for (int i = 0; i < xDataDim; i++) {
			    
			    vol[i][j][k][n] = ds.nextVoxel()[0] * localScaleSlope + localScaleInt; 
			}
		    }
		}
	    }
	    
            return vol;
            
        } catch (DataSourceException e) {
            throw new LoggedException(e);
        }
        
    }

   
    /**
     * Get the header file name associated with the header object.
     *
     */
    public String getHeaderFilename() {
        return hdrFile;
    }


    /**
     * Get the data file name (.img or .img.gz) associated with the header object. 
     *
     */
    public String getDataFilename() {
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
        
        if (data.length != width || data[0].length != height || data[0][0].length != depth) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + width +                                      " " + height + " " + depth);
        }

        AnalyzeHeader ah = new AnalyzeHeader(this, fileRoot);

        ah.nImages = 1;

        ah.setDataTypeToSigned();
        
        double[][][][] fourD = new double[width][height][depth][1];

        double min = 0.0;
        double max = 0.0;
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int k = 0; k < depth; k++) { 
                    fourD[i][j][k][0] = data[i][j][k];

                    if (data[i][j][k] < min) {
                        min = data[i][j][k];
                    }
                    if (data[i][j][k] > max) {
                        max = data[i][j][k];
                    }
                }
            }
        }

        ah.glmin = (int)Math.round(min);
        ah.glmax = (int)Math.round(max);

        ah.scaleSlope = 1.0f;
        ah.scaleInter = 0.0f;

        
        try {
            ah.writeImage(fourD);
        }
        catch(IOException e) {
            throw new LoggedException(e.getMessage());
        }

        return ah;

        
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

        if (data.length != width || data[0].length != height || data[0][0].length != depth) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + width +                                      " " + height + " " + depth);
        }

        AnalyzeHeader ah = new AnalyzeHeader(this, fileRoot);

        ah.nImages = (short)data[0][0][0].length;

        ah.setDataTypeToSigned();
        
        ah.glmin = 0;
        ah.glmax = 0;

        ah.scaleSlope = 1.0f;
        ah.scaleInter = 0.0f;

        try {
            ah.writeImage(data);
        }
        catch(IOException e) {
            throw new LoggedException(e.getMessage());
        }

        return ah;

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

        if (red.length != width || red[0].length != height || red[0][0].length != depth) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + width +                                      " " + height + " " + depth);
        }


        AnalyzeHeader ah = new AnalyzeHeader(this, fileRoot);

	FileOutputStream fout = null;
	DataOutputStream dout = null;
        
        try {
            
            fout = new FileOutputStream(ah.dataFile);
            
            if (ah.gzip) {
                dout = new DataOutputStream(new BufferedOutputStream(new GZIPOutputStream(fout), BUFFERSIZE));
            }
            else {
                dout = new DataOutputStream(new BufferedOutputStream(fout, BUFFERSIZE));
            }
            
            ah.datatype = DT_RGB;
            
            ah.bitpix = 24;
            
            ah.scaleSlope = 1.0f;
            ah.scaleInter = 0.0f;
            
            ah.glmin = 0;
            ah.glmax = 0;
            
            ah.writeHeader();
            
            for (int i = 0; i < width; i++) {
                for (int j = 0; j < height; j++) {
                    for (int k = 0; k < depth; k++) { 
                        dout.write(red[i][j][k]);
                        dout.write(green[i][j][k]);
                        dout.write(blue[i][j][k]);
                    }
                }
            }

            
            dout.close();

        }
        catch(IOException e) {
            throw new LoggedException(e.getMessage());
        }
        

        return ah;

    }



    /**
     * 
     * Use the current header to write a tensor image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, with fields altered where necessary 
     * for writing tensor data. 
     *
     * @param data will be checked against header data dimensions. The data array should be in the 
     * same voxel space as this image, to ensure a correct definition of physical space. data[i][j][k] 
     * contains six components in upper-triangular order. The data will be written to disk in  
     * lower triangular order.
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings. 
     *
     * @return the image header associated with the output data set
     *
     */
    public ImageHeader writeTensorImage(double[][][][] data, String fileRoot) {

        if (data.length != width || data[0].length != height || data[0][0].length != depth) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + width +                                      " " + height + " " + depth);
        }
        
        AnalyzeHeader ah = new AnalyzeHeader(this, fileRoot);

        ah.nImages = 6;

        ah.setDataTypeToSigned();
        
        ah.glmin = 0;
        ah.glmax = 0;

        ah.scaleSlope = 1.0f;
        ah.scaleInter = 0.0f;

        double[][][][] lowerTriangular = new double[width][height][depth][6];
        

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                for (int k = 0; k < depth; k++) { 
                    
                    lowerTriangular[i][j][k][0] = data[i][j][k][0];
                    lowerTriangular[i][j][k][1] = data[i][j][k][1];
                    lowerTriangular[i][j][k][2] = data[i][j][k][3];
                    lowerTriangular[i][j][k][3] = data[i][j][k][2];
                    lowerTriangular[i][j][k][4] = data[i][j][k][4];
                    lowerTriangular[i][j][k][5] = data[i][j][k][5];
                    
                }
            }
        }

        try {
            ah.writeImage(lowerTriangular);
        }
        catch(IOException e) {
            throw new LoggedException(e.getMessage());
        }

        return ah;

    }


  

    /**
     * Set the data type of the image.
     *
     * @param a Camino data type string.
     *
     */
    public void setDataType(String type) {
   	
	if (type.equals("ubyte")) { 
            datatype = DT_UNSIGNED_CHAR;
            bitpix = 8;
	}
	else if (type.equals("byte")) { 
            datatype = DT_SIGNED_CHAR;
            bitpix = 8;
	}
	else if (type.equals("short")) { 
            datatype = DT_SIGNED_SHORT;
            bitpix = 16;
	}
	else if (type.equals("int")) {  
            datatype = DT_SIGNED_INT;
            bitpix = 16;
	}
	else if (type.equals("ushort")) {  
            datatype = DT_UNSIGNED_SHORT;
            bitpix = 16;
	}
	else if (type.equals("uint")) { 
            datatype = DT_UNSIGNED_INT;
            bitpix = 32;
	}
	else if (type.equals("float")) { 
            datatype = DT_FLOAT;
            bitpix = 32;
	}
	else if (type.equals("double")) { 
            datatype = DT_DOUBLE;
            bitpix = 64;
	}
        else {
            throw new LoggedException("Header does not have a supported Camino data type");
        }
    }





    /**
     * Sets gzip compression. This may alter the file name(s) of the image.
     *
     */
    public void setGzip(boolean gz) {
        if (gz) {
            if (gzip) {
                // nothing to do
            }
            else {
                dataFile = dataFile + ".gz";
                gzip = true;
            }
        }
        else {
            
            if (!gzip) {

            }
            else {
                dataFile = dataFile.substring(0, dataFile.length() - 4);
                gzip = false;
            }

        }
    }



    public boolean lowerTriangularSymmMatrix() { 
        return true;
    }
    
}


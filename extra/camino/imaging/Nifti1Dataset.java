package imaging;

import data.*;
import misc.*;
import numerics.*;
import tools.*;

import java.io.*;
import java.nio.*;
import java.text.*;
import java.util.*;
import java.util.logging.*;
import java.util.zip.*;



/** 
 * This class is adapted from the NIFTI Java code provided by Kate Fissell at
 * the University of Pittsburgh. It supports scalar and vector (meaning any 
 * multi-valued data set) images. 
 *
 *
 *
 * @author Philip Cook
 * @version $Id$
 */
public final class Nifti1Dataset extends ImageHeader {

    private static Logger logger = Logger.getLogger("camino.imaging.Nifti1Dataset");

    //////////////////////////////////////////////////////////////////
    //
    //		Nifti-1 Defines        
    //
    //////////////////////////////////////////////////////////////////
    public static final String ANZ_HDR_EXT = ".hdr";
    public static final String ANZ_DAT_EXT = ".img";
    public static final String NI1_EXT = ".nii";
    public static final String GZIP_EXT = ".gz";
    public static final int ANZ_HDR_SIZE = 348;
    public static final long NII_HDR_SIZE = 352;
    public static final int  EXT_KEY_SIZE = 8;   // esize+ecode
    public static final String NII_MAGIC_STRING = "n+1";
    public static final String ANZ_MAGIC_STRING = "ni1";


    //////////////////////////////////////////////////////////////////
    //
    //		Nifti-1 Codes from C nifti1.h
    //
    //////////////////////////////////////////////////////////////////

    // intent codes for field intent_code
    public static final short NIFTI_INTENT_NONE        = 0;
    public static final short NIFTI_INTENT_CORREL      = 2;
    public static final short NIFTI_INTENT_TTEST       = 3;
    public static final short NIFTI_INTENT_FTEST       = 4;
    public static final short NIFTI_INTENT_ZSCORE      = 5;
    public static final short NIFTI_INTENT_CHISQ       = 6;
    public static final short NIFTI_INTENT_BETA        = 7;
    public static final short NIFTI_INTENT_BINOM       = 8;
    public static final short NIFTI_INTENT_GAMMA       = 9;
    public static final short NIFTI_INTENT_POISSON    = 10;
    public static final short NIFTI_INTENT_NORMAL     = 11;
    public static final short NIFTI_INTENT_FTEST_NONC = 12;
    public static final short NIFTI_INTENT_CHISQ_NONC = 13;
    public static final short NIFTI_INTENT_LOGISTIC   = 14;
    public static final short NIFTI_INTENT_LAPLACE    = 15;
    public static final short NIFTI_INTENT_UNIFORM    = 16;
    public static final short NIFTI_INTENT_TTEST_NONC = 17;
    public static final short NIFTI_INTENT_WEIBULL    = 18;
    public static final short NIFTI_INTENT_CHI        = 19;
    public static final short NIFTI_INTENT_INVGAUSS   = 20;
    public static final short NIFTI_INTENT_EXTVAL     = 21;
    public static final short NIFTI_INTENT_PVAL       = 22;
    public static final short NIFTI_INTENT_ESTIMATE  = 1001;
    public static final short NIFTI_INTENT_LABEL     = 1002;
    public static final short NIFTI_INTENT_NEURONAME = 1003;
    public static final short NIFTI_INTENT_GENMATRIX = 1004;
    public static final short NIFTI_INTENT_SYMMATRIX = 1005;
    public static final short NIFTI_INTENT_DISPVECT  = 1006;
    public static final short NIFTI_INTENT_VECTOR    = 1007;
    public static final short NIFTI_INTENT_POINTSET  = 1008;
    public static final short NIFTI_INTENT_TRIANGLE  = 1009;
    public static final short NIFTI_INTENT_QUATERNION =  1010;
    public static final short NIFTI_FIRST_STATCODE	= 2;
    public static final short NIFTI_LAST_STATCODE		= 22;


    // datatype codes for field datatype
    public static final short DT_NONE                    = 0;
    public static final short DT_BINARY                  = 1;
    public static final short NIFTI_TYPE_UINT8           = 2;
    public static final short NIFTI_TYPE_INT16           = 4;
    public static final short NIFTI_TYPE_INT32           = 8;
    public static final short NIFTI_TYPE_FLOAT32        = 16;
    public static final short NIFTI_TYPE_COMPLEX64      = 32;
    public static final short NIFTI_TYPE_FLOAT64        = 64;
    public static final short NIFTI_TYPE_RGB24         = 128;
    public static final short DT_ALL                   = 255;
    public static final short NIFTI_TYPE_INT8          = 256;
    public static final short NIFTI_TYPE_UINT16        = 512;
    public static final short NIFTI_TYPE_UINT32        = 768;
    public static final short NIFTI_TYPE_INT64        = 1024;
    public static final short NIFTI_TYPE_UINT64       = 1280;
    public static final short NIFTI_TYPE_FLOAT128     = 1536;
    public static final short NIFTI_TYPE_COMPLEX128   = 1792;
    public static final short NIFTI_TYPE_RGBA32       = 2304;
    public static final short NIFTI_TYPE_COMPLEX256   = 2048;

    // units codes for xyzt_units
    public static final short NIFTI_UNITS_UNKNOWN = 0;
    public static final short NIFTI_UNITS_METER   = 1;
    public static final short NIFTI_UNITS_MM      = 2;
    public static final short NIFTI_UNITS_MICRON  = 3;
    public static final short NIFTI_UNITS_SEC     = 8;
    public static final short NIFTI_UNITS_MSEC   = 16;
    public static final short NIFTI_UNITS_USEC   = 24;
    public static final short NIFTI_UNITS_HZ     = 32;
    public static final short NIFTI_UNITS_PPM    = 40;

    // slice order codes for slice_code
    public static final short NIFTI_SLICE_SEQ_INC =  1;
    public static final short NIFTI_SLICE_SEQ_DEC =  2;
    public static final short NIFTI_SLICE_ALT_INC =  3;
    public static final short NIFTI_SLICE_ALT_DEC =  4;

    // codes for qform_code sform_code
    public static final short NIFTI_XFORM_UNKNOWN      = 0;
    public static final short NIFTI_XFORM_SCANNER_ANAT = 1;
    public static final short NIFTI_XFORM_ALIGNED_ANAT = 2;
    public static final short NIFTI_XFORM_TALAIRACH    = 3;
    public static final short NIFTI_XFORM_MNI_152      = 4;


    // variables derived from header info and filenames
    private String ds_hdrname;	// file name for header
    private String ds_datname;	// file name for data
    private boolean ds_is_nii;	// does dataset use single file .nii 
    private boolean gzip;       // does dataset use gzip
    private boolean big_endian;	// does hdr appear to have BE format
    private short XDIM,YDIM,ZDIM,TDIM,DIM5,DIM6,DIM7;	// from dim[] field
    private short freq_dim,phase_dim,slice_dim;  // unpack dim_info
    private short xyz_unit_code, t_unit_code;	// unpack xyzt_units;
    private short qfac;				// unpack pixdim[0]
    private Vector<Object> extensions_list;		// vector of size/code pairs for ext.
    private Vector<Object> extension_blobs;		// vector of extension data


    // variables for fields in the nifti header 
    private int		sizeof_hdr;	// must be 348 bytes
    private StringBuffer	data_type_string;	// 10 char UNUSED
    private StringBuffer	db_name;	// 18 char UNUSED
    private int		extents;	// UNUSED
    private short		session_error;	// UNUSED
    private StringBuffer	regular;	// 1 char UNUSED
    private StringBuffer	dim_info;	// 1 char MRI slice ordering
    private short		dim[];	// data array dimensions (8 shorts)
    private float		intent[];	// intents p1 p2 p3
    private short		intent_code;	// nifti intent code for dataset
    private short		datatype;	// datatype of image blob
    private short		bitpix;		// #bits per voxel
    private short		slice_start;	// first slice index
    private float		pixdim[];	// grid spacings
    private float		vox_offset;	// offset to data blob in .nii file
    private float		scl_slope;	// data scaling: slope
    private float		scl_inter;	// data scaling: intercept
    private short		slice_end;	// last slice index
    private byte		slice_code;	// slice timing order
    private byte		xyzt_units;	// units of pixdim[1-4]
    private float		cal_max;	// max display intensity
    private float		cal_min;	// min display intensity
    private float		slice_duration;	// time to acq. 1 slice
    private float		toffset;	// time axis shift
    private int		glmax;		// UNUSED
    private int		glmin;		// UNUSED
    private StringBuffer	descrip;	// 80 char any text you'd like (comment)
    private StringBuffer	aux_file;	// 24 char auxiliary file name
    private short		qform_code;	// code for quat. transform
    private short		sform_code;	// code for affine transform
    private float		quatern[];	// 3 floats Quaternion b,c,d params
    private float		qoffset[];	// 3 floats Quaternion x,y,z shift
    private float		srow_x[];	// 4 floats 1st row affine xform
    private float		srow_y[];	// 4 floats 2nd row affine xform
    private float		srow_z[];	// 4 floats 3rd row affine xform
    private StringBuffer	intent_name;	// 16 char name/meaning/id for data
    private StringBuffer	magic;		// 4 char id must be "ni1\0" or "n+1\0"
    private byte		extension[];	// 4 byte array, byte 0 is 0/1 indicating extensions or not




    //////////////////////////////////////////////////////////////////
    /**
     * Constructor for a dataset existing on disk.
     *
     * @param name - name of the nifti-1 dataset on disk.
     *
     */
    public Nifti1Dataset(String name) {

	setDefaults();

        setFilename(name);

        try { 
            readHeader();
        }
        catch (IOException e) { 
            throw new LoggedException(e);
        }

        // Convert a 2D image to 3D with dimension 1
        if (dim[0] == 2) {
            
            logger.warning("Attempting to read 2D image as 3D image with ZDIM == 1");
            
            dim[0] = 3;
            dim[2] = 1;
            
            // probably wasn't a good idea to have duplicate dimensions, but it was like that when I got here
            ZDIM = 1;

            // Attempt to make a reasonable estimate
            pixdim[2] = (pixdim[0] + pixdim[1]) / 2.0f;
            
        }

	return;
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Constructor for creation of a new dataset. 
     *
     */
    public Nifti1Dataset() {

	setDefaults();

	return;
    }



    /**
     * Copy constructor for creation of a new dataset.
     *
     * @param nds a source header. All data will be copied from here, except
     * for the file name.
     *
     * @param fileRoot the root name of the new file. Extension is determined from the settings in <code>nds</code>.
     *
     */
    public Nifti1Dataset(Nifti1Dataset nds, String fileRoot) {

        setDefaults();

        copyHeader(nds);

        setFilename(fileRoot, ds_is_nii, gzip);

	return;
    }




    //////////////////////////////////////////////////////////////////
    /**
     * Read header information into memory
     * @exception IOException 
     * @exception FileNotFoundException
     */
    private void readHeader() throws IOException, FileNotFoundException {

	DataInputStream dis;
	EndianCorrectInputStream ecs;
	short s, ss[];
	byte bb[];
	int i;

	if (ds_hdrname.endsWith(".gz")) {
	    dis = new DataInputStream(new GZIPInputStream(new FileInputStream(ds_hdrname)));
            gzip = true;
        }
        else {
	    dis = new DataInputStream(new FileInputStream(ds_hdrname));
            gzip = false;
        }


	try {

	    ///// first, read dim[0] to get endian-ness
	    dis.skipBytes(40);  
	    s = dis.readShort();
	    dis.close();
	    if ((s < 1) || (s > 7))
		big_endian = false;
	    else
		big_endian = true;


	    ///// get input stream that will flip bytes if necessary 
	    if (gzip)
		ecs = new EndianCorrectInputStream(new GZIPInputStream(new FileInputStream(ds_hdrname)),big_endian);
	    else
		ecs = new EndianCorrectInputStream(ds_hdrname,big_endian);

	    sizeof_hdr = ecs.readIntCorrect();

	    bb = new byte[10];
	    ecs.readFully(bb,0,10);
	    data_type_string = new StringBuffer(new String(bb));

	    bb = new byte[18];
	    ecs.readFully(bb,0,18);
	    db_name = new StringBuffer(new String(bb));

	    extents = ecs.readIntCorrect();

	    session_error = ecs.readShortCorrect();
			
	    regular = new StringBuffer();
	    regular.append((char)(ecs.readUnsignedByte()));
		
	    dim_info = new StringBuffer();
	    dim_info.append((char)(ecs.readUnsignedByte()));
	    ss = unpackDimInfo((int)dim_info.charAt(0));
	    freq_dim = ss[0];
	    phase_dim = ss[1];
	    slice_dim = ss[2];

	    for (i=0; i<8; i++)
		dim[i] = ecs.readShortCorrect();
	    if (dim[0] > 0)
		XDIM = dim[1];
	    if (dim[0] > 1)
		YDIM = dim[2];
	    if (dim[0] > 2)
		ZDIM = dim[3];
	    if (dim[0] > 3)
		TDIM = dim[4];

	    for (i=0; i<3; i++)
		intent[i] = ecs.readFloatCorrect();

	    intent_code = ecs.readShortCorrect();

	    datatype = ecs.readShortCorrect();

	    bitpix = ecs.readShortCorrect();

	    slice_start = ecs.readShortCorrect();

	    for (i=0; i<8; i++)
		pixdim[i] = ecs.readFloatCorrect();
	    qfac = (short) Math.floor((double)(pixdim[0]));
		
	    vox_offset = ecs.readFloatCorrect();

	    scl_slope = ecs.readFloatCorrect();
	    scl_inter = ecs.readFloatCorrect();

	    slice_end = ecs.readShortCorrect();

	    slice_code = (byte) ecs.readUnsignedByte();

	    xyzt_units = (byte) ecs.readUnsignedByte();
	    ss = unpackUnits((int)xyzt_units);
	    xyz_unit_code =  ss[0];
	    t_unit_code =  ss[1];

	    cal_max = ecs.readFloatCorrect();
	    cal_min = ecs.readFloatCorrect();

	    slice_duration = ecs.readFloatCorrect();

	    toffset = ecs.readFloatCorrect();
		
	    glmax = ecs.readIntCorrect();
	    glmin = ecs.readIntCorrect();

	    bb = new byte[80];
	    ecs.readFully(bb,0,80);
	    descrip = new StringBuffer(new String(bb));

	    bb = new byte[24];
	    ecs.readFully(bb,0,24);
	    aux_file = new StringBuffer(new String(bb));

	    qform_code = ecs.readShortCorrect();
	    sform_code = ecs.readShortCorrect();

	    for (i=0; i<3; i++)
		quatern[i] = ecs.readFloatCorrect();
	    for (i=0; i<3; i++)
		qoffset[i] = ecs.readFloatCorrect();

	    for (i=0; i<4; i++)
		srow_x[i] = ecs.readFloatCorrect();
	    for (i=0; i<4; i++)
		srow_y[i] = ecs.readFloatCorrect();
	    for (i=0; i<4; i++)
		srow_z[i] = ecs.readFloatCorrect();


	    bb = new byte[16];
	    ecs.readFully(bb,0,16);
	    intent_name = new StringBuffer(new String(bb));

	    bb = new byte[4];
	    ecs.readFully(bb,0,4);
	    magic = new StringBuffer(new String(bb));

	}
	catch (IOException ex) {
	    throw new IOException("Error: unable to read header file "+ds_hdrname+": "+ex.getMessage());
	}
		

	/////// Read possible extensions
	if (ds_is_nii) 
	    readNiiExt(ecs);
	else
	    readNp1Ext(ecs);


	ecs.close();

	return;	
    }

    ////////////////////////////////////////////////////////////////////
    //
    // Copy all in memory header field settings from datset A to this dataset
    // Extension data not set, fields set to no extension
    //
    ////////////////////////////////////////////////////////////////////
    private void copyHeader(Nifti1Dataset A) {

	int i;

        ds_hdrname = 	new String(A.ds_hdrname);  
	ds_datname = 	new String(A.ds_datname);
	ds_is_nii = 	A.ds_is_nii;
	big_endian = 	A.big_endian;
	sizeof_hdr = 	A.sizeof_hdr;
	data_type_string = new StringBuffer(A.data_type_string.toString());
	db_name = new StringBuffer(A.db_name.toString());
	extents = 	A.extents;
	session_error =	A.session_error;
	regular = new StringBuffer(A.regular.toString());
	dim_info = new StringBuffer(A.dim_info.toString());
	freq_dim=A.freq_dim; 
	phase_dim=A.phase_dim; 
	slice_dim=A.slice_dim;
	for (i=0; i<8; i++)
	    dim[i] = A.dim[i];
	XDIM=A.XDIM; YDIM=A.YDIM; ZDIM=A.ZDIM; TDIM=A.TDIM;
	DIM5=A.DIM5; DIM6=A.DIM6; DIM7=A.DIM7;
	for (i=0; i<3; i++)
	    intent[i] = A.intent[i];
	intent_code = A.intent_code;
	datatype = A.datatype;
	bitpix = A.bitpix;	
	slice_start = A.slice_start;
	qfac = 1;
	for (i=0; i<8; i++)
	    pixdim[i] = A.pixdim[i];

        qfac = (short) Math.floor((double)(pixdim[0]));
	
	vox_offset = A.vox_offset;
	scl_slope = A.scl_slope;
	scl_inter = A.scl_inter;
	slice_end = A.slice_end;
	slice_code = A.slice_code;
	xyzt_units = A.xyzt_units;
	xyz_unit_code = A.xyz_unit_code;
	t_unit_code = A.t_unit_code;

	cal_max = A.cal_max;
	cal_min = A.cal_min;
	slice_duration = A.slice_duration;
	toffset = A.toffset;
	glmax = A.glmax;
	glmin = A.glmin;

	descrip = new StringBuffer(A.descrip.toString());
	aux_file = new StringBuffer(A.aux_file.toString());

	qform_code = A.qform_code;
	sform_code = A.sform_code;

	for (i=0; i<3; i++) {
	    quatern[i] = A.quatern[i];
	    qoffset[i] = A.qoffset[i];
	}

	for (i=0; i<4; i++) {
	    srow_x[i] = A.srow_x[i];
	    srow_y[i] = A.srow_y[i];
	    srow_z[i] = A.srow_z[i];
	}

	intent_name = new StringBuffer(A.intent_name.toString());

	magic = new StringBuffer(A.magic.toString());

	for (i=0; i<4; i++)
	    extension[i] = (byte)0;

        gzip = A.gzip;

	return;
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Read extension data from nii (1 file) dataset
     * @param ecs is an InputStream open and pointing to begining of 
     * 	extensions array
     * @exception IOException 
     */
    private void readNiiExt(EndianCorrectInputStream ecs) throws IOException {

	int[] size_code;
	int start_addr;
	byte[] eblob;
        byte[] extension = new byte[4];

	// read 4 ext bytes if it is nii, bytes 348-351 must be there
	try {
	    int bytesRead = ecs.read(extension);
            
            if (bytesRead < 4) {
                throw new IOException("Error: i/o error reading extension bytes on header file " + ds_hdrname);
            }
	}
	catch (IOException ex) {
	    throw new IOException("Error: i/o error reading extension bytes on header file "+ds_hdrname+": "+ex.getMessage());
	}
		

	/// jump thru extensions getting sizes and codes
	size_code = new int[2];

	if ( extension[0] == (byte) 0 ) {
            return;
        }

        start_addr=ANZ_HDR_SIZE+4;
        
        size_code[0] = 0;
        size_code[1] = 0;
        
	    /// in nii files, vox_offset is end of hdr/ext,
	    /// beginning of data.  
        while (start_addr < (int) vox_offset) {
            try {
                size_code = new int[2];
                size_code[0] = ecs.readIntCorrect();
                size_code[1] = ecs.readIntCorrect();
                eblob = new byte[size_code[0]-EXT_KEY_SIZE];
                ecs.readFully(eblob,0,size_code[0]-EXT_KEY_SIZE);
                extension_blobs.add(eblob);
            }
            catch (IOException ex) {
                printHeader();
                throw new EOFException("Error: i/o error reading extension data for extension "+(extensions_list.size()+1)+" on header file "+ds_hdrname+": "+ex.getMessage());
            }
            
            extensions_list.add(size_code);
            start_addr += (size_code[0]);
            
            // check if extensions appeared to overrun data blob
            // when extensions are done, start_addr should == vox_offset
            if (start_addr > (int) vox_offset) {
                printHeader();
                throw new IOException("Error: Data  for extension "+(extensions_list.size())+" on header file "+ds_hdrname+" appears to overrun start of image data.");
            }
        } // while not yet at data blob
    
    }

    

    //////////////////////////////////////////////////////////////////
    ////
    /*
     * Read extension data from n+1 (2 file) dataset
     * @param ecs is an InputStream open and pointing to begining of 
     * 	extensions array
     * @exception IOException 
     */
    private void readNp1Ext(EndianCorrectInputStream ecs) throws IOException, EOFException {

	int size_code[];
	byte eblob[];


	// read 4 ext bytes if it is n+1, bytes 348-351 do NOT
	// need to be there
	try {
	    ecs.readFully(extension,0,4);
	}
	catch (EOFException ex) {
	    return;
	}
	catch (IOException ex) {
	    throw new IOException("Error: i/o error reading extension bytes on header file "+ds_hdrname+": "+ex.getMessage());
	}
		

	/// jump thru extensions getting sizes and codes
	size_code = new int[2];
	if ( extension[0] != (byte) 0 ) {

	    size_code[0] = 0;
	    size_code[1] = 0;

	    /// in np1 files, read to end of hdr file looking
	    /// for extensions
	    while (true) {
		try {
		    size_code = new int[2];
		    size_code[0] = ecs.readIntCorrect();
		    size_code[1] = ecs.readIntCorrect();
		    eblob = new byte[size_code[0]-EXT_KEY_SIZE];
		    ecs.readFully(eblob,0,size_code[0]-EXT_KEY_SIZE);
		    extension_blobs.add(eblob);
		}
		catch (EOFException ex) {
		    return;
		}
		catch (IOException ex) {
		    throw new EOFException("Error: i/o error reading extension data for extension "+(extensions_list.size()+1)+" on header file "+ds_hdrname+": "+ex.getMessage());
		}

		extensions_list.add(size_code);
	    } // while not yet at EOF

	}	// if there are extensions


	return;
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Get list of extensions and return it as nx2 array
     * @return nx2 array where n = # of extensions, array elem 0
     * is the size in bytes of that extension and array elem 1 is
     * the extension code.
     */
    public int[][] getExtensionsList() {

	int n,i;
	int size_code[];
	int extlist[][];

	size_code = new int[2];
	n = extensions_list.size();
	extlist = new int[n][2];
		

	for (i=0; i<n; i++) {
	    size_code = (int[]) extensions_list.get(i);
	    extlist[i][0] = size_code[0];
	    extlist[i][1] = size_code[1];
	}

	return(extlist);
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Remove an extension from a header
     * @param index number of the extension to remove (0 based)
     */
    public void removeExtension(int index) {

	int n;
	int size_code[] = new int[2];

	n = extensions_list.size();

	if (index >= n) {
	    System.out.println("\nERROR: could not remove extension "+index+1+" from "+ds_hdrname+". It only has "+n+" extensions.");
	    return;
	}

	// remove extension from lists
	size_code = (int[]) extensions_list.get(index);
	extensions_list.remove(index);
	extension_blobs.remove(index);

	// readjust vox_offset if necessary
	if (ds_is_nii)
	    vox_offset -= size_code[0];

	return;
    }

    //////////////////////////////////////////////////////////////////
    /**
     * Add an extension stored in a file to a header
     * @param code -- code identifying the extension
     * @param filename -- filename containing the extension.  The entire
     * file will be added as an extension
     */
    public void addExtension(int code, String filename) throws IOException {

	File f;
	long l;
	int size_code[] = new int[2];
	DataInputStream dis;
	byte b[];
	int i, il, pad;

	f = new File(filename);
	l = f.length();

	//// check length, compute padding
	//// 8bytes of size+code plus ext data must be mult. of 16
	if (l > Integer.MAX_VALUE) {
	    throw new IOException("Error: maximum extension size is "+Integer.MAX_VALUE+"bytes. "+filename+" is "+l+" bytes.");
	}
	il = (int)l;
	pad =  (il+EXT_KEY_SIZE)%16;
	if (pad != 0)
	    pad = 16-pad;
	///System.out.println("\next file size is "+l+", padding with "+pad);

	/// read the extension data from the file
	b = new byte[il+pad];
	try {
	    dis = new DataInputStream(new FileInputStream(filename));
	    dis.readFully(b,0,il);
	    dis.close();
	}
	catch(IOException ex) {
	    throw new IOException("Error reading extension data for "+ds_hdrname+" from file "+filename+". :"+ex.getMessage());
	}
	for(i=il; i<il+pad; i++)
	    b[i] = 0;


	/// well, if we got this far, I guess we really have to add it
	size_code[0] = il+pad+EXT_KEY_SIZE;
	size_code[1] = code;
	extensions_list.add(size_code);
	extension_blobs.add(b);
	extension[0] = 1;

	// update vox_offset for nii files
	if (ds_is_nii)
	    vox_offset += size_code[0];

	return;
    }



    /**
     * Write header information to disk.
     *
     * @throws IOException 
     *
     */
    public void writeHeader() throws IOException {

        byte[] header = headerToBytes();

        FileOutputStream fos = new FileOutputStream(ds_hdrname);

        OutputStream out;

        if (gzip) {
            out = new GZIPOutputStream(fos, 1024);
        }
        else {
            out = new BufferedOutputStream(fos, 1024);
        }

        out.write(header);

        out.close();
    }



    /**
     * @return the header as a byte array, with endianness set by the header
     *
     */
    private byte[] headerToBytes() throws IOException {

	EndianCorrectOutputStream ecs;
	ByteArrayOutputStream baos;
	short s, ss[];
	byte b, bb[], ext_blob[];
	int hsize;
	int i,n;
	int extlist[][];


	// header is 348 except nii and anz/hdr w/ extensions is 352
	hsize = ANZ_HDR_SIZE;
	if ( (ds_is_nii) || (extension[0] != 0) )
	    hsize += 4;

        baos = new ByteArrayOutputStream(hsize + EXT_KEY_SIZE);
	
        ecs = new EndianCorrectOutputStream(baos,big_endian);


        ecs.writeIntCorrect(sizeof_hdr);

        if (data_type_string.length() >= 10) {
            ecs.writeBytes(data_type_string.substring(0,10));
        }
        else {
            ecs.writeBytes(data_type_string.toString());
            for (i=0; i<(10-data_type_string.length()); i++)
                ecs.writeByte(0);
        }

        if (db_name.length() >= 18) {
            ecs.writeBytes(db_name.substring(0,18));
        }
        else {
            ecs.writeBytes(db_name.toString());
            for (i=0; i<(18-db_name.length()); i++)
                ecs.writeByte(0);
        }

        ecs.writeIntCorrect(extents);

        ecs.writeShortCorrect(session_error);

        ecs.writeByte((int) regular.charAt(0));

        b = packDimInfo(freq_dim, phase_dim, slice_dim);
        ecs.writeByte((int) b);

        for (i=0; i<8; i++)
            ecs.writeShortCorrect(dim[i]);

        for (i=0; i<3; i++)
            ecs.writeFloatCorrect(intent[i]);
			
        ecs.writeShortCorrect(intent_code);

        ecs.writeShortCorrect(datatype);

        ecs.writeShortCorrect(bitpix);

        ecs.writeShortCorrect(slice_start);
		
        for (i=0; i<8; i++)
            ecs.writeFloatCorrect(pixdim[i]);


        ecs.writeFloatCorrect(vox_offset);

        ecs.writeFloatCorrect(scl_slope);
        ecs.writeFloatCorrect(scl_inter);

        ecs.writeShortCorrect(slice_end);

        ecs.writeByte((int)slice_code);

        ecs.writeByte((int)packUnits(xyz_unit_code,t_unit_code));


        ecs.writeFloatCorrect(cal_max);
        ecs.writeFloatCorrect(cal_min);

        ecs.writeFloatCorrect(slice_duration);

        ecs.writeFloatCorrect(toffset);
		
        ecs.writeIntCorrect(glmax);
        ecs.writeIntCorrect(glmin);

        ecs.write(setStringSize(descrip,80),0,80);
        ecs.write(setStringSize(aux_file,24),0,24);


        ecs.writeShortCorrect(qform_code);
        ecs.writeShortCorrect(sform_code);

        for (i=0; i<3; i++)
            ecs.writeFloatCorrect(quatern[i]);
        for (i=0; i<3; i++)
            ecs.writeFloatCorrect(qoffset[i]);

        for (i=0; i<4; i++)
            ecs.writeFloatCorrect(srow_x[i]);
        for (i=0; i<4; i++)
            ecs.writeFloatCorrect(srow_y[i]);
        for (i=0; i<4; i++)
            ecs.writeFloatCorrect(srow_z[i]);


        ecs.write(setStringSize(intent_name,16),0,16);
        ecs.write(setStringSize(magic,4),0,4);


        // nii or anz/hdr w/ ext. gets 4 more
        if ( (ds_is_nii) || (extension[0] != 0) ) {
            for (i=0; i<4; i++)
                ecs.writeByte((int)extension[i]);
        }
            
        ////// extensions
        if (extension[0] != 0) {
                
            ecs = new EndianCorrectOutputStream(baos,big_endian);
            extlist = getExtensionsList();
            n = extlist.length;
            for(i=0; i<n; i++) {
                // write size, code
                ecs.writeIntCorrect(extlist[i][0]);
                ecs.writeIntCorrect(extlist[i][1]);

                ext_blob = (byte[]) extension_blobs.get(i);
                ecs.write(ext_blob,0,extlist[i][0]-EXT_KEY_SIZE);
            }
        }

        ecs.close();
        return baos.toByteArray();
            
    
            
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Print header information to standard out.
     */
    public void printHeader() {
	System.out.print(toString());
    }


    public String toString() {
	int i;
	int extlist[][], n;

	StringBuffer buffer = new StringBuffer();
	
	buffer.append("\n");
	buffer.append("Dataset header file:\t\t\t\t"+ds_hdrname+"\n");
	buffer.append("Dataset data file:\t\t\t\t"+ds_datname+"\n");
	buffer.append("Size of header:\t\t\t\t\t"+sizeof_hdr+"\n");
	buffer.append("File offset to data blob:\t\t\t"+vox_offset+"\n");

	buffer.append("Endianness:\t\t\t\t\t");
	if (big_endian)
	    buffer.append("big"+"\n");
	else
	    buffer.append("little"+"\n");

	buffer.append("Magic filetype string:\t\t\t\t"+magic+"\n");



	///// Dataset datatype, size, units
	buffer.append("Datatype:\t\t\t\t\t"+datatype+" ("+decodeDatatype(datatype)+")\n");
	buffer.append("Bits per voxel:\t\t\t\t\t"+bitpix+"\n");
	buffer.append("Scaling slope and intercept:\t\t\t"+scl_slope+" "+scl_inter+"\n");


	buffer.append("Dataset dimensions (Count, X,Y,Z,T...):\t\t");
	for (i=0; i<=dim[0]; i++)
	    buffer.append(dim[i]+" ");
	buffer.append("\n");

	buffer.append("Grid spacings (X,Y,Z,T,...):\t\t\t");
	for (i=1; i<=dim[0]; i++)
	    buffer.append(pixdim[i]+" ");
	buffer.append("\n");

	buffer.append("XYZ  units:\t\t\t\t\t"+xyz_unit_code+" ("+decodeUnits(xyz_unit_code)+")\n");
	buffer.append("T units:\t\t\t\t\t"+t_unit_code+" ("+decodeUnits(t_unit_code)+")\n");
	buffer.append("T offset:\t\t\t\t\t"+toffset+"\n");


	buffer.append("Intent parameters:\t\t\t\t");
	for (i=0; i<3; i++)
	    buffer.append(intent[i]+" ");
	buffer.append("\n");
	buffer.append("Intent code:\t\t\t\t\t"+intent_code+" ("+decodeIntent(intent_code)+")\n");

	buffer.append("Cal. (display) max/min:\t\t\t\t"+cal_max+" "+cal_min+"\n");


	///// Slice order/timing stuff
	buffer.append("Slice timing code:\t\t\t\t"+slice_code+" ("+decodeSliceOrder((short)slice_code)+")\n");
	buffer.append("MRI slice ordering (freq, phase, slice index):\t"+freq_dim+" "+phase_dim+" "+slice_dim+"\n");

	buffer.append("Start/end slice:\t\t\t\t"+slice_start+" "+slice_end+"\n");
	buffer.append("Slice duration:\t\t\t\t\t"+slice_duration+"\n");

	///// Orientation stuff
	buffer.append("Q factor:\t\t\t\t\t"+qfac+"\n");
	buffer.append("Qform transform code:\t\t\t\t"+qform_code+" ("+decodeXform(qform_code)+")\n");
	buffer.append("Quaternion b,c,d params:\t\t\t"+quatern[0]+" "+quatern[1]+" "+quatern[2]+"\n");
	buffer.append("Quaternion x,y,z shifts:\t\t\t"+qoffset[0]+" "+qoffset[1]+" "+qoffset[2]+"\n");

	buffer.append("Affine transform code:\t\t\t\t"+sform_code+" ("+decodeXform(sform_code)+")\n");
	buffer.append("1st row affine transform:\t\t\t");
	for (i=0; i<4; i++)
	    buffer.append(srow_x[i]+" ");
	buffer.append("\n");
	buffer.append("2nd row affine transform:\t\t\t");
	for (i=0; i<4; i++)
	    buffer.append(srow_y[i]+" ");
	buffer.append("\n");
	buffer.append("3rd row affine transform:\t\t\t");
	for (i=0; i<4; i++)
	    buffer.append(srow_z[i]+" ");
	buffer.append("\n");


	///// comment stuff
	buffer.append("Description:\t\t\t\t\t"+descrip+"\n");
	buffer.append("Intent name:\t\t\t\t\t"+intent_name+"\n");
	buffer.append("Auxiliary file:\t\t\t\t\t"+aux_file+"\n");
	buffer.append("Extension byte 1:\t\t\t\t\t"+(int)extension[0]+"\n");


	///// unused stuff
	buffer.append("\n\nUnused Fields\n");
	buffer.append("----------------------------------------------------------------------\n");
	buffer.append("Data type string:\t\t\t"+data_type_string+"\n");
	buffer.append("db_name:\t\t\t\t\t"+db_name+"\n");
	buffer.append("extents:\t\t\t\t\t"+extents+"\n");
	buffer.append("session_error:\t\t\t\t\t"+session_error+"\n");
	buffer.append("regular:\t\t\t\t\t"+regular+"\n");
	buffer.append("glmax/glmin:\t\t\t\t\t"+glmax+" "+glmin+"\n");
	buffer.append("Extension bytes 2-4:\t\t\t\t"+(int)extension[1]+" "
			   +(int)extension[2]+" "+(int)extension[3]+"\n");
		
		
	////// extensions
	if (extension[0] != 0) {
	    extlist = getExtensionsList();
	    n = extlist.length;
	    buffer.append("\n\nExtensions\n");
	    buffer.append("----------------------------------------------------------------------\n");
	    buffer.append("#\tCode\tSize\n");
	    for(i=0; i<n; i++)
		buffer.append((i+1)+"\t"+extlist[i][1]+"\t"+extlist[i][0]+"\n");
	    buffer.append("\n\n");
	}

	return buffer.toString();
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Print a voxel timecourse to standard out.
     * @param d 1D double array of timecourse values, length TDIM
     */
    public void printDoubleTmcrs(double d[]) {

	short i;
	NumberFormat nf;

	nf = NumberFormat.getInstance();
	nf.setMaximumFractionDigits(6);
	nf.setGroupingUsed(false);

	for (i=0; i<TDIM; i++)
	    System.out.println(nf.format(d[i]));

	return;
    }



    /**
     * Sets the dataset file name(s) and additional options inferred from the file name. Use this to
     * make an image ready to load an existing header and image.
     *
     * @param name either a .nii[.gz] file or a .hdr file.
     */
    public void setFilename(String name) {
        
        if (name.endsWith(".nii") || name.endsWith(".nii.gz")) {
            setToNii();

            ds_hdrname = name;
            ds_datname = name;

            if (name.endsWith(".nii.gz")) {
                gzip = true;
            }
        }
        
        if (name.endsWith(".hdr")) {
            setToNi1();
            
            ds_hdrname = name;

            String root = name.substring(0, name.length() - 4);

            File f = new File(root + ".img");
            
            if (f.exists()) {
                ds_datname = root + ".img";
            }
            else {
                f = new File(root + ".img.gz");
                
                if (f.exists()) {
                    ds_datname = root + ".img.gz";
                    gzip = true;
                }
            }
                
        }

    }


     /**
     * Sets the dataset file name. 
     *
     * @param root, path and filename for the dataset, without extension.
     * @param nii if true, set the image file name to root.nii[.gz], if false, set separate root.hdr / img[.gz] names.
     * @param gz if true, output will be compressed, image extension will change accordingly.
     *
     */
    public void setFilename(String root, boolean nii, boolean gz) {

        String hdrExt = null;
        String datExt = null;

        gzip = gz;

        if (nii) {

            setToNii();
            
            if (gzip) {
                hdrExt = NI1_EXT + GZIP_EXT;
            }
            else {
                hdrExt = NI1_EXT;
            }
    
            datExt = hdrExt;
        }
        else {

            setToNi1();
            
            hdrExt = ANZ_HDR_EXT;
            
            if (gzip) {
                datExt = ANZ_DAT_EXT + GZIP_EXT;
            }
            else {
                datExt = ANZ_DAT_EXT;
            }
        }

        ds_hdrname = root + hdrExt;
        ds_datname = root + datExt;
        
    }

    
    protected void setScale(double slope, double inter) {

        scl_slope = slope == 0.0 ? 1.0f : (float)slope;
        
        scl_inter = (float)(inter);
    }
   


    //////////////////////////////////////////////////////////////////
    /*
     * Set fields to make this a nii (n+1) (single file) dataset
     * switching from nii to anz/hdr affects
     * -- magic field "n+1\0" not "ni1\0"
     * -- vox_offset must include 352+extensions
     * -- internal ds_is_nii flag
     *
     * NOTE: all changes are in memory, app still needs to
     * write header and data for change to occur on disk
     *
     * maybe add auto set of dat name to hdr name, strip img/hdr ??
     */
    private void setToNii() {
	int i,n;
	int extlist[][];

	ds_is_nii = true;
	magic = new StringBuffer(NII_MAGIC_STRING); 

	vox_offset = NII_HDR_SIZE;
	if (extension[0] != 0) {
	    extlist = getExtensionsList();
	    n = extlist.length;
	    for(i=0; i<n; i++)
		vox_offset += extlist[i][0];
	}
	return;
    }


    //////////////////////////////////////////////////////////////////
    /*
     * Set fields to make this a ni1 (2 file img/hdr) dataset
     * switching from nii to anz/hdr affects
     * -- magic field "n+1\0" vs "ni1\0"
     * -- vox_offset does notinclude 352+extensions
     * -- internal ds_is_nii flag
     *
     * NOTE: all changes are in memory, app still needs to
     * write header and data for change to occur on disk
     *
     * maybe add auto set of dat name to hdr name, strip img/hdr ??
     */
    private void setToNi1() {
	int n;
	int extlist[][];

	ds_is_nii = false;
	magic = new StringBuffer(ANZ_MAGIC_STRING); 

	// if there was stuff after header before data it cannot
	// survive transition from nii to ni1
	vox_offset = 0;

	return;
    }



    //////////////////////////////////////////////////////////////////
    /**
     * Set the dataset dimensions
     */
    public void setDims(short a, short x, short y, short z, short t, short d5, short d6, short d7) {

	dim[0] = a;
	dim[1] = x;
	dim[2] = y;
	dim[3] = z;
	dim[4] = t;
	dim[5] = d5;
	dim[6] = d6;
	dim[7] = d7;
		
	XDIM = x;
	YDIM = y;
	ZDIM = z;
	TDIM = t;

	return;
    }


    /**
     * Casts args to short and calls setDims.
     */
    public void setDims(int a, int x, int y, int z, int t, int d5, int d6, int d7) {
	setDims((short)a, (short)x, (short)y, (short)z, (short)t, (short)d5, (short)d6, (short)d7);
    }

    
    public void setDims(short[] dims) {
	setDims(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], dims[6], dims[7]);
    }
    

    /**
     * Sets the pixel dims in order. The first parameter is the qfac, which controls the orientation
     * of the image. The next three are spatial dimensions, then the time spacing, then the remaining
     * dimensions.
     *
     * @see http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/pixdim.html
     */
    public void setPixDims(float q, float x, float y, float z, float t, float d5, float d6, float d7) {
	pixdim[0] = q < 0.0f ? -1.0f : 1.0f;
	pixdim[1] = x;
	pixdim[2] = y;
	pixdim[3] = z;
	pixdim[4] = t;
	pixdim[5] = d5;
	pixdim[6] = d6;
	pixdim[7] = d7;

    }


    /**
     * Sets the pixel dims in order. The first parameter is the qfac, which controls the orientation
     * of the image. The next three are spatial dimensions, then the time spacing, then the remaining
     * dimensions.
     *
     * @see http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/pixdim.html
     */
    public void setPixDims(float[] dims) {
	setPixDims(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], dims[6], dims[7]);
    }


    /**
     *  Casts args to float and calls setPixDims.
     */
    public void setPixDims(double q, double x, double y, double z, double t, double d5, double d6, 
			     double d7) {
	setPixDims((float)q, (float)x, (float)y, (float)z, (float)t, (float)d5, (float)d6, (float)d7);
    }

    //////////////////////////////////////////////////////////////////
    /**
     * Set the dataset datatype.  (bitpix will also be set accordingly.)
     * @param code nifti-1 datatype code
     */
    private void setDataType(short code) {
	datatype = code;
	bitpix = (short)(bytesPerVoxel(code)*8);
	return;
    }

    

    //////////////////////////////////////////////////////////////////
    /**
     * Get the dataset datatype.  
     * @return  datatype code (note: it is not guaranteed to be a valid
     * code, what is there is what you get...)
     */
    public short getDatatype() {
	return(datatype);
    }
	
    //////////////////////////////////////////////////////////////////
    /**
     * Get the bitpix field
     * @return  bitpix: number of bits per pixel
     */
    public short getBitpix() {
	return(bitpix);
    }
	

    //////////////////////////////////////////////////////////////////
    /**
     * Decode the nifti intent codes                            
     * @param icode nifti intent code
     * @return a terse string describing the intent
     */
    public String decodeIntent(short icode) {

	switch(icode) {

	case NIFTI_INTENT_NONE:
	    return("NIFTI_INTENT_NONE");
	case NIFTI_INTENT_CORREL:
	    return("NIFTI_INTENT_CORREL");
	case NIFTI_INTENT_TTEST:
	    return("NIFTI_INTENT_TTEST");
	case NIFTI_INTENT_FTEST:
	    return("NIFTI_INTENT_FTEST");
	case NIFTI_INTENT_ZSCORE:
	    return("NIFTI_INTENT_ZSCORE");
	case NIFTI_INTENT_CHISQ:
	    return("NIFTI_INTENT_CHISQ");
	case NIFTI_INTENT_BETA:
	    return("NIFTI_INTENT_BETA");
	case NIFTI_INTENT_BINOM:
	    return("NIFTI_INTENT_BINOM");
	case NIFTI_INTENT_GAMMA:
	    return("NIFTI_INTENT_GAMMA");
	case NIFTI_INTENT_POISSON:
	    return("NIFTI_INTENT_POISSON");
	case NIFTI_INTENT_NORMAL:
	    return("NIFTI_INTENT_NORMAL");
	case NIFTI_INTENT_FTEST_NONC:
	    return("NIFTI_INTENT_FTEST_NONC");
	case NIFTI_INTENT_CHISQ_NONC:
	    return("NIFTI_INTENT_CHISQ_NONC");
	case NIFTI_INTENT_LOGISTIC:
	    return("NIFTI_INTENT_LOGISTIC");
	case NIFTI_INTENT_LAPLACE:
	    return("NIFTI_INTENT_LAPLACE");
	case NIFTI_INTENT_UNIFORM:
	    return("NIFTI_INTENT_UNIFORM");
	case NIFTI_INTENT_TTEST_NONC:
	    return("NIFTI_INTENT_TTEST_NONC");
	case NIFTI_INTENT_WEIBULL:
	    return("NIFTI_INTENT_WEIBULL");
	case NIFTI_INTENT_CHI:
	    return("NIFTI_INTENT_CHI");
	case NIFTI_INTENT_INVGAUSS:
	    return("NIFTI_INTENT_INVGAUSS");
	case NIFTI_INTENT_EXTVAL:
	    return("NIFTI_INTENT_EXTVAL");
	case NIFTI_INTENT_PVAL:
	    return("NIFTI_INTENT_PVAL");
	case NIFTI_INTENT_ESTIMATE:
	    return("NIFTI_INTENT_ESTIMATE");
	case NIFTI_INTENT_LABEL:
	    return("NIFTI_INTENT_LABEL");
	case NIFTI_INTENT_NEURONAME:
	    return("NIFTI_INTENT_NEURONAME");
	case NIFTI_INTENT_GENMATRIX:
	    return("NIFTI_INTENT_GENMATRIX");
	case NIFTI_INTENT_SYMMATRIX:
	    return("NIFTI_INTENT_SYMMATRIX");
	case NIFTI_INTENT_DISPVECT:
	    return("NIFTI_INTENT_DISPVECT");
	case NIFTI_INTENT_VECTOR:
	    return("NIFTI_INTENT_VECTOR");
	case NIFTI_INTENT_POINTSET:
	    return("NIFTI_INTENT_POINTSET");
	case NIFTI_INTENT_TRIANGLE:
	    return("NIFTI_INTENT_TRIANGLE");
	case NIFTI_INTENT_QUATERNION:
	    return("NIFTI_INTENT_QUATERNION");
	default:
	    return("INVALID_NIFTI_INTENT_CODE");
	}
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Decode the nifti datatype codes                            
     * @param dcode nifti datatype code
     * @return a terse string describing the datatype
     */
    public String decodeDatatype(short dcode) {

	switch(dcode) {

	case DT_NONE:
	    return("DT_NONE");
	case DT_BINARY:
	    return("DT_BINARY");
	case NIFTI_TYPE_UINT8:
	    return("NIFTI_TYPE_UINT8");
	case NIFTI_TYPE_INT16:
	    return("NIFTI_TYPE_INT16");
	case NIFTI_TYPE_INT32:
	    return("NIFTI_TYPE_INT32");
	case NIFTI_TYPE_FLOAT32:
	    return("NIFTI_TYPE_FLOAT32");
	case NIFTI_TYPE_COMPLEX64:
	    return("NIFTI_TYPE_COMPLEX64");
	case NIFTI_TYPE_FLOAT64:
	    return("NIFTI_TYPE_FLOAT64");
	case NIFTI_TYPE_RGB24:
	    return("NIFTI_TYPE_RGB24");
	case DT_ALL:
	    return("DT_ALL");
	case NIFTI_TYPE_INT8:
	    return("NIFTI_TYPE_INT8");
	case NIFTI_TYPE_UINT16:
	    return("NIFTI_TYPE_UINT16");
	case NIFTI_TYPE_UINT32:
	    return("NIFTI_TYPE_UINT32");
	case NIFTI_TYPE_INT64:
	    return("NIFTI_TYPE_INT64");
	case NIFTI_TYPE_UINT64:
	    return("NIFTI_TYPE_UINT64");
	case NIFTI_TYPE_FLOAT128:
	    return("NIFTI_TYPE_FLOAT128");
	case NIFTI_TYPE_COMPLEX128:
	    return("NIFTI_TYPE_COMPLEX128");
	case NIFTI_TYPE_COMPLEX256:
	    return("NIFTI_TYPE_COMPLEX256");
	default:
	    return("INVALID_NIFTI_DATATYPE_CODE");
	}
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Return bytes per voxel for each nifti-1 datatype           
     * @param dcode nifti datatype code
     * @return a short with number of bytes per voxel, 0 for unknown,
     *  -1 for 1 bit
     */
    public short bytesPerVoxel(short dcode) {

	switch(dcode) {

	case DT_NONE:
	    return(0);
	case DT_BINARY:
	    return(-1);
	case NIFTI_TYPE_UINT8:
	    return(1);
	case NIFTI_TYPE_INT16:
	    return(2);
	case NIFTI_TYPE_INT32:
	    return(4);
	case NIFTI_TYPE_FLOAT32:
	    return(4);
	case NIFTI_TYPE_COMPLEX64:
	    return(8);
	case NIFTI_TYPE_FLOAT64:
	    return(8);
	case NIFTI_TYPE_RGB24:
	    return(3);
	case NIFTI_TYPE_RGBA32:
	    return(4);
	case DT_ALL:
	    return(0);
	case NIFTI_TYPE_INT8:
	    return(1);
	case NIFTI_TYPE_UINT16:
	    return(2);
	case NIFTI_TYPE_UINT32:
	    return(4);
	case NIFTI_TYPE_INT64:
	    return(8);
	case NIFTI_TYPE_UINT64:
	    return(8);
	case NIFTI_TYPE_FLOAT128:
	    return(16);
	case NIFTI_TYPE_COMPLEX128:
	    return(16);
	case NIFTI_TYPE_COMPLEX256:
	    return(32);
	default:
	    return(0);
	}
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Decode the nifti slice order codes                            
     * @param code nifti slice order code
     * @return a terse string describing the slice order
     */
    public String decodeSliceOrder(short code) {

	switch(code) {

	case NIFTI_SLICE_SEQ_INC:
	    return("NIFTI_SLICE_SEQ_INC");
	case NIFTI_SLICE_SEQ_DEC:
	    return("NIFTI_SLICE_SEQ_DEC");
	case NIFTI_SLICE_ALT_INC:
	    return("NIFTI_SLICE_ALT_INC");
	case NIFTI_SLICE_ALT_DEC:
	    return("NIFTI_SLICE_ALT_DEC");
	default:
	    return("INVALID_NIFTI_SLICE_SEQ_CODE");
	}
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Decode the nifti xform codes                            
     * @param code nifti xform code
     * @return a terse string describing the coord. system
     */
    public String decodeXform(short code) {

	switch(code) {
	case NIFTI_XFORM_UNKNOWN:
	    return("NIFTI_XFORM_UNKNOWN");
	case NIFTI_XFORM_SCANNER_ANAT:
	    return("NIFTI_XFORM_SCANNER_ANAT");
	case NIFTI_XFORM_ALIGNED_ANAT:
	    return("NIFTI_XFORM_ALIGNED_ANAT");
	case NIFTI_XFORM_TALAIRACH:
	    return("NIFTI_XFORM_TALAIRACH");
	case NIFTI_XFORM_MNI_152:
	    return("NIFTI_XFORM_MNI_152");
	default:
	    return("INVALID_NIFTI_XFORM_CODE");
	}
    }


    //////////////////////////////////////////////////////////////////
    /**
     * Decode the nifti unit codes                            
     * @param code nifti units code
     * @return a terse string describing the unit
     */
    public String decodeUnits(short code) {

	switch(code) {
	case NIFTI_UNITS_UNKNOWN:
	    return("NIFTI_UNITS_UNKNOWN");
	case NIFTI_UNITS_METER:
	    return("NIFTI_UNITS_METER");
	case NIFTI_UNITS_MM:
	    return("NIFTI_UNITS_MM");
	case NIFTI_UNITS_MICRON:
	    return("NIFTI_UNITS_MICRON");
	case NIFTI_UNITS_SEC:
	    return("NIFTI_UNITS_SEC");
	case NIFTI_UNITS_MSEC:
	    return("NIFTI_UNITS_MSEC");
	case NIFTI_UNITS_USEC:
	    return("NIFTI_UNITS_USEC");
	case NIFTI_UNITS_HZ:
	    return("NIFTI_UNITS_HZ");
	case NIFTI_UNITS_PPM:
	    return("NIFTI_UNITS_PPM");
	default:
	    return("INVALID_NIFTI_UNITS_CODE");
	}
    }


    // //////////////////////////////////////////////////////////////////
    // /**
    //  * Check the header fields for valid settings
    //  * @return 0 if all checks are passed else, error code
    //  */
    // public short checkHeader() {
	

    //     // check that each dim claimed to be in use is not 0
    //     // check for negative voxel sizes
    //     // check for 348
    //     // check for magic code n+1 or ni1
    //     // check bitpix divisible by 8 and in sync with datatype

    //     return 0;
    // }

    
    /**
     * Checks whether an .hdr header is Analyze or NIfTI. It is assumed to not be NIfTI unless
     * it ends in .hdr, and the "magic" field is ANZ_MAGIC_STRING.
     * 
     */
    public static boolean hdrIsNifti(String file) throws IOException {
        if (! file.endsWith(".hdr") ) {
            return false;
        }

	Nifti1Dataset nds = new Nifti1Dataset(file);
	
	nds.readHeader();

        if (nds.magic.toString().equals(ANZ_MAGIC_STRING)) {
            return true;
        }
        else {
            return false;
        }
    }




    ////////////////////////////////////////////////////////////////////
    //
    // Set default settings for all nifti1 header fields and some
    // other fields for this class
    //
    ////////////////////////////////////////////////////////////////////
    private void setDefaults() {

	int i;
	ByteOrder bo;

	ds_hdrname = 	new String("");  
	ds_datname = 	new String("");
	ds_is_nii = 	true;		// use single file .nii format or not
	bo = ByteOrder.nativeOrder();
	if (bo == ByteOrder.BIG_ENDIAN)
	    big_endian = 	true;	// endian flag, need to flip little end.
	else
	    big_endian = 	false;


        gzip = false; // whether to use GZIP I/O
	
	sizeof_hdr = 	ANZ_HDR_SIZE;		// must be 348 bytes
	data_type_string = new StringBuffer(); // 10 char UNUSED
	for (i=0; i<10; i++)
	    data_type_string.append("\0");
	db_name = new StringBuffer();	// 18 char UNUSED
	for (i=0; i<18; i++)
	    db_name.append("\0");
	extents = 	0;		// UNUSED
	session_error =	0;		// UNUSED
	regular = new StringBuffer("\0");	// UNUSED
	dim_info = new StringBuffer("\0");	// MRI slice ordering
	freq_dim=0; phase_dim=0; slice_dim=0;
	dim = new short[8];		// data array dimensions (8 shorts)
	for (i=0; i<8; i++)
	    dim[i] = 0;
	XDIM=0; YDIM=0; ZDIM=0; TDIM=0;
	intent = new float[3];		// intents p1 p2 p3
	for (i=0; i<3; i++)
	    intent[i] = (float)0.0;
	intent_code = NIFTI_INTENT_NONE;
	datatype = DT_NONE;		// datatype of image blob
	bitpix = 0;			// #bits per voxel
	slice_start = 0;		// first slice index
	pixdim = new float[8];		// grid spacings
	pixdim[0] = 1; qfac = 1;
	for (i=1; i<8; i++)
	    pixdim[i] = (float)0.0;

	vox_offset = NII_HDR_SIZE;	// offset to data blob in .nii file
					// for .nii files default is 352
					
	scl_slope = (float)0.0;		// data scaling: slope
	scl_inter = (float)0.0;		// data scaling: intercept
	slice_end = 0;			// last slice index
	slice_code = (byte) 0;		// slice timing order
	xyzt_units = (byte) 0;		// units of pixdim[1-4]
	xyz_unit_code = NIFTI_UNITS_UNKNOWN;
	t_unit_code = NIFTI_UNITS_UNKNOWN;

	cal_max = (float) 0.0;		// max display intensity
	cal_min = (float) 0.0;		// min display intensity
	slice_duration = (float) 0.0;	// time to acq. 1 slice
	toffset = (float) 0.0;		// time axis shift
	glmax = 0;			// UNUSED
	glmin = 0;			// UNUSED

	descrip = new StringBuffer();	// 80 char any text you'd like (comment)
	for (i=0; i<80; i++)
	    descrip.append("\0");
	aux_file = new StringBuffer();	// 24 char auxiliary filename
	for (i=0; i<24; i++)
	    aux_file.append("\0");


	qform_code = NIFTI_XFORM_UNKNOWN;	// code for quat. xform
	sform_code = NIFTI_XFORM_UNKNOWN;	// code for affine xform

	quatern = new float[3];		// 3 floats Quaternion b,c,d params
	qoffset = new float[3];		// 3 floats Quaternion x,y,z shift
	for (i=0; i<3; i++) {
	    quatern[i] = (float)0.0;
	    qoffset[i] = (float)0.0;
	}

	srow_x = new float[4];		// 4 floats 1st row affine xform
	srow_y = new float[4];		// 4 floats 2nd row affine xform
	srow_z = new float[4];		// 4 floats 3rd row affine xform
	for (i=0; i<4; i++) {
	    srow_x[i] = (float)0.0;
	    srow_y[i] = (float)0.0;
	    srow_z[i] = (float)0.0;
	}

	intent_name = new StringBuffer(); // 16 char name/meaning/id for data
	for (i=0; i<16; i++)
	    intent_name.append("\0");

        
	magic = new StringBuffer(NII_MAGIC_STRING); // 4 char id must be "ni1\0" or "n+1\0"

	extension = new byte[4];	// 4 byte array, [0] is for extensions
	for (i=0; i<4; i++)
	    extension[i] = (byte)0;

	extensions_list = new Vector<Object>(5);	// list of int[2] size/code pairs for exts.
	extension_blobs = new Vector<Object>(5);	// vector to hold data in each ext.

	return;
    }

    ////////////////////////////////////////////////////////////////////
    //
    // Unpack/pack the 3 bitfields in the dim_info char                 
    //	bits 0,1 freq_dim
    //	bits 2,3 phase_dim
    //	bits 4,5 slice_dim
    //
    ////////////////////////////////////////////////////////////////////
    private short[] unpackDimInfo(int b) {
	short s[];

	s = new short[3];
	s[0] = (short) (b & 3);
	s[1] = (short) ((b >>> 2) & 3);
	s[2] = (short) ((b >>> 4) & 3);
	return s;
    }
    private byte packDimInfo(short freq, short phase, short slice) {
		
	int i = 0;

	i = (i & ((int)(slice)&3)) << 2;
	i = (i & ((int)(phase)&3)) << 2;
	i = (i & ((int)(freq)&3));
	return ((byte)i);
    }


    ////////////////////////////////////////////////////////////////////
    //
    // Unpack/pack the 2 bitfields in the xyzt_units field
    //	bits 0,1,2 are the code for units for xyz
    //	bits 3,4,5 are the code for units for t, no need to shift
    //	bits for t code, the code is a multiple of 8.
    //
    ////////////////////////////////////////////////////////////////////
    private short[] unpackUnits(int b) {
	short s[];

	s = new short[2];
	s[0] = (short) (b & 007);
	s[1] = (short) (b & 070);
	return s;
    }

    private byte packUnits(short space, short time) {

		
	return( (byte) (((int) (space) & 007) | ((int)(time) & 070)) );
    }







    /**
     * Write data to an output stream using current header settings for data type, scaling etc.
     *
     * @param data the 4D data array.
     * @param out an OutputStream to write the data.
     * @exception IOException
     * 
     */
    private void writeData(double data[][][][], OutputStream out) throws IOException {
	
	short i,j,k,n;
	short ZZZ;
	EndianCorrectOutputStream ecs;

	// for 2D volumes, zdim may be 0
	ZZZ = ZDIM;
	if (dim[0] == 2)
	    ZZZ = 1;

	ecs = new EndianCorrectOutputStream(out,big_endian);

        
        int comps = components();

	switch (datatype) {

	case NIFTI_TYPE_INT8:
            for (n=0; n<comps; n++)
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            if (scl_slope == 0)
                                ecs.write((byte)Math.round(data[i][j][k][n]));
                            else
                                ecs.write((byte)Math.round((data[i][j][k][n] - scl_inter) / scl_slope));
                        }
	    break;
            

	case NIFTI_TYPE_INT16:
            for (n=0; n<comps; n++)
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            if (scl_slope == 0)
                                ecs.writeShortCorrect((short)Math.round(data[i][j][k][n]));
                            else
                                ecs.writeShortCorrect((short)Math.round((data[i][j][k][n] - scl_inter) / scl_slope));
                        }
	    break;

            
	case NIFTI_TYPE_INT32:
            for (n=0; n<comps; n++)
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            if (scl_slope == 0)
                                ecs.writeIntCorrect((int)Math.round(data[i][j][k][n]));
                            else
                                ecs.writeIntCorrect((int)Math.round((data[i][j][k][n] - scl_inter) / scl_slope));
                        }
	    break;
            
	case NIFTI_TYPE_INT64:
            for (n=0; n<comps; n++)
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            if (scl_slope == 0)
                                ecs.writeLongCorrect(Math.round(data[i][j][k][n]));
                            else
                                ecs.writeLongCorrect(Math.round((data[i][j][k][n] - scl_inter) / scl_slope));
                        }
	    break;
	case NIFTI_TYPE_FLOAT32:
            for (n=0; n<comps; n++)
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            if (scl_slope == 0) 
                                ecs.writeFloatCorrect((float)(data[i][j][k][n]));
                            
                            else
                                ecs.writeFloatCorrect((float)((data[i][j][k][n] - scl_inter) / scl_slope));
                        }
	    break;
	case NIFTI_TYPE_FLOAT64:
            for (n=0; n<comps; n++)
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            if (scl_slope == 0)
                                ecs.writeDoubleCorrect(data[i][j][k][n]);
                            else
                                ecs.writeDoubleCorrect((data[i][j][k][n] - scl_inter) / scl_slope);
                        }
	    break;


	case NIFTI_TYPE_RGB24:
            for (n=0; n<comps; n++)
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            int rgb = (int)(data[i][j][k][n]);

                            ecs.write( rgb >> 16 );
                            ecs.write( (rgb >> 8) & 0x00ff );
                            ecs.write( rgb & 0x00ff );
                        }
	    break;

	case NIFTI_TYPE_RGBA32:
            for (n=0; n<comps; n++)
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            int rgba = (int)(data[i][j][k][n]);

                            // don't want this to ever get byte swapped
                            ecs.write( rgba >> 24 );
                            ecs.write( (rgba >> 16) & 0x00ff );
                            ecs.write( (rgba >> 8) & 0x00ff );
                            ecs.write( rgba & 0x00ff );
                        }
	    break;
            

	case DT_NONE:
	case DT_BINARY:
	case NIFTI_TYPE_UINT8:
	case NIFTI_TYPE_UINT16:
	case NIFTI_TYPE_UINT64:
	case NIFTI_TYPE_COMPLEX64:
	case NIFTI_TYPE_FLOAT128:
	case NIFTI_TYPE_COMPLEX128:
	case NIFTI_TYPE_COMPLEX256:
	case DT_ALL:
	default:
	    throw new IOException("Sorry, cannot yet write nifti-1 datatype "+decodeDatatype(datatype));

	}
        
        ecs.close();

	return;
    }

 



    /**
     * Writes the header and image data to disk.
     * 
     *
     */
    protected void writeImage(double[][][][] data) {
        
        try {

            byte[] hdrBytes = headerToBytes();

            FileOutputStream hdrFos = new FileOutputStream(ds_hdrname);

            // may be the same object, or not
            OutputStream hdrOut, dsOut;

            int BUFFERSIZE = 1024*1024*32;

            if (gzip) {
                hdrOut = new BufferedOutputStream(new GZIPOutputStream(hdrFos, BUFFERSIZE / 2), BUFFERSIZE / 2);
            }
            else {
                hdrOut = new BufferedOutputStream(hdrFos, BUFFERSIZE);
            }

            if (ds_is_nii) {
                dsOut = hdrOut;
            }
            else {
            
                FileOutputStream dataFos = new FileOutputStream(ds_datname); 

                if (gzip) {
                    dsOut = new BufferedOutputStream(new GZIPOutputStream(dataFos, BUFFERSIZE / 2), BUFFERSIZE / 2);
                }
                else {
                    dsOut = new BufferedOutputStream(dataFos, BUFFERSIZE);
                }
            
            }
       
            hdrOut.write(hdrBytes);

            writeData(data, dsOut);

            hdrOut.close();
            dsOut.close();
        }
        catch(IOException e) {
            throw new LoggedException(e.getMessage());
        }
    }


    /**
     * truncate or pad a string to make it the needed length
     * @param s input string
     * @param n desired length
     * @return s String of length n with as much of s as possible, padded
     * with 0 if necessary
     */
    private byte[] setStringSize(StringBuffer s, int n) {

	byte b[];
	int i,slen;

	slen = s.length();

	if (slen >= n)
	    return(s.toString().substring(0,n).getBytes());

	b = new byte[n];
	for (i=0; i<slen; i++)
	    b[i] = (byte)s.charAt(i);
	for (i=slen; i<n; i++)
	    b[i] = 0;

	return(b);
    }
		

    /**
     * Gets the Camino data type string, for passing to the data.DataSource classes. 
     *
     * @throws LoggedException if no Camino equivalent exists for the NIFTI datatype.
     *
     */
    public String caminoDataTypeString() {
	switch (datatype) {

	case NIFTI_TYPE_UINT8: return "ubyte";

	case NIFTI_TYPE_INT8: return "byte";

	case NIFTI_TYPE_INT16: return "short";

	case NIFTI_TYPE_UINT16: return "ushort";

	case NIFTI_TYPE_INT32: return "int";

	case NIFTI_TYPE_UINT32: return "uint";

	case NIFTI_TYPE_INT64: return "long";

	case NIFTI_TYPE_FLOAT32: return "float";

	case NIFTI_TYPE_FLOAT64: return "double";

	default: throw new LoggedException("Camino does not support this data type: " + decodeDatatype(datatype));

	}
    }

    
    // BEGIN ImageHeader methods

    public int xDataDim() {
	return XDIM;
    }
    
    public int yDataDim() {
	return YDIM; 
    }
    
    public int zDataDim() {
	return ZDIM;
    }
    
    public int[] getDataDims() {
	return new int[] {XDIM, YDIM, ZDIM};
    }
    
    public double xVoxelDim() {
	return pixdim[1];
    }
    
    public double yVoxelDim() {
	return pixdim[2];
    }
    
    public double zVoxelDim() {
	return pixdim[3];
    }

    public double[] getVoxelDims() {
	return new double[] {pixdim[1], pixdim[2], pixdim[3]};
    }


    /** 
     * @return dim[4] if dim[4] > 1, or dim[5] if dim[5] > 0. If dim[5] == 0 or dim[0] == 3 (3D image), 
     * return 1.
     *
     */
    public int components() {

	if (dim[0] < 3) {
	    throw new LoggedException("Can't handle " + dim[0] + "D images");
	}
	if (dim[0] == 3) {
	    return 1;
	}
	if (dim[0] == 4) {
	    // treat as 4D image
	    return dim[4] > 0 ? dim[4] : 1;
	}

	if (dim[4] > 1 && dim[5] > 1) {
	    throw new LoggedException("Can't handle multi-component images with multiple time points");
	}
	
	return dim[5] > 0 ? dim[5] : 1;
	
    }


    /**
     * Gets the transformation matrix R that provides a transformation from a = (i, j, k, 1) in 
     * voxel space to b = (x, y, z, 1) in physical space. Following the NIfTI standard sform is used
     * if sform_code is in the range [1,2], else qform is used if qform_code is in that range,
     * else we default to Analyze. Note that ITK will ignore sform unless qform == 0, and will override
     * sform to match qform. DTI-TK will also change output to be LPI.
     * <p>
     * We don't use sform or qform if they are set to a value greater than 2, denoting alignment of the image
     * to a standard space.  
     *
     */
    public RealMatrix getVoxelToPhysicalTransform() {

        RealMatrix R = new RealMatrix(4,4);
        
        // always
        R.entries[3][0] = 0.0;
        R.entries[3][1] = 0.0;
        R.entries[3][2] = 0.0;
        R.entries[3][3] = 1.0;
        

        if (sform_code > 0 && sform_code < 3) {
            
            // developer, Y U NO store sform as matrix?
            
            // anyhow this is the easy way
            
            for (int i = 0; i < 4; i++) {
                R.entries[0][i] = srow_x[i]; 
                R.entries[1][i] = srow_y[i]; 
                R.entries[2][i] = srow_z[i]; 
            }
            
        }
        else if (qform_code > 0) {

            float b = quatern[0];
            float c = quatern[1];
            float d = quatern[2];
            
            float a = (float)Math.sqrt(1.0 - (b*b + c*c + d*d) );
            
            RealMatrix rotation = new RealMatrix(3,3);
            
            rotation.entries[0][0] = a*a+b*b-c*c-d*d;
            rotation.entries[1][1] = a*a+c*c-b*b-d*d;
            rotation.entries[2][2] = a*a+d*d-c*c-b*b;
            
            rotation.entries[1][0] = 2.0*b*c+2*a*d;
            rotation.entries[0][1] = 2.0*b*c-2*a*d;

            rotation.entries[0][2] = 2.0*b*d+2.0*a*c;
            rotation.entries[2][0] = 2.0*b*d-2.0*a*c;

            rotation.entries[1][2] = 2.0*c*d-2*a*b;
            rotation.entries[2][1] = 2.0*c*d+2*a*b;

            RealMatrix scaling = new RealMatrix(3,3);

            scaling.entries[0][0] = pixdim[1];
            scaling.entries[1][1] = pixdim[2];
            scaling.entries[2][2] = qfac * pixdim[3];

            RealMatrix product = rotation.product(scaling);
            
            
            for (int i = 0; i < 3; i++) { 
                for (int j = 0; j < 3; j++) { 
                    R.entries[i][j] = product.entries[i][j];
                }
            }
            
            R.entries[0][3] = qoffset[0];
            R.entries[1][3] = qoffset[1];
            R.entries[2][3] = qoffset[2];

        }
        else {

            // Last ditch attempt 

            // do it the old fashioned way, which is equivalent to Analyze, assuming LPI and no negative
            // pixdim shenanigans
            R.entries[0][0] = pixdim[1];
            R.entries[1][1] = pixdim[2];
            R.entries[2][2] = pixdim[3];

            logger.warning("Neither qform nor sform set for anatomical alignment, assuming identity rotation");
            
        }

        return R;

    }
       
       
    /**
     * Get the data source for this header. For multi-component data, scanner order is assumed.
     *
     */
    public DataSource getImageDataSource() {


	double scaleSlope = scl_slope == 0.0 ? 1.0 : scl_slope;
	double scaleInter = scl_slope == 0.0 ? 0.0 : scl_inter;

	boolean intelByteOrder = !big_endian;
	String dataType = caminoDataTypeString();

	// voxel offset
	int offset = (int)vox_offset;
	
	String fileName = ds_datname;
	
	// read data, scale it, and then wrap in a data source
	if (!(dim[0] > 1 && dim[0] < 6)) {
	    throw new LoggedException("Can't handle this " + dim[0] + "-dimensional dataset. Camino " +
				      "supports 3D scalar volumes or 4D " + "or 5D multivariate datasets");
	}

	int voxels = dim[1] * dim[2] * dim[3];

	int components = components();

	if (components == 1) {
	    // try to save some memory
	    // return voxel order data source if we have a 3D image 
	    return new VoxelOrderScaledDataSource(ds_datname, 1, dataType, intelByteOrder, 
						  offset, scaleSlope, scaleInter);
	}
	
        // Much faster than trying to do random access to a compressed file
        // if you read each volume independently, you need to decompress the file over and over
        // to get to the right index for each 3D volume. This takes a long time for large files
        DataSource source = new ScannerOrderScaledDataSource(ds_datname, voxels, components, dataType, 
                                                             intelByteOrder, offset, scaleSlope, scaleInter);
        
	return source;
	
    }


    /**
     * Gets the underlying input stream for the image data, after skipping the header
     *
     */
    private EndianNeutralDataInputStream getImageDataInputStream() {

        InputStream is;

        // Note different input stream; to header methods, because
        // EndianCorrectInputStream doesn't support unsigned types
	EndianNeutralDataInputStream ecs;

        int comps = components();

        try {
            
            if (ds_datname.endsWith(".gz")) {
                // Buffer more on the decompression side than the disk side, since that's where a 
                // lot of time is wasted unzipping small chunks
                
                // buffer this much on the (compressed) disk side, more on the uncompressed side
                int bufferSize = 1024*1024*8;
                
                // BufferedInputStream around FileInputStream here? Need to test
                is = new BufferedInflaterInputStream(new GZIPInputStream(new FileInputStream(ds_datname), bufferSize), bufferSize * 3);
                
                // Note: it may be possible to improve performance here by reading the entire file in one sweep into a byte array,
                // at the cost of using more memory
                
            }
            else {
                
                int bufferSize = 1024*1024*24;
                
                is = new BufferedInputStream(new FileInputStream(ds_datname), bufferSize);
            }
            
            is.skip((int)vox_offset);
            
            // read the correct endian datatype from the byte array
            // add scaling if necessary
            
            ecs = new EndianNeutralDataInputStream(is, !big_endian);
            
        }
        catch (IOException e) {
            throw new LoggedException("Error reading image data: " + e.getMessage());
        }

        return ecs;

    }


    /**
     * Reads image volume data into a 4D double array.
     *
     * @param ecs an input stresm object that is positioned to read the next volume of data from the file.
     */
    private double[][][][] readNextVolume(EndianNeutralDataInputStream ecs) throws IOException {

        
        short ZZZ;
        int i,j,k,n;


        // for 2D volumes, zdim may be 0
        ZZZ = ZDIM;
        if (dim[0] == 2)
            ZZZ = 1;

        double[][][][] data = new double[XDIM][YDIM][ZZZ][1];
        

        switch (datatype) {
            
        case NIFTI_TYPE_INT8:
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            data[i][j][k][0] = (double) (ecs.readByte());
                            if (scl_slope != 0.0)
                                data[i][j][k][0] = data[i][j][k][0] * scl_slope + scl_inter;
                        }
            break;

        case NIFTI_TYPE_UINT8:
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            data[i][j][k][0] = (double) (ecs.readUnsignedByte());
                            if (scl_slope != 0.0)
                                data[i][j][k][0] = data[i][j][k][0] * scl_slope + scl_inter;
                        }
            break;

        case NIFTI_TYPE_INT16:
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            data[i][j][k][0] = (double) (ecs.readShort());
                            if (scl_slope != 0.0)
                                data[i][j][k][0] = data[i][j][k][0] * scl_slope + scl_inter;
                        }
            break;

        case NIFTI_TYPE_UINT16:
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            data[i][j][k][0] = (double) (ecs.readUnsignedShort());
                            if (scl_slope != 0.0)
                                data[i][j][k][0] = data[i][j][k][0] * scl_slope + scl_inter;
                        }
            break;

        case NIFTI_TYPE_INT32:
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            data[i][j][k][0] = (double) (ecs.readInt());
                            if (scl_slope != 0.0)
                                data[i][j][k][0] = data[i][j][k][0] * scl_slope + scl_inter;
                        }
            break;


        case NIFTI_TYPE_UINT32:
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            data[i][j][k][0] = (double) (ecs.readUnsignedInt());
                            if (scl_slope != 0.0)
                                data[i][j][k][0] = data[i][j][k][0] * scl_slope + scl_inter;
                        }
            break;


        case NIFTI_TYPE_INT64:
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            data[i][j][k][0] = (double) (ecs.readLong());
                            if (scl_slope != 0.0)
                                data[i][j][k][0] = data[i][j][k][0] * scl_slope + scl_inter;
                        }
            break;


        case NIFTI_TYPE_FLOAT32:
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            data[i][j][k][0] = (double) (ecs.readFloat());
                            if (scl_slope != 0.0)
                                data[i][j][k][0] = data[i][j][k][0] * scl_slope + scl_inter;
                        }
            break;


        case NIFTI_TYPE_FLOAT64:
                for (k=0; k<ZZZ; k++)
                    for (j=0; j<YDIM; j++)
                        for (i=0; i<XDIM; i++) {
                            data[i][j][k][0] = ecs.readDouble();
                            if (scl_slope != 0.0)
                                data[i][j][k][0] = data[i][j][k][0] * scl_slope + scl_inter;
                        }
            break;


        case DT_NONE:
        case DT_BINARY:
        case NIFTI_TYPE_COMPLEX64:
        case NIFTI_TYPE_UINT64:
        case NIFTI_TYPE_FLOAT128:
        case NIFTI_TYPE_RGB24:
        case NIFTI_TYPE_COMPLEX128:
        case NIFTI_TYPE_COMPLEX256:
        case DT_ALL:
                
        default:
            throw new LoggedException("Sorry, cannot yet read nifti-1 datatype "+decodeDatatype(datatype));
        }


        return data;
        
        
    }
    

    /**
     * Reads image volume data into a 4D double array. Slower than getImageDataSource (read one vol at a time then copies) 
     * but uses less memory.
     *
     */
    public double[][][][] readVolumeData() {


        double data[][][][];

        int comps = components();

        short ZZZ;
        
        // for 2D volumes, zdim may be 0
        ZZZ = ZDIM;
        if (dim[0] == 2)
            ZZZ = 1;
        
        // the data matrix to be populated
        data = new double[XDIM][YDIM][ZZZ][comps];
        
        // Note different input stream to header methods, because
        // EndianCorrectInputStream doesn't support unsigned types
        EndianNeutralDataInputStream ecs = getImageDataInputStream();

        try {
            for (int n = 0; n < comps; n++) {
                double[][][][] vol = readNextVolume(ecs);
                
                for (int k=0; k<ZZZ; k++) {
                    for (int j = 0; j < YDIM; j++) {
                        for (int i = 0; i < XDIM; i++) {
                            data[i][j][k][n] = vol[i][j][k][0];
                        }
                    }
                }
                
            }
            
            ecs.close();
        }
        catch(IOException e) {
            throw new LoggedException("Error reading volume data: " + e.getMessage());
        }

        return data;

    }


    /**
     * Get the header file name associated with the header object.
     *
     */
    public String getHeaderFilename() {
        return ds_hdrname;
    }


    /**
     * Get the data file name associated with the header object. May be the same as the header
     * file name.
     *
     */
    public String getDataFilename() {
        return ds_datname;
    }

    
    /**
     * Use the current header to write a scalar image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, with fields altered where necessary 
     * for writing the data.
     *
     * @param data will be checked against header data dimensions. The data array should be in the 
     * same voxel space as this image, to ensure a correct definition of physical space.
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings. 
     *
     *
     * @return the image header associated with the output data set
     */
    public ImageHeader writeScalarImage(double[][][] data, String fileRoot) {
        
        if (data.length != XDIM || data[0].length != YDIM || data[0][0].length != ZDIM) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + XDIM + 
                                      " " + YDIM + " " + ZDIM);
        }

        Nifti1Dataset nds = new Nifti1Dataset(this, fileRoot);

        nds.setDataTypeToSigned();

        double[][][][] fourD = new double[XDIM][YDIM][ZDIM][1];

        double min = 0.0;
        double max = 0.0;
        
        for (int i = 0; i < XDIM; i++) {
            for (int j = 0; j < YDIM; j++) {
                for (int k = 0; k < ZDIM; k++) { 
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

        nds.glmin = (int)Math.round(min);
        nds.glmax = (int)Math.round(max);

        nds.scl_slope = 1.0f;
        nds.scl_inter = 0.0f;

        // set fields appropriate for this data type
        nds.dim[0] = 3;
        nds.dim[4] = 0;
        nds.dim[5] = 0;
        nds.dim[6] = 0;
        nds.dim[7] = 0;

        nds.setPixDims(pixdim[0], pixdim[1], pixdim[2], pixdim[3], 0.0f, 0.0f, 0.0f, 0.0f);
        

        nds.intent_code = NIFTI_INTENT_NONE;
        nds.intent_name = new StringBuffer();
        
        for (int i = 0; i < 3; i++) {
            nds.dim[i+5] = 0;
        }

        nds.writeImage(fourD);

        return nds;

    }


    /**
     * Use the current header to write a vector image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, The new header is a copy of this one, 
     * with fields altered where necessary for writing vector data. 
     * <p>
     * If the number of components is 1, a 3D image will be written, otherwise a 5D image is written with
     * the fifth dimension set to the number of components.
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
   
        if (data.length != XDIM || data[0].length != YDIM || data[0][0].length != ZDIM) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + XDIM +                                      " " + YDIM + " " + ZDIM);
        }

        Nifti1Dataset nds = new Nifti1Dataset(this, fileRoot);

        nds.setDataTypeToSigned();

        // set fields appropriate for this data type

        short comps = (short)data[0][0][0].length;

        if (comps == 1) {
            // write image as 3D
            nds.dim[0] = 3;
            nds.dim[4] = 0;
            nds.dim[5] = 0;
            nds.dim[6] = 0;
            nds.dim[7] = 0;
            
        } 
        else {
            nds.dim[0] = 5;
            nds.dim[4] = 1;
            nds.dim[5] = comps;
            nds.dim[6] = 0;
            nds.dim[7] = 0;
        }

        nds.setPixDims(pixdim[0], pixdim[1], pixdim[2], pixdim[3], 0.0f, 0.0f, 0.0f, 0.0f);
        
        // could set intent to vector here if number of components is 3
        // AFAIK only the DT intent matters, because it affects how the data is read by ITK

        // 1007 seems to be a generic intent code for an N-vector, trouble is we have no way of knowing if 
        // a true vector is being written or just a composite of stuff, so we leave intent alone

        nds.intent_code = NIFTI_INTENT_NONE;
        nds.intent_name = new StringBuffer();

        // set gray levels to zero
        nds.glmin = 0;
        nds.glmax = 0;

        nds.scl_slope = 1.0f;
        nds.scl_inter = 0.0f;
        
        nds.writeImage(data);

        return nds;
    }


    /**
     * Use the current header to write an RGB image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, with fields altered where necessary 
     * for writing RGB data. 
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings. 
     *
     * @param red the red channel intensity, should be normalized between 0 and 255.
     * @param green the green channel intensity, should be normalized between 0 and 255.
     * @param blue the blue channel intensity, should be normalized between 0 and 255.
     *
     * @return the image header associated with the output data set
     */
    public ImageHeader writeRGB_Image(int[][][] red, int[][][] green, int[][][] blue, String fileRoot) {
        
        if (red.length != XDIM || red[0].length != YDIM || red[0][0].length != ZDIM) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + XDIM +                                      " " + YDIM + " " + ZDIM);
        }

        Nifti1Dataset nds = new Nifti1Dataset(this, fileRoot);

        nds.setDataTypeToSigned();

        // set fields appropriate for this data type
        nds.dim[0] = 3;
        nds.dim[4] = 0;
        nds.dim[5] = 0;
        nds.dim[6] = 0;
        nds.dim[7] = 0;

        nds.setPixDims(pixdim[0], pixdim[1], pixdim[2], pixdim[3], 0.0f, 0.0f, 0.0f, 0.0f);

        nds.datatype = NIFTI_TYPE_RGB24;
        
        nds.intent_code = NIFTI_INTENT_NONE;
        nds.intent_name = new StringBuffer();

        // set gray levels to zero
        nds.glmin = 0;
        nds.glmax = 0;

        nds.scl_slope = 1.0f;
        nds.scl_inter = 0.0f;


        double[][][][] data = new double[XDIM][YDIM][ZDIM][1];

        // put RGB bits into double form
        for (int i = 0; i < XDIM; i++) {
            for (int j = 0; j < YDIM; j++) {
                for (int k = 0; k < ZDIM; k++) { 
                    int rgb = blue[i][j][k] + 256 * (green[i][j][k] + 256 * red[i][j][k]);

                    data[i][j][k][0] = (double)rgb;
                }
            }
        }

        nds.writeImage(data);
        
        return nds;
    }



    /**
     * Use the current header to write an RGBA image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, with fields altered where necessary 
     * for writing RGBA data. 
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings. 
     *
     * @param red the red channel intensity, should be normalized between 0 and 255.
     * @param green the green channel intensity, should be normalized between 0 and 255.
     * @param blue the blue channel intensity, should be normalized between 0 and 255.
     * @param alpha the alph channel intensity, should be normalized between 0 and 1.
     *
     * @return the image header associated with the output data set
     */
    public ImageHeader writeRGBA_Image(int[][][] red, int[][][] green, int[][][] blue, double[][][] alpha, String fileRoot) {
        
        if (red.length != XDIM || red[0].length != YDIM || red[0][0].length != ZDIM) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + XDIM + " " + YDIM + " " + ZDIM);
        }

        Nifti1Dataset nds = new Nifti1Dataset(this, fileRoot);

        nds.setDataTypeToSigned();

        // set fields appropriate for this data type
        nds.dim[0] = 3;
        nds.dim[4] = 0;
        nds.dim[5] = 0;
        nds.dim[6] = 0;
        nds.dim[7] = 0;

        nds.setPixDims(pixdim[0], pixdim[1], pixdim[2], pixdim[3], 0.0f, 0.0f, 0.0f, 0.0f);

        nds.datatype = NIFTI_TYPE_RGBA32;
        
        nds.intent_code = NIFTI_INTENT_NONE;
        nds.intent_name = new StringBuffer();

        // set gray levels to zero
        nds.glmin = 0;
        nds.glmax = 0;

        nds.scl_slope = 1.0f;
        nds.scl_inter = 0.0f;


        double[][][][] data = new double[XDIM][YDIM][ZDIM][1];

        // put RGB bits into double form
        for (int i = 0; i < XDIM; i++) {
            for (int j = 0; j < YDIM; j++) {
                for (int k = 0; k < ZDIM; k++) { 
                    
                    int a = (int)(255.0 * alpha[i][j][k]);

                    if (a < 0) {
                        a = 0;
                    }
                    if (a > 255) {
                        a = 255;
                    }
                    
                    int rgb = a + 256 * (blue[i][j][k] + 256 * green[i][j][k] + 256 * 256 * red[i][j][k]);
                    
                    data[i][j][k][0] = (double)rgb;
                }
            }
        }

        nds.writeImage(data);
        
        return nds;
    }



    /**
     * 
     * Use the current header to write a tensor image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, with fields altered where necessary 
     * for writing tensor data. 
     *
     * @param data will be checked against header data dimensions. The data array should be in the 
     * same voxel space as this image, to ensure a correct definition of physical space. data[i][j][k] 
     * contains six components in upper-triangular order. The data will be written to disk in lower triangular order.
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings. 
     *
     * @return the image header associated with the output data set
     */
    public ImageHeader writeTensorImage(double[][][][] data, String fileRoot) {

        if (data.length != XDIM || data[0].length != YDIM || data[0][0].length != ZDIM) {
            throw new LoggedException("Attempted to write image inconsistent with header dimensions " + XDIM +                                      " " + YDIM + " " + ZDIM);
        }
        
        if (data[0][0][0].length != 6) {
            throw new LoggedException("DT output data must have six components");
        }
        

        Nifti1Dataset nds = new Nifti1Dataset(this, fileRoot);

        nds.setDataTypeToSigned();

        // set fields appropriate for this data type
        nds.dim[0] = 5;
        nds.dim[4] = 1;
        nds.dim[5] = 6;
        nds.dim[6] = 0;
        nds.dim[7] = 0;

        nds.setPixDims(pixdim[0], pixdim[1], pixdim[2], pixdim[3], 0.0f, 0.0f, 0.0f, 0.0f);

        nds.intent_code = NIFTI_INTENT_SYMMATRIX;
        nds.intent_name = new StringBuffer();

        // set gray levels to zero
        nds.glmin = 0;
        nds.glmax = 0;
        
        nds.scl_slope = 1.0f;
        nds.scl_inter = 0.0f;

        double[][][][] lowerTriangular = new double[XDIM][YDIM][ZDIM][6];
        
        for (int i = 0; i < XDIM; i++) {
            for (int j = 0; j < YDIM; j++) {
                for (int k = 0; k < ZDIM; k++) { 
                    
                    lowerTriangular[i][j][k][0] = data[i][j][k][0];
                    lowerTriangular[i][j][k][1] = data[i][j][k][1];
                    lowerTriangular[i][j][k][2] = data[i][j][k][3];
                    lowerTriangular[i][j][k][3] = data[i][j][k][2];
                    lowerTriangular[i][j][k][4] = data[i][j][k][4];
                    lowerTriangular[i][j][k][5] = data[i][j][k][5];
                    
                }
            }
        }
          
        nds.writeImage(lowerTriangular);

        return nds;

    }



    /**
     * Sets the data type to a supported Camino type.
     *
     * @throws LoggedException if no NIFTI datatype exists for the Camino equivalent.
     *
     */
    public void setDataType(String type) {
	if (type.equals("ubyte")) {
            datatype = NIFTI_TYPE_UINT8;
        }
	else if (type.equals("byte")) {
            datatype = NIFTI_TYPE_INT8;
        }
	else if (type.equals("short")) {
            datatype = NIFTI_TYPE_INT16;
        }
	else if (type.equals("ushort")) {
            datatype = NIFTI_TYPE_UINT16;
        }
	else if (type.equals("int")) {
            datatype = NIFTI_TYPE_INT32;
        }
	else if (type.equals("uint")) {
            datatype = NIFTI_TYPE_UINT32;
        }
	else if (type.equals("long")) {
            datatype = NIFTI_TYPE_INT64;
        }
	else if (type.equals("float")) {
            datatype = NIFTI_TYPE_FLOAT32;
        }
	else if (type.equals("double")) {
            datatype = NIFTI_TYPE_FLOAT64;
        }
        else {
            throw new LoggedException("Camino does not support this data type: " + decodeDatatype(datatype));
	}

        // set bitpix to standard values
	bitpix = (short)(bytesPerVoxel(datatype)*8);
    }


    public void setGzip(boolean gz) {
        
        if (gz == gzip) {
            return;
        }

        String fileRoot = null;

        if (gzip) {
            fileRoot = ds_hdrname.substring(0, ds_hdrname.length() - 7);
        }
        else { 
            fileRoot = ds_hdrname.substring(0, ds_hdrname.length() - 4);
        }

        setFilename(fileRoot, ds_is_nii, gz);
        
    }


    // END ImageHeader methods



    /**
     * Gets the Camino data type string.
     *
     * @throws LoggedException if no Camino equivalent exists for the NIFTI datatype.
     *
     */
    public String getDataType() {
	switch (datatype) {

	case NIFTI_TYPE_UINT8: return "char";

	case NIFTI_TYPE_INT8: return "byte";

	case NIFTI_TYPE_INT16: return "short";

	case NIFTI_TYPE_UINT16: return "ushort";

	case NIFTI_TYPE_INT32: return "int";

	case NIFTI_TYPE_UINT32: return "uint";

	case NIFTI_TYPE_INT64: return "long";

	case NIFTI_TYPE_FLOAT32: return "float";

	case NIFTI_TYPE_FLOAT64: return "double";

	default: throw new LoggedException("Camino does not support this data type: " + decodeDatatype(datatype));

	}
    }


    public static ImageHeader readHeader(String filename) throws IOException {
	Nifti1Dataset nds = new Nifti1Dataset(filename);

	return nds;
    }


    /**
     * Convert a 4D nii file to a series of 3D NIFTI volumes. 
     * <p>
     * Volumes will be written as <code>outputRoot</code>0000.ext, where 0000 is a fixed width 
     * volume index starting from 1, and ext is the same as the extension of the input data.
     */
    public static void convertTo3D(String filename, String outputRoot) throws IOException {
	Nifti1Dataset input = new Nifti1Dataset(filename);

	input.readHeader();

	int numVols = input.components();

	DecimalFormat df = new DecimalFormat("0000");

        EndianNeutralDataInputStream ecs = input.getImageDataInputStream();

	for (int n = 0; n < numVols; n++) {

            double[][][][] vol = input.readNextVolume(ecs);

            Nifti1Dataset singleVolume = new Nifti1Dataset();

	    singleVolume.copyHeader(input);

	    singleVolume.setDims(3, singleVolume.xDataDim(), singleVolume.yDataDim(), 
				 singleVolume.zDataDim(), 1, 0, 0, 0);
    
	    singleVolume.intent_code = NIFTI_INTENT_NONE;

	    singleVolume.intent_name = new StringBuffer("comp " + df.format(n+1));
	    
	    singleVolume.setFilename(outputRoot + df.format(n+1), input.ds_is_nii, input.gzip);
	   
	    singleVolume.writeImage(vol);

	}

    }




    /**
     * Check the qfac is consistent with the quaternion
     *
     * @return true if everything is OK, false otherwise. The quaternion will not be modified
     */
    public boolean checkQuaternion() {
        
        float b = quatern[0];
        float c = quatern[1];
        float d = quatern[2];

        if ( (b*b + c*c + d*d) > 1.00000001 ) {
            logger.warning("Quaternian does not have unit magnitude [b,c,d] = [" + b + ", " + c + ", " + 
                           d + "]");
            return false;
        }

        float aSq = 1.0f - (b*b + c*c + d*d);

        float a = 0.0f;

        if (aSq > 0.0) {
            a = (float)Math.sqrt( aSq );
        }
        
        RealMatrix R = new RealMatrix(3,3);

        R.entries[0][0] = a*a+b*b-c*c-d*d;
        R.entries[1][1] = a*a+c*c-b*b-d*d;
        R.entries[2][2] = a*a+d*d-c*c-b*b;

        R.entries[0][1] = 2.0*b*c-2*a*d;
        R.entries[0][2] = 2.0*b*d+2.0*a*c;
        R.entries[1][2] = 2.0*c*d-2*a*b;

        R.entries[1][0] = 2.0*b*c+2*a*d;
        R.entries[2][0] = 2.0*b*d-2.0*a*c;
        R.entries[2][1] = 2.0*c*d+2*a*b;
        
        double det = R.det();

        if (det < 0.0) {
            logger.warning("Determinant of quaternian is negative");
            return false;
        }

        //        logger.info("Quaternion rotation:\n" + R.toString());

        //        logger.info("qfac: " + qfac);

        return true;

    }


    /**
     * Sets the quaternion transform for this header.
     *
     * @param qformCode one of the NIFTI_XFORM codes. To use qform as the definition of physical space in Camino, this should 
     * be either NIFTI_XFORM_SCANNER_ANAT or NIFTI_XFORM_ALIGNED_ANAT. Even with these codes set, the sform takes precedence 
     * if it is defined with the appropriate code.
     *
     * @param qfactor either 1 or -1. Used to ensure that the  determinant of the quaternion rotation is positive.
     *
     * @param quaternBCD the b, c, d elements of the quaternion, as defined in the NIfTI-1 standard. The sum of squares
     * of  these must be less than 1.0. The complete quaternion has unit magnitude, ie a*a + b*b + c*c + d*d = 1.0.
     *
     * @param offset a 3D translation specifying the origin of physical space.
     *
     @see http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html
     */
    public void setQuaternion(short qformCode, short qfactor, float[] quaternBCD, float[] offset) {
        qform_code = qformCode;

        qfac = qfactor;

        pixdim[0] = (float)qfac;
        
        for (int i = 0; i < 3; i++) {
            quatern[i] = quaternBCD[i];
            qoffset[i] = offset[i];
        }

        if (!checkQuaternion()) {
            throw new LoggedException("Invalid quarternion parameters set");
        }
        
    }


    /**
     * Sets the sform transformation.
     *
     * @param sformCode one of the NIFTI_XFORM codes. If this is anything other than NIFTI_XFORM_UNKNOWN,
     * then sform is used in preference to qform in setting the physical space transformation.
     *
     * @param sformTrans a 4D affine matrix representing the sform
     */
    public void setSform(short sformCode, RealMatrix sformTrans) {
        sform_code = sformCode;

        for (int i = 0; i < 4; i++) {
            srow_x[i] = (float)sformTrans.entries[0][i]; 
            srow_y[i] = (float)sformTrans.entries[1][i]; 
            srow_z[i] = (float)sformTrans.entries[2][i]; 
        }
    }


    /**
     * Change the data type of an image from unsigned to signed
     *
     */
    private void setDataTypeToSigned() {
       
        switch (datatype) {

	case NIFTI_TYPE_UINT8: datatype = NIFTI_TYPE_INT8; break;

	case NIFTI_TYPE_UINT16: datatype = NIFTI_TYPE_INT16; break;

	case NIFTI_TYPE_UINT32: datatype = NIFTI_TYPE_INT32; break;

	case NIFTI_TYPE_UINT64: datatype = NIFTI_TYPE_INT64; break;


	}

    }
   


    /**
     * This main method would be much more useful if you could use it to edit a header, rather than just write
     * a new one. It's not clear how to do this with gzipped images. There is code on the web for
     * random access to a gzip file - maybe this could be used to allow users to edit header fields.
     *
     */
    public static void main(String[] args) {
		
	Nifti1Dataset nds = null;

	OutputManager.outputFile = "niftiheader.nii";
	CL_Initializer.CL_init(args);

	// get initial values from header first, then process other args
	for (int i = 0; i < args.length; i++) {
	    if (args[i].equals("-initfromheader")) {
		nds = new Nifti1Dataset(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
		try {
		    nds.readHeader();
		}
		catch (IOException ex) {
		    System.out.println("\nCould not read header file for "+args[i+1]+": "+ex.getMessage());
		}
		
	    }
	}

	if (nds == null) {
            nds = new Nifti1Dataset();
        }

	for (int i = 0; i < args.length; i++) {
	    
	    if (args[i].equals("-readheader")) {
		
		nds = new Nifti1Dataset(args[i+1]);

		CL_Initializer.markAsParsed(i,2);

		try {
		    nds.readHeader();
		}
		catch (IOException ex) {
		    System.out.println("\nCould not read header file for "+args[i+1]+": "+ex.getMessage());
		}

                nds.checkQuaternion();

		nds.printHeader();
		return;
		
	    }
	    if (args[i].equals("-datatype")) {
		nds.setDataType(args[i+1]);
		CL_Initializer.markAsParsed(i, 2);
	    }
	    if (args[i].equals("-intent")) {
		nds.intent_code = Short.parseShort(args[i+1]);
		CL_Initializer.markAsParsed(i, 2);
	    }
	    if (args[i].equals("-networkbyteorder")) {
		nds.big_endian = true;
		CL_Initializer.markAsParsed(i);
	    }
	    if (args[i].equals("-intelbyteorder")) {
		nds.big_endian = false;
		CL_Initializer.markAsParsed(i);
	    }
	    if (args[i].equals("-scaleslope")) {
		nds.scl_slope = Float.parseFloat(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
	    if (args[i].equals("-scaleinter")) {
		nds.scl_inter = Float.parseFloat(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
	    if (args[i].equals("-dims")) {
		
		// number of dimensions
		short numDims = Short.parseShort(args[i+1]);

		short[] dims = new short[8];

                dims[0] = numDims;
                
		for (int d = 0; d < numDims; d++) {
		    dims[d+1] = Short.parseShort(args[i+2+d]);
		} 

		nds.setDims(dims);

		CL_Initializer.markAsParsed(i, numDims + 2);

	    }
	    if (args[i].equals("-pixdims")) {
                float[] pixdims = new float[8];

		// 1 extra arg for qfac
		int numDims = Integer.parseInt(args[i+1]) + 1;

		for (int d = 0; d < numDims; d++) {
                    pixdims[d] = Float.parseFloat(args[i+2+d]);
                } 

                nds.setPixDims(pixdims);

		CL_Initializer.markAsParsed(i, numDims + 2);
            }
            if (args[i].equals("-qform")) {

                nds.qform_code = Short.parseShort(args[i+1]);
                
                // qfac should be 1 or -1 to give rotation matrix
                // a positive determinant
                float qfac = Float.parseFloat(args[i+2]);

                if (qfac < 0.0f) { 
                    qfac = -1.0f;
                }
                else {
                    qfac = 1.0f;
                }

                nds.pixdim[0] = qfac;
                
                // internal variable qfac is set when header is read
                nds.qfac = (short)qfac;
                                    
                // quatern b,c,d params
                nds.quatern[0] = Float.parseFloat(args[i+3]);
                nds.quatern[1] = Float.parseFloat(args[i+4]);
                nds.quatern[2] = Float.parseFloat(args[i+5]);
                
                // translation params
                nds.qoffset[0] = Float.parseFloat(args[i+6]);
                nds.qoffset[1] = Float.parseFloat(args[i+7]);
                nds.qoffset[2] = Float.parseFloat(args[i+8]);

                CL_Initializer.markAsParsed(i,9);

                nds.checkQuaternion();

            }
            if (args[i].equals("-sform")) {
                // sform_code, then row-wise elements of srow
                nds.sform_code = Short.parseShort(args[i+1]);
                
                for (int j = 0; j < 4; j++) {
                    nds.srow_x[j] = Float.parseFloat(args[i+2+j]);
                }
                for (int j = 0; j < 4; j++) {
                    nds.srow_y[j] = Float.parseFloat(args[i+6+j]);
                }
                for (int j = 0; j < 4; j++) {
                    nds.srow_z[j] = Float.parseFloat(args[i+10+j]);
                }
                
                CL_Initializer.markAsParsed(i,14);
            }
            if (args[i].equals("-dt2symmatrix")) {
                nds.intent_code = NIFTI_INTENT_SYMMATRIX;
                nds.dim[0] = 5;
                nds.dim[4] = 1;
	        nds.dim[5] = 6;
                CL_Initializer.markAsParsed(i,1);
            }    
	    
	}

	CL_Initializer.checkParsing(args);

        boolean gzip = false;
        boolean nii = false;

        if (OutputManager.outputFile.endsWith(".gz")) {
            gzip = true;
            OutputManager.outputFile = OutputManager.outputFile.substring(0, OutputManager.outputFile.length() - 3);
        }

        if (OutputManager.outputFile.endsWith(".nii")) {
            nii = true;
        }

        OutputManager.outputFile = OutputManager.outputFile.substring(0, OutputManager.outputFile.length() - 4);

	try {
	    nds.setFilename(OutputManager.outputFile, nii, gzip);

            if (new File(nds.getHeaderFilename()).exists()) {
                // This capability could be added with a new method to overwrite the existing header - barring any extensions
                throw new LoggedException("This program cannot modify existing headers (image data would be overwritten)");
            }

	    nds.writeHeader();
	}
	catch(IOException e) {
	    throw new LoggedException(e);
	}
    }


    public boolean lowerTriangularSymmMatrix() {
        return true;
    }


    // END Camino extensions

    

}



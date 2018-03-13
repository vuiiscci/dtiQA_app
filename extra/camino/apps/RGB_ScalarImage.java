package apps;

import java.io.*;

import java.util.*;
import java.util.logging.Logger;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;
import sphfunc.*;
import fitters.*;



/**
 *
 * Produces a scalar-modulated RGB image, colored by the tensor eigenvector or ODF peak direction.
 * The application element reads command line options and outputs the RGB image.
 * <p>
 * This class also serves as an image object for PD_OrientationViewer.
 * 
 * @author Philip Cook
 * @version $Id$
 *  
 */
public class RGB_ScalarImage {

    // 16 Mb output buffer
    public static final int BUFFERSIZE = 16777216;

    // array of vectors. Allow direct access for speed
    protected Vector3D[][][][] vectors = null;

    protected double[][][] scalarVol = null;

    protected double[][][] normScalarVol = null;

    protected ImageHeader scalarHeader;

    // between 0.0 and 1.0
    private double[][][] red = null;
    private double[][][] green = null;
    private double[][][] blue = null;
    private double[][][] alpha = null;


    private boolean haveAlpha = false;

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.RGB_ScalarImage");

    
    protected int xDataDim = 0;
    protected int yDataDim = 0;
    protected int zDataDim = 0;

    protected double xVoxelDim = 1.0;
    protected double yVoxelDim = 1.0;
    protected double zVoxelDim = 1.0;
     

    // range of scalar values for normScalarVol calculation. 
    // Anything > max set to 1, < min set to 0.
    protected double minScalarValue = 0.0;
    protected double maxScalarValue = 0.0;
    
    // gamma correction factor for the background scalar image
    private double gsGamma = 1.0;

    // gamma correction factor applied independently to R, G, B channels
    private double rgbGamma = 1.0;

   
    /**
     * @param vecs an array of vectors, where the last index is the PDs in each voxel.
     * @param hdr an image header used for writing the RGB image, if required.
     * @param scalars scalar values used as the grayscale contrast
     * @param minScalar minimum scalar for greyscale contrast, values below this clamped to black
     * @param maxScalar maximum scalar for greyscale contrast, values above this clamped to white
     *
     */    
    public RGB_ScalarImage(Vector3D[][][][] vecs, ImageHeader hdr, double[][][] scalars,
			   double minScalar, double maxScalar) {
	
	vectors = vecs;
        
       	scalarVol = scalars;
        
        scalarHeader = hdr;
	
	xDataDim = vecs.length;
	yDataDim = vecs[0].length;
	zDataDim = vecs[0][0].length;

	xVoxelDim = hdr.xVoxelDim();
	yVoxelDim = hdr.yVoxelDim();
	zVoxelDim = hdr.zVoxelDim();

	setNormalizedScalars(minScalar, maxScalar);		

	calculateRGB();
	
    }




    /**
     * Writes this image, determines output type from extension (either .mha, .mhd, .nii, .nii.gz, .vtk).
     * <p>
     * Because of ongoing problems reading RGB nii images in Paraview, we allow the user to write a different 
     * file format to the input data. This may cause issues with image orientation.
     *
     */
    public void writeImage(String file) throws IOException {

        logger.info("Saving RGB image as " + file);

	if (file.endsWith(".mha") || file.endsWith(".mhd")) {
            
            String fileRoot = file.substring(0, file.length() - 4);
            
	    writeMeta(fileRoot);
	}
	else if (file.endsWith(".vtk")) {

            String fileRoot = file.substring(0, file.length() - 4);

	    writeVTK(fileRoot);
	}
	else if (file.endsWith(".nii") || file.endsWith(".nii.gz")) {

            String fileRoot = null;

            if (file.endsWith(".gz")) {
                fileRoot = file.substring(0, file.length() - 7);
            }
            else {
                fileRoot = file.substring(0, file.length() - 4);
            }

	    writeNifti(fileRoot);
	}

	else {
	    logger.warning("No recognized extension to file " + file + " - image not written");
	}

        logger.info("RGB image saved");
    }


    /**
     * Writes the RGB image in Meta I/O format. Data is binary rgb byte triplets for each voxel.
     *
     *
     */
    public void writeNifti(String fileRoot) throws IOException {


        // Hack, eventually we will have one method, and will write output in the format
        // of the scalar header, and we'll deprecate VTK / Analyze and support alpha in all types

        Nifti1Dataset hdr = null;

        if (scalarHeader != null && scalarHeader instanceof Nifti1Dataset) {
            hdr = (Nifti1Dataset)scalarHeader;
        }
        else {

            hdr = new Nifti1Dataset();

            hdr.setDims(3, xDataDim, yDataDim, zDataDim, 1, 0, 0, 0);

            hdr.setPixDims(1.0, xVoxelDim, yVoxelDim, zVoxelDim, 0.0, 0.0, 0.0, 0.0);

            logger.warning("Header template is not NIfTI, physical space information may be lost");
        }
        

        int[][][] red_int = new int[xDataDim][yDataDim][zDataDim];
        int[][][] green_int = new int[xDataDim][yDataDim][zDataDim];
        int[][][] blue_int = new int[xDataDim][yDataDim][zDataDim];
        
        for (int k = 0; k < zDataDim; k++) { 
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
                    int[] rgb = rgbTriplet(i,j,k);
                    
                    red_int[i][j][k] = rgb[0];
                    green_int[i][j][k] = rgb[1];
                    blue_int[i][j][k] = rgb[2];
                    
                }
            }
        }
        
    
        if (haveAlpha) {
            
            logger.info("Writing RGBA image, alpha is 1 in foreground and 0 in background");

            Nifti1Dataset niiHdr = (Nifti1Dataset)hdr;
            niiHdr.writeRGBA_Image(red_int, green_int, blue_int, alpha, fileRoot);

        }
        else {
            hdr.writeRGB_Image(red_int, green_int, blue_int, fileRoot);
        }
            



    }


    /**
     * Writes the RGB image in Meta I/O format. Data is binary rgb byte triplets for each voxel.
     *
     *
     */
    public void writeMeta(String fileRoot) throws IOException {


        MetaImageHeader hdr = null;

        if (scalarHeader != null && scalarHeader instanceof MetaImageHeader) {
            hdr = (MetaImageHeader)scalarHeader;
        }
        else {

            hdr = new MetaImageHeader();

            hdr.setDataDims(new int[] {xDataDim, yDataDim, zDataDim});

            hdr.setVoxelDims(new double[] {xVoxelDim, yVoxelDim, zVoxelDim});

            logger.warning("Header template is not Meta, physical space information may be lost");

        }
            
        int[][][] red_int = new int[xDataDim][yDataDim][zDataDim];
        int[][][] green_int = new int[xDataDim][yDataDim][zDataDim];
        int[][][] blue_int = new int[xDataDim][yDataDim][zDataDim];
        
        // use meta header to write image
        for (int k = 0; k < zDataDim; k++) { 
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
                    int[] rgb = rgbTriplet(i,j,k);
                    
                    red_int[i][j][k] = rgb[0];
                    green_int[i][j][k] = rgb[1];
                    blue_int[i][j][k] = rgb[2];
                    
                }
            }
        }
        

        hdr.writeRGB_Image(red_int, green_int, blue_int, fileRoot);

    }


    /**
     * Writes a VTK 2.0 (ie, not the new XML) file. The RGB triplets are stored as scalar values.
     *
     */
    public void writeVTK(String fileRoot) throws IOException {

	String fileName = null;

	if (fileRoot.endsWith(".vtk")) {
	    fileName = fileRoot;
	}
	else {
	    fileName = fileRoot + ".vtk"; 
	}

	DataOutputStream dout = 
	    new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fileName), BUFFERSIZE));


	
	dout.write(new String("# vtk DataFile Version 2.0\n").getBytes("US-ASCII"));
	dout.write(new String("Camino RGB-scalar image\n").getBytes("US-ASCII"));
	dout.write(new String("BINARY\n").getBytes("US-ASCII"));
	dout.write(new String("DATASET STRUCTURED_POINTS\n").getBytes("US-ASCII"));

	String dataDimString = "DIMENSIONS " + xDataDim + " " + yDataDim + " " + zDataDim + "\n"; 

	dout.write(new String(dataDimString).getBytes("US-ASCII"));

	String voxelDimString = "SPACING " + xVoxelDim + " " + yVoxelDim + " " + zVoxelDim + "\n"; 

	dout.write(voxelDimString.getBytes("US-ASCII"));

	dout.write(new String("ORIGIN 0 0 0\n").getBytes("US-ASCII"));
	dout.write(new String("POINT_DATA " + xDataDim * yDataDim * zDataDim + "\n").getBytes("US-ASCII"));

	dout.write(new String("COLOR_SCALARS scalars 3\n").getBytes("US-ASCII"));

    	// now write data
	for (int k = 0; k < zDataDim; k++) { 
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    int rgb = rgbIndex(i,j,k);
		    
		    dout.writeByte( (byte)(rgb >> 16) );
		    
		    dout.writeByte( (byte)( (rgb >> 8) & 0xff ) );

		    dout.writeByte( (byte)(rgb & 0xff) );
		    
		}
	    }
	}
	
	dout.close();


    }


    /**
     * @return the 24-bit rgb triplet for the specified volumetric index.
     *
     */
    protected final int rgbIndex(int i, int j, int k) {


	// gamma first, then scale
	double r = Math.pow(red[i][j][k], rgbGamma);
	double g = Math.pow(green[i][j][k], rgbGamma);
	double b = Math.pow(blue[i][j][k], rgbGamma);
	
	double gsScalar = Math.pow(normScalarVol[i][j][k], gsGamma);
	
	if (!(gsScalar >= 0.0)) {
	    gsScalar = 0.0;
	}
	if (gsScalar > 1.0) {
	    gsScalar = 1.0;
	}

	int ri = 0;
	int gi = 0;
	int bi = 0;
	
	ri = (int)(255.0 * gsScalar * r);
	gi = (int)(255.0 * gsScalar * g);
	bi = (int)(255.0 * gsScalar * b);
    
	return bi + 256 * (gi + 256 * ri);

    }


    /**
     * @return the 24-bit rgb triplet for the specified volumetric index.
     *
     */
    protected final int[] rgbTriplet(int i, int j, int k) {


	// gamma first, then scale
	double r = Math.pow(red[i][j][k], rgbGamma);
	double g = Math.pow(green[i][j][k], rgbGamma);
	double b = Math.pow(blue[i][j][k], rgbGamma);
	
	double gsScalar = Math.pow(normScalarVol[i][j][k], gsGamma);
	
	if (!(gsScalar >= 0.0)) {
	    gsScalar = 0.0;
	}
	if (gsScalar > 1.0) {
	    gsScalar = 1.0;
	}
	
	return new int[] {(int)(255.0 * gsScalar * r), (int)(255.0 * gsScalar * g), (int)(255.0 * gsScalar * b)};
    
    }



    /**
     * Calculates scalar range over which grey levels are mapped. Anything less than the minimum
     * is set to black, anything greater than the maximum is set to white. 
     *
     * @param chop the fraction of voxels to exclude. The smallest and largest (chop * 100) % of
     * the scalar values will be mapped to black and white respectively.
     */
    protected void calculateScalarRange(double chop) {
	int counter = 0;

	double[] scalarFlat = new double[xDataDim * yDataDim * zDataDim];

	for (int k = 0; k < zDataDim; k++) { 
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
                   
		    scalarFlat[counter++] = scalarVol[i][j][k];
		    
		}
	    }
	}

	// calculate scalar range from data
	Arrays.sort(scalarFlat);
        
	int minIndex = (int)(chop * scalarFlat.length);

	int maxIndex = scalarFlat.length - minIndex - 1;
        
	minScalarValue = scalarFlat[minIndex];
	maxScalarValue = scalarFlat[maxIndex];

	if (maxScalarValue < Double.MAX_VALUE) {
	    // we're OK
	}
	else {
	    // NaN or Infinity
	    while (maxIndex >= 0 && !(maxScalarValue <= Double.MAX_VALUE)) {
		maxIndex--;
		maxScalarValue = scalarFlat[maxIndex];
	    }
	}


	if (minScalarValue > -1.0 * Double.MAX_VALUE) {
	    // we're OK
	}
	else {
	    // -Infinity
	    while (minIndex < (scalarFlat.length  - 1) && !(minScalarValue >= -1.0 * Double.MAX_VALUE)) {
		minIndex++;
		minScalarValue = scalarFlat[minIndex];
	    }
	}
	

    }


    /**
     * Sets scalar range from the range present in the data, after excluding outliers.
     * Anything less than the minimum
     * is set to black, anything greater than the maximum is set to white. 
     *
     * @param chop the fraction of voxels to exclude. The smallest and largest (chop * 100) % of
     * the scalar values will be mapped to black and white respectively.
     */
    public void setNormalizedScalars(double chop) {

	if ( !( (chop >= 0.0) && (chop < 1.0) ) ) {
	    throw new LoggedException("Fraction of scalar values to exclude was : " + chop);
	}

	calculateScalarRange(chop);
	setNormalizedScalars();
    }

    
    /**
     * Sets scalar range between the specified minimum (maps to black) and maximum 
     * (maps to white).
     * Values outside the range are clamped. If the two parameters are equal, the scalar range
     * is computed from the data.
     */
    public void setNormalizedScalars(double minScalar, double maxScalar) {

	// if the values are the same, we calculate from the data 
	if (minScalar == maxScalar) {
	    calculateScalarRange(0.005);
	}
	else {
	    minScalarValue = minScalar;
	    maxScalarValue = maxScalar;
	}

	setNormalizedScalars();
    }


    /**
     * Sets normalized scalars according to the current scalar range. 
     */
    private void setNormalizedScalars() {

	normScalarVol = new double[xDataDim][yDataDim][zDataDim];

	for (int k = 0; k < zDataDim; k++) { 
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    
		    if (scalarVol[i][j][k] < minScalarValue) {
			normScalarVol[i][j][k] = 0.0;
		    }
		    else if (scalarVol[i][j][k] > maxScalarValue) {
			normScalarVol[i][j][k] = 1.0;
		    }
		    else {
 			if ((maxScalarValue - minScalarValue) != 0) {
			    normScalarVol[i][j][k] = 
				(scalarVol[i][j][k] - minScalarValue) / (maxScalarValue - minScalarValue);   
 			}
 			else {
 			    normScalarVol[i][j][k] = 1.0;
 			}
 		    }
		    
		    
		}
	    }
	}
    }


    public static void main(String[] args) {

	CL_Initializer.inputDataType = "double";
        CL_Initializer.maxTensorComponents = 1;

	CL_Initializer.inputModel = "dteig";

	CL_Initializer.numPDsIO = -1;

        CL_Initializer.CL_init(args);
        
	String picoPDF = "bingham";

        int xDataDim = 0;
        int yDataDim = 0;
        int zDataDim = 0;

        double xVoxelDim = 0.0;
        double yVoxelDim = 0.0;
        double zVoxelDim = 0.0;

        double minScalar = 0.0;
        double maxScalar = 0.0;

        String scalarFile = null;


	double gsGamma = 1.0;
	double rgbGamma = 1.0;

	int eigIndex = 0;

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-scalarrange")) { 
                minScalar = Double.parseDouble(args[i+1]);
                maxScalar = Double.parseDouble(args[i+2]);
                CL_Initializer.markAsParsed(i, 3);
            }
            if (args[i].equals("-scalarfile")) {
                scalarFile = args[i+1];
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-gsgamma")) {
		gsGamma = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-rgbgamma")) {
		rgbGamma = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
	    if (args[i].equals("-e1")) {
                eigIndex = 0;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-e2")) {
                eigIndex = 1;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-e3")) {
                eigIndex = 2;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-pdf")) {
                picoPDF = args[i+1];
                CL_Initializer.markAsParsed(i, 2);
            }

        }


        CL_Initializer.checkParsing(args);


	if (CL_Initializer.numPDsIO < 0) {
	    if (CL_Initializer.inputModel.equals("pds")) {
		// default is 3 for sfpeaks
		CL_Initializer.numPDsIO = 3;
	    }
	    else {
		// default is 1 for everything else
		CL_Initializer.numPDsIO = 1;
	    }
	}

        if (CL_Initializer.headerTemplateFile == null) {
            if (scalarFile != null) {
                CL_Initializer.headerTemplateFile = scalarFile;
            }
            else if (CL_Initializer.bgMaskFile != null) {
                CL_Initializer.headerTemplateFile = CL_Initializer.bgMaskFile;
            }
        }        

        CL_Initializer.initInputSpaceAndHeaderOptions();
        
	xDataDim = CL_Initializer.dataDims[0];
	yDataDim = CL_Initializer.dataDims[1];
	zDataDim = CL_Initializer.dataDims[2];
	
	xVoxelDim = CL_Initializer.voxelDims[0];
	yVoxelDim = CL_Initializer.voxelDims[1];
	zVoxelDim = CL_Initializer.voxelDims[2];


	double[][][] scalarVol = null;

        ImageHeader ih = CL_Initializer.headerTemplate;
        
        CL_Initializer.initMaskSource();
        
	if (scalarFile != null) {
            
	    DataSource scalars = null;
            
	    if (ImageHeader.imageExists(scalarFile)) {
                
                ImageHeader scalarHdr = null;
                
		try {
		    scalarHdr = ImageHeader.readHeader(scalarFile);
		}
		catch (IOException e) {
		    throw new LoggedException(e);
                    
		}
                
                if (!CL_Initializer.headerTemplate.sameSpace(scalarHdr)) {
                    throw new LoggedException("Scalar image space does not match input space");
                }
                
		scalars = scalarHdr.getImageDataSource();
	    }
	    else {
		scalars = new VoxelOrderDataSource(scalarFile, 1, CL_Initializer.inputDataType);
	    }
            
	    scalarVol = new double[xDataDim][yDataDim][zDataDim];
            
	    
	    for (int k = 0; k < zDataDim; k++) { 
		for (int j = 0; j < yDataDim; j++) {
		    for (int i = 0; i < xDataDim; i++) {
			scalarVol[i][j][k] = scalars.nextVoxel()[0];
		    }
		}
	    }
	}

	RGB_ScalarImage image = null;

	if (CL_Initializer.inputModel.equals("dteig")) {
            DataSource data = 
		ExternalDataSource.getDataSource(CL_Initializer.inputFile, 12 * CL_Initializer.maxTensorComponents,
                                                 CL_Initializer.inputDataType);
	    
	    image = imageFromTensorEigenSys(data, CL_Initializer.headerTemplate, scalarVol, 
					    minScalar, maxScalar, eigIndex);
	}
	else if (CL_Initializer.inputModel.equals("pds")) {

            DataSource data = ExternalDataSource.getDataSource(CL_Initializer.inputFile, 
                                                               6 + 8 * CL_Initializer.numPDsIO, 
                                                               CL_Initializer.inputDataType);

	    image = imageFromSphFuncPDs(data, CL_Initializer.headerTemplate, scalarVol, minScalar, maxScalar);
	}
	else if (CL_Initializer.inputModel.equals("pico")) {

	    int paramsPerPD = 1;

	    if (picoPDF.equals("bingham")) {
		paramsPerPD = 2;
	    }
	    if (picoPDF.equals("acg")) {
		paramsPerPD = 3;
	    }

            
            DataSource data = 
		ExternalDataSource.getDataSource(CL_Initializer.inputFile, 
                                                 1 + (10 + paramsPerPD) * CL_Initializer.numPDsIO, 
                                                 CL_Initializer.inputDataType);


	    image = RGB_ScalarImage.imageFromPICoPDFs(data, CL_Initializer.numPDsIO, paramsPerPD,
                                                      CL_Initializer.headerTemplate,
                                                      scalarVol, minScalar, maxScalar);

        }
	else if (CL_Initializer.inputModel.equals("ballstick")) {
            
            DataSource data = 
               ExternalDataSource.getDataSource(CL_Initializer.inputFile, 7, CL_Initializer.inputDataType);
            

	    image = RGB_ScalarImage.imageFromBallStick(data, CL_Initializer.headerTemplate,
						       scalarVol, minScalar, maxScalar);
            
	}

        else {
            throw new misc.LoggedException("Unrecognized input model " + CL_Initializer.inputModel);
        }
	
	// apply gamma correction
	image.setScalarGamma(gsGamma);
	image.setRGB_Gamma(rgbGamma);


        // set up mask if we have it
        if (CL_Initializer.bgMask != null) {
            image.haveAlpha = true;

            image.alpha = new double[xDataDim][yDataDim][zDataDim];

            for (int k = 0; k < zDataDim; k++) { 
                for (int j = 0; j < yDataDim; j++) {
                    for (int i = 0; i < xDataDim; i++) {
                        if (CL_Initializer.bgMask.nextVoxel()[0] > 0.0) {
                            image.alpha[i][j][k] = 1.0;
                        }
                        else {
                            image.alpha[i][j][k] = 0.0;
                        }
                    }
                }
            }
        }



	
	try {
            image.writeImage(OutputManager.outputFile);
        }
	catch (IOException e) {
	    throw new LoggedException(e);
	}
	
    }



    /**
     * Sets RGB values for display. The RGB value for each voxel is taken from the vectors.
     *
     */
    protected void calculateRGB() {
	
	if (red == null) {
	    red = new double[xDataDim][yDataDim][zDataDim];
	    green = new double[xDataDim][yDataDim][zDataDim];
	    blue = new double[xDataDim][yDataDim][zDataDim];
	}
	
	for (int k = 0; k < zDataDim; k++) { 
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    double r = Math.abs(vectors[i][j][k][0].x);
		    double g = Math.abs(vectors[i][j][k][0].y);
		    double b = Math.abs(vectors[i][j][k][0].z);

		    if (vectors[i][j][k].length > 1) {
			Vector3D average = 
			    new Vector3D(Math.abs(vectors[i][j][k][0].x) + Math.abs(vectors[i][j][k][1].x),
					 Math.abs(vectors[i][j][k][0].y) + Math.abs(vectors[i][j][k][1].y),
					 Math.abs(vectors[i][j][k][0].z) + 
					 Math.abs(vectors[i][j][k][1].z)).normalized();
			
			r = average.x;
			g = average.y;
			b = average.z;

		    }

		    red[i][j][k] = r;
		    green[i][j][k] = g;
		    blue[i][j][k] = b;
		    
		}
	    }
	}
	
    }


    /**
     * Read image from ODF peak input. The image has a maximum of two directions per voxel.
     *
     * @param scalars if null, then the Hessian trace of the first peak is used.
     */
    public static RGB_ScalarImage imageFromSphFuncPDs(DataSource data, ImageHeader ih, double[][][] scalars, 
						      double minScalar, double maxScalar) {

	int xDataDim = ih.xDataDim();
	int yDataDim = ih.yDataDim();
	int zDataDim = ih.zDataDim();
	
	double[][][] hessianTr = new double[xDataDim][yDataDim][zDataDim];
	
	Vector3D[][][][] vecs = new Vector3D[xDataDim][yDataDim][zDataDim][];

	for (int k = 0; k < zDataDim; k++) { 
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    double[] voxel = data.nextVoxel();
                    
		    Vector3D v1 = null;
		    Vector3D v2 = null;

		    v1 = new Vector3D(voxel[6], voxel[7], voxel[8]);

                    hessianTr[i][j][k] = Math.abs(voxel[10] + voxel[13]);

		    if (hessianTr[i][j][k] > 0.0) {
			// it's OK
		    }
		    else {
			// negative Hessian, or NaN
			hessianTr[i][j][k] = 0.0;
		    }

                   
                    if (voxel[2] >= 2.0) {
			
			// get second peak
                        v2 = new Vector3D(voxel[14], voxel[15], voxel[16]);

                    }
		    
		    if (v2 == null) {
			vecs[i][j][k] = new Vector3D[] {v1};
		    }
		    else {
			vecs[i][j][k] = new Vector3D[] {v1, v2};
		    }

                }
	    }
	}


	if (scalars == null) {
	    scalars = hessianTr;
	}
	
	return new RGB_ScalarImage(vecs, ih, scalars, minScalar, maxScalar);

    }


    /**
     * Read an image from PICoPDFs. 
     * The image has a maximum of two directions per voxel.
     *
     * @param maxPDs maximum number of PDs in the voxel. Only two will be displayed.
     * @param paramsPerPD scalar PICo parameters per PD.
     * @param scalars if null, then the sum of the concentration parameters is used. 
     */
    public static RGB_ScalarImage imageFromPICoPDFs(DataSource data, int maxPDs, 
						    int paramsPerPD, ImageHeader ih, double[][][] scalars, double minScalar,
						   double maxScalar) {

	return imageFromPICoPDFs(data, maxPDs, paramsPerPD, ih, scalars, minScalar, maxScalar, 0);
    }

    /**
     * Read an image from PICoPDFs. 
     * The image has a maximum of two directions per voxel.
     *
     * @param maxPDs maximum number of PDs in the voxel. Only two will be displayed.
     * @param paramsPerPD scalar PICo parameters per PD.
     * @param scalars if null, then the sum of the concentration parameters is used. 
     * @param eigIndex which of the scatter matrix eigenvectors to use (0-2).
     */
    public static RGB_ScalarImage imageFromPICoPDFs(DataSource data, int maxPDs, 
						    int paramsPerPD, ImageHeader ih, double[][][] scalars, double minScalar,
						   double maxScalar, int eigIndex) {
	
	int xDataDim = ih.xDataDim();
	int yDataDim = ih.yDataDim();
	int zDataDim = ih.zDataDim();

	Vector3D[][][][] vecs = new Vector3D[xDataDim][yDataDim][zDataDim][];
	
	double[][][] sumK = new double[xDataDim][yDataDim][zDataDim];
	
	
	for (int k = 0; k < zDataDim; k++) { 
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    
		    double[] voxel = data.nextVoxel();

		    int numPDs = (int)voxel[0];
		    
                    Vector3D v1 = new Vector3D(voxel[2 + 3*eigIndex], voxel[3 + 3*eigIndex], voxel[4 + 3*eigIndex]);
                    Vector3D v2 = null;

		    for (int p = 0; p < paramsPerPD; p++) {
			sumK[i][j][k] += Math.abs(voxel[11+p]);
		    }

		    if (numPDs >= 2) {		
			int start = 11 + paramsPerPD;
			
			v2 = new Vector3D(voxel[start + 1 + 3*eigIndex], voxel[start + 2 + 3*eigIndex], 
					  voxel[start + 3 + 3*eigIndex]);

			double tmp = 0.0;

			for (int p = 0; p < paramsPerPD; p++) {
			    tmp += Math.abs(voxel[start + 10 + p]);
			}
			
			sumK[i][j][k] = (sumK[i][j][k] + tmp) / 2.0;
                    }

		    if (v2 == null) {
			vecs[i][j][k] = new Vector3D[] {v1};
		    }
		    else {
			vecs[i][j][k] = new Vector3D[] {v1, v2};
		    }
		    
		}
	    }
	}
	
	if (scalars == null) {
	    scalars = sumK;
	}
	
	return new RGB_ScalarImage(vecs, ih, scalars, minScalar, maxScalar);
	
    }


    /**
     * Read image from Ball and stick input. The image has a maximum of one direction per voxel.
     *
     * @param scalars if null, then the mixing fraction f is used.
     */
    public static RGB_ScalarImage imageFromBallStick(DataSource data, ImageHeader ih, double[][][] scalars, 
                                                     double minScalar, double maxScalar) {

	int xDataDim = ih.xDataDim();
	int yDataDim = ih.yDataDim();
	int zDataDim = ih.zDataDim();
	
	double[][][] f = new double[xDataDim][yDataDim][zDataDim];
	
	Vector3D[][][][] vecs = new Vector3D[xDataDim][yDataDim][zDataDim][];

	for (int k = 0; k < zDataDim; k++) { 
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    double[] voxel = data.nextVoxel();
                    
		    Vector3D v1 = null;

		   vecs[i][j][k] = new Vector3D[] {new Vector3D(voxel[4], voxel[5], voxel[6])};

                    f[i][j][k] = Math.abs(voxel[3]);

		    if (f[i][j][k] >= 0.0) {
			// it's OK
		    }
		    else {
			f[i][j][k] = 0.0;
		    }

                }
	    }
	}


	if (scalars == null) {
	    scalars = f;
	}
	
	return new RGB_ScalarImage(vecs, ih, scalars, minScalar, maxScalar);

    }


    /**
     * Read image from fitting a compartment model. The image has a
     * maximum of one direction per voxel.  The method assumes the
     * convention that element three of the set of values per voxel is
     * a sensible parameter to use as the scalar background, eg
     * intra-cellular volume fraction, if parameter scalars is not
     * specified.  It also assumes that elements N-2 and N-1 of the N
     * values encode the fibre orientation as spherical polar
     * coordinates.
     *
     * @param if scalars is null, then the mixing fraction f is used.
     */
    public static RGB_ScalarImage imageFromCompartmentModel(ImageHeader ih, 
                                                            double[][][] scalars, double minScalar, double maxScalar) {

        // This is just to figure out the number of values in each
        // voxel.  Recreate the fitter used to generate the data:
        // FitAlgorithm.LM just says there is one set of parameters
        // per voxel; the choice of noise model is irrelevant but
        // required by the constructor.
        FitModel fm = FitModel.getFitModel(CL_Initializer.fitModel);
        // Constructing the fitter requires a scheme object to be
        // defined, although this program doesn't need it so we 
        // create a dummy
        if(CL_Initializer.imPars == null)
            CL_Initializer.imPars = DW_Scheme.nullScheme();
        // Now we can make the fitter object...
        Fitter fitter = ModelFit.getFitter(fm, NoiseModel.GAUSSIAN, FitAlgorithm.LM);
        // ...and figure out how many values per voxel.
        int valsPerVoxel = fitter.getNumValuesPerRun();

        // Also need to indentify the index of the 
        // orientation parameter for the model.
        int dirInd = fitter.getDirectionIndex();
        if(dirInd == 0)
            logger.warning("Direction index is zero: suggests the model contains no orientational information.");

        // For this particular method, we need to create the data
        // source object inside, because only then do we know the
        // number of values per voxel.
        DataSource data = ExternalDataSource.getDataSource(CL_Initializer.inputFile, valsPerVoxel, CL_Initializer.inputDataType);

	int xDataDim = ih.xDataDim();
	int yDataDim = ih.yDataDim();
	int zDataDim = ih.zDataDim();

	double[][][] f = new double[xDataDim][yDataDim][zDataDim];
	
	Vector3D[][][][] vecs = new Vector3D[xDataDim][yDataDim][zDataDim][];

	for (int k = 0; k < zDataDim; k++) { 
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    double[] voxel = data.nextVoxel();
                    
		    Vector3D v1 = null;

                    // By convention, the orientation is
                    // encoded as spherical polar coordinates.
                    PD pd = new PD(voxel[dirInd], voxel[dirInd+1], 1);
                    // Note the y-component needs negating for consistency
                    // with the other input models.
                    vecs[i][j][k] = new Vector3D[] {new Vector3D(pd.getPDX(), pd.getPDY(), pd.getPDZ())};

                    // Also by convention, the third element is the
                    // intra-axonal volume fraction.
                    f[i][j][k] = Math.abs(voxel[2]);

		    if (f[i][j][k] >= 0.0) {
			// it's OK
		    }
		    else {
			f[i][j][k] = 0.0;
		    }
                }
	    }
	}


	if (scalars == null) {
	    scalars = f;
	}
	
	return new RGB_ScalarImage(vecs, ih, scalars, minScalar, maxScalar);

    }


    /**
     * Read an image from tensor eigen system input. 
     * The image has a maximum of two directions per voxel.
     *
     * @param scalars if null, then the FA of the tensor (or mean FA of the first two tensors)
     * is used.
     */
    public static RGB_ScalarImage imageFromTensorEigenSys(DataSource data, ImageHeader ih, 
							  double[][][] scalars, double minScalar,
							  double maxScalar) {
	return imageFromTensorEigenSys(data, ih, scalars, minScalar, maxScalar, 0);
    }



    /**
     * Read an image from tensor eigen system input. 
     * The image has a maximum of two directions per voxel.
     *
     * @param scalars if null, then the FA of the tensor (or mean FA of the first two tensors)
     * is used.
     * @param eigIndex 
     */
    public static RGB_ScalarImage imageFromTensorEigenSys(DataSource data, ImageHeader ih,
							  double[][][] scalars, double minScalar,
							  double maxScalar, int eigIndex) {
	
	int xDataDim = ih.xDataDim();
	int yDataDim = ih.yDataDim();
	int zDataDim = ih.zDataDim();

	Vector3D[][][][] vecs = new Vector3D[xDataDim][yDataDim][zDataDim][];
	
	double[][][] fa = new double[xDataDim][yDataDim][zDataDim];
	
	for (int k = 0; k < zDataDim; k++) { 
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    
		    double[] voxel = data.nextVoxel();
		    
                    Vector3D v1 = null;

		    switch (eigIndex) {

		    case 0 : v1 = new Vector3D(voxel[1], voxel[2], voxel[3]); break;
		    case 1 : v1 = new Vector3D(voxel[5], voxel[6], voxel[7]); break;
		    case 2 : v1 = new Vector3D(voxel[9], voxel[10], voxel[11]); break;
		    default : throw new LoggedException("Invalid eigenvector index " + eigIndex);
		    }
			
			

                    Vector3D v2 = null;

		    double lmean = (voxel[0] + voxel[4] + voxel[8]) / 3.0;

		    if (lmean > 0.0) {

			fa[i][j][k] = 
			    Math.sqrt(1.5 * ((voxel[0] - lmean) * (voxel[0] - lmean) + 
					     (voxel[4] - lmean) * (voxel[4] - lmean) + 
					     (voxel[8] - lmean) * (voxel[8] - lmean)) / 
				      (voxel[0] * voxel[0] + voxel[4] * voxel[4] + voxel[8] * voxel[8]));
		    }
		    
                    // positive scalars only
		    if ( fa[i][j][k] >= 0.0 && fa[i][j][k] < 1.0) {
			// Sensible value
		    }
		    else {
			if (fa[i][j][k] < 0.0) {
			    fa[i][j][k] = 0.0;
			}
			else if (fa[i][j][k] > 1.0) {
			    fa[i][j][k] = 1.0;
			}
		    }

                    if (voxel.length > 12) {
    
                        lmean = (voxel[12] + voxel[16] + voxel[20]) / 3.0;
			
                        if (lmean > 0.0) {
			    switch (eigIndex) {
				
			    case 0 : v2 = new Vector3D(voxel[13], voxel[14], voxel[15]); break;
			    case 1 : v2 = new Vector3D(voxel[17], voxel[18], voxel[19]); break;
			    case 2 : v2 = new Vector3D(voxel[21], voxel[22], voxel[23]); break;
			    default : throw new LoggedException("Invalid eigenvector index " + eigIndex);
			    }
                        
			    double tmp = 
				Math.sqrt(1.5 * ((voxel[12] - lmean) * (voxel[12] - lmean) + 
						 (voxel[16] - lmean) * (voxel[16] - lmean) + 
						 (voxel[20] - lmean) * (voxel[20] - lmean)) / 
					  (voxel[12] * voxel[12] + voxel[16] * voxel[16] + 
					   voxel[20] * voxel[20]));
                            
                            if (tmp < 0.0) {
                                tmp = 0.0;
                            }
                            else if (tmp > 1.0) {
                                tmp = 1.0;
                            }
			    
			    fa[i][j][k] = (fa[i][j][k] + tmp) / 2.0;

			}
			    
                    }

		    if (v2 == null) {
			vecs[i][j][k] = new Vector3D[] {v1};
		    }
		    else {
			vecs[i][j][k] = new Vector3D[] {v1, v2};
		    }
		    
		}
	    }
	}
	
	if (scalars == null) {
	    scalars = fa;
	}
	
	return new RGB_ScalarImage(vecs, ih, scalars, minScalar, maxScalar);
	
    }



    /**
     * Sets the gamma constant <code>g</code> for the rgb component. For each color channel c, the 
     * image value is 255 * c^g, where c is between 0 and 1.
     *
     * If <code>g == 0.0</code>, the image is greyscale.
     *
     */
    public void setRGB_Gamma(double g) {
	
	if (g >= 0.0) {
	    rgbGamma = g;
	}
	else {
	    throw new IllegalArgumentException("Cannot use gamma value " + g);
	}
    }

    
    public double rgbGamma() {
	return rgbGamma;
    }


    /**
     * Sets the gamma constant <code>g</code> for the scalar component. For each scalar value s, the 
     * image value is 255 * s^g, where s is between 0 and 1.
     *
     *
     */
    public void setScalarGamma(double g) {
	
	if (g >= 0.0) {
	    gsGamma = g;
	}
	else {
	    throw new IllegalArgumentException("Cannot use gamma value " + g);
	}
    }

    public double scalarGamma() {
	return gsGamma;
    }



    public double[] getVoxelDims() {
	return new double[] {xVoxelDim, yVoxelDim, zVoxelDim};
    }

    
}

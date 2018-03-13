package imaging;

import java.io.*;
import java.util.logging.*;

import data.*;
import misc.LoggedException;
import numerics.RealMatrix;

/**
 *
 * Superclass for medical image headers, provides an interface for 
 * basic header information and a method to get the underlying data in voxel order.
 * <p>
 * Constructors create header objects, they do not read and write headers. 
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public abstract class ImageHeader {

	 
    private static Logger logger = Logger.getLogger("camino.imaging.ImageHeader");


    /** 
     * Analyze : width<br>
     * MetaIO  : DimSize[0]<br>
     * NIFTI-1 : XDIM<br>
     */
    public abstract int xDataDim();


    /** 
     * Analyze : height<br>
     * Meta    : DimSize[1]<br>
     * NIFTI-1 : YDIM<br>
     */
    public abstract int yDataDim();

   
    /** 
     * Analyze : depth<br>
     * Meta    : DimSize[2]<br>
     * NIFTI-1 : ZDIM<br>
     */
    public abstract int zDataDim();


    /** 
     * @return {xDataDim(), yDataDim(), zDataDim()}.
     */
    public abstract int[] getDataDims();


    /** 
     * Analyze : pixelwidth<br> (may be negative).
     * Meta    : ElementSpacing[0]<br>
     * NIFTI-1 : pixdim[0]<br>
     */
    public abstract double xVoxelDim();


    /** 
     * Analyze : pixelheight<br> (may be negative).
     * Meta    : ElementSpacing[1]<br>
     * NIFTI-1 : pixdim[1]<br>
     */
    public abstract double yVoxelDim();


    /** 
     * Analyze : pixeldepth<br> (may be negative).
     * Meta    : ElementSpacing[2]<br>
     * NIFTI-1 : pixdim[2]<br>
     */
    public abstract double zVoxelDim();


    /** 
     *
     * @return {xVoxelDim(), yVoxelDim(), zVoxelDim()}.
     */
    public abstract double[] getVoxelDims();


    /** 
     * Analyze : nImages<br>
     * Meta    : channels<br>
     * NIFTI-1 : DIM5<br>
     *
     * @return the number of components in the image, or 1 if this field is zero.
     */
    public abstract int components();

    

    /**
     * Gets the transformation matrix that provides a transformation from a = (i, j, k, 1) in 
     * voxel space to b = (x, y, z, 1) in physical space. The origin of voxel space is the centre 
     * of voxel (0,0,0). If your coordinates are in Camino space, you will need to translate 
     * by half a voxel before applying this transformation.
     *
     */
    public abstract RealMatrix getVoxelToPhysicalTransform();



    /**
     * Gets the image data source, scaling applied if applicable. This may
     * involve reading the entire data file from disk if the data is in scanner order.
     */
    public abstract DataSource getImageDataSource();


    /**
     * Gets all the data in the image. Images that store data in voxel order
     * will incur an overhead here while they are converted to scanner order.
     *
     * @return an array of image data, in the order [x][y][z][volume].
     */
    public abstract double[][][][] readVolumeData();


    /**
     * Convenience method for reading scalar 3D images (or the first volume of a 4D image.
     * 
     * @return an array of image data, in the order [x][y][z].
     */
    public double[][][] readSingleVolumeData() {
        
        // it may be tempting to add an integer parameter here and offer the capability
        // to get any volume, but this gets slow for large, compressed images.
        // Better for code that wants to iterate over all volumes to just do the 4D I/O once

        double[][][][] fourD = readVolumeData();

        int xDim = xDataDim();
        int yDim = yDataDim();
        int zDim = zDataDim();
        
        double[][][] threeD = new double[xDim][yDim][zDim];
        
        for (int i = 0; i < xDim; i++) {
            for (int j = 0; j < yDim; j++) {
                for (int k = 0; k < zDim; k++) { 
                    threeD[i][j][k] = fourD[i][j][k][0];
                }
            }
        }
        
        return threeD;
    }
    

    /**
     * Get the header file name associated with the header object.
     *
     */
    public abstract String getHeaderFilename();

    /**
     * Get the data file name associated with the header object. May be the same as the header
     * file name.
     *
     */
    public abstract String getDataFilename();


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
     *
     * @return the image header associated with the output data set
     */
    public abstract ImageHeader writeScalarImage(double[][][] data, String fileRoot);


    /**
     * Use the current header to write a vector image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, The new header is a copy of this one, 
     * with fields altered where necessary for writing Vector data. 
     *
     * @param data will be checked against header data dimensions. The data array should be in the 
     * same voxel space as this image, to ensure a correct definition of physical space.
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
N     * based upon the header settings. 
     *
     * @return the image header associated with the output data set
     */
    public abstract ImageHeader writeVectorImage(double[][][][] data, String fileRoot);


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
    public abstract ImageHeader writeRGB_Image(int[][][] red, int[][][] green, int[][][] blue, String fileRoot);


    /**
     * 
     * Use the current header to write a tensor image to a file, and return the header 
     * associated with that file. The new header is a copy of this one, with fields altered where necessary 
     * for writing tensor data. 
     *
     * @param data will be checked against header data dimensions. The data array should be in the 
     * same voxel space as this image, to ensure a correct definition of physical space. data[i][j][k] 
     * contains six components in upper-triangular order. The data will be written to disk in upper or 
     * lower triangular order, as appropriate for the image type.
     *
     * @param fileRoot the path to the output file name, minus any the extension. The extension may vary
     * based upon the header settings.
     *
     * @return the image header associated with the output data set
     */
    public abstract ImageHeader writeTensorImage(double[][][][] data, String fileRoot);

  

    /**
     * Set the data type of the image.
     *
     * @param a Camino data type string.
     *
     */
    public abstract void setDataType(String datatype);



    /**
     * Sets gzip compression. This may alter the file name(s) of the image.
     *
     */
    public abstract void setGzip(boolean gz);



    /**
     * Gets the transformation matrix that is the inverse of the voxel to physical space transformation. 
     * <p>
     * Assumes the transform is composed of a 3x3 invertible transform A, and a translation vector B. 
     * In the inverse, A' = A^{i-1} and B' = -A^{i-1}B.
     * <p>
     * Note that this converts to NIfTI voxel space, the minimum voxel coordinate is (-0.5, -0.5, -0.5)
     * because the coordinate system begins at the centre of voxel (0,0,0).
     *
     */
    public RealMatrix getPhysicalToVoxelTransform() {
        
        RealMatrix voxelToPhysical = getVoxelToPhysicalTransform();

        RealMatrix physicalToVoxel = new RealMatrix(4,4);
        
        // voxel to physical trans is 
        // [A,       B]
        // [0, 0, 0, 1]
        //
        // where A is 3x3 transform, and B is 3x1 translation
        //
        RealMatrix A = new RealMatrix(3,3);

        RealMatrix B = new RealMatrix(3,1);

        for (int i = 0; i < 3; i++) {
            
            B.entries[i][0] = voxelToPhysical.entries[i][3];

            for (int j = 0; j < 3; j++) {
                A.entries[i][j] = voxelToPhysical.entries[i][j];
            }
        }

        RealMatrix Ainv = A.inverse();

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                physicalToVoxel.entries[i][j] = Ainv.entries[i][j];
            }
        }
        
        Ainv.scale(-1.0);

        RealMatrix trans = Ainv.product(B);

        for (int i = 0; i < 3; i++) {
            physicalToVoxel.entries[i][3] = trans.entries[i][0];
        }

        physicalToVoxel.entries[3][3] = 1.0;


        return physicalToVoxel;

    }



    /**
     * Compares the image space of another header to this one. Checks that images have the same dimensions,
     * but doesn't care about the physical space. Checks first three dimensions only.
     * 
     *
     * @return true if the images have the same dimensions, zero otherwise.
     */
    public final boolean sameDimensions(ImageHeader header) {
        return sameDimensions(header.getDataDims(), header.getVoxelDims());
    }



   /**
    * For legacy compatibility, check data and voxel dims. First three dimensions only.
    * 
    *
    * @return true if the dimensions match, false otherwise.
    */
    public final boolean sameDimensions(int[] dataDims, double[] voxelDims) {


        double epsilon = 1E-5;

        if (xDataDim() != dataDims[0] || yDataDim() != dataDims[1] || zDataDim() != dataDims[2]) {
            return false;
        }
        
        double[] myVoxDims = getVoxelDims();

        for (int i = 0; i < 3; i++) {
            
            if (Math.abs(myVoxDims[i] - voxelDims[i]) > epsilon) {
                return false;
            }
        }
        
        return true;
    }


    /**
     * Compares the image space of another header to this one. Images are considered to be 
     * in the same space if they have the same dimensions, and the voxel to physical transformation 
     * is the same (within some numerical epsilon). If images are in the same space, then it means
     * that the voxel to physical transformation is identical between the two. 
     * <p>
     * Images might be aligned anatomically but still not be in the same space according to this
     * method. This method just tells you whether the images are aligned in both the voxel and physical space.
     *
     * @return true if the images are in the same space, false otherwise.
     */
    public final boolean sameSpace(ImageHeader header) {

        // Can't be in the same space if the image is a different size
        if (!sameDimensions(header)) {
            return false;
        }


        RealMatrix vox2Phys = getVoxelToPhysicalTransform();

        RealMatrix headerVox2Phys = header.getVoxelToPhysicalTransform();
	
        // don't complain if differences are less than this
        double epsilon = 1E-5;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if ( Math.abs(vox2Phys.entries[i][j] - headerVox2Phys.entries[i][j]) > epsilon ) {
                    logger.info("Different voxel to physical space transforms between " + getHeaderFilename() + 
                                   " and " + header.getHeaderFilename());

                    return false;
                }
            }
        }

	return true;
        
    }




    /**
     * Compares the image space of another header to this one. Images are considered to be 
     * in the same space if they have the same dimensions, and the voxel to physical transformation 
     * is the same (within some numerical epsilon). It's possible that this could be
     * misleading because the data might be stored differently on disk, but there's no way to tell this
     * from the header.
     * 
     *
     * @return true if the images are in the same space, false otherwise. Returns false if the image
     * cannot be read.
     */
    public final boolean sameSpace(String file) {
        
        try {
            return sameSpace(ImageHeader.readHeader(file));
        }
        catch (IOException e) {
            logger.warning("Cannot read image " + file);
            return false;
        }

    }



    /**
     * Returns the appropriate header based on the given file name.
     *
     * @return the header object for the image.
     *
     */
    public static ImageHeader readHeader(String hdrFile) throws IOException {
	
	if (hdrFile == null) {
	    throw new LoggedException("File name required to read image header " + 
				      "(format determined by extension)");
	}


	if (AnalyzeHeader.getImageRoot(hdrFile) != null) {
            
            // Check if this is a NIfTI hdr / img pair
            if (Nifti1Dataset.hdrIsNifti(hdrFile)) {
                return Nifti1Dataset.readHeader(hdrFile);
            }
            
	    String root = AnalyzeHeader.getImageRoot(hdrFile);

	    File f = new File(root + ".hdr");

	    long length = f.length();

            return AnalyzeHeader.readHeader(root + ".hdr");
            
	}
	if (hdrFile.endsWith(".nii") || hdrFile.endsWith(".nii.gz")) {
	    // NIFTI
	    return Nifti1Dataset.readHeader(hdrFile);
	}
	if (hdrFile.endsWith(".mhd") || hdrFile.endsWith(".mha")) {
	    // meta
	    return MetaImageHeader.readHeader(hdrFile);
	}


	throw new LoggedException("Can't find image for input file name " + hdrFile);
    }


    /**
     * @param fileName the full path to a file, including the extension (works without extension for Analyze, backwards compatibility thing).
     *
     * @return true if the file exists and ends with a recognized extension: .hdr, .nii[.gz], .mha, .mhd
     *
     */
    public static boolean imageExists(String fileName) {
        
        if (fileName == null) {
            return false;
        }

	// first check if the file ends in a recognized extension, if so check for
	// that image type only
	if (fileName.endsWith(".nii") || fileName.endsWith(".nii.gz") || fileName.endsWith(".mha") ||
	    fileName.endsWith(".mhd") || fileName.endsWith(".hdr")) {
	    
	    File f = new File(fileName);
	    return f.exists();
	    
	}

        return false;

        
    }


    /**
     * Strips away any recognized extensions. If the file doesn't end with a recognized extension, this method just 
     * returns the original string.
     *
     */
    public static final String getFileRoot(String fileName) {
        
        String root = fileName;

        if (fileName.endsWith(".nii") || fileName.endsWith(".mha") || fileName.endsWith(".mhd") || fileName.endsWith(".hdr")) {

            root = fileName.substring(0, fileName.length() - 4);

        }
        if (fileName.endsWith(".nii.gz")) {
            root = fileName.substring(0, fileName.length() - 7);
        }

        return root;

    }

    /**
     * Returns the extension of the image, if it's a recognized format. Otherwise an empty string. In all cases,
     * getFileRoot() + getFileExtension() == the original fileName.
     *
     */
    public static final String getFileExtension(String fileName) {
        
        String root = fileName;

        if (fileName.endsWith(".nii") || fileName.endsWith(".mha") || fileName.endsWith(".mhd") || fileName.endsWith(".hdr")) {

            return fileName.substring(fileName.length() - 4, fileName.length());

        }
        if (fileName.endsWith(".nii.gz")) {
            return fileName.substring(fileName.length() - 7, fileName.length());
        }

        return "";

    }



    /**
     * @return true if a symmetric matrix in this format will be lower-triangular.
     */
    public abstract boolean lowerTriangularSymmMatrix();

}

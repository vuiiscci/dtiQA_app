package tractography;

import data.*;
import misc.*;
import numerics.*;


/**
 * Tractography image for spherical function peak data (output of sfpeaks). Aside from offering methods to 
 * construct an image from such data files, this class is identical to PD_TractographyImage.
 * 
 *
 *
 * @version $Id: SF_TractographyImage.java 1147 2012-11-23 23:58:59Z ucacpco $
 * @author  Philip Cook
 * @see apps.SphFuncPD_Stats
 */
public class SF_TractographyImage extends PD_TractographyImage {


    
    /**
     * Constructs an image from a data source. The source should provide data from <code>SphFuncPD_Stats</code>.
     *
     * @param pdSource the input data, in the output format of <code>SphFuncPD_Stats</code>.
     * @param maxPDs the maximum number of PDs in a voxel.
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     *
     */
    public SF_TractographyImage(DataSource pdSource, int maxPDs, int[] dataDims, 
                                double[] voxelDims) {
	
	super(dataDims, voxelDims);


        for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {

                    double[] voxel = pdSource.nextVoxel();

                    if (voxel[0] < 0.0) {
                        numPDs[i][j][k] = 0;
                    }
                    else {
                        numPDs[i][j][k] = (int)voxel[2];

                        if (numPDs[i][j][k] > maxPDs) {
                            numPDs[i][j][k] = maxPDs;
                        } 
                    }
                    
                    if (numPDs[i][j][k] < 0) {
                        // can't check for number of PDs because SphFuncPD_Stats can exceed
                        // maxPDs (but the extra ones are not in the file)
                        throw new LoggedException("Invalid number of components in input data. " + 
                                                  "Check -inputmodel and -numpds options");
                    }
                    // might be zero length
                    vectors[i][j][k] = new Vector3D[numPDs[i][j][k]];

                    for (int p = 0; p < numPDs[i][j][k]; p++) {
                        int start = 6 + p * 8;
                        vectors[i][j][k][p] = new Vector3D(voxel[start],voxel[start+1],voxel[start+2]);
                    }
                }
            }
        }

        // default mask just segments background
	computeIsotropicMask();
    }


    /**
     * @param vectors vector data.
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public SF_TractographyImage(Vector3D[][][][] vectors, double[] voxelDims) {
	super(vectors, voxelDims);
    }


    /**
     * Copy constructor
     * 
     */
    protected SF_TractographyImage(SF_TractographyImage im) {
	super(im.vectors, im.getVoxelDims());
    }


  
    /**
     * Gets an image from the data file. 
     * If <code>anisMapFile</code> is not <code>null</code>, it is read and used 
     * for isotropic masking.
     * 
     *
     * @param inputFile the data file.
     * @param dataType the data type of the data file and <code>anisMapFile</code>.
     * @param maxPDs the maximum number of PDs in a voxel.
     * @param anisMap the anisotropy map, which is used to create the tract mask.
     * May be <code>null</code> if not required.
     * @param anisThresh threshold for the anisotropy in the computation of the tract mask.
     * 
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public static final SF_TractographyImage getTractographyImage(String inputFile, 
                                                                     String dataType, int maxPDs, 
                                                                     double[][][] anisMap, 
                                                                     double anisThresh, 
                                                                     int[] dataDims, 
                                                                     double[] voxelDims) {

        // components per voxel in the data
        int numComponents = 6 + 8 * maxPDs;


	DataSource dataSource = ExternalDataSource.getDataSource(inputFile, numComponents, dataType);

        SF_TractographyImage image = new SF_TractographyImage(dataSource, maxPDs, dataDims, voxelDims);

        if (anisMap != null) {
            image.computeIsotropicMask(anisMap, anisThresh);
        }

        return image;	    

    }


    
}

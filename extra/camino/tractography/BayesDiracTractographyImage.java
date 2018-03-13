package tractography;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;

import java.util.Random;


/**
 * Bayesian tractography image using Dirac priors for "nuisance parameters", 
 * see Friman et al, TMI 25:965-978 (2006). 
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class BayesDiracTractographyImage extends ProbabilisticTractographyImage {

   
    protected final float[][][][] data;
    protected final DW_Scheme ip;

    protected final BayesDiracRandomizer randomizer;


    /**
     * Matlab constructor.
     *
     * @param data DWI data
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param nPDs PDs per voxel (currently must be 1 or 0).
     * @param scheme imaging scheme for data
     * @param model Choice of Bayesian model
     * @param pointSetInd index of ISCodes point set to use. Larger point sets slow computation, but 
     * have finer angular resolution. Index 1 (1922 directions) should be sufficient, if memory is a 
     * concern you might also do as well with 0 (1082) but this is untested.
     *
     * @param ran a random number generator. tools.MTRandom preferred.
     *
     * @see sphfunc.ISCodes
     * @see tractography.BayesDiracRandomizer
     */
    public BayesDiracTractographyImage(float[][][][] data, double[] voxelDims, int[][][] nPDs, 
					  DW_Scheme scheme, BayesDataModel model, int pointSetInd, Random ran ) { 
	
	super(new int[] {data.length, data[0].length, data[0][0].length}, voxelDims, ran);


	ip = scheme;
	this.data = data;

	switch (model) {
	    
	case CYL_SYMM_DT: 
	    randomizer = new BayesDiracCylSymmTensorRandomizer(data, ip, pointSetInd, ran);
	    break;

	case BALL_STICK: 
	    randomizer = new BayesDiracBallStickRandomizer(data, ip, pointSetInd, ran);
	    break;

	default:
	    throw new LoggedException("Unsupported data model " + model);
	}

        for (int j = 0; j < yDataDim; j++) {
            for (int i = 0; i < xDataDim; i++) {
                System.arraycopy(nPDs[i][j], 0, numPDs[i][j], 0, zDataDim);
            }
        }

	
	computeIsotropicMask();
    }



    protected void setVectors(int i, int j, int k) {
        vectors[i][j][k] = randomizer.getRandomizedPDs(i,j,k);
    }
    
    
    protected void setVectors(int i, int j, int k, Vector3D fibreOrientation) {
        vectors[i][j][k] = randomizer.getRandomizedPDs(i,j,k, fibreOrientation);
    }


    public void setExternalPriors(PICoTractographyImage priorImage) {

	// check that image dimensions are the same

	int[] dataDims = priorImage.getDataDims();
	double[] voxelDims = priorImage.getVoxelDims();

	if (dataDims[0] != xDataDim || dataDims[1] != yDataDim || dataDims[2] != zDataDim) {
	    throw new LoggedException("External prior image has different data dimensions");
	}

	if (voxelDims[0] != xVoxelDim || voxelDims[1] != yVoxelDim || voxelDims[2] != zVoxelDim) {
	    throw new LoggedException("External prior image has different voxel dimensions");
	}

	randomizer.setExternalPriors(priorImage.randomizer);
    }



    /**
     * Gets an image from the data file. 
     * If <code>anisMapFile</code> is not <code>null</code>, it is read and used 
     * for isotropic masking.
     * 
     *
     * @param inputFile the data file.
     * @param dataType the data type of the data file.
     * @param numPDs array containing number of PDs in each voxel.
     * @param anisMap the anisotropy map, which is used to create the tract mask.
     * May be <code>null</code> if not required.
     * @param anisThresh threshold for the anisotropy in the computation of the tract mask.
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public static final BayesDiracTractographyImage getTractographyImage(String inputFile,
									 String dataType,
									 DW_Scheme scheme,
									 BayesDataModel dataModel,
									 int[][][] numPDs,
									 double[][][] anisMap,
									 double anisThresh,
									 int[] dataDims,
									 double[] voxelDims,
									 int pointSetInd,
									 Random ran) {
	

	int numMeas = scheme.numMeasurements();
	
        DataSource dataSource = ExternalDataSource.getDataSource(inputFile, numMeas, dataType);

	float[][][][] data = new float[dataDims[0]][dataDims[1]][dataDims[2]][numMeas];

	
	for (int k = 0; k < dataDims[2]; k++) {
	    for (int j = 0; j < dataDims[1]; j++) {
		for (int i = 0; i < dataDims[0]; i++) {

		    double[] voxel = dataSource.nextVoxel();
		    
		    for (int n = 0; n < numMeas; n++) {
			data[i][j][k][n] = (float)voxel[n];
		    }
		}
	    }
	}
				       
	BayesDiracTractographyImage image = new BayesDiracTractographyImage(data, voxelDims, numPDs, 
									    scheme, dataModel, pointSetInd, ran);
	
        if (anisMap != null) {
            image.computeIsotropicMask(anisMap, anisThresh);
        }
	
	return image;
    }

    
    /**
     * Sets the curvature prior concentration using a Watson distribution.
     *
     */
    public void setCurvePriorKappa(double k) {
	randomizer.setCurvatureK(k);
    }

    /**
     * Sets the curvature prior concentration to |v_i \dot v_{i-1}|^\gamma given the current
     * orientation v_i and the previous one v_{i-1}.
     *
     */
    public void setCurvePriorGamma(double g) {
	randomizer.setCurvatureG(g);
    }


}

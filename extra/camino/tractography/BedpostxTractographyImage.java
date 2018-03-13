package tractography;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.CL_Initializer;

import java.io.*;
import java.util.Random;


/**
 * Does bootstrap sampling of data produced by FSL's bedpostx program. See FSL website for bedpost details
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class BedpostxTractographyImage extends ProbabilisticTractographyImage {

   
    private final float[][][][][] theta;
    private final float[][][][][] phi;

    private final int numSamples;
    
    /**
     * 
     * @param th theta samples in order [i][j][k][p][n], for voxel (i,j,k), compartment p ranges from 0 to numPDs[i][j][k] - 1, 
     * and n ranges from 0 to numSamples - 1.
     *
     * @param ph phi samples in the same order as theta samples.
     *
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     *
     * @param nPDs number of PDs per voxel.
     *
     * @param nSamples number of samples per voxel, usually 50.
     *
     * @param ran a random number generator. tools.MTRandom preferred.
     *
     * 
     */
    public BedpostxTractographyImage(float[][][][][] th, float[][][][][] ph, double[] voxelDims,  int[][][] nPDs, int nSamples,
				     Random ran ) { 
	
	super(new int[] {th.length, th[0].length, th[0][0].length}, voxelDims, ran);

	theta = th;
	phi = ph;

	numSamples = nSamples;

	for (int j = 0; j < yDataDim; j++) {
            for (int i = 0; i < xDataDim; i++) {
                System.arraycopy(nPDs[i][j], 0, numPDs[i][j], 0, zDataDim);
            }
        }
	
	computeIsotropicMask();
	
   }


    protected void setVectors(int i, int j, int k) {

	int pdsVoxel = numPDs[i][j][k];
	
	vectors[i][j][k] = new Vector3D[pdsVoxel];

	int sample = ran.nextInt(numSamples);

	for(int n = 0; n < pdsVoxel; n++) {
	    vectors[i][j][k][n] = Vector3D.vectorFromSPC(1.0, theta[i][j][k][n][sample], phi[i][j][k][n][sample]);
	}
	
    }
    

    /**
     * Gets an image from the specified bedpostx data. 
     * 
     *
     * @param inputRoot directory containing bedpost output.
     * @param probabilistic if true, load the MCMC samples of the principal directions and return a probabilistic image.
     * Else, load the mean dyads and return a deterministic image.
     *
     * @param minCompartmentSize the minimum mean volume fraction for a principal direction. The number of PDs in a voxel is determined
     * by this, if the mean is below the minimum then the PD is ignored. 
     *
     * @param anisMap the anisotropy map, which is used to create the tracking mask. 
     * @param anisThresh anisotropy threshold corresponding to the anisMap. If anisMap is null, the threshold is applied based on the 
     * mean anisotropic compartment fractions. The total anisotropy in a voxel is the sum of all fractions that are "sticks".
     * May be <code>null</code> if not required.
     *
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param ran a random number generator for probabilistic sampling. Can be <code>null</code> if not required.
     * 
     */
    public static final TractographyImage getTractographyImage(String bedpostxDir,
							       boolean probabilistic,
							       double minCompartmentFraction,
							       double[][][] anisMap,
							       double anisThresh,
							       int[] dataDims,
							       double[] voxelDims,
							       Random ran) {
	

	String inputRoot = getBedpostxInputRoot(bedpostxDir);

	String fileExt = getBedpostxImageExtension(inputRoot);
	
	int maxPDs = 0;
	
	String[] dyadFiles = new String[100];
	
	String nextDyad = inputRoot + "dyads" + (maxPDs + 1) + fileExt;
	
	File f = new File(nextDyad);
	
	if (! f.exists()) {
	    throw new LoggedException("Can't find any bedpost dyads matching " + inputRoot);
	}

	// quick sanity check that physical space is set correctly
	try {
	    ImageHeader ih = ImageHeader.readHeader(nextDyad);
	    
	    if (!CL_Initializer.headerTemplate.sameSpace(ih)) {
		throw new LoggedException("Bedpost data does not agree with definition of input physical space");
	    }
	}
	catch (IOException e) {
	    throw new LoggedException(e);
	}
		
	
	while (f.exists()) {
	    dyadFiles[maxPDs] = nextDyad;
	    
	    maxPDs++;
	    
	    nextDyad = inputRoot + "dyads" + (maxPDs + 1) + fileExt;

	    f = new File(nextDyad);
	}
	
	String[] tmpDyads = new String[maxPDs];
	
	System.arraycopy(dyadFiles, 0, tmpDyads, 0, maxPDs);
	
	dyadFiles = tmpDyads;

	DataSource[] fmSources = new DataSource[maxPDs];
	
	try {
	    for (int i = 0; i < maxPDs; i++) {
		// file names count from 1
		ImageHeader header = ImageHeader.readHeader(inputRoot + "mean_f" + (i+1) + "samples" + fileExt);
		
		fmSources[i] = header.getImageDataSource();
	    }   
	}
	catch(IOException e) {
	    throw new LoggedException(e);
	}

	
	double[][][][] fm = new double[dataDims[0]][dataDims[1]][dataDims[2]][maxPDs];

	// use as an anisotropy map if one is not provided
	double[][][] fSum = new double[dataDims[0]][dataDims[1]][dataDims[2]];
	
	for (int k = 0; k < dataDims[2]; k++) {
	    for (int j = 0; j < dataDims[1]; j++) {
		for (int i = 0; i < dataDims[0]; i++) {
		    for (int p = 0; p < maxPDs; p++) {

			fm[i][j][k][p] = fmSources[p].nextVoxel()[0];

			fSum[i][j][k] += fm[i][j][k][p];

		    }
		}
	    }
	}

	if (anisMap == null) {
	    anisMap = fSum;
	}


	if ( !probabilistic ) {   

	    int numComponents = 3;
	    
	    // Read vectors into an array and then construct image
	    
	    Vector3D[][][][] vectors = new Vector3D[dataDims[0]][dataDims[1]][dataDims[2]][];
	    
	    DataSource[] dyadSources = new DataSource[maxPDs];
	    
	    for (int i = 0; i < maxPDs; i++) {
		dyadSources[i] = ExternalDataSource.getDataSource(dyadFiles[i], numComponents, CL_Initializer.inputDataType);
	    }
	    
	    for (int k = 0; k < dataDims[2]; k++) {
		for (int j = 0; j < dataDims[1]; j++) {
		    for (int i = 0; i < dataDims[0]; i++) {
			
			// determine how many PDs exist for this voxel
			Vector3D[] tmpVec = new Vector3D[maxPDs];
			
			int numPDs = 0;
		
			for (int n = 0; n < maxPDs; n++) {
			    double[] components = dyadSources[n].nextVoxel();
			    
			    if (fm[i][j][k][n] >= minCompartmentFraction) {
				tmpVec[numPDs] = new Vector3D(components);
				numPDs++;
			    }
			}
			
			vectors[i][j][k] = new Vector3D[numPDs];
			
			System.arraycopy(tmpVec, 0, vectors[i][j][k], 0, numPDs);
			
		    }
		}
	    }

	    
	    PD_TractographyImage image = new PD_TractographyImage(vectors, voxelDims);
	
	    image.computeIsotropicMask(anisMap, anisThresh);
	    
	    return image;
	}

	
	// Probabilistic image; read theta and phi samples
	
	DataSource[] thSources = new DataSource[maxPDs];
	DataSource[] phSources = new DataSource[maxPDs];

	// default samples is 50 but check from file 
	int numSamples = 0;
	
	try {
	    ImageHeader header = ImageHeader.readHeader(inputRoot + "merged_th1samples" + fileExt);

	    numSamples = header.components();

	    for (int i = 0; i < maxPDs; i++) {
		// file names count from 1
		header = ImageHeader.readHeader(inputRoot + "merged_th" + (i+1) + "samples" + fileExt);
		
		thSources[i] = header.getImageDataSource();

		header = ImageHeader.readHeader(inputRoot + "merged_ph" + (i+1) + "samples" + fileExt);

		phSources[i] = header.getImageDataSource();
	    }  
	}
	catch(IOException e) {
	    throw new LoggedException(e);
	}

	float[][][][][] theta = new float[dataDims[0]][dataDims[1]][dataDims[2]][maxPDs][];
	float[][][][][] phi = new float [dataDims[0]][dataDims[1]][dataDims[2]][maxPDs][];

	// Cumbersome procedure for copying because we store things as float but read as double, thus no System.arraycopy

	// We set numPDs for each voxel as the number of compartments in which the mean f is above the minimum 
	int[][][] numPDs = new int[dataDims[0]][dataDims[1]][dataDims[2]];
	
	for (int k = 0; k < dataDims[2]; k++) {
	    for (int j = 0; j < dataDims[1]; j++) {
		for (int i = 0; i < dataDims[0]; i++) {

		    int numPDsVoxel = 0;
		    
		    for (int p = 0; p < maxPDs; p++) {
			
			double[] thSamples = thSources[p].nextVoxel();

			double[] phSamples = phSources[p].nextVoxel();

			if (fm[i][j][k][numPDsVoxel] >= minCompartmentFraction) {
			    
			    theta[i][j][k][numPDsVoxel] = new float[numSamples];
			    phi[i][j][k][numPDsVoxel] = new float[numSamples];
		    
			    for (int n = 0; n < numSamples; n++) {
				theta[i][j][k][numPDsVoxel][n] = (float)thSamples[n];
				phi[i][j][k][numPDsVoxel][n] = (float)phSamples[n];
			    }

			    numPDsVoxel++; 
			}
			
		    }
		    
		    numPDs[i][j][k] = numPDsVoxel;
		}
	    }
	}
	
	TractographyImage image = new BedpostxTractographyImage(theta, phi, voxelDims, numPDs, numSamples, ran);			       

	image.computeIsotropicMask(anisMap, anisThresh);

	return image;
    }



     /**
     * Default FSL behavior is to have files inside the directory with no prefix. Optionally, a prefix may be specified
     * so there might be dyads1.nii.gz or prefix_dyads1.nii.gz for example. This method figures out the complete input root
     * - which by default is just the directory.
     *
     * @return the input root to the bedpost files if they exist, or throws an exception if the dyads1 file cannot be found.
     */
    public static String getBedpostxInputRoot(String bedpostxDir) {

	if (bedpostxDir.endsWith("/")) {
	    bedpostxDir = bedpostxDir.substring(0, bedpostxDir.length());
	}

	File dir = new File(bedpostxDir);
	
	File[] files = dir.listFiles();

	for (int i = 0; i < files.length; i++) {
	    
	    if (files[i].getName().matches(".*dyads1.*")) {
		String[] tokens = files[i].getName().split("dyads1");
		   
		return bedpostxDir + "/" + tokens[0];
	    }

	}

	// no dyads
	throw new LoggedException("Can't locate dyads1 in bedpostxdir " + bedpostxDir);

    }
    

    /**
     * Looks for .nii.gz, .nii, .hdr
     */
    public static String getBedpostxImageExtension(String bedpostxInputRoot) {

	String[] extensions = {".nii.gz", ".nii", ".hdr"};

	for (int i = 0; i < extensions.length; i++) {
	    
	    File f = new File(bedpostxInputRoot + "dyads1" + extensions[i]);
	    
	    if (f.exists()) {
		return extensions[i];
	    }
	}

	throw new LoggedException("Can't find any supported image format matching " + bedpostxInputRoot);

    }

}

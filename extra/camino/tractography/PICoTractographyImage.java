package tractography;

import data.*;
import misc.LoggedException;
import numerics.*;

import java.util.Random;

/**
 * Tractography image with all three eigenvectors of spherical PDF(s) specified in each voxel. 
 * These can be from the diffusion tensor or from other models. 
 *
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class PICoTractographyImage extends ProbabilisticTractographyImage implements TractographyImage {

  
    protected final double[][][][] pdfParams;

    protected final PICoRandomizer randomizer;

    protected final PICoPDF pdf;

    /**
     * Array of eigenvectors for each PDF. 
     * 
     */
    protected Vector3D[][][][] eigenvectors;

 
    /**
     * Constructs an image from the data sources. The data source should be the output of PICoPDFs.
     *
     * @param dataSource source for the data, in the format of the picopdfs program.
     * @param maxPDs the maximum number of PDs in a voxel.
     * @param pdf PICo PDF type.
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * @param ran a source of random numbers.
     *
     */
    public PICoTractographyImage(DataSource dataSource, int maxPDs, PICoPDF pdf,
                                 int[] dataDims, double[] voxelDims, Random ran) {
	
        super(dataDims, voxelDims, ran);

	this.pdf = pdf;
      
        pdfParams = new double[xDataDim][yDataDim][zDataDim][];

	int paramsPerPD = pdf.numParams;
        
        eigenvectors = new Vector3D[xDataDim][yDataDim][zDataDim][];

   	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
                    double[] voxel = dataSource.nextVoxel();

                    numPDs[i][j][k] = (int)voxel[0];

                    if (numPDs[i][j][k] < 0 || numPDs[i][j][k] > maxPDs) {
                        throw new LoggedException("Invalid number of components, " + numPDs[i][j][k] + 
                                                  ", in input data. " + 
                                                  "Check -inputmodel, PICo PDF and -numpds options.");
                    }

                    eigenvectors[i][j][k] = new Vector3D[3 * numPDs[i][j][k]];

                    pdfParams[i][j][k] = new double[numPDs[i][j][k]*paramsPerPD];

                    for (int p = 0; p < numPDs[i][j][k]; p++) {
                        int start = 2 + p * (10 + paramsPerPD);

                        eigenvectors[i][j][k][3*p] = 
                            new Vector3D(voxel[start],voxel[start+1],voxel[start+2]);

                        eigenvectors[i][j][k][3*p+1] = 
                            new Vector3D(voxel[start+3],voxel[start+4],voxel[start+5]);

                        eigenvectors[i][j][k][3*p+2] = 
                            new Vector3D(voxel[start+6],voxel[start+7],voxel[start+8]);
                        
                        for (int a = 0; a < paramsPerPD; a++) {
                            pdfParams[i][j][k][paramsPerPD*p + a] = voxel[start + 9 + a];
                        }
                    }
                    
                    
                }
            }
        }

        // default mask just segments background
        computeIsotropicMask();

	randomizer = pdf.getRandomizer(this, ran);

    }


    /**
     *
     * @param eVectors vector data in the order {e1, e2, e3}.
     *
     * @param pdfParams the array of PICo parameters.
     *
     * @param pdf the PICoPDF corresponding to the pdfParams.
     *
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public PICoTractographyImage(Vector3D[][][][] eVectors, double[][][][] pdfParams, PICoPDF pdf,
                                 double[] voxelDims, Random ran) {

        // note do not pass eVectors to superclass constructor, because it will try to calculate
        // numPDs from that
	super(new int[] {eVectors.length, eVectors[0].length, eVectors[0][0].length}, voxelDims, ran);

        this.pdfParams = pdfParams;

	this.pdf = pdf;

        for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
                    numPDs[i][j][k] = eVectors[i][j][k].length / 3;
                }
            }
        }
        
        eigenvectors = eVectors;
        
	randomizer = pdf.getRandomizer(this, ran);

        computeIsotropicMask();
    }


    /**
     * Copy constructor for testing. Arrays in this object will not be copied.
     *
     */
    protected PICoTractographyImage(PICoTractographyImage im, Random r) {
        this(im.eigenvectors, im.pdfParams, im.pdf, im.getVoxelDims(), r);
    }


 
    /**
     * Gets the eigenvectors of the PDF for the voxel in the order  <code>{e1, e2, e3, e1...}</code>. 
     *
     */
    protected Vector3D[] getEigenvectors(int i, int j, int k) {
	return eigenvectors[i][j][k];
    }


    /**
     * @return the PICo PDF params for this voxel.
     */
    protected double[] getPICoPDFParams(int i, int j, int k) {
        return pdfParams[i][j][k];
    }


    protected void setVectors(int i, int j, int k) {
        vectors[i][j][k] = randomizer.getRandomizedPDs(i,j,k);
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
     * @param pdf PICo PDF type.
     * @param anisMap the anisotropy map, which is used to create the tract mask.
     * May be <code>null</code> if not required.
     * @param anisThresh threshold for the anisotropy in the computation of the tract mask.
     * 
     * @param dataDims array of data dimensions {xDataDim, yDataDim, zDataDim}.
     * @param voxelDims array of voxel dimensions (in mm) {xVoxelDim, yVoxelDim, zVoxelDim}.
     * 
     */
    public static final PICoTractographyImage getTractographyImage(String inputFile, 
                                                                   String dataType, int maxPDs, 
                                                                   PICoPDF pdf,
                                                                   double[][][] anisMap, 
                                                                   double anisThresh, 
                                                                   int[] dataDims, 
                                                                   double[] voxelDims,
                                                                   Random ran) {

	int paramsPerPD = pdf.numParams;
        
        // components per voxel in the data
        int numComponents = 1 + maxPDs * (10 + paramsPerPD);

        DataSource dataSource = ExternalDataSource.getDataSource(inputFile, numComponents, dataType);
	       
        PICoTractographyImage image = new PICoTractographyImage(dataSource, maxPDs, pdf, 
                                                                dataDims, voxelDims, ran);

        if (anisMap != null) {
            image.computeIsotropicMask(anisMap, anisThresh);
        }

        return image;	    

    }


   
}

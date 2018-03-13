package tractography;

import misc.*;
import numerics.*;

import java.util.*;
import java.util.logging.*;

/**
 * Labels a seed ROI with the target that has the maximum connection probability to each seed.
 * The output is a labeled copy of the seed ROI image. It is assumed that there is one seed
 * point per voxel of this image. If there are multiple seeds per voxel, the label is determined
 * by the last seed to be processed for that voxel.
 * <p>
 * This class is required by procstreamlines but is deprecated, users should call cbsmat which uses
 * apps.ConnectivityMatrix. This class doesn't work with labels that are not numbered 1,2,...,N.
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 * 
 */
public class ConnectivitySegmentedImage {

    /**
     *
     * Logging object
     *
     */
    private static Logger logger = Logger.getLogger("camino.tractography.ConnectivitySegmentedImage");

    // voxel dims of seed space
    private final double xVoxelDim;
    private final double yVoxelDim;
    private final double zVoxelDim;


    // dimensions of seed space
    private final int xDataDim;
    private final int yDataDim;
    private final int zDataDim;

    private RegionOfInterest roi = null;

    private Voxel[] seedVoxels = null;

    private int numSeeds = 0;

    private int[][][] targets;

    private boolean directional = false;

    // if directional, vector tells us forwards / backwards (or left / right, or whatever)
    private Vector3D forwards = null;

    private boolean countFirstEntry = true;

    // lazy initialized
    private double[][][] labelledSeedImage = null;

    // streamline counts to associated label
    private double[][][] targetSCImage; 

    private double[][][] targetCPImage;

    TargetCP_Image[] seedTargetCP;

    /**
     * Initializes the image with the dimensions of the seed space.
     *
     *
     */
    public ConnectivitySegmentedImage(RegionOfInterest seeds, int[][][] targets,
				      double xVoxelDim, double yVoxelDim, double zVoxelDim) {

	roi = seeds;

	seedVoxels = roi.getSeedVoxels(xVoxelDim, yVoxelDim, zVoxelDim);

	numSeeds = seedVoxels.length;

	seedTargetCP = new TargetCP_Image[numSeeds];

	for (int s = 0; s < numSeeds; s++) {
	    seedTargetCP[s] = new TargetCP_Image(targets, xVoxelDim, yVoxelDim, zVoxelDim);
	}

	this.targets = targets;


	xDataDim = targets.length;
	yDataDim = targets[0].length;
	zDataDim = targets[0][0].length;

	this.xVoxelDim = xVoxelDim;
	this.yVoxelDim = yVoxelDim;
	this.zVoxelDim = zVoxelDim;

    }




    /**
     * Add some tracts to the segmentation, relating to a particular seed point. 
     *
     * @param seedIndex ranging from 0 to the maximum seed index, where index i means the voxel
     * containing seed point i in the RegionOfInterest to be segmented.
     * @param tracts 
     */
    public void processTracts(int seedIndex, TractCollection tracts) {

	// reset these
	labelledSeedImage = null;
	targetSCImage = null;
	targetCPImage = null;

        // This implementation allows you to add tracts one at a time, which means
        // it keeps all connectivity information in memory until all seed points
        // are done. For very large ROIs with many targets, it might be necessary to 
        // process all tracts for a particular seed together

	seedTargetCP[seedIndex].processTracts(tracts);


    }

    


    /**
     * @return an image <code>im</code> where the seed voxels are labelled with the target 
     * ROI to which they have maximum connectivity.
     *
     */
    public double[][][] getSegmentedSeeds() {


	if (labelledSeedImage == null) {

	    double[][][] image = new double[xDataDim][yDataDim][zDataDim];
	    double[][][] sc =  new double[xDataDim][yDataDim][zDataDim];
	    double[][][] cp =  new double[xDataDim][yDataDim][zDataDim];

	    for (int s = 0; s < numSeeds; s++) {
		
		double norm = (double)seedTargetCP[s].totalStreamlines();
		
		int x = seedVoxels[s].x;
		int y = seedVoxels[s].y;
		int z = seedVoxels[s].z;
		
                
                int[][] nonZeroTargetSC = seedTargetCP[s].getNonZeroStreamlineCounts();
       
                int maxSC = 0;
                int maxLabel = 0;

                for (int i = 0; i < nonZeroTargetSC.length; i++) {
                    if (nonZeroTargetSC[i][1] > maxSC) {
                        maxLabel = nonZeroTargetSC[i][0];
                        maxSC = nonZeroTargetSC[i][1];
                    }
                }
                
		    
                image[x][y][z] = maxLabel;
                sc[x][y][z] = maxSC;
                cp[x][y][z] = norm > 0.0 ? maxSC / norm : 0.0;
		    
            }
	    
	    
	    targetCPImage = cp;
	    targetSCImage = sc;
	    labelledSeedImage = image;
	    
	}
	
	double[][][] defCopy = new double[xDataDim][yDataDim][zDataDim];
	    
        for (int i = 0; i < xDataDim; i++) {
            for (int j = 0; j < yDataDim; j++) {
                System.arraycopy(labelledSeedImage[i][j], 0, defCopy[i][j], 0, zDataDim);
            }
        }

	return defCopy;

    }


    /**
     * @return the seeds labelled with the connection probability to the target with the highest
     * connectivity to the seed. If a seed voxel s is labelled with target N, then the value of 
     * s returned from this image is the connection probability to target N.
     */
    public double[][][] getMaxTargetCP() {

	if (targetCPImage == null) {
	    getSegmentedSeeds();
	}

	double[][][] defCopy = new double[xDataDim][yDataDim][zDataDim];
	    
        for (int i = 0; i < xDataDim; i++) {
            for (int j = 0; j < yDataDim; j++) {
                System.arraycopy(targetCPImage[i][j], 0, defCopy[i][j], 0, zDataDim);
            }
        }
        

	return defCopy;

	
    }

    /**
     * @return the seeds labelled with the streamline count to the target with the highest
     * connectivity to the seed. 
     *
     */
    public double[][][] getMaxTargetSC() {

	if (targetSCImage == null) {
	    getSegmentedSeeds();
	}

	double[][][] defCopy = new double[xDataDim][yDataDim][zDataDim];
	    
        for (int i = 0; i < xDataDim; i++) {
            for (int j = 0; j < yDataDim; j++) {
                System.arraycopy(targetSCImage[i][j], 0, defCopy[i][j], 0, zDataDim);
            }
        }
    

	return defCopy;

	
    }





    /**
     * The default behaviour is to count the first entry of a streamline to a target zone. 
     * If this method is called with a <code>false</code> parameter, then streamlines are
     * allowed to connect to any number of target zones.
     *
     * @see #setDirectional(numerics.Vector3D)
     */
    public void setCountFirstEntry(boolean b) {
	countFirstEntry = b;
	
	for (int s = 0; s < numSeeds; s++) {
	    seedTargetCP[s].setCountFirstEntry(countFirstEntry);
	}
    }


}

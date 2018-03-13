package tractography;

import numerics.*;

import java.util.*;

/**
 * Connection probability image with targets. Connection probabilities are defined to targets 
 * and not to voxels.
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public class TargetCP_Image extends ConnectivityMatrix {


    // total number of streamlines used to make this image
    private int totalStreamlines = 0;
   

    // update only the target that the streamline hits first
    private boolean countFirstEntry = true;


    // Streamline counts
    private HashMap<Integer, int[]> targetSC;


    /**
     * Initializes the image with the dimensions of the target space.
     *
     */
    public TargetCP_Image(int[][][] targets, double[] voxelDims) {

	super(targets, voxelDims);

        targetSC = new HashMap<Integer, int[]>(20, 0.75f);

    }


    /**
     * Initializes the image with optional target names.
     *
     * @param labelFile should contain one label ID and name per line, separated by comma or 
     * white space, eg "1 some_brain_area"
     *
     */
    public TargetCP_Image(int[][][] targets, double[] voxelDims, String labelFile) {

	super(targets, voxelDims, labelFile);

        targetSC = new HashMap<Integer, int[]>(20, 0.75f);

    }


    /**
     * Initializes the image with the dimensions of the target space.
     *
     */
    public TargetCP_Image(int[][][] targets, double xVoxelDim, double yVoxelDim, double zVoxelDim) {

	this(targets, new double[] {xVoxelDim, yVoxelDim, zVoxelDim});

    }


    /**
     * Add some streamlines to the image.
     *
     */
    public final void processTracts(TractCollection tc) {

	for (int tCounter = 0; tCounter < tc.numberOfTracts(); tCounter++) {

	    processTract(tc.getTract(tCounter));
	}
    }


    /**
     * Add a single streamline to the image.
     *
     */
    private final void processTract(Tract t) {
	
	totalStreamlines++;
	
	VoxelList voxelList = t.toVoxelList(xVoxelDim, yVoxelDim, zVoxelDim);
	
	Voxel[] voxels = voxelList.getVoxels();
	
	int voxelSeedIndex = voxelList.seedPointIndex();
	int tractSeedIndex = t.seedPointIndex();
	
	int numVoxels = voxels.length;
	
	int numPoints = t.numberOfPoints();
	    

		
        if (countFirstEntry) {
            
            double forwardDistance = -1.0;
            double backwardDistance = -1.0;
            
            int forwardTarget = 0;
            int backwardTarget = 0;
            
            for (int p = tractSeedIndex; p < numPoints; p++) {
                
                Point3D point = t.getPoint(p);
                
                int i = (int)(point.x / xVoxelDim);
                int j = (int)(point.y / yVoxelDim);
                int k = (int)(point.z / zVoxelDim);

                int voxelTarget = targets[i][j][k];
		
                if (voxelTarget > 0) {
                    forwardDistance = t.pathLengthFromSeed(p);
                    forwardTarget = voxelTarget;
                    break;
                }
            }
            
            // count down to 0 from seed
            for (int p = tractSeedIndex; p >= 0; p--) {

                Point3D point = t.getPoint(p);
                
                int i = (int)(point.x / xVoxelDim);
                int j = (int)(point.y / yVoxelDim);
                int k = (int)(point.z / zVoxelDim);
                
                int voxelTarget = targets[i][j][k];
		
                if (voxelTarget > 0) {
                    backwardDistance = t.pathLengthFromSeed(p);
                    backwardTarget = voxelTarget;
                    break;
                }
            }

            // lazy initialize target SC, since there might be many targets
            
            int targetToUpdate = -1;

            if (forwardDistance < 0.0) {
                if (backwardDistance > 0.0) {
                    targetToUpdate = backwardTarget;
                }
            }
            else if (backwardDistance < 0.0) {
                if (forwardDistance > 0.0) {
                    targetToUpdate = forwardTarget;
                }
            }
            else if (forwardDistance < backwardDistance) {
                targetToUpdate = forwardTarget;
            }
            else { // happens if distances are precisely equal
                targetToUpdate = backwardTarget;
            }
            
            if (targetToUpdate > 0) {
                int[] sc = targetSC.get(new Integer(targetToUpdate));
                
                if (sc == null) {
                    targetSC.put(new Integer(targetToUpdate), new int[] {1});
                }
                else {
                    sc[0] = sc[0] + 1;
                }
            }
            
        }
        else {
            // increase cp to all connected targets

            // Collect all targets intersected, this insures we only update each one
            // once
            HashMap<Integer, Boolean> hitTargets = new HashMap<Integer, Boolean>();
            
            int vl = voxels.length;
            
            int v = 0;
            
           
            while (v < vl) {
                
                int voxelTarget = targets[voxels[v].x][voxels[v].y][voxels[v].z];
                
                if (voxelTarget > 0) {
                    
                    Integer key = new Integer(voxelTarget);
                
                    if ( !hitTargets.containsKey(key) ) {
                        hitTargets.put(new Integer(voxelTarget), new Boolean(true));
                    }
                }
                
                v++;
		
            }

            // Now increment count for relevant targets
            
            for (Integer label : hitTargets.keySet()) {
                int[] sc = targetSC.get(label);

                if (sc == null) {
                    targetSC.put(label, new int[] {1});
                }
                else {
                    sc[0] = sc[0] + 1;
                }
            }
        }


    }
    




    /**
     * Get the streamline counts as an image. Each voxel of a target has the same value,
     * which is the number of streamlines that enter that target.
     *
     */
    public double[][][] getStreamlineCounts() {

	double[][][] scImage = null;

        scImage = new double[xDataDim][yDataDim][zDataDim];
	
	for (int i = 0; i < xDataDim; i++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int k = 0; k < zDataDim; k++) {
		    if (targets[i][j][k] > 0) {

                        int[] sc = targetSC.get( new Integer(targets[i][j][k]) );

                        if (sc == null) {
                            scImage[i][j][k] = 0;
                        }
                        else {
                            scImage[i][j][k] = sc[0];
                        }
		    }
		}
	    }
	}

	return scImage;
    }



    /**
     * Get the connection probabilities as an image, which are the streamline counts
     * divided by the total number of streamlines.
     *
     */
    public double[][][] getConnectionProbabilities() {
	
        double[][][] scImage = getStreamlineCounts();

        double[][][] cpImage =  new double[xDataDim][yDataDim][zDataDim];

	double norm = (double)(totalStreamlines);

	for (int i = 0; i < xDataDim; i++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int k = 0; k < zDataDim; k++) {
                    if (scImage[i][j][k] > 0) {
                        cpImage[i][j][k] = scImage[i][j][k] / norm; 
                    }
		    
		}
	    }
	}

	return cpImage;

    }


    /**
     * @return the total number of streamlines processed by <code>processTracts</code>.
     *
     */
    public int totalStreamlines() {
	return totalStreamlines;
    }


    /**
     * Gets the streamline count for a given target.
     *
     */
    public int getStreamlineCount(int label) {

	// get the index of the label in the re-mapped image
	Integer i = targetMap.get(new Integer(label));

        return targetSC.get(i)[0];
    }


    /**
     * Gets all streamline counts for targets with non-zero counts
     *
     * @return an array of dimensions [N][2], where N is the number of targets with non-zero
     * counts. array[n][0] is the target label and array[n][1] is the number of streamlines 
     * connecting to that target. The list is sorted by target label intensity.
     */
    public int[][] getNonZeroStreamlineCounts() {
        Integer[] keys = targetSC.keySet().toArray(new Integer[0]); 

	Arrays.sort(keys);

        int numKeys = keys.length;
        
        int[][] sc = new int[numKeys][2];

        for (int i = 0; i < numKeys; i++) {
            sc[i][0] = targetUnMap.get(keys[i]);
            sc[i][1] = targetSC.get(keys[i])[0];
        }

        return sc;

    }




    /**
     * Get the streamline counts for all targets.
     * 
     * @return the streamline counts in order of label intensity.  
     * 
     */
    public int[] getStreamlineCountList() {

	int[] scArray = new int[numNodes];

	for (int i = 0; i < numNodes; i++) {

	    int[] sc = targetSC.get(new Integer(i+1));

	    if (sc != null) {
		scArray[i] = sc[0];
	    }
	}

	return scArray;
    }



    /**
     * Get the list of target intensities, in order. This gives you
     * the original target labels in the same order that the streamline counts would appear
     * if you called getStreamlineCountList().
     * 
     * @return the target labels in order of label intensity.  
     * 
     */
    public int[] getTargetLabels() {

	int[] labels = new int[numNodes];

	for (int i = 0; i < numNodes; i++) {
	    labels[i] = targetUnMap.get(new Integer(i+1)).intValue();
	}

	return labels;
    }


    /**
     * The default behaviour is to count the first entry of a streamline to a target zone. 
     * If this method is called with a <code>false</code> parameter, then streamlines are
     * allowed to connect to any number of target zones.
     *
     */
    public void setCountFirstEntry(boolean b) {
	countFirstEntry = b;
    }


    /**
     * Resets all streamline counts to zero
     *
     */
    public void reset() {
	targetSC = new HashMap<Integer, int[]>(20, 0.75f);
	totalStreamlines = 0;
    }

}

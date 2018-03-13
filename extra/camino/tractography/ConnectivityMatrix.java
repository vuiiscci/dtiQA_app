package tractography;

import misc.*;
import numerics.*;

import java.io.*;
import java.util.*;
import java.text.*;


/**
 * Connectivity matrix superclass. Contains basic methods to map target ROIs in arbitrary order
 * to the range 1...N, give them names, or read a list of names from disk.
 *
 * @author Philip Cook
 */
public abstract class ConnectivityMatrix  {

    /** 
     * Image of target regions, these are nodes of the connectivity graph. 
     * These are re-mapped to the range 1...N for internal use 
     */
    protected int[][][] targets;


    /** Target names in same order as re-mapped target labels. */
    protected String[] targetNames;


    /** Maps target integers from their values on disk to the range 1...N. */
    protected HashMap<Integer, Integer> targetMap;
    

    /** Maps target integers back to their original values */
    protected HashMap<Integer, Integer> targetUnMap;
    
    /** Number of nodes in the graph, ie the number of target regions. */
    protected int numNodes;


    // voxel dims of seed space
    protected final double xVoxelDim;
    protected final double yVoxelDim;
    protected final double zVoxelDim;


    // dimensions of seed space
    protected final int xDataDim;
    protected final int yDataDim;
    protected final int zDataDim;


   
    /**
     * Unlabeled targets will be given default names for output
     *
     */
    public ConnectivityMatrix(int[][][] targetROIs, double[] voxelDims) {
        xDataDim = targetROIs.length;
	yDataDim = targetROIs[0].length;
	zDataDim = targetROIs[0][0].length;

        xVoxelDim = voxelDims[0];
        yVoxelDim = voxelDims[1];
        zVoxelDim = voxelDims[2];

        setNodes(targetROIs);
        
    }


    /**
     * Set up a connectivity matrix with label names
     * 
     * @param labelNameFile should contain one label ID and name per line, separated by white space, eg "1 some_brain_area"
     */
    public ConnectivityMatrix(int[][][] targetROIs, double[] voxelDims, String labelNameFile) {
       
        xDataDim = targetROIs.length;
	yDataDim = targetROIs[0].length;
	zDataDim = targetROIs[0][0].length;

        xVoxelDim = voxelDims[0];
        yVoxelDim = voxelDims[1];
        zVoxelDim = voxelDims[2];
        
        setNodes(targetROIs, labelNameFile);   
    }


    /**
     * Set the targets to use as nodes of the graph. Each positive integer-valued region 
     * is one node. They need not be contiguous.
     */
    private void setNodes(int[][][] targetROIs) {

        // get list of targets that exist
        
        HashMap<Integer, String> targetsInImage = new HashMap<Integer, String>();

        int numberOfTargets = 0;

	DecimalFormat dfInt = new DecimalFormat("0000");

	for (int i = 0; i < xDataDim; i++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int k = 0; k < zDataDim; k++) {
		    if (targetROIs[i][j][k] > 0) {
                        
                        Integer label = new Integer(targetROIs[i][j][k]);

                        if (targetsInImage.get(label) == null) {
                            targetsInImage.put(label, "Label_" + dfInt.format(label));
                            
                            numberOfTargets++;
                        }
                    }
                }
            }
        }


        Integer[] keys = targetsInImage.keySet().toArray(new Integer[0]);

        Arrays.sort(keys);

        targetNames = new String[keys.length];
        
	targetMap = new HashMap<Integer, Integer>(keys.length * 2);

	targetUnMap = new HashMap<Integer, Integer>(keys.length * 2);


        for (int i = 0; i < keys.length; i++) {

            targetMap.put(keys[i], new Integer(i+1));
	    targetUnMap.put(new Integer(i+1), keys[i]);

            targetNames[i] = targetsInImage.get(keys[i]);
        }
     

	remapTargets(targetROIs);

    }


    /**
     *
     * Remap targets to a continuous series of integers. Targets that have no name will be discarded.
     *
     * @param targetROIs target image, these will be mapped to the range 1...N where there
     * are N label names.
     *
     */
    private void remapTargets(int[][][] targetROIs) {
	
        targets = new int[xDataDim][yDataDim][zDataDim];
	
	for (int i = 0; i < xDataDim; i++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int k = 0; k < zDataDim; k++) {
		    if (targetROIs[i][j][k] > 0) {
                        Integer index = targetMap.get(new Integer(targetROIs[i][j][k]));
                        
			if (index == null) {
                            targets[i][j][k] = 0;
			}
			else {
                            targets[i][j][k] = index;
                        }
                    }
                }
            }
        }
	
        numNodes = targetNames.length;

   }

    
    /**
     * Set targets and a list of names, these will be used as column names for the matrix. 
     *
     * @param targetListFile a CSV file containing a list of targets
     *
     */
    private void setNodes(int[][][] targetROIs, String targetListFile) {
        
       
        Scanner scanner;

	try {
	    scanner = new Scanner(new File(targetListFile));
	}        
	catch (FileNotFoundException e) {
	    throw new LoggedException("Can't find target definition file " + targetListFile);
	}
	
        // Get each line of the file
        scanner.useDelimiter("\r\n|\n");

        // because the label file might not be in order
        HashMap<Integer, String> labelNames = new HashMap<Integer, String>();

        while(scanner.hasNext()) {
            // ignore empty lines, trim leading / trailing space
            
            String next = scanner.next().trim();
            
            if (next.length() > 0) {
                String[] tmp = next.split("[,\\s]+", 2);
                Integer key = Integer.parseInt(tmp[0]);
                
                labelNames.put(key, tmp[1]);
            }
        }

        scanner.close();
        
        Integer[] keys = labelNames.keySet().toArray(new Integer[0]);

        Arrays.sort(keys);

        targetNames = new String[keys.length];
        
        targetMap = new HashMap<Integer, Integer>(keys.length * 2);

	targetUnMap = new HashMap<Integer, Integer>(keys.length * 2);

        for (int i = 0; i < keys.length; i++) {

            targetMap.put(keys[i], new Integer(i+1));
	    targetUnMap.put(new Integer(i+1), keys[i]);

            targetNames[i] = labelNames.get(keys[i]);
        }
        
        remapTargets(targetROIs);

    }



    public final String[] getTargetNames() {
        String[] copy = new String[numNodes];

	System.arraycopy(targetNames, 0, copy, 0, numNodes);

	return copy;
    }


    public final int numNodes() {
	return numNodes;
    }

        
}

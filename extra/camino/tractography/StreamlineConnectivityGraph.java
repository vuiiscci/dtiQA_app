package tractography;

import misc.*;
import numerics.*;

import java.io.*;
import java.util.*;
import java.text.*;


/**
 * Connectivity matrix derived from streamlines. The entries for the connectivity can be 
 * based on streamline counts, or derived from tract statistics. The graph is defined by labeled
 * ROIs, connected by edges (streamlines). Tract statistics are always averaged by the streamline count,
 * so if you set the tract statistic to length, the output matrix will give you the average length of 
 * streamlines forming each edge. Likewise you can ask for the average tract-mean FA, tract-max FA, and so on.
 *
 * @author Philip Cook
 */
public class StreamlineConnectivityGraph extends ConnectivityMatrix {


    private TractStatisticFilter statFilter;

    private boolean useStatFilter = false;

    // total streamlines added
    private int totalStreamlines = 0;


    private final DecimalFormat df = new DecimalFormat("0.00000E00");

    
    // Use matrices for speed - acceptable memory use up to a few thousand ROIs
    //
    // Reduce memory by using a symmetric matrix - also avoids having to keep track of [i][j] vs [j][i]
    
    // number of streamlines connecting an edge
    private SymmetricMatrix streamlineCounts;


    // optional, used for tract stats
    // keeps a running sum of tract FA, etc. Use streamline counts to normalize
    private SymmetricMatrix tractStats = null;


   
    /**
     * Unlabeled targets will be given default names for output
     *
     */
    public StreamlineConnectivityGraph(int[][][] targetROIs, double[] voxelDims) {
        super(targetROIs, voxelDims);
	streamlineCounts = new SymmetricMatrix(numNodes);
        
    }

    /**
     * Set up a connectivity graph with label names
     * 
     * @param labelNameFile should contain one label ID and name per line, separated by white space, eg "1 some_brain_area"
     */
    public StreamlineConnectivityGraph(int[][][] targetROIs, double[] voxelDims, String labelNameFile) {
	super(targetROIs, voxelDims, labelNameFile);  
	streamlineCounts = new SymmetricMatrix(numNodes);
 
    }


    /**
     * Gets the matrix of connectivity values as a RealMatrix
     *
     */
    public RealMatrix getStreamlineCountMatrix() {
        return streamlineCounts.toRealMatrix();
    }


    /**
     * Gets the matrix of tract statistics as a RealMatrix
     *
     */
    public RealMatrix getTractStatisticMatrix() {

        if (!useStatFilter) {
            return null;
        }
                                    
        SymmetricMatrix norm = new SymmetricMatrix(numNodes);
       
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j <= i; j++) {

		double sc = streamlineCounts.get(i,j);

                norm.set(i,j, sc > 0.0 ? tractStats.get(i,j) / sc : 0.0);
            }
        }
        
        return norm.toRealMatrix();
    }



    /**
     * Gets a matrix of connectivity values as a CSV string. The matrix is symmetric.
     *
     */
    public String getStreamlineCountMatrixCSV() {
        return getMatrixCSV(streamlineCounts);
    }

    /**
     * Gets the normalized tract statistic matrix; divides each edge value by the number of streamlines
     *
     */
    public String getTractStatisticMatrixCSV() {

        if (!useStatFilter) {
            return null;
        }

        SymmetricMatrix norm = new SymmetricMatrix(numNodes);
       
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j <= i; j++) {
                
                double normC = streamlineCounts.get(i,j);

                if (normC == 0.0) {
                    // no connections
                    norm.set(i,j,0.0);
                }
                else {
                    norm.set(i,j, tractStats.get(i,j) / streamlineCounts.get(i,j));
                }
            }
        }
        
        return getMatrixCSV(norm);
    }


    /**
     * Gets the Matrix in CSV format, with column names - so m's columns must match
     */
    private String getMatrixCSV(SymmetricMatrix m) {

        StringBuilder b = new StringBuilder();
        
        int size = m.size();

        b.append(targetNames[0]);

        for (int i = 1; i < size; i++) {
            b.append(",");
            b.append(targetNames[i]);
        }

        b.append("\n");

        for (int i = 0; i < size; i++) {

            b.append(df.format(m.get(i,0)));

            for (int j = 1; j < size; j++) {
                b.append(",");
                b.append(df.format(m.get(i,j)));
            }
            
            b.append("\n");
        }
        

        return b.toString();
        
    }


    /**
     * Set the tract stat filter used to generate the tract stat matrix. Only univariate tract stats are 
     * supported.
     *
     */
    public void setTractStatFilter(TractStatisticFilter filter) {

        statFilter = filter;
        useStatFilter = true;

        tractStats = new SymmetricMatrix(numNodes);
    }


    /**
     *
     * Process a bunch of tracts
     *
     */
    public void processTracts(TractCollection tc) {
        for (int i = 0; i < tc.numberOfTracts(); i++) {
            processTract(tc.getTract(i));
        }
    }

    


    /**
     * Add a single streamline to the image. The streamline is added to the graph if it connects two nodes. The first pair of connections
     * is used, one in each direction from the seed point.
     * <p>
     * If a seed point is within a target region, then that counts as one node, but not both. The tract should connect to exactly one target
     * outside the seed region. If this isn't the case then the tract is not counted. 
     * <p>
     * For seed points outside a target region, the same target may be connected at both ends, ie, the tract will update a diagonal entry on the 
     * connectivity matrix.
     *
     */
    public final void processTract(Tract t) {
        
        totalStreamlines++;
	
	int tractSeedIndex = t.seedPointIndex();
	
	int numPoints = t.numberOfPoints();
	    
        int forwardTarget = 0;
        int backwardTarget = 0;

        int forwardTargetIndex = -1;
        int backwardTargetIndex = -1;
      
        Point3D point = t.getPoint(tractSeedIndex);
            
        int i = (int)(point.x / xVoxelDim);
        int j = (int)(point.y / yVoxelDim);
        int k = (int)(point.z / zVoxelDim);
        
        int seedTarget = targets[i][j][k];
        
        
        for (int p = tractSeedIndex; p < numPoints; p++) {
            
            point = t.getPoint(p);
            
            i = (int)(point.x / xVoxelDim);
            j = (int)(point.y / yVoxelDim);
            k = (int)(point.z / zVoxelDim);
            
            int voxelTarget = targets[i][j][k];
            
            if (voxelTarget > 0 && voxelTarget != seedTarget) {
                forwardTarget = voxelTarget;
                forwardTargetIndex = p;
                break;
            }
        }
        
        // count down to 0 from seed
        for (int p = tractSeedIndex; p >= 0; p--) {
            
            point = t.getPoint(p);
            
            i = (int)(point.x / xVoxelDim);
            j = (int)(point.y / yVoxelDim);
            k = (int)(point.z / zVoxelDim);
            
            int voxelTarget = targets[i][j][k];
            
            if (voxelTarget > 0 && voxelTarget != seedTarget) {
                backwardTarget = voxelTarget;
                backwardTargetIndex = p;
                break;
            }
        }

        if (seedTarget > 0) {
            
            if (forwardTarget > 0 && backwardTarget > 0) {
                
                // connects three regions, choose closest one to seed
                
                if (t.pathLengthFromSeed(forwardTargetIndex) < t.pathLengthFromSeed(backwardTargetIndex)) {
                    backwardTargetIndex = tractSeedIndex;
                    backwardTarget = seedTarget;
                }
                else { 
                    forwardTargetIndex = tractSeedIndex;
                    forwardTarget = seedTarget;
                }

            }
            else if (forwardTarget == 0) {
                forwardTarget = seedTarget;
            }
            else {
                backwardTarget = seedTarget;
            }
            
        } 

     
        if (forwardTarget > 0 && backwardTarget > 0) {
            streamlineCounts.add(forwardTarget - 1, backwardTarget - 1, 1.0);

            if (useStatFilter) {
                tractStats.add(forwardTarget - 1, backwardTarget - 1, statFilter.processTract(t)[0]);
            }
        } 

    }

        
}

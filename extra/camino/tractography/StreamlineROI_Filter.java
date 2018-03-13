package tractography;

import data.*;
import imaging.*;
import inverters.*;
import misc.*;
import numerics.*;

import java.util.logging.*;
import java.util.*;



/**
 * Takes as input the streamline output from a <code>FibreTracker</code>. Reads streamlines
 * in Camino space, and performs ROI operations using images in that same space.
 *
 *
 *
 * @version $Id$
 * @author  Philip Cook
 */
public class StreamlineROI_Filter {

    // voxel dims of seed space
    private final double xVoxelDim;
    private final double yVoxelDim;
    private final double zVoxelDim;


    // dimensions of seed space
    private final int xDataDim;
    private final int yDataDim;
    private final int zDataDim;

    private int[][][] waypoints;

    private int numWaypoints = 0;

    private int maxWaypointIndex = 0;


    // min and max tract length. If zero, they are ignored, otherwise each tract is evaluated
   
    private int minTractPoints = 0;

    private double minTractLength = 0.0;

    private int maxTractPoints = 0;

    private double maxTractLength = 0.0;

    private int[][][] exclusionROI;

    private boolean haveExclusionROIs = false;

    private int[][][] endZones;

    private boolean haveEndZones = false;
    
    /** 
     * Discard entire streamline that enters an exclusion ROI (true)
     * terminate tracking in an exclusion ROI (false) 
     */
    private boolean discardInExclusionROI = false;


    /**
     * Check streamlines for multiple entry to waypoints, and truncate if found.
     * 
     */
    private boolean truncateLoops = false;


    /**
     * Check streamlines for multiple entry to waypoints, and discard whole streamline if found.
     * 
     */
    private boolean discardLoops = false;


    /**
     * Resample tracts to high resolution. If the step size of the tracts is not small compared
     * to the size of the seed space voxels, then successive points might cross multiple voxels.
     * 
     */
    private boolean resample = false;


    /**
     * Maximum resampled step size.
     * 
     */
    private double resampleStepSize = 0.5;


    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.tractography.StreamlineROI_Filter");



    /**
     * Initializes the filter, with the header corresponding to the ROI images. They must
     * all share the same voxel and physical space with the streamlines.
     *
     */
    public StreamlineROI_Filter(int[] dataDims, double[] voxelDims) {
        
        xDataDim = dataDims[0];
        yDataDim = dataDims[1];
        zDataDim = dataDims[2];

        xVoxelDim = voxelDims[0];
        yVoxelDim = voxelDims[1];
        zVoxelDim = voxelDims[2];

	exclusionROI = new int[xDataDim][yDataDim][zDataDim];
    }


    /**
     * Set the waypoint volume. All positive values in the volume constitute a separate waypoint.
     *
     * @param waypoints a volume in the same voxel space as the streamlines
     */
    public void setWaypoints(int[][][] waypoints) {
	// sanity check
	if (waypoints.length != xDataDim || waypoints[0].length != yDataDim || 
	    waypoints[0][0].length != zDataDim) {
	    throw new LoggedException("Waypoint dimensions " + waypoints.length + " " + 
				      waypoints[0].length +
				      " " + waypoints[0][0].length + 
				      " do not match seed space dimensions " + 
				      xDataDim + " " + yDataDim + " " + zDataDim);
	}

	this.waypoints = waypoints;

	numWaypoints = 0;

        HashMap<Integer, Integer> wpMap = new HashMap<Integer, Integer>(100, 0.75f);

        // index of each waypoint in the order in which we encounter it
        int nextWaypointIndex = 1;

	for (int i = 0; i < xDataDim; i++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int k = 0; k < zDataDim; k++) {

		    if (waypoints[i][j][k] > 0) {		    
                        Integer wp = new Integer(waypoints[i][j][k]);
                        
                        if (!wpMap.containsKey(wp)) {
                            wpMap.put(wp, nextWaypointIndex);
                            nextWaypointIndex++;
                        }
                    }
		}
	    }
	}
        
        Integer[] waypointList = wpMap.keySet().toArray(new Integer[wpMap.size()]);
        
        numWaypoints = waypointList.length;
        
        Arrays.sort(waypointList);
        
        maxWaypointIndex = waypointList[numWaypoints - 1].intValue();
        
        // map waypoints 1-N, since their value doesn't matter but we don't want the maximum value to be too big, 
        // because we want to allocate a boolean array of dimension [maxWaypointIndex + 1] so that we can quickly 
        // mark which waypoints a streamline has hit
        
        if (maxWaypointIndex > Short.MAX_VALUE) {
            if (numWaypoints > Short.MAX_VALUE) { 
                throw new LoggedException("Too many waypoints " + numWaypoints);
            }
            else {
                
                logger.info("Re-mapping " + numWaypoints + " waypoints to range 1..." + numWaypoints + " for processing");
                
                for (int i = 0; i < xDataDim; i++) {
                    for (int j = 0; j < yDataDim; j++) {
                        for (int k = 0; k < zDataDim; k++) {
                            
                            if (waypoints[i][j][k] > 0) {
                                
                                int remapped = wpMap.get(new Integer(waypoints[i][j][k])).intValue();
                                
                                waypoints[i][j][k] = remapped;
                            }
                        }
                    }
                }
                
            }
            
        }


    }
    

    /**
     * Set the end zones volume. All positive values in the volume constitute an end zone. No distinction
     * is made between end zones.
     *
     * @param endZones a volume in the same voxel space as the streamlines
     */
    public void setEndZones(int[][][] endZones) {
	// sanity check
	if (endZones.length != xDataDim || endZones[0].length != yDataDim || 
	    endZones[0][0].length != zDataDim) {
	    throw new LoggedException("EndZone dimensions " + endZones.length + " " + 
				      endZones[0].length + " " + endZones[0][0].length + 
				      " do not match seed space dimensions " + 
				      xDataDim + " " + yDataDim + " " + zDataDim);
	}

	this.endZones = endZones;

	haveEndZones = true;

    }
    


    /**
     * Set the exclusion volume. Any positive value is an exclusion ROI, no distinction is 
     * made between different values.
     *
     * @param exclusionROI a volume in the same voxel space as the streamlines
     */
    public void setExclusionROIs(int[][][] exclusionROI) {
	// sanity check
	if (exclusionROI.length != xDataDim || exclusionROI[0].length != yDataDim || 
	    exclusionROI[0][0].length != zDataDim) {
	    throw new LoggedException("Exclusion ROI dimensions " + exclusionROI.length + " " + 
				      exclusionROI[0].length + " " + exclusionROI[0][0].length + 
				      " do not match seed space dimensions " + 
				      xDataDim + " " + yDataDim + " " + zDataDim);
	}

	this.exclusionROI = exclusionROI;

	haveExclusionROIs = true;

    }


    /**
     *
     * @return a TractCollection that either contains <code>t</code> (truncated if necessary) or
     * is empty.
     *
     */
    public TractCollection processTract(Tract t) {
	TractCollection container = new TractCollection(2, 100.0);
	container.addTract(t);
	return processTracts(container);
    }


    /**
     *
     * @return all <code>Tract</code>'s that pass through all waypoints, with excluded tracts 
     * trimmed or discarded.
     */
    public TractCollection processTracts(TractCollection tc) {

	TractCollection processed = new TractCollection(tc.numberOfTracts() + 1, 100.0);

	for (int fibre = 0; fibre < tc.numberOfTracts(); fibre++) {
	    
	    Tract tract = tc.getTract(fibre);

	    if (tract.numberOfPoints() < minTractPoints) {
		continue;
	    }
	    if (maxTractPoints > 0 && tract.numberOfPoints() > maxTractPoints) {
		// chop tract down to maximum length
		tract.truncateToMaxPoints(maxTractPoints);
	    }

	    // first condition avoids computation of length if it is not needed

	    if (minTractLength > 0.0 && tract.length() < minTractLength) {
		continue;
	    }
	    if (maxTractLength > 0.0 && tract.length() > maxTractLength) {
		tract.truncateToMaxLength(maxTractLength);
	    }

            // Gets rid of any points outside the image volume. This operation removes all points outside
            // the boundary of the image, whereas truncate to exclusion / end points will leave one point 
            // inside those regions (because they might also be waypoints)
            truncateToBounds(tract);

            if (tract.numberOfPoints() == 0) {
                // nothing left of this tract
                continue;
            }

            // Now resample if required
	    if (resample) {
		tract = tract.resample(resampleStepSize);
	    }

	    boolean entersExclusion = false;
            
            // truncates upon entry to exclusion ROIs
            if (haveExclusionROIs) {
                entersExclusion = truncateToExclusion(tract);
            }

            
	    boolean passesWaypoints = true;

	    if (numWaypoints > 0) {
		passesWaypoints = passesWaypoints(tract, truncateLoops, discardLoops);
	    }
	    

	    // always true if there are no end zones
	    boolean connectsEndZones = true;
	    
	    // if we go through the waypoints, check for exclusion
	    if (passesWaypoints) {
		
		if (haveEndZones) {
		    connectsEndZones = truncateToEndZones(tract);
		}
	    }

	    if (passesWaypoints && connectsEndZones && !(entersExclusion && discardInExclusionROI)) {
               
                processed.addTract(tract);
	    }

	}

	return processed;

    }

    
    /**
     * @return true if the streamline enters all waypoints.
     * 
     * @param t the Tract to check
     * @param truncate if true, truncate streamlines that loop, ie those that pass through
     * a particular waypoint more than once.
     *
     * @param discard if true, discard streamlines that loop, ie those that pass through
     * a particular waypoint more than once.
     *
     * 
     *
     */
    private final boolean passesWaypoints(Tract t, boolean truncate, boolean discard) {

        if (numWaypoints == 0) {
            return true;
        }


        boolean[] hitWaypoint = new boolean[maxWaypointIndex + 1];
                    
        int waypointsHit = 0;
           
	if (truncateLoops || discardLoops) {

	    Point3D[] points = t.getPoints();
	    
	    int seedIndex = t.seedPointIndex();
	    
	    int minPointIndex = 0;
	    int maxPointIndex = points.length - 1;
	    
	    int currentPoint = seedIndex;


	    // We want to catch streamlines that enter a waypoint multiple times
	    // set this to the waypoint index upon entry, so that we continue to 
	    // count points inside the waypoint. Reset to 0 upon leaving the waypoint.
	    // This tells us if we are traversing contiguous waypoint voxels
	    // or have left and returned to the waypoint.
	    int inWaypoint = 0;
	 
	    while (currentPoint >= 0) {
		
		int x = (int)(points[currentPoint].x / xVoxelDim);
		int y = (int)(points[currentPoint].y / yVoxelDim);
		int z = (int)(points[currentPoint].z / zVoxelDim);
		
		int wp = waypoints[x][y][z];

		countdown: 
		if (wp > 0) {
		    
		    if (hitWaypoint[wp]) {
			if (inWaypoint != wp) {
			    // back to a place we hit before, terminate here
			    
			    if (discardLoops) {
				return false;
			    }
			    
			    minPointIndex = currentPoint;
			    break countdown;
			}
		    }
		    else {
			hitWaypoint[wp] = true;
			waypointsHit++;
			inWaypoint = wp;
		    }
		    
		}
		else {
		    inWaypoint = 0;
		}

		currentPoint--;
		
	    }
	    
	    currentPoint = seedIndex;

	    inWaypoint = 0;

	    
	    while (currentPoint < points.length) {
		
		int x = (int)(points[currentPoint].x / xVoxelDim);
		int y = (int)(points[currentPoint].y / yVoxelDim);
		int z = (int)(points[currentPoint].z / zVoxelDim);

		int wp = waypoints[x][y][z];

		// If the seed is in a waypoint, waypoint intersectoin is recorded already
		if (currentPoint == seedIndex) {
		    inWaypoint = wp;
		}
		
		countup:
		if (wp > 0) {

		    if (hitWaypoint[wp]) {

			if (inWaypoint != wp) {
			    // back to a place we hit before, terminate here

			    if (discardLoops) {
				return false;
			    }
			    
			    maxPointIndex = currentPoint;
			    break countup;
			}
		    }
		    else {
			hitWaypoint[wp] = true;
			waypointsHit++;
			inWaypoint = wp;
		    }
		}
		else {
		    inWaypoint = 0;
		}
		
		currentPoint++;
	    }
	    
	    if (minPointIndex > 0 || maxPointIndex < points.length - 1) {
		t.chop(minPointIndex, maxPointIndex);
	    }
	    
	}
	else {

	    int voxelCounter = 0;

	    VoxelList voxels = t.toVoxelList(xVoxelDim, yVoxelDim, zVoxelDim);

	    int numVoxels = voxels.size();
	    
	    while (waypointsHit < numWaypoints && voxelCounter < numVoxels) {

		Voxel v = voxels.getVoxel(voxelCounter);
		
		int point = waypoints[ v.x ][ v.y ][ v.z ];
		
		if (point > 0 && !hitWaypoint[point]) {
		    hitWaypoint[point] = true;
		    waypointsHit++;
		}
		
		voxelCounter++;
		
	    }
	 
	}   

	if (waypointsHit != numWaypoints) {
	    return false;
	}
	
	return true;
        
    }

    

    /**
     * @return true if the streamline enters an exclusion ROI.
     *
     */
    private final boolean entersExclusion(Tract t) {

	if (!haveExclusionROIs) {
	    return false;
	}

	// check for seeds placed outside volume
	if (t.numberOfPoints() == 1) {
	    Point3D p = t.getPoint(0);

	    int x = (int)(p.x / xVoxelDim);
	    int y = (int)(p.y / yVoxelDim);
	    int z = (int)(p.z / zVoxelDim);
	    
	    if (x < 0 || x >= xDataDim || y < 0 || y >= yDataDim || z < 0 || z >= zDataDim) {
		return false;
	    }
	    
	}

	VoxelList voxels = t.toVoxelList(xVoxelDim, yVoxelDim, zVoxelDim);

	int length = voxels.size();

	for (int i = 0; i < length; i++) {

	    Voxel v = voxels.getVoxel(i);

	    if (exclusionROI[ v.x ][ v.y ][ v.z ] > 0) {
		return true;
	    }
	}
	
	return false;
    }


    /**
     * Cuts out sections that are outside the dimensions of the image space.
     *
     *
     */
    private final void truncateToBounds(Tract t) {

	Point3D[] points = t.getPoints();

	int seedIndex = t.seedPointIndex();

	int minPointIndex = 0;
	int maxPointIndex = points.length - 1;

	int currentPoint = seedIndex;

	while (currentPoint >= 0) {
	    
	    int x = (int)(points[currentPoint].x / xVoxelDim);
	    int y = (int)(points[currentPoint].y / yVoxelDim);
	    int z = (int)(points[currentPoint].z / zVoxelDim);

	    if (x < 0 || x >= xDataDim || y < 0 || y >= yDataDim || z < 0 || z >= zDataDim) {
		minPointIndex = currentPoint + 1;
		break;
	    }

	    currentPoint--;
	}

	currentPoint = seedIndex;

	while (currentPoint < points.length) {
	    
	    int x = (int)(points[currentPoint].x / xVoxelDim);
	    int y = (int)(points[currentPoint].y / yVoxelDim);
	    int z = (int)(points[currentPoint].z / zVoxelDim);

	    if (x < 0 || x >= xDataDim || y < 0 || y >= yDataDim || z < 0 || z >= zDataDim) {
		maxPointIndex = currentPoint - 1;
		break;
	    }
	   
	    currentPoint++;
	}

	if (minPointIndex >= maxPointIndex) {
	    // seed point outside volume, cut this tract to nothing
            t.truncateToMaxPoints(0);
	}
	else if (minPointIndex > 0 || maxPointIndex < points.length - 1) {
            t.chop(minPointIndex, maxPointIndex);
	}

        

    }



    /**
     * Cut away sections of a tract that pass through an exclusion ROI. Also cuts out sections
     * that are outside the dimensions of the seed space.
     *
     * @return true if the tract enters an exclusion ROI, false otherwise. 
     */
    private final boolean truncateToExclusion(Tract t) {

	Point3D[] points = t.getPoints();

	int seedIndex = t.seedPointIndex();

	int minPointIndex = 0;
	int maxPointIndex = points.length - 1;

	int currentPoint = seedIndex;

	boolean entersExclusion = false;

	while (currentPoint >= 0) {
	    
	    int x = (int)(points[currentPoint].x / xVoxelDim);
	    int y = (int)(points[currentPoint].y / yVoxelDim);
	    int z = (int)(points[currentPoint].z / zVoxelDim);

	    if (exclusionROI[x][y][z] > 0) {
		minPointIndex = currentPoint;
		entersExclusion = true;
		break;
	    }
	   
	    currentPoint--;
	}

	currentPoint = seedIndex;

	while (currentPoint < points.length) {
	    
	    int x = (int)(points[currentPoint].x / xVoxelDim);
	    int y = (int)(points[currentPoint].y / yVoxelDim);
	    int z = (int)(points[currentPoint].z / zVoxelDim);

	    if (exclusionROI[x][y][z] > 0) {
		maxPointIndex = currentPoint;
		entersExclusion = true;
		break;
	    }
	   
	    currentPoint++;
	}

	if (minPointIndex >= maxPointIndex) {
	    // seed point in exclusion
	    t.chop(seedIndex, seedIndex);
	}
	else if (minPointIndex > 0 || maxPointIndex < points.length - 1) {
	    t.chop(minPointIndex, maxPointIndex);
	}

	return entersExclusion;
	
    }


    /**
     * Truncates a tract such that it does not extend past the end zones. Returns true if the Tract 
     * (with truncation, if necessary) connects exactly two end zones and does not extend beyond them, false otherwise.
     * <p>
     * If the seed point is in a labeled region, the method finds the shortest path from the seed to another labeled
     * region. If the streamline does not enter another labeled region, the method returns false.
     * <p>
     * If the seed point is not in a labeled region, the method searches in both directions from the seed 
     * until a labeled region has been encountered in both directions. If the same region is encountered
     * in both directions, the tract is considered a loop and the method returns false. If a different region
     * is encountered in each direction, then the tract is truncated. 
     * <p>
     * When a Tract is truncated, the point of truncation is the last point before the Tract leaves the end
     * zone. 
     *
     * @return true if the tract connects two distinct end zones, false otherwise.
     */
    private final boolean truncateToEndZones(Tract t) {

	Point3D[] points = t.getPoints();

	if (points.length == 1) {
	    return false;
	}
	
	int seedIndex = t.seedPointIndex();

	int minPointIndex = seedIndex;
	int maxPointIndex = seedIndex;

	int currentPoint = seedIndex;


	int seedZone = 
	    endZones[(int)(points[seedIndex].x / xVoxelDim)][(int)(points[seedIndex].y / yVoxelDim)][(int)(points[seedIndex].z / zVoxelDim)];

	int lowerZone = 0; // found at some point between the seed and 0
	int upperZone = 0; // found at some point between the seed and numberOfPoints - 1

	while (currentPoint >= 0) {
	    
	    int x = (int)(points[currentPoint].x / xVoxelDim);
	    int y = (int)(points[currentPoint].y / yVoxelDim);
	    int z = (int)(points[currentPoint].z / zVoxelDim);

	    if (endZones[x][y][z] > 0 && endZones[x][y][z] != seedZone) {
		minPointIndex = currentPoint;
		lowerZone = endZones[x][y][z];
		break;
	    }
	   
	    currentPoint--;
	}

	currentPoint = seedIndex;

	while (currentPoint < points.length) {
	    
	    int x = (int)(points[currentPoint].x / xVoxelDim);
	    int y = (int)(points[currentPoint].y / yVoxelDim);
	    int z = (int)(points[currentPoint].z / zVoxelDim);

	    if (endZones[x][y][z] > 0 && endZones[x][y][z] != seedZone) {
		maxPointIndex = currentPoint;
		upperZone = endZones[x][y][z];
		break;
	    }
	   
	    currentPoint++;
	}

	
	if (seedZone > 0) {
	    
	    // seed point inside end zone

	    if (lowerZone == 0 && upperZone == 0) {
		// tract does not connect two end zones
		return false;
	    }
	    
	    // find shortest path to an end zone
	    double pathLengthToLowerZone = 
		minPointIndex < seedIndex ? t.pathLengthFromSeed(minPointIndex) : Double.MAX_VALUE;

	    double pathLengthToUpperZone = 
		maxPointIndex > seedIndex ? t.pathLengthFromSeed(maxPointIndex) : Double.MAX_VALUE;
	    
	    if (pathLengthToLowerZone < pathLengthToUpperZone) {
		t.chop(minPointIndex, seedIndex);
	    }
	    else {
		t.chop(seedIndex, maxPointIndex);
	    }

	    return true;

	}
	else {

	    if (lowerZone == 0 || upperZone == 0) {
		// tract does not connect two end zones
		return false;
	    }
	    if (lowerZone == upperZone) {
		// tract loops into the same zone
		return false;
	    }
	    else {
		if (minPointIndex > 0 || maxPointIndex < points.length - 1) {
		    t.chop(minPointIndex, maxPointIndex);
		}
		return true;
	    }
	}
    }

    public void setDiscardOnExclusionEntry(boolean discard) {
	discardInExclusionROI = discard;
    }


    /**
     * Determines whether tracts will be resampled.
     * 
     */
    public void setResampleTracts(boolean resampleTracts) {
	resample = resampleTracts;
    }


    /**
     * Filtered streamlines will be resampled to the specified step size.
     *
     */
    public void setResampleStepSize(double size) {
	if (size > 0.0) {
	    resampleStepSize = size;
            resample = true;
	}
	else {
            resample = false;
	}
    }

    /**
     * Automatically remove tracts if they have less than <code>p</code> points
     * before resampling and before any exclusion ROIs are applied.
     *
     */
    public void setMinTractPoints(int p) { 
	minTractPoints = p;
    }


    /**
     * Automatically remove tracts if their length is less than <code>length</code> mm
     * before any exclusion ROIs are applied.
     *
     */
    public void setMinTractLength(double length) { 
	minTractLength = length;
    }


    /**
     * Truncates tracts if they have more than <code>p</code> points
     * before resampling and before any exclusion ROIs are applied.
     *
     */
    public void setMaxTractPoints(int p) { 
	maxTractPoints = p;
    }


    /**
     * Truncates tracts if their length is more than <code>length</code> mm
     * before any exclusion ROIs are applied. 
     *
     */
     public void setMaxTractLength(double length) { 
	 maxTractLength = length;
     }

    

    /**
     * @param trunc if true, truncate loops that enter waypoint regions more than once.
     *
     */
    public void setTruncateLoops(boolean trunc) {
	truncateLoops = trunc;
    }

    /**
     * @param disc if true, discard streamlines with loops that enter waypoint regions more than once.
     *
     */
    public void setDiscardLoops(boolean disc) {
	discardLoops = disc;
	
	if (discardLoops) {
	    truncateLoops = false;
	}
    }

    
}

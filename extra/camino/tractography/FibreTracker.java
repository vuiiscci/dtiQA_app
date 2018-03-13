package tractography;

import numerics.*;

/**
 *  Abstract superclass for classes that perform streamline tractography. 
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public abstract class FibreTracker {

    protected short[][][] visitedVoxel = null; 

    public static final boolean FORWARD = true;
    public static final boolean BACKWARD = false;
    
    protected final double xVoxelDim;
    protected final double yVoxelDim;
    protected final double zVoxelDim; // mm per voxel
    
    protected final int xDataDim;
    protected final int yDataDim;
    protected final int zDataDim;
    
    // Upper boundaries for tracking, in mm
    private final double xUBound;
    private final double yUBound;
    private final double zUBound;

    protected final boolean[][][] isotropic; 
         
    protected final TractographyImage image; 

    private short[] zeros;
    

    /** Check that we don't curve more than ipThreshold over some scale (mm) */
    private double checkCurveInterval = 10.0;


    /** Maximum curvature between successive steps */
    private double ipThreshold = 0.0;



    // /** 
    //  * The loop check stops tracking if we previously visited a voxel, unless we visited it less than
    //  * LOOP_ALLOWANCE steps ago. This is steps, not length, so its physical length will vary a bit depending on 
    //  * step size. It's designed to stop loops, but allow tracking that is tangential to a voxel boundary.
    //  */
    public static final short LOOP_ALLOWANCE = 5;


    /** 
     * Maximum number of points in a streamline. This is an internal limit that cannot be changed. At 0.1mm step size,
     * this equates to 3.27 m tracts, longer than anything in the body.   
     *
     */
    public static final short MAX_TRACT_POINTS=32767;


    /**
     * Sets an upper bound on the length of fiber segments in mm.
     * 
     */
    private double maxSegmentLength = 1000.0;
    


    /** 
     * Construct a FibreTracker.
     * @param image the image within which the tracking will take place.
     *
     */
    protected FibreTracker(TractographyImage image) {
        
	xDataDim = image.xDataDim();
	yDataDim = image.yDataDim();
	zDataDim = image.zDataDim();

	xVoxelDim = image.xVoxelDim();
	yVoxelDim = image.yVoxelDim();
	zVoxelDim = image.zVoxelDim();

	// Upper bounds
	xUBound = xDataDim * xVoxelDim;
	yUBound = yDataDim * yVoxelDim;
	zUBound = zDataDim * zVoxelDim;

	visitedVoxel = new short[xDataDim][yDataDim][zDataDim];

	isotropic = image.getIsotropicMask();

	this.image = image;

        zeros = new short[zDataDim];

    }


    
    /**
     * Set how far to track between curvature checks, in mm. At the specified intervals, the current fiber orientation
     * will be compared to the one from the last interval, and tracking will terminate if the angle between them exceeds
     * some threshold. Since we compare only the end points, it is not a good idea to set this value too large. The default
     * is 10mm. The value of this impacts how you should set your curvature threshold.
     *
     * @param length a distance to track before checking curvature, in mm. The curvature is checked at the first
     * step that where the path length to the last curvature check exceeds this value.
     */
    public void setCurveCheckInterval(double length) {
        checkCurveInterval = length;
    }


    /**
     * Set an upper bound on the length of streamlines. Use this to prevent erroneous lines from growing indefinitely, or to limit the 
     * length to a particular value. The default value is 1000 mm.
     * 
     * @param maxLengthMM the maximum length of each fibre segment, The total tract length may be twice this parameter since 
     * we may track up to <code>maxLengthMM</code> in both directions from the seed.
     *

     */
    public void setMaxSegmentLength(double maxLengthMM) {
        maxSegmentLength = maxLengthMM;
    }


    /**
     * Sets the minimum cosine between two vectors describing the current fibre orientation at two points whose separation is 
     * defined by calling <code>setCheckCurveDistance</code>. If the curvature is evaluated every 1mm, then the value passed to this
     * method controls the maximum curvature per mm. 
     *
     * @param ipThresh the threshold value, in the range of -1 (no threshold) to 1 (complete collinearity). A more realistic value would 
     * be 1.0 / Math.sqrt(2.0), which would allow a curvature of up to 45 degrees.
     *
     */
    public void setIP_Threshold(double ipThresh) {
        ipThreshold = ipThresh;
    }

 
    /** 
     * Track paths from seed points placed at the centre of all voxels within the ROI. 
     * @param roi Voxel region within which to sow seeds. Checks bounds and truncates ROI
     * if necessary.
     * @return a <code>TractCollection</code> containing the results of tracking from all seed 
     * points in the ROI.
     * @see tractography.TractCollection
     */
    public final TractCollection trackPaths(RegionOfInterest roi) {

	Point3D[] points = roi.getSeedPoints();

	TractCollection paths = new TractCollection(points.length + 1, 100.0);

	for (int i = 0; i < points.length; i++) {
	    // don't call trackFromSeed(Point, pdIndex, direction)
	    // directly or you will mess up probabilistic methods
	    paths.addTractCollection(trackFromSeed(points[i]));
	    
	}
  

	return paths;
	
  
    }


  
    /** 
     * Track paths from a single seed point within the ROI. 
     *
     * @param point the point in mm to track from
     * @return a <code>TractCollection</code> containing the results of tracking from this 
     * seed point. For single fibre trackers, this <code>TractCollection</code> will contain 
     * a single tract; multi-fibre trackers may return more than one tract.  
     * 
     */
    public final TractCollection trackFromSeed(Point3D seedPoint) {

	TractCollection collection = new TractCollection(4, 100.0);	

	if(!inBounds(seedPoint)) {
	    Tract t = new Tract(2, 100.0);
	    t.addPoint(seedPoint);
	    collection.addTract(t);
	    return collection;
	}
	
	int i = (int)(seedPoint.x / xVoxelDim);
	int j = (int)(seedPoint.y / yVoxelDim);
	int k = (int)(seedPoint.z / zVoxelDim);

	int numPDs = image.numberOfPDs(i,j,k);

	for (int p = 0; p < numPDs; p++) {
	    
	    Tract t1 = trackFromSeed(seedPoint, p);

	    collection.addTract(t1);

	}

	return collection;

    }



    /** 
     * Track paths from a single seed point within the ROI, using a particular principal direction as the initial direction.
     *
     * @param point the point in mm to track from.
     * @param pdIndex the principal direction to follow.
     * 
     */
    public final Tract trackFromSeed(Point3D seedPoint, int pdIndex) {

	if(!inBounds(seedPoint)) {
	    Tract t = new Tract(2, 100.0);
	    t.addPoint(seedPoint);
	    return t;
	}

        Tract t1 = trackFromSeed(seedPoint, pdIndex, FORWARD);
        Tract t2 = trackFromSeed(seedPoint, pdIndex, BACKWARD);
            
        t1.joinTract(t2);
            
        for (int x = 0; x < xDataDim; x++) {
            for (int y = 0; y < yDataDim; y++) {
                System.arraycopy(zeros, 0, visitedVoxel[x][y], 0, zDataDim);
            }
        }

        image.resetRandomization();
        
 	return t1;

    }

    
    /** 
     * Track paths from a single seed point within the ROI. This method should only be called from trackFromSeed(Point3D),
     * because that also resets the randomization and stuff appropriately after tracking.
     * @param point the point in mm to track from.
     * @param pdIndex the principal direction to follow.
     * @param direction if {@link tractography.FibreTracker#FORWARD FORWARD}, start tracking along 
     * the PD, otherwise track in the opposite direction.
     *
     * @return a <code>Tract</code> containing the results of tracking. 
     *
     * @see #trackFromSeed(numerics.Point3D)
     * 
     */
    private final Tract trackFromSeed(Point3D seedPoint, int pdIndex, boolean direction) {

        // The next step displacement
        Vector3D step = null;

        short pointsAdded = 1;

        // unit vector of current fibre orientation, equals step.normalized()
        Vector3D currentBearing = null;

        // Previous step orientation
        Vector3D previousBearing = null;

        // A previous step orientation from some time ago, that we use to assess curvature
        // We stop tracking if the fibre has curved by more than some threshold amount
	Vector3D checkCurveBearing = null;


        // total tract length
        double tractLength = 0.0;

        // length since last curve check
        double distanceSinceCurveCheck = 0.0;

	Tract path = new Tract(400, 100.0);
	
	Point3D currentPos = seedPoint;

	path.addPoint(seedPoint, 0.0);

	if (!inBounds(seedPoint)) {
	    return path;
	}

	int i = (int)(currentPos.x / xVoxelDim);
	int j = (int)(currentPos.y / yVoxelDim);
	int k = (int)(currentPos.z / zVoxelDim);
        
	if (isotropic[i][j][k]) {
	    return path;
	}
        
	step = getFirstStep(currentPos, pdIndex, direction);
        
        // Note that this may vary with each step. Subclasses are not constrained to a fixed step size
        double stepSize = step.mod();

        double prevStepSize = stepSize;

        currentBearing = step.scaled(1.0 / stepSize);
        
        checkCurveBearing = currentBearing;
        
	// Now track the rest of the way, always moving in direction that has the largest dot product
	// with the vector currentPos - previousPos

	visitedVoxel[i][j][k] = pointsAdded;
	
	while (true) {

	    currentPos = currentPos.displace( step );   

	    // Quit if we are out of bounds
	    if (!inBounds(currentPos)) {
		return path;
	    }

	    i = (int)(currentPos.x / xVoxelDim);
	    j = (int)(currentPos.y / yVoxelDim);
	    k = (int)(currentPos.z / zVoxelDim);
	    
	    if (isotropic[i][j][k]) {
		return path;
	    }

            if (visitedVoxel[i][j][k] > 0 && pointsAdded - visitedVoxel[i][j][k] > LOOP_ALLOWANCE) {
                return path;
            }
            else {
                visitedVoxel[i][j][k] = pointsAdded;
            }
            
            tractLength += stepSize;
            
            if (tractLength > maxSegmentLength) {
                return path;
            }

            distanceSinceCurveCheck += stepSize;

            // Now get the next step and check curvature if necessary - still haven't added point yet

            step = getNextStep(currentPos, currentBearing);
            
            stepSize = step.mod();
    
            currentBearing = step.scaled(1.0 / stepSize);

            
	    if ( distanceSinceCurveCheck >= checkCurveInterval) {
		// no Math.abs here, since we want to know if we're going backwards
		if ( checkCurveBearing.dot(currentBearing) < ipThreshold) {
		    return path;
		}
		checkCurveBearing = currentBearing;
                distanceSinceCurveCheck = 0.0;
	    }


            // Passed curvature test, record position and displacement since previous step
	    path.addPoint(currentPos, prevStepSize);
            
            prevStepSize = stepSize;

	    pointsAdded++;

	    if (pointsAdded == MAX_TRACT_POINTS) {
                return path;
	    }

          
	}
        
        
    }


    /**
     * Get the initial step by following a particular principal direction at the seed point.
     */
    protected abstract Vector3D getFirstStep(Point3D currentPos, int pdIndex, boolean direction);


    /**
     * Get the next step given the previous direction.
     */
    protected abstract Vector3D getNextStep(Point3D currentPos, Vector3D previousBearing);


    /** 
     * @param position a position to test.
     * @return <code>false</code> if the position is outside the dataset, 
     * <code>true</code> otherwise.
     */
    public final boolean inBounds(Point3D point) {
	double xPos = point.x;
	double yPos = point.y;
	double zPos = point.z;
	return (xPos >= 0.0 && yPos >= 0.0 && zPos >= 0.0 
		&& xPos < xUBound && yPos < yUBound
		&& zPos < zUBound);
    }


  
    public double getIP_Threshold() {
	return ipThreshold;
    }



}







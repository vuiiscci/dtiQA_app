package tractography;

import imaging.*;
import misc.*;
import numerics.*;

import java.io.*;


/**
 * <dl>
 * <dt>Purpose: To contain a list of points describing one streamline.
 * <BR><BR>
 * <dt>Description:
 * <dd> Objects of this class represent one tract. They hold a list of 3D points in mm space. 
 * <code>Tract</code>'s perform bounds checking as points are added, and increase their 
 * capacity when full.
 *
 * </dl>
 *
 * @version $Id$
 * @author  Philip Cook
 * @see tractography.TractCollection
 * 
 */
public final class Tract {

    /** Number of points Tract can hold. */
    private int capacity; 
    
    /** Percentage to increase capacity by when Tract is full. */
    private double growBy; 

    /** Number of points currently in the Tract. */
    private int numberOfPoints;

    /** Holds the point data. */
    private Point3D[] points; 

    /** Displacement of each point from the last. */
    private double[] displacements;

    /** Path length in each voxel. */
    private double[] voxelPathLengths;

    /** Which point is the seed. */
    private int seedPointIndex;
    
    /** Which voxel of the voxelList is the seed. */
    private int seedPointVoxelIndex;

    /** List of voxels in this streamline. Lazy initialized. */
    private VoxelList voxels = null;
    
    private boolean displacementsCalculated = false;

    private boolean voxelPathLengthsCalculated = false;





    /** 
     * Construct a Tract with a specified capacity and growth factor.
     * @param initialCapacity the initial number of 3D points the Tract can hold.
     * @param growFactor the percentage to increase capacity by when this Tract is full.
     */
    public Tract(int initialCapacity, double growFactor) {
	numberOfPoints = 0;
	seedPointIndex = 0;
	capacity = initialCapacity;
	growBy = growFactor;
	initialise();
    }


    /** Copy constructor
     *
     */
    public Tract(Tract t) {
	numberOfPoints = t.numberOfPoints;
	capacity = t.capacity + 1;
	growBy = t.growBy;
	seedPointIndex = t.seedPointIndex;
	initialise();

	displacementsCalculated = t.displacementsCalculated;
	
        System.arraycopy(t.points, 0, points, 0, numberOfPoints);
        System.arraycopy(t.displacements, 0, displacements, 0, numberOfPoints);
    }


    /** Construct a Tract with a capacity of 100 and a growth factor of 100%. */
    public Tract() {
	numberOfPoints = 0;
	seedPointIndex = 0;
	capacity = 100;
	growBy = 100.0;
	initialise();
    }
   

    /** Allocate memory for points. */
    private void initialise() {
	points = new Point3D[capacity];
        displacements = new double[capacity];
    }
    

    /** Increase size of array when full. */
    private void growArray() {
    
	int newCapacity = (int)((double)capacity + (double)capacity * growBy / 100.0);
	if (newCapacity == capacity) {
	    newCapacity = capacity + 1;
	}
	
	Point3D[] newArray = new Point3D[newCapacity];
	
        double[] newDisplacements = new double[newCapacity];

        System.arraycopy(points, 0, newArray, 0, numberOfPoints);
        System.arraycopy(displacements, 0, newDisplacements, 0, numberOfPoints);
	

	capacity = newCapacity;
        displacements = newDisplacements;
	points = newArray;
    }
    

    /** 
     * Join another <code>Tract</code> to this one. This method assumes the two 
     * <code>Tract's</code> begin at the same point. It is used to join <code>Tract's</code> 
     * tracked in opposite directions from the same seed point. 
     * @param t the Tract to be joined.
     */
    protected void joinTract(Tract t) {

	if (t.numberOfPoints == 0) {
	    return;
	}

	voxels = null;

	int tPoints = t.numberOfPoints;
	int newCapacity = numberOfPoints + tPoints;
	Point3D[] newArray = new Point3D[newCapacity];

        double[] newDisplacements = new double[newCapacity];
	
	int i;
	
	// Reverse order of tract to be joined on left
	for (i = 0; i < (tPoints - 1); i++) {
	    newArray[i] = t.points[(tPoints - 1) - i];
            newDisplacements[i] = t.displacements[(tPoints - 1) - i];
	}

	seedPointIndex = (tPoints - 1);

	// Note we have not included seed point from t.points
	// This ensures the seed point is not repeated
	for (i = 0; i < numberOfPoints; i++) {
	    newArray[i + (tPoints - 1)] = points[i];
            newDisplacements[i + (tPoints - 1)] = displacements[i];
	}
	
	// Set number of points of joined tract
	numberOfPoints += (tPoints - 1);
	capacity = newCapacity;
	points = newArray;
        displacements = newDisplacements;

	displacementsCalculated = displacementsCalculated && t.displacementsCalculated;
	
    }


    
    
    /**
     * Add a point.
     * @param point the point to add.
     * 
     */
    public void addPoint(Point3D point) {

	voxels = null;
	voxelPathLengthsCalculated = false;

	displacementsCalculated = false;

	points[numberOfPoints] = point;

	numberOfPoints++;

	if (numberOfPoints == capacity) {
	    growArray();
	}

    }

    /**
     * Adds a point.
     * @param point the point to add.
     * @param displacement the displacement of this point from the last.
     * 
     */
    public void addPoint(Point3D point, double displacement) {

	voxels = null;
	voxelPathLengthsCalculated = false;

	points[numberOfPoints] = point;
        displacements[numberOfPoints] = displacement;

	numberOfPoints++;

	if (numberOfPoints == capacity) {
	    growArray();
	}

    }


    
    private void calculateDisplacements() {

	// nothing to do
	if (numberOfPoints < 2) {
	    return;
	}
	
	int sp = seedPointIndex;

	for (int i = sp - 1; i >= 0; i--) {

	    Point3D point = points[i];
	    Point3D last = points[i+1];

	    displacements[i] = Math.sqrt( (point.x - last.x) * (point.x - last.x) + 
					  (point.y - last.y) * (point.y - last.y) + 
					  (point.z - last.z) * (point.z - last.z) );
	}
	for (int i = sp + 1; i < numberOfPoints; i++) {

	    Point3D point = points[i];
	    Point3D last = points[i-1];

	    displacements[i] = Math.sqrt( (point.x - last.x) * (point.x - last.x) + 
					  (point.y - last.y) * (point.y - last.y) + 
					  (point.z - last.z) * (point.z - last.z) );
	}

        displacementsCalculated = true;
	
    }



    /**
     * Gets the length of the streamline from the seed to point <code>index</code>.
     *
     */
    public double pathLengthFromSeed(int index) {

	if (!displacementsCalculated) {
	    calculateDisplacements();
	}

        if (index == seedPointIndex) {
            return 0.0;
        }
        else if (index < 0 || index >= numberOfPoints) {
            throw new IllegalArgumentException("Illegal index " + index + 
                                               " is < 0 or exceeds index of last point");
        }
        else if (index < seedPointIndex) {
            
            double distance = 0.0;

            for (int i = seedPointIndex - 1; i >= index; i--) {
                distance += displacements[i];
            }

            return distance;

        }
        else { // index > seedPointIndex
            
            double distance = 0.0;

            for (int i = seedPointIndex + 1; i <= index; i++) {
                distance += displacements[i];
            }

            return distance;

        }
    }

    
  
    /**
     * @return the total length of the streamline.
     *
     */
    public double length() {

	if (!displacementsCalculated) {
	    calculateDisplacements();
	}

	double length = 0.0;
	
	for (int i = 0; i < numberOfPoints; i++) {
	    length += displacements[i];
	}
	
	return length;
    }


    /**
     *
     * @return the distance between the end points, which is equivalent to length() if the tract is a straight line.
     *
     */
    public double endpointSeparation() {
        return points[0].distance(points[numberOfPoints - 1]);
    }

    
    /** @return the number of points in the Tract. */
    public int numberOfPoints() {
	return numberOfPoints;
    }
    

    /**
     * List of voxel coordinates of the tract. If the tract has multiple steps in the same voxel, 
     * the redundancies will be removed, though no check is made for looping fibres 
     * (these should be prevented by the tracker) or fibres that cross back and forth along voxel 
     * boundaries.
     * 
     * @return List of voxels.
     */
    public VoxelList toVoxelList(double[] voxelDims) {
	return toVoxelList(voxelDims[0], voxelDims[1], voxelDims[2]);
    }
    
    /**
     * List of voxel coordinates of the tract. If the tract has multiple steps in the same voxel, 
     * the redundancies will be removed, though no check is made for looping fibres 
     * (these should be prevented by the tracker) or fibres that cross back and forth along voxel 
     * boundaries.
     * <p>
     * This method only returns meaningful results if the Tract is in Camino space.
     * 
     * @return List of voxels.
     */
    public VoxelList toVoxelList(double xVoxelDim, double yVoxelDim, double zVoxelDim) {

	if (voxels != null) {
	    
	    double[] voxelDims = voxels.getVoxelDims();

	    if (xVoxelDim == voxelDims[0] && yVoxelDim == voxelDims[1] && zVoxelDim == voxelDims[2]) {
		return voxels;
	    }
	    else {
		voxels = null;
		voxelPathLengthsCalculated = false;
	    }
	}

	Voxel[] tmp = new Voxel[numberOfPoints];

        int x = (int)(points[0].x / xVoxelDim);
        int y = (int)(points[0].y / yVoxelDim);
        int z = (int)(points[0].z / zVoxelDim);

        tmp[0] = new Voxel(x,y,z);

	int voxelsAdded = 1;	

	int seedPointVoxelIndex = 0;

	for (int i = 1; i < numberOfPoints; i++) {

	    x = (int)(points[i].x / xVoxelDim);
	    y = (int)(points[i].y / yVoxelDim);
	    z = (int)(points[i].z / zVoxelDim);

	    if (x != tmp[voxelsAdded-1].x || y != tmp[voxelsAdded-1].y || z != tmp[voxelsAdded-1].z) {

		tmp[voxelsAdded] = new Voxel(x,y,z);
		voxelsAdded++;
	    } 

	    if (i == seedPointIndex) {
		seedPointVoxelIndex = voxelsAdded - 1;
	    }

	    
	}

	Voxel[] voxelList = new Voxel[voxelsAdded];
	
        System.arraycopy(tmp, 0, voxelList, 0, voxelsAdded);

	Vector3D tangent = null;
	Vector3D negTangent = null;

	if (seedPointIndex < numberOfPoints - 1) {
	    tangent = new Vector3D(points[seedPointIndex + 1], points[seedPointIndex]);
	}
	else {
	    if (numberOfPoints > 1) {
		tangent = new Vector3D(points[seedPointIndex], points[seedPointIndex - 1]);
	    }
	    else {
		// one point, write a default
		tangent = new Vector3D(1.0, 0.0, 0.0);
	    }
	}

	if (seedPointIndex > 0) {
	    negTangent = new Vector3D(points[seedPointIndex - 1], points[seedPointIndex]);
	}
	else {
	    negTangent = tangent.negated();
	}

	VoxelList tractVoxels = new VoxelList(voxelList, seedPointVoxelIndex, xVoxelDim, yVoxelDim, zVoxelDim,
					      tangent, negTangent);

	voxels = tractVoxels;

	return voxels;

    }




    /** 
     * @return tab separated list of (x,y,z) points in this tract.
     */
    public String toString() {

	StringBuffer buff = new StringBuffer();

	buff.append("tractography.Tract\n" + numberOfPoints + " points:\nseed point index:" + 
		    seedPointIndex + "\n");

	for (int i = 0; i < numberOfPoints; i++) {
	    buff.append(points[i].x);
	    buff.append("\t");
	    buff.append(points[i].y);
	    buff.append("\t");
	    buff.append(points[i].z);
	    buff.append("\n");
	}

	return buff.toString();
    }

    /** 
     * Get a point.
     * @param pointNum the point index, 0 to (numberOfPoints-1).
     * @return the point <code>pointNum</code>.
     */
    public Point3D getPoint(int pointNum) {
	
	if (pointNum >= 0 && pointNum < numberOfPoints) {
	    // Point3D is immutable so no need for defensive copy
	    return points[pointNum];
	}
	else {
	    throw new IndexOutOfBoundsException("Point " + pointNum + " does not exist");
	}
    }


    /**
     * Get the seed point
     */
    public Point3D getSeedPoint() {
        return points[seedPointIndex];
    }

    /** 
     * @return all points.
     *
     */
    public Point3D[] getPoints() {
	
	Point3D[] defCopy = new Point3D[numberOfPoints];

	for (int i = 0; i < numberOfPoints; i++) {
	    // safe to do this, points are immutable
	    defCopy[i] = points[i];
	}

	return defCopy;
    }


    /**
     * @return the index of the seed point.
     */
    public int seedPointIndex() {
	return seedPointIndex;
    }


    /**
     * Resample the points to a minimum resolution, such that the distance between
     * points is no more than <code>stepSize</code>. 
     *
     * @return a copy of this tract, with points resampled to the minimum resolution.
     */
    public Tract resample(double stepSize) {

	if (!(stepSize > 0.0)) {
	    throw new LoggedException("Can't use step size of " + stepSize);
	}
	
	Tract resampled = new Tract(2 * numberOfPoints, 100.0);

	resampled.addPoint(points[0]);
		
	for (int i = 1; i < numberOfPoints; i++) {
	    Vector3D segment = new Vector3D(points[i], points[i-1]);

	    double mod = segment.mod();

	    int steps = (int)(mod / stepSize);

	    for (int s = 0; s < steps; s++) {
		resampled.addPoint( 
				   points[i-1].displace( segment.scaled((s+1) * stepSize / mod) ) 
				   );
	    }
	    
	    // unless we put a point very close to the end, place the end point in the
	    // resampled streamline
	    if (mod - steps * stepSize > 1E-3 || steps == 0) {
		resampled.addPoint(points[i]);
	    }

	    if (i == seedPointIndex) {
		resampled.seedPointIndex = resampled.numberOfPoints - 1;
	    }
	}

	return resampled;

    }


    /**
     * Reduce the length of a tract to the specified maximum number of points.
     *
     * @param maxPoints must be >= 1. If it is greater than the number of points
     * in this tract, nothing is done.
     */
    public void truncateToMaxPoints(int maxPoints) {
	
	if (!(numberOfPoints > maxPoints)) {
	    return;
	}

        if (maxPoints == 0) {
            // empty this tract
            numberOfPoints = 0;
            return;
        }

	if (!(maxPoints >= 0)) {
	    throw new IllegalArgumentException("can't truncate tract to " + maxPoints + " points");
	}

	int maxSegmentLength = maxPoints / 2;

	int minPointIndex = seedPointIndex - maxSegmentLength;
	int maxPointIndex = seedPointIndex + maxSegmentLength;

	if (maxPoints % 2 == 0) {
	    // need one point less if maxPoints is even, since 
	    // (2 * maxSegmentLength) + (seed) == maxPoints + 1 in that case
	    // if maxPoints is odd, (2 * maxSegmentLength) + (seed) == maxPoints
	    minPointIndex += 1;
	}

	// initial segment to chop out is [minPointIndex, maxPointIndex], which is centered 
	// on the seed point (with an extra point added to the min if needed)
	// if this range does not overlap the tract, we shift it

	if (minPointIndex < 0) {
	    maxPointIndex += -minPointIndex;
	    
	    minPointIndex = 0;
	}
	
	if (maxPointIndex > numberOfPoints - 1) {
	    minPointIndex -= maxPointIndex - numberOfPoints + 1;
	    maxPointIndex = numberOfPoints - 1;
	}

	

	if (maxPointIndex - minPointIndex + 1 != maxPoints) {
	    // should never get here
	    throw new LoggedException("Could not truncate tract to " + maxPoints + 
				      " points. Indices are " + minPointIndex + " " + 
				      maxPointIndex + ". Tract is\n" + toString());
	}
	
	chop(minPointIndex, maxPointIndex);

    }


    /**
     * Reduce the length of a Tract to the specified maximum length in mm. 
     *
     */
    public void truncateToMaxLength(double maxLength) {

	if (!(maxLength > 0.0)) {
	    throw new IllegalArgumentException("can't truncate tract to " + maxLength + " mm");
	}
	if (!(length() > maxLength)) {
	    return;
	}
	
	int minIndex = seedPointIndex;
	int maxIndex = seedPointIndex;

	double choppedLength = 0.0;

	double lowerSegmentLength = 0.0;
	double upperSegmentLength = 0.0;

	if (minIndex > 0 && choppedLength + displacements[minIndex - 1] < maxLength) {
	    minIndex--;
	    choppedLength += displacements[minIndex];
	    lowerSegmentLength += displacements[minIndex];
	}
	if ((maxIndex < numberOfPoints - 1) && (choppedLength + displacements[maxIndex + 1] < maxLength)) {
	    maxIndex++;
	    choppedLength += displacements[maxIndex];
	    upperSegmentLength += displacements[maxIndex];
	}

	boolean increasedLength = true;

	while (choppedLength < maxLength && increasedLength) {

	    double chopBeforeExpand = choppedLength;

	    if (lowerSegmentLength < upperSegmentLength) {
		
		// try to do upper segment first
		
		if ((maxIndex < numberOfPoints - 1) && (choppedLength + displacements[maxIndex + 1] < maxLength)) {
		    maxIndex++;
		    choppedLength += displacements[maxIndex];
		    upperSegmentLength += displacements[maxIndex];
		}
		if (minIndex > 0 && choppedLength + displacements[minIndex - 1] < maxLength) {
		    minIndex--;
		    choppedLength += displacements[minIndex];
		    lowerSegmentLength += displacements[minIndex];
		}
		
	    }
	    else {

		if (minIndex > 0 && choppedLength + displacements[minIndex - 1] < maxLength) {
		    minIndex--;
		    choppedLength += displacements[minIndex];
		    lowerSegmentLength += displacements[minIndex];
		}
		if ((maxIndex < numberOfPoints - 1) && (choppedLength + displacements[maxIndex + 1] < maxLength)) {
		    maxIndex++;
		    choppedLength += displacements[maxIndex];
		    upperSegmentLength += displacements[maxIndex];
		}
		
	    }

	    increasedLength = choppedLength > chopBeforeExpand;
	}

	chop(minIndex, maxIndex);
	
    }


    /**
     * Chop off one or both ends of a tract.
     *  
     * @param minPointIndex an index between 0 and seedPointIndex().
     * @param maxPointIndex an index between seedPointIndex() and numberOfPoints() - 1.
     */
    public void chop(int minPointIndex, int maxPointIndex) {
	if (!(minPointIndex <= seedPointIndex && maxPointIndex < numberOfPoints 
	      && maxPointIndex >= seedPointIndex)) {
	    throw new LoggedException("invalid indices " + minPointIndex + " " + maxPointIndex);
	}


	// got to null out voxels if this method is void
	voxels = null;
	voxelPathLengthsCalculated = false;

	// number of points in chopped streamline
	int choppedPoints = maxPointIndex - minPointIndex + 1;

	for (int i = 0; i < choppedPoints; i++) {
	    points[i] = points[i + minPointIndex];
	    displacements[i] = displacements[i + minPointIndex];
	}

 	for (int i = choppedPoints; i < numberOfPoints; i++) {
 	    points[i] = null;
 	    displacements[i] = 0.0;
 	}

	numberOfPoints = choppedPoints;

	seedPointIndex = seedPointIndex - minPointIndex;

    }

    

    /**
     * Applies a transformation matrix directly.
     *
     */
    public void transform(RealMatrix trans) { 
	
	voxels = null;
	voxelPathLengthsCalculated = false;
	displacementsCalculated = false;

	Point3D[] newPoints = new Point3D[numberOfPoints + 1];
        
	for (int p = 0; p < numberOfPoints; p++) {
	    newPoints[p] = points[p].transform(trans);
	}
        
	points = newPoints;
    }


    /**
     * Writes tract so that it can be read later by a TractSource.
     *
     */
    public void writeRaw(DataOutputStream dout) throws IOException {
		
	dout.writeFloat((float)numberOfPoints);
	dout.writeFloat((float)seedPointIndex);

	for (int p = 0; p < numberOfPoints; p++) {
	    dout.writeFloat((float)points[p].x);
	    dout.writeFloat((float)points[p].y);
	    dout.writeFloat((float)points[p].z);
	}
		
    }


    /**
     * Flattens tract into an array. This method is to facilitate output from 
     * data.OutputManager. Note that all raw tracts should be written as floats.
     *
     * @return {numberOfPoints, seedPointIndex, p0.x, p0.y, p0.z, p1.x...}.
     */
    public double[] toArray() {

	double[] array = new double[2 + 3 * numberOfPoints];
		
	array[0] = (double)numberOfPoints;
	array[1] = (double)seedPointIndex;

	for (int p = 0; p < numberOfPoints; p++) {
	    
	    array[2 + p * 3] = points[p].x;
	    array[3 + p * 3] = points[p].y;
	    array[4 + p * 3] = points[p].z;
	}

	return array;
		
    }


    public boolean equals(Object o) {
	
	if (!(o instanceof Tract)) {
	    return false;
	}
	if (o == null) {
	    return false;
	}
	
	Tract t = (Tract)o;

	if (t == this) {
	    return true;
	}

	if (t.numberOfPoints != numberOfPoints) {
	    return false;
	}

	if (t.seedPointIndex != seedPointIndex) {
	    return false;
	}

	for (int i = 0; i < numberOfPoints; i++) {
	    if (!points[i].equals(t.points[i])) {
		return false;
	    }
	}

	return true;

    }

    
    public int hashCode() {
	return 13 * (int)(points[seedPointIndex].x) + 29 * numberOfPoints + 37 * seedPointIndex;
    }



    /**
     * Divides the streamline into line segments. Each segment connects two voxel boundaries. The
     * method computes the length of the line segment in each voxel.
     *
     * @return a list of path lengths that corresponds to the ordering of voxels in <code>t.toVoxelList()</code>.
     */
    public double[] getVoxelPathLengths(double[] voxelDims) {
	return getVoxelPathLengths(voxelDims[0], voxelDims[1], voxelDims[2]);
    }



    /**
     * Divides the streamline into line segments. Each segment connects two voxel boundaries. The
     * method computes the length of the line segment in each voxel.
     *
     * @return a list of path lengths that corresponds to the ordering of voxels in <code>t.toVoxelList()</code>.
     */
    public double[] getVoxelPathLengths(double xVoxelDim, double yVoxelDim, double zVoxelDim) {

	// check for cached value
	// any time voxels is set to null, the voxel path lengths need re-calculating
	if (voxelPathLengthsCalculated && voxels != null) {
	    
	    double[] voxelDims = voxels.getVoxelDims();
	    
	    if (xVoxelDim == voxelDims[0] && yVoxelDim == voxelDims[1] && zVoxelDim == voxelDims[2]) {
		// can use cached version
		double[] defCopy = new double[voxelPathLengths.length];

		System.arraycopy(voxelPathLengths, 0, defCopy, 0, voxelPathLengths.length);
		
		return defCopy;
	    }
	    else {
		toVoxelList(xVoxelDim, yVoxelDim, zVoxelDim);
	    }
	    
	}
	else if (voxels != null) {
	    
	    double[] voxelDims = voxels.getVoxelDims();
	    
	    if (xVoxelDim == voxelDims[0] && yVoxelDim == voxelDims[1] && zVoxelDim == voxelDims[2]) {
		// can use voxels
	    }
	    else {
		toVoxelList(xVoxelDim, yVoxelDim, zVoxelDim);
	    }
	}
	else {
	    toVoxelList(xVoxelDim, yVoxelDim, zVoxelDim);

	}
	

	double[] segmentLengths = new double[voxels.size()];

	int voxelCounter = 0;
	
	Voxel firstVox = voxels.getVoxel(0);
	
	// voxel indices
	int x = firstVox.x;
	int y = firstVox.y;
	int z = firstVox.z;

	for (int i = 0; i < numberOfPoints - 1; i++) {

	    int xp = (int)(points[i+1].x / xVoxelDim);
	    int yp = (int)(points[i+1].y / yVoxelDim);
	    int zp = (int)(points[i+1].z / zVoxelDim);


	    // displacement vector from point i to point i+1
	    Vector3D displacement = new Vector3D(points[i+1], points[i]);

	    double[] planeIntersections = new double[6];

	    double mod = displacement.mod();

	    if (x != xp || y != yp || z != zp) {
		// changed voxel
		
		displacement = displacement.scaled(1.0 / mod);

		double rightPlane = xVoxelDim * (1 + x);
		double leftPlane = xVoxelDim * x;

		double frontPlane = yVoxelDim * (1 + y);
		double backPlane = yVoxelDim * y;

		double topPlane = zVoxelDim * (1 + z);
		double bottomPlane = zVoxelDim * z;


		if (displacement.x != 0.0 ) {
		    planeIntersections[0] = (leftPlane - points[i].x) / displacement.x;
		    planeIntersections[1] = (rightPlane - points[i].x) / displacement.x;
		}
		else {
		    planeIntersections[0] = Double.MAX_VALUE;
		    planeIntersections[1] = Double.MAX_VALUE;
		}
		
		if (displacement.y != 0.0) {
		    
		    planeIntersections[2] = (frontPlane - points[i].y) / displacement.y;
		    planeIntersections[3] = (backPlane - points[i].y) / displacement.y;
		    
		}
		else {
		    // meets at infinity...
		    planeIntersections[2] = Double.MAX_VALUE;
		    planeIntersections[3] = Double.MAX_VALUE;
		}		

		if (displacement.z != 0.0) {
		    planeIntersections[4] = (topPlane - points[i].z) / displacement.z;
		    planeIntersections[5] = (bottomPlane - points[i].z) / displacement.z;
		    
		}
		else {
		    planeIntersections[4] = Double.MAX_VALUE;
		    planeIntersections[5] = Double.MAX_VALUE;
		}		

		
		// now sort
		java.util.Arrays.sort(planeIntersections);
		
		// scale by smallest positive element
		// allow value of 0.0; occurs only if point sits precisely on a voxel boundary
		int planeCounter = 0;
		while (planeIntersections[planeCounter] < 0.0) {
		    planeCounter++;
		    
		}

		segmentLengths[voxelCounter] += planeIntersections[planeCounter];
		
		// There could be a small intersection with a third voxel in between these two. The tract
		// resolution determines the potential magnitude of this error. This method assumes that 
		// the tract has been sampled / resampled to an acceptable resolution.

		double remainder = mod - planeIntersections[planeCounter];

		segmentLengths[voxelCounter+1] += remainder > 0.0 ? remainder : 0.0;

		voxelCounter++;

		x = xp;
		y = yp;
		z = zp;
		
	    }
	    else {
		
		segmentLengths[voxelCounter] += mod;
	    }


	}

	
	voxelPathLengths = segmentLengths;
	voxelPathLengthsCalculated = true;

	double[] defCopy = new double[voxelPathLengths.length];
	
	System.arraycopy(voxelPathLengths, 0, defCopy, 0, voxelPathLengths.length);
	
	return defCopy;

    }


    /**
     * Does a point-by-point transformation from Camino space to physical space 
     * as defined by the transformation.
     *
     */
    public final void transformToPhysicalSpace(RealMatrix trans, double xVoxelDim, double yVoxelDim, double zVoxelDim) {

        // clear all cached information
	voxels = null;
	voxelPathLengthsCalculated = false;
	displacementsCalculated = false;

        for (int i = 0; i < numberOfPoints; i++) {

            Point3D voxel = new Point3D( points[i].x / xVoxelDim - 0.5, points[i].y / yVoxelDim - 0.5, 
                                         points[i].z / zVoxelDim - 0.5 );

            points[i] = voxel.transform(trans);
            
        }
        
    }
    

    /**
     * Does a point-by-point transformation from Camino space to physical space 
     * as defined by the image header.
     *
     */
    public final void transformToCaminoSpace(RealMatrix trans, double xVoxelDim, double yVoxelDim, double zVoxelDim) {


        // clear all cached information
	voxels = null;
	voxelPathLengthsCalculated = false;
	displacementsCalculated = false;


        for (int i = 0; i < numberOfPoints; i++) {

            Point3D voxel = points[i].transform(trans);

            // Back to Camino space, remembering the half voxel translation

             points[i] = new Point3D( (voxel.x + 0.5) * xVoxelDim, (voxel.y + 0.5) * yVoxelDim, 
                                      (voxel.z + 0.5) * zVoxelDim );
            
        }


    }

   

}

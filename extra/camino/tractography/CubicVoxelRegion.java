package tractography;

import numerics.Point3D;

/**
 * <dl>
 * <dt>Purpose: to encapsulate a cuboid region.
 * <dt> Description
 * <dd> This class is useful for describing a region of an image in voxel coordinates. It operates entirely
 * within Camino space, without reference to physical space. Use for testing only
 *
 * </dl>
 *
 * @author Philip Cook
 * @version $Id$
 */
final class CubicVoxelRegion implements RegionOfInterest {

    // Region in Voxels
    protected final int iminX;
    protected final int iminY;
    protected final int iminZ;
    protected final int imaxX; 
    protected final int imaxY; 
    protected final int imaxZ;

    // And also the bounds in MM coordinates
    protected final double minXMM;
    protected final double minYMM;
    protected final double minZMM;
    protected final double maxXMM;
    protected final double maxYMM;
    protected final double maxZMM;

    private final double xVoxelDim;
    private final double yVoxelDim;
    private final double zVoxelDim;


    private final int label;

    
    /**
     * Create a region of size 1 positioned at (0,0,0).
     * Voxel scale is 1.0 mm cubed
     */
    private CubicVoxelRegion() {
        
        label = 1;

	iminX = 0;
	iminY = 0;
	iminZ = 0;
	imaxX = 0;
	imaxY = 0;
	imaxZ = 0;

	xVoxelDim = 1.0;
	yVoxelDim = 1.0;
	zVoxelDim = 1.0;

	minXMM = (double)iminX * xVoxelDim;
	minYMM = (double)iminY * yVoxelDim;
	minZMM = (double)iminZ * zVoxelDim;
	maxXMM = (double)(imaxX + 1) * xVoxelDim; 
	maxYMM = (double)(imaxY + 1) * yVoxelDim;
	maxZMM = (double)(imaxZ + 1) * zVoxelDim; 

    }


    /**
     * Creates a region with the given minimum and maximum values of x, y
     * and z, in VOXEL coordinates. Boundaries are inclusive.
     *  
     * @param bounds {xMin, xMax, yMin, yMax, zMin, zMax}
     * @param voxelDims {x, y, z} voxel dimensions, in mm.
     */
    public CubicVoxelRegion(int[] bounds, double[] voxelDims) {
	this(bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5], voxelDims[0], 
	     voxelDims[1], voxelDims[2]);
    }

    
    
    /**
     * Creates a region with the given minimum and maximum values of x, y
     * and z, in VOXEL coordinates, and label 1.
     *  
     */
    public CubicVoxelRegion(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, 
			    double xVoxDim, double yVoxDim, double zVoxDim) {

        this(xmin, xmax, ymin, ymax, zmin, zmax, xVoxDim, yVoxDim, zVoxDim, 1);

    }



    /**
     * Creates a region with the given minimum and maximum values of x, y
     * and z, in VOXEL coordinates, and a specified region label.
     *  
     */
    public CubicVoxelRegion(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, 
			    double xVoxDim, double yVoxDim, double zVoxDim, int regionLabel) {


        label = regionLabel;

	if (xmin <= xmax) {
	    iminX = xmin;
	    imaxX = xmax;
	}
	else {
	    iminX = xmax;
	    imaxX = xmin;
	}
	if (ymin <= ymax) {
	    iminY = ymin;
	    imaxY = ymax;
	}
	else {
	    iminY = ymax;
	    imaxY = ymin;
	}
	if (zmin <= zmax) {
	    iminZ = zmin;
	    imaxZ = zmax;
	}
	else {
	    iminZ = zmax;
	    imaxZ = zmin;
	}

	xVoxelDim = xVoxDim;
	yVoxelDim = yVoxDim;
	zVoxelDim = zVoxDim;

	minXMM = (double)iminX * xVoxelDim;
	minYMM = (double)iminY * yVoxelDim;
	minZMM = (double)iminZ * zVoxelDim;
	maxXMM = (double)(imaxX + 1) * xVoxelDim; 
	maxYMM = (double)(imaxY + 1) * yVoxelDim;
	maxZMM = (double)(imaxZ + 1) * zVoxelDim; 


    }

    /**
     * Returns the lower x VOXEL bound of the region.
     */
    public int getMinXVoxel() {
	return iminX;
    }

    /**
     * Returns the lower y VOXEL bound of the region.
     */
    public int getMinYVoxel() {
	return iminY;
    }

    /**
     * Returns the lower z VOXEL bound of the region.
     */
    public int getMinZVoxel() {
	return iminZ;
    }

    /**
     * Returns the upper x VOXEL bound of the region.
     */
    public int getMaxXVoxel() {
	return imaxX;
    }

    /**
     * Returns the upper y VOXEL bound of the region.
     */
    public int getMaxYVoxel() {
	return imaxY;
    }

    /**
     * Returns the upper z VOXEL bound of the region.
     */
    public int getMaxZVoxel() {
	return imaxZ;
    }



    /**
     * Tests to see if a specific voxel is in this region.
     *
     * @param x the x index of the voxel
     * @param y the y index of the voxel
     * @param z the z index of the voxel
     * @return true if the voxel is within the region bounds; false otherwise
     */
    public boolean containsVoxel(int x, int y, int z) {
	return (x >= iminX && x <= imaxX && y >= iminY && y <= imaxY && z >= iminZ && z <= imaxZ);
    }


    /**
     * Tests to see if a point measured in VOXEL coordinates is in this region.
     *
     * @param point the point to test
     * @return true if the point is within the region bounds; false otherwise
     */
    public boolean containsMMPoint(Point3D point) {

	return ( point.x >= minXMM && point.x < maxXMM && 
		 point.y >= minYMM && point.y < maxYMM && 
		 point.z >= minZMM && point.z < maxZMM );
    }



    /**
     * @return Seed points as voxels in some resampled (but aligned) space
     */
    public Voxel[] getSeedVoxels(double xVoxelDim, double yVoxelDim, double zVoxelDim) {
        PointListROI plROI = new PointListROI(getSeedPoints(), label);

        return plROI.getSeedVoxels(xVoxelDim, yVoxelDim, zVoxelDim);

    }


    public Voxel[] getSeedVoxels() {
        Voxel[] voxels = new Voxel[(imaxX - iminX + 1) * (imaxY - iminY + 1) * (imaxZ - iminZ + 1)];
        
        int voxelCounter = 0;
        
        // must match order of seed points
        for (int k = iminZ; k <= imaxZ; k++) {
            for (int j = iminY; j <= imaxY; j++) {
                for (int i = iminX; i <= imaxX; i++) {
                    voxels[voxelCounter++] = new Voxel(i,j,k);
                }
	    }
	}
        
	return voxels;
    }


    /**
     * @return a point at the centre of every voxel in this region.
     */
    public Point3D[] getSeedPoints() {

	Point3D[] centres = new Point3D[(imaxX - iminX + 1) * (imaxY - iminY + 1) * (imaxZ - iminZ + 1)];

	int counter = 0;

	for (int k = iminZ; k <= imaxZ; k++) {
	    for (int j = iminY; j <= imaxY; j++) {
		for (int i = iminX; i <= imaxX; i++) {
		    centres[counter++] = new Point3D((i + 0.5) * xVoxelDim,
						     (j + 0.5) * yVoxelDim,
						     (k + 0.5) * zVoxelDim);
		}
	    }
	}

	return centres;

    }


    /**
     * Returns true if the region has identical bounds and voxel dimensions as this one.
     */
    public boolean equals(Object o) {

	if ( (o instanceof CubicVoxelRegion) == false ) {
	    return false;
	}
	else if (o == this) {
	    return true;
	}
	else {

	    CubicVoxelRegion oRegion = (CubicVoxelRegion)o;
	    
	    return (iminX == oRegion.iminX) && (imaxX == oRegion.imaxX) && (iminY == oRegion.iminY) && 
		(imaxY == oRegion.imaxY) && (iminZ == oRegion.iminZ) && (imaxZ == oRegion.imaxZ) && 
		(xVoxelDim == oRegion.xVoxelDim) && (yVoxelDim == oRegion.yVoxelDim) && 
		(zVoxelDim == oRegion.zVoxelDim);
	}

    }

    
    


    /**
     * Returns the CubicRegion in String format.
     */
    public String toString() {
	String s = "CubicRegion (VOXEL BOUNDS):\n";
	s += "xMin: " + iminX + "\txMax: " + imaxX + "\n";
	s += "yMin: " + iminY + "\tyMax: " + imaxY + "\n";
	s += "zMin: " + iminZ + "\tzMax: " + imaxZ + "\n";
	return s;
    }


    /**
     * @return 1.
     */
    public int numberOfRegions() {
	return 1;
    }


    public int getRegionLabel() {
        return label;
    }


}

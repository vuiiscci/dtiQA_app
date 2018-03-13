package tractography;

import numerics.*;

import java.io.*;
import java.util.*;


/**
 * Holds a streamline as a list of voxels.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class VoxelList {


    private final Voxel[] voxels;

    private final int seedPointIndex;

    private final int size;

    private final short[] tangentAtSeed;
    private final short[] negTangentAtSeed;

    private double xVoxelDim;
    private double yVoxelDim;
    private double zVoxelDim;


    /**
     *
     * @param voxels list of voxels in this streamline. 
     *
     * @param seedPointIndex the voxel <code>voxels[seedPointIndex]</code> contains the 
     * seed point of this streamline.
     *
     * @param xVoxelDim the x voxel dimension in mm.
     * @param yVoxelDim the y voxel dimension in mm.
     * @param zVoxelDim the z voxel dimension in mm.
     *
     * @param tangentAtSeed the vector that is tangential to the seed point trajectory at the seed point,
     * pointing upwards (from point (seedPointIndex) towards point (seedPointIndex+1).
     * @param negativeTangentAtSeed the vector that is tangential to the seed point trajectory at 
     * the seed point,
     * pointing downwards (from point (seedPointIndex) towards point (seedPointIndex-1).
     *
     */
    public VoxelList(Voxel[] voxels, int seedPointIndex, double xVoxelDim, double yVoxelDim,
		     double zVoxelDim, Vector3D tanAtSeed, Vector3D negativeTanAtSeed) {

	this.voxels = voxels;
	
	size = voxels.length;

	this.seedPointIndex = seedPointIndex;

	this.xVoxelDim = xVoxelDim;
	this.yVoxelDim = yVoxelDim;
	this.zVoxelDim = zVoxelDim;

	tangentAtSeed = vectorToShort(tanAtSeed.normalized());
	negTangentAtSeed = vectorToShort(negativeTanAtSeed.normalized());

	// represent as short


    }


    public Voxel getVoxel(int i) {
	return voxels[i];
    }


    /**
     * @return the index <code>i</code> of the seed, such that calling <code>getVoxel(i)</code> 
     * returns the voxel containing the seed point.
     * 
     */
    public int seedPointIndex() {
	return seedPointIndex;
    }

    /**
     * @return the number of voxels in this list.
     *
     */
    public int size() {
	return size;
    }


    public Voxel[] getVoxels() {

	Voxel[] copy = new Voxel[size];

	System.arraycopy(voxels, 0, copy, 0, size);

	return copy;
    }


    public String toString() {
	
	StringBuffer buffer = new StringBuffer();

	buffer.append("Number of voxels: " + size + "\n");
	buffer.append("seed voxel: " + seedPointIndex + "\n");
	buffer.append("Voxels:\n");
	
	for (int i = 0; i < size; i++) {
	    buffer.append(voxels[i].x);
	    buffer.append(" ");
	    buffer.append(voxels[i].y);
	    buffer.append(" ");
	    buffer.append(voxels[i].z);
	    buffer.append("\n");
	}

	return buffer.toString();
    }



    
    /**
     * Makes a Tract with a point at the centre of each voxel in the list.
     *
     */
    public Tract toTract() {
	Tract t1 = new Tract(size / 2 + 1, 100.0);
	Tract t2 = new Tract(size / 2 + 1, 100.0);

	Point3D seedPoint = new Point3D((voxels[seedPointIndex].x + 0.5) * xVoxelDim, 
					(voxels[seedPointIndex].y + 0.5) * yVoxelDim, 
					(voxels[seedPointIndex].z + 0.5) * zVoxelDim);

	Vector3D tan = new Vector3D(tangentAtSeed[0], tangentAtSeed[1], tangentAtSeed[2]).normalized();

	Vector3D negTan = new Vector3D(negTangentAtSeed[0], negTangentAtSeed[1], 
				       negTangentAtSeed[2]).normalized();

	t1.addPoint(seedPoint);

	// encode tangent by putting in extra point, 0.01 mm from seed
	t1.addPoint(seedPoint.displace(negTan.scaled(0.01)));

	t2.addPoint(seedPoint);

	// encode tangent by putting in extra point, 0.01 mm from seed
	t2.addPoint(seedPoint.displace(tan.scaled(0.01)));


	for (int i = seedPointIndex - 1; i >= 0; i--) {

	    Point3D point = new Point3D((voxels[i].x + 0.5) * xVoxelDim, 
					(voxels[i].y + 0.5) * yVoxelDim, 
					(voxels[i].z + 0.5) * zVoxelDim);
	    t1.addPoint(point);
	}
	for (int i = seedPointIndex + 1; i < size; i++) {
	    
	    Point3D point = new Point3D((voxels[i].x + 0.5) * xVoxelDim, 
					(voxels[i].y + 0.5) * yVoxelDim, 
					(voxels[i].z + 0.5) * zVoxelDim);

	    t2.addPoint(point);
	}

	t2.joinTract(t1);
	
	return t2;

    }



    public void writeVoxelList(DataOutputStream dout) throws IOException {
	
	dout.writeShort((short)voxels.length);
	dout.writeShort((short)seedPointIndex);

	dout.writeShort(tangentAtSeed[0]);
	dout.writeShort(tangentAtSeed[1]);
	dout.writeShort(tangentAtSeed[2]);

	dout.writeShort(negTangentAtSeed[0]);
	dout.writeShort(negTangentAtSeed[1]);
	dout.writeShort(negTangentAtSeed[2]);
	
	for (int point = 0; point < voxels.length; point++) {
	    dout.writeShort((short)voxels[point].x);
	    dout.writeShort((short)voxels[point].y);
	    dout.writeShort((short)voxels[point].z);
	}
    }
    

    private short[] vectorToShort(Vector3D v) {

	double x = v.x;
	double y = v.y;
	double z = v.z;

	// correct for errors
	if (x > 1.0) {
	    x = 1.0;
	}
	else if (x < -1.0) {
	    x = -1.0;
	}


	if (y > 1.0) {
	    y = 1.0;
	    
	}
	else if (y < -1.0) {
	    y = -1.0;
	}


	if (z > 1.0) {
	    z = 1.0;
	}
	else if (z < -1.0) {
	    z = -1.0;
	}


	return new short[] {(short)(Short.MAX_VALUE * x), (short)(Short.MAX_VALUE * y),
			    (short)(Short.MAX_VALUE * z)};

    }
    
    
    public double[] getVoxelDims() {

	return new double[] {xVoxelDim, yVoxelDim, zVoxelDim};
    } 

    
}

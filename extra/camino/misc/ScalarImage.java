package misc;

import java.io.*;

import numerics.*;
import imaging.*;
import tractography.*;

/**
 * Interpolated scalar image. Provides a scalar value at any real-valued point (in mm space) within the
 * image.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class ScalarImage {


    /** Nearest neighbour interpolation  */
    private boolean interpNN = true;

    /** Linear (eight neighbour) interpolation  */
    private boolean interpLinear = false;


    // Public for fast access to data array.
    public double[][][] data;
    
    
     // voxel dims in mm
    private int[] dataDims;

    private double[] voxelDims;


    private EightNeighbourInterpolator linearInterpolator;


    /**
     * Construct an image, default interpolation is nearest neighbour.
     *
     * @param dataDims image dimensions.
     * @param voxelDims voxel dimensions, in mm.
     */
    public ScalarImage(double[][][] data, double[] voxelDims) {
	
	this.data = data;
	
	dataDims = new int[3];

	dataDims[0] = data.length;
	dataDims[1] = data[0].length;
	dataDims[2] = data[0][0].length;

	this.voxelDims = voxelDims;

	linearInterpolator = new EightNeighbourInterpolator(dataDims, voxelDims);

    }


    
    /**
     * Construct an image from a file.
     *
     * @param fileName the name of the image, format is determined by the extension.
     */
    public ScalarImage(String fileName) throws IOException {
	loadData(fileName, 0);	
    }


    /**
     * Construct an image from one component of a 4D file, default
     * interpolation is nearest neighbour.
     *
     * @param fileName File name 
     * @param index The index of the component volume to read.
     */
    public ScalarImage(String fileName, int index) throws IOException {
	loadData(fileName, index);
    }


    /**
     * Load the image from one component of an image file, default
     * interpolation is nearest neighbour.
     *
     * @param fileName File name with extension
     *
     * @param index The index of the component volume to read.
     */
    private void loadData(String fileName, int index) throws IOException {

	ImageHeader header = ImageHeader.readHeader(fileName);

	double[][][][] imageData = header.readVolumeData();
        
	dataDims = new int[3];
	dataDims[0] = header.xDataDim();
	dataDims[1] = header.yDataDim();
	dataDims[2] = header.zDataDim();

	voxelDims = new double[3];
	voxelDims[0] = Math.abs(header.xVoxelDim());
	voxelDims[1] = Math.abs(header.yVoxelDim());
	voxelDims[2] = Math.abs(header.zVoxelDim());

        data = new double[dataDims[0]][dataDims[1]][dataDims[2]];

        for (int i = 0; i < dataDims[0]; i++) {
            for (int j = 0; j < dataDims[1]; j++ ) {
                for (int k = 0; k < dataDims[2]; k++) {
                    data[i][j][k] = imageData[i][j][k][index];
                }
            }
        }


	linearInterpolator = new EightNeighbourInterpolator(dataDims, voxelDims);

    }

    

    /**
     * If the point is outside the image array, the function returns
     * zero.
     *
     * @return the value of the image at a point p with mm coordinates.
     *
     */
    public double valueAt(Point3D p) {

	double value = 0.0;
	
        if(p.x<0 || p.y<0 || p.z<0 || p.x>dataDims[0]*voxelDims[0] || p.y>dataDims[1]*voxelDims[1] || p.z>dataDims[2]*voxelDims[2])
            return 0.0;
   

	if (interpNN) {
            int xInd = (int)(p.x / voxelDims[0]);
            int yInd = (int)(p.y / voxelDims[1]);
            int zInd = (int)(p.z / voxelDims[2]);
            xInd = (xInd<0)?0:xInd;
            yInd = (yInd<0)?0:yInd;
            zInd = (zInd<0)?0:zInd;
            xInd = (xInd>=dataDims[0])?(dataDims[0]-1):xInd;
            yInd = (yInd>=dataDims[1])?(dataDims[1]-1):yInd;
            zInd = (zInd>=dataDims[2])?(dataDims[2]-1):zInd;
	    value = data[xInd][yInd][zInd];
	}
	else if (interpLinear) {
        
            double[] interpFraction = new double[8];
            int[] dims = new int[6];

	    // get the interpolation parameters
	    int inVoxel = linearInterpolator.setInterpolationVoxels(p, interpFraction, dims);

	    for (int i = 0; i < 8; i++) {

		int x = dims[i / 4];
		int y = dims[2 + ((i / 2) % 2)];
		int z = dims[4 + (i % 2)];
		
		value += interpFraction[i] * data[x][y][z];
	    }
        }
	else {
	    throw new LoggedException("Internal error: unknown interpolation method");
	}

	return value;
    }


    /**
     *
     * @return the values at each of the points.
     *
     */
    public double[] valuesAt(Point3D[] points) {

	int p = points.length;

	double[] values = new double[p];

	for (int i = 0; i < p; i++) {
	    values[i] = valueAt(points[i]);
	}

	return values;
    }


    /**
     * @return the values at the centre of the voxels. The voxels are assumed to be defined in the
     * same space as the image.
     *
     */
    public double[] valuesAt(Voxel[] voxels) {

	int v = voxels.length;

	double[] values = new double[v];

	for (int i = 0; i < v; i++) {
	    values[i] = data[voxels[i].x][voxels[i].y][voxels[i].z];
	}

	return values;
	
    }


    /**
     *
     * @return the derivative of the image at a point p with mm coordinates.
     *
     */
    public double[] derivAt(Point3D p) {

	double[] deriv = {0.0, 0.0, 0.0};
	
        if(p.x<0 || p.y<0 || p.z<0 || p.x>dataDims[0]*voxelDims[0] || p.y>dataDims[1]*voxelDims[1] || p.z>dataDims[2]*voxelDims[2])
            return deriv;
   
	if (interpNN) {
            int iX = (int)(p.x / voxelDims[0]);
            int iY = (int)(p.y / voxelDims[1]);
            int iZ = (int)(p.z / voxelDims[2]);

            // Check they are inside the image.
            iX = (iX>=0)?iX:0;
            iX = (iX<dataDims[0])?iX:(dataDims[0]-1);
            iY = (iY>=0)?iY:0;
            iY = (iY<dataDims[1])?iY:(dataDims[1]-1);
            iZ = (iZ>=0)?iZ:0;
            iZ = (iZ<dataDims[2])?iZ:(dataDims[2]-1);

            // Perturbed indices
            int iXp = iX+1;
            iXp = (iXp<dataDims[0])?iXp:(dataDims[0]-1);
            int iYp = iY+1;
            iYp = (iYp<dataDims[1])?iYp:(dataDims[1]-1);
            int iZp = iZ+1;
            iZp = (iZp<dataDims[2])?iZp:(dataDims[2]-1);

	    deriv[0] = (data[iXp][iY][iZ] - data[iX][iY][iZ])/voxelDims[0];
	    deriv[1] = (data[iX][iYp][iZ] - data[iX][iY][iZ])/voxelDims[1];
	    deriv[2] = (data[iX][iY][iZp] - data[iX][iY][iZ])/voxelDims[2];
	}
	else if (interpLinear) {
        
	    // get the interpolation parameters
            int[] dims = new int[6];

            double[] interpFractionDX = new double[8];
	    int inVoxelX = linearInterpolator.setInterpolationVoxelsDX(p, interpFractionDX, dims);

            double[] interpFractionDY = new double[8];
	    int inVoxelY = linearInterpolator.setInterpolationVoxelsDY(p, interpFractionDY, dims);

            double[] interpFractionDZ = new double[8];
	    int inVoxelZ = linearInterpolator.setInterpolationVoxelsDZ(p, interpFractionDZ, dims);

	    for (int i = 0; i < 8; i++) {

		int x = dims[i / 4];
		int y = dims[2 + ((i / 2) % 2)];
		int z = dims[4 + (i % 2)];
		
		deriv[0] += interpFractionDX[i] * data[x][y][z];
		deriv[1] += interpFractionDY[i] * data[x][y][z];
		deriv[2] += interpFractionDZ[i] * data[x][y][z];
	    }
        }
	else {
	    throw new LoggedException("Internal error: unknown interpolation method");
	}

	return deriv;
    }


    /**
     * Sets the interpolation method for this image. 
     * @param interpMethod either "nearestneighbour" or "linear" for nearest neighbor or linear interpolation
     * respectively.
     */
    public void setInterpolation(String interpMethod) {
	if (interpMethod.equals("linear")) {
	    interpLinear = true;
	    interpNN = false;
	}
	else if (interpMethod.equals("nearestneighbour")) {
	    interpNN = true;
	    interpLinear = false;
	}
	else {
	    throw new LoggedException("Unknown interpolation method: " + interpMethod);
	}

    }


    /**
     *
     * @return the data dimensions of the image.
     *
     */
    public int[] getDataDims() {
	return new int[] {dataDims[0], dataDims[1], dataDims[2]};
    }
 

    /**
     *
     * @return the voxel dimensions of the image.
     *
     */
    public double[] getVoxelDims() {
	return new double[] {voxelDims[0], voxelDims[1], voxelDims[2]};
    }

 
}

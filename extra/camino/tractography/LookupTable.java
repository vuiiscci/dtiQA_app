package tractography;

import java.io.*;
import misc.LoggedException;

/**
 * Lookup table for PICo parameters.
 * 
 * @author Philip Cook
 * @version $Id$
 */
public class LookupTable {

    // number of values at each point in the LUT
    private final int valuesPerPosition;
    
    private final double[][][] lut;

    private final int dimension;
   
    private final double xMin;
    private final double xMax;   
    private final double yMin;
    private final double yMax;
    private final double zMin;
    private final double zMax;

    private final double xStep;
    private final double yStep;
    private final double zStep;

    /**
     * 
     * @param lut the actual data. The lut value at a given x,y,z is lut[(int)(x - xMin) / step][(int)(y - yMin) / step][values * (int)(z - zMin) / step)] for a 3D lut, lut[0][(int)(y - yMin) / step][values * (int)(z - zMin) / step] for 2D, and lut[0][0][values * (int)(z - zMin) / step] for 1D.
     * @param dims {xMin xMax yMin yMax zMin zMax}. Set unused dimensions to zero, eg for a 1D LUT pass {0.0, 0.0, 0.0, 0.0, zMin, zMax}.
     * @param step {xStep yStep zStep}. Follow the same procedure as above for lower dimensions.
     *
     */
    public LookupTable(double[][][] lut, double[] dims, double[] step, int values) {

	// check that params make sense
	this.lut = lut;
	

	xMin = dims[0];
	xMax = dims[1];
	yMin = dims[2];
	yMax = dims[3];
	zMin = dims[4];
	zMax = dims[5];
	
	xStep = step[0];
	yStep = step[1];
	zStep = step[2];

	valuesPerPosition = values;
	
	if (lut.length == 1) {
	    if (lut[0].length == 1) {
		dimension = 1;

		if ( ( 1 + (int)Math.round((zMax - zMin) / zStep) ) * values != lut[0][0].length ) {
		    throw new java.lang.IllegalArgumentException("LUT has the wrong dimensions for given range and step size: zRange " + zMin + " " + zMax + " step " + zStep + " values per position " + values);
		}


	    }
	    else {
		dimension = 2;
		if ( 1 + (int)( Math.round((yMax - yMin) / yStep) ) != lut[0].length || 
		    values * (1 + (int)( Math.round((zMax - zMin) / zStep) )) != lut[0][0].length ) {
		    throw new java.lang.IllegalArgumentException("LUT has the wrong dimensions for given range and step size: yRange " + yMin + " " + yMax + " zRange " + zMin + " " + zMax + " step " + yStep + " " + zStep + " values per position " + values);
		}

	    }
	}
	else {
	    dimension = 3;
	    if ( 1 + (int)( Math.round((xMax - xMin) / xStep) ) != lut.length || 
		 1 + (int)( Math.round((yMax - yMin) / yStep) ) != lut[0].length || 
		 values * (1 + (int)( Math.round( (zMax - zMin) / zStep ) )) != lut[0][0].length ) {
		    throw new java.lang.IllegalArgumentException("LUT has the wrong dimensions for given range and step size: xRange " + xMin + " " + xMax + " yRange " + yMin + " " + yMax +  " zRange " + zMin + " " + zMax + " step " + xStep + " " + yStep + " " + zStep + " values per position " + values + ". LUT dimensions: " + lut.length + " " + lut[0].length + " " + lut[0][0].length);
		}

	}
    }
    

    /**
     * Get the value(s) of the LUT at the point x,y.
     * @param x the x index of the LUT: zero for 2D LUTs.
     * @param y the y index of the LUT: zero for 1D LUTs.
     * @param z the z index of the LUT.
     * @param interpolate if true, interpolate between positions.
     * @param clamp if true, return max or min value instead of throwing an exception.
     * @throws OutsideLUTRangeException if x and y are outside the LUT range, and clamp is false.
     */
    public double[] getValues(double x, double y, double z, boolean interpolate, boolean clamp) throws OutsideLUTRangeException {

	
	double xPos = 0.0;
	double yPos = 0.0;
	double zPos = 0.0;


	if (clamp) {
	    if (z < zMin) {
		z = zMin;
	    }
	    else if (z > zMax) {
		z = zMax;
	    }
	    
	}
	else {
	    if (z < zMin || z > zMax) {
		throw new OutsideLUTRangeException("LUT index " + x + " " + y + " " + z + " is outside LUT range.");
	    }
	}

	
	zPos = (z - zMin) / zStep;

	if (dimension > 1) {
	    if (clamp) {
		if (y < yMin) {
		    y = yMin;
		}
		else if (y > yMax) {
		    y = yMax;
		}
		
	    }
	    else {
		if (y < yMin || y > yMax) {
		    throw new OutsideLUTRangeException("LUT index " + x + " " + y + " is outside LUT range.");
		}
	    }
	    
	    yPos = (y - yMin) / yStep;
	}

	if (dimension > 2) {
	    if (clamp) {
		if (x < xMin) {
		    x = xMin;
		}
		else if (x > xMax) {
		    x = xMax;
		}
		
	    }
	    else {
		if (x < xMin || x > xMax) {
		    throw new OutsideLUTRangeException("LUT index " + x + " " + y + " is outside LUT range.");
		}
	    }

	    xPos = (x - xMin) / xStep;
	    
	}



	double[] values = new double[valuesPerPosition];

	if (dimension == 1) {
	    
	    int zIndex = (int)(zPos);

	    double relZ = zPos - zIndex;
	    
	    if (z == zMax) {
		zIndex = (int)Math.round(zPos);
                relZ = 0.0;
	    }

	    if (interpolate) {
		for (int v = 0; v < valuesPerPosition; v++) {
		    
		    double lut0 = lut[0][0][valuesPerPosition * zIndex + v];
		    double lut1 = 0.0;
		    
		    if (z < zMax) {
			lut1 = lut[0][0][valuesPerPosition * (zIndex + 1) + v];
		    }
		    else {
			lut1 = lut0;
		    }
		    
		    values[v] = (1.0 - relZ) * lut0 + relZ * lut1;
		    
		}
	    }
	    else {
		// nearest neighbour
		int zNearest = relZ > 0.5 ? zIndex + 1 : zIndex;
		
		for (int v = 0; v < valuesPerPosition; v++) {
		    values[v] = lut[0][0][valuesPerPosition * zNearest + v];
		}
		
	    }
		
	
	}
	else if (dimension == 2) { // 2D
	    
	    int yIndex = (int)(yPos);
	    int zIndex = (int)(zPos);
	    
	    double relY = yPos - yIndex;
	    double relZ = zPos - zIndex;

	    if (y== yMax) {
		yIndex = (int)Math.round(yPos);
                relY = 0.0;
	    }
	    if (z == zMax) {
		zIndex = (int)Math.round(zPos);
                relZ = 0.0;
	    }

	    if (interpolate) {
		for (int v = 0; v < valuesPerPosition; v++) {
		    
		    int y0 = yIndex;
		    int y1;

		    int z0 = zIndex;
		    int z1;
		    
		    y1 = y == yMax ? y0 : y0 + 1;
		    z1 = z == zMax ? z0 : z0 + 1;
		    
		    double lut00, lut01, lut10, lut11;

		    lut00 = lut[0][y0][valuesPerPosition * z0 + v];
		    lut01 = lut[0][y0][valuesPerPosition * z1 + v];
		    lut10 = lut[0][y1][valuesPerPosition * z0 + v];
		    lut11 = lut[0][y1][valuesPerPosition * z1 + v];
		    		    
		    double[] interpComponents = new double[4];

		    values[v] = 
			(1.0 - relY) * (1.0 - relZ) * lut00 + 
			relY * (1.0 - relZ) * lut10 + 
			(1.0 - relY) * relZ * lut01 + relY * relZ * lut11;
		    
		}
	    }
	    else {
		int yNearest = relY < 0.5 ? yIndex : yIndex + 1;
		int zNearest = relZ < 0.5 ? zIndex : zIndex + 1;
		for (int v = 0; v < valuesPerPosition; v++) {
		    values[v] = lut[0][yNearest][valuesPerPosition * zNearest + v];	
		}
	    }
	    
	} // end 2D
	else { // dimension == 3

	    int xIndex = (int)(xPos);
	    int yIndex = (int)(yPos);
	    int zIndex = (int)(zPos);

	    // xPos, yPos, is bounded by a cube, which contain eight LUT values
	    // position relative to lower left corner of square 
	    double relX = xPos - xIndex;
	    double relY = yPos - yIndex;
	    double relZ = zPos - zIndex;

	    if (x == xMax) {
                xIndex = (int)Math.round(xPos);
                relX = 0.0;
            }
	    if (y == yMax) {
		yIndex = (int)Math.round(yPos); 
                relY = 0.0;

	    }
	    if (z == zMax) {
		zIndex = (int)Math.round(zPos);
                relZ = 0.0;
	    }
	    

	    if (interpolate) {

		for (int v = 0; v < valuesPerPosition; v++) {
		    double lut000, lut001, lut010, lut011, lut100, lut101, lut110, lut111;
		    
		    int x0 = xIndex;
		    int x1;

		    int y0 = yIndex;
		    int y1;

		    int z0 = zIndex;
		    int z1;
		    
		    x1 = x == xMax ? x0 : x0 + 1;
		    y1 = y == yMax ? y0 : y0 + 1;
		    z1 = z == zMax ? z0 : z0 + 1;
		    
		    lut000 = lut[x0][y0][valuesPerPosition * z0 + v];
		    lut001 = lut[x0][y0][valuesPerPosition * z1 + v];
		    lut010 = lut[x0][y1][valuesPerPosition * z0 + v];
		    lut011 = lut[x0][y1][valuesPerPosition * z1 + v];
		    lut100 = lut[x1][y0][valuesPerPosition * z0 + v];
		    lut101 = lut[x1][y0][valuesPerPosition * z1 + v];
		    lut110 = lut[x1][y1][valuesPerPosition * z0 + v];
		    lut111 = lut[x1][y1][valuesPerPosition * z1 + v];
		    
		    
		    double[] interpComponents = new double[8];

		    interpComponents[0] = (1.0 - relX) * (1.0 - relY) * (1.0 - relZ);
		    interpComponents[1] = (1.0 - relX) * (1.0 - relY) * relZ; 
		    interpComponents[2] = (1.0 - relX) * relY * (1.0 - relZ);
		    interpComponents[3] = (1.0 - relX) * relY * relZ;
		    interpComponents[4] = relX * (1.0 - relY) * (1.0 - relZ);
		    interpComponents[5] = relX *  (1.0 - relY) * relZ;
		    interpComponents[6] = relX * relY * (1.0 - relZ);
		    interpComponents[7] = relX * relY * relZ;


		    values[v] = lut000 * interpComponents[0] +
			lut001 * interpComponents[1] +
			lut010 * interpComponents[2] +
			lut011 * interpComponents[3] +
			lut100 * interpComponents[4] +
			lut101 * interpComponents[5] +
			lut110 * interpComponents[6] +
			lut111 * interpComponents[7]; 
		    
		}

	    }
	    else {
		int xNearest = relX < 0.5 ? xIndex : xIndex + 1;
		int yNearest = relY < 0.5 ? yIndex : yIndex + 1;
		int zNearest = relZ < 0.5 ? zIndex : zIndex + 1;
		for (int v = 0; v < valuesPerPosition; v++) {
		    values[v] = lut[xNearest][yNearest][valuesPerPosition * zNearest + v];	
		}
		
	    }

	}

	return values;

    }


    /**
     * Get the value(s) of the LUT at the point z. Should only be called for 1D LUTs.
     * 
     * @param interpolate if true, interpolate between positions.
     * @param clamp if true, return max or min value instead of throwing an exception.
     * @throws OutsideLUTRangeException if the LUT is not 1D, or y is outside of the LUT range.
     */
    public double[] getValues(double z, boolean interpolate, boolean clamp) throws OutsideLUTRangeException {
	
       if (dimension == 1) {
	   return getValues(0.0, 0.0, z, interpolate, clamp);
       }
       else {
	   throw new OutsideLUTRangeException("LUT requires " + dimension + " indices");
       }
    }



  /**
     * Get the value(s) of the LUT at the point y,z. Should only be called for 2D LUTs.
     * 
     * @param interpolate if true, interpolate between positions.
     * @param clamp if true, return max or min value instead of throwing an exception.
     * @throws OutsideLUTRangeException if the LUT is not 2D, or y or z are outside of the LUT range.
     */
    public double[] getValues(double y, double z, boolean interpolate, boolean clamp) throws OutsideLUTRangeException {
	
       if (dimension == 2) {
	   return getValues(0.0, y, z, interpolate, clamp);
       }
       else {
	   throw new OutsideLUTRangeException("LUT requires " + dimension + " indices");
       }
    }



    public double xMin() {
	return xMin;
    }

    public double yMin() {
	return yMin;
    }


    public double xMax() {
	return xMax;
    }

    public double yMax() {
	return yMax;
    }


    public double zMin() {
	return zMin;
    }

    public double zMax() {
	return zMax;
    }

    public double xStep() {
	return xStep;
    }

    public double yStep() {
	return yStep;
    }

    public double zStep() {
	return zStep;
    }


    /**
     * Dimensionality of the LUT.
     *
     */
    public int dimension() {
	return dimension;
    }

    /**
     * The length of the array returned by 
     * <code>getValues</code>.
     */
    public int valuesPerPosition() {
	return valuesPerPosition;
    }


    /**
     * LUT Format: xMin xMax yMin yMax zMin zMax xStep yStep zStep valuesPerPosition <data>.
     * 
     *
     * @return null if the file cannot be opened.
     *
     */
    public static LookupTable readLUT(String fileName) {

	// 
	// lut is 1D if xMin == xMax == 0

	FileInputStream fit = null;

        if (fileName == null) {
            throw new LoggedException("Tried to read LUT, but file name is null");
        }

	try {
   
	    fit = new FileInputStream(fileName);
	    
	    DataInputStream dit = new DataInputStream(new BufferedInputStream(fit, 4194304));

	    double xMin = dit.readDouble();
	    double xMax = dit.readDouble();
	    
	    double yMin = dit.readDouble();
	    double yMax = dit.readDouble();

	    double zMin = dit.readDouble();
	    double zMax = dit.readDouble();
	    
	    double xStep = dit.readDouble();
	    double yStep = dit.readDouble();
	    double zStep = dit.readDouble();
	
	    int valuesPerPosition = (int)dit.readDouble();

	    double[][][] lut = new double[1 + (int)(Math.round( (xMax - xMin ) / xStep ) )]
		[1 + (int)(Math.round((yMax - yMin) / yStep)) ]
		[valuesPerPosition * ( 1 + (int)(Math.round((zMax - zMin) / zStep)) ) ];

	    for (int v = 0; v < valuesPerPosition; v++) {
		for (int k = 0; k < lut[0][0].length / valuesPerPosition; k++) {
		    for (int j = 0; j < lut[0].length; j++)  {
			for (int i = 0; i < lut.length; i++) {
			    lut[i][j][k * valuesPerPosition + v] = dit.readDouble();
			}
		    }
		}
	    }
	    
	    dit.close();

	    LookupTable l = new LookupTable(lut, new double[] {xMin, xMax, yMin, yMax, zMin, zMax}, 
					    new double[] {xStep, yStep, zStep}, valuesPerPosition);
	    
	    return l;

	}
	catch (IOException e) {
            throw new LoggedException(Thread.currentThread().getName(), e);
	}
    }

    
}

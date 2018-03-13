package data;

import misc.*;
import tools.*;

import java.io.*;



/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Source of data from a scanner-order data file or standard input stream.
 * 
 * <dt>Description:
 * 
 * <dd>Reads the data into a 4D array and returns 1D arrays of data from each
 * individual voxel in turn.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: ScannerOrderDataSource.java,v 1.4 2005/08/18 10:59:37 ucacmgh
 *          Exp $
 *  
 */
public class ScannerOrderDataSource extends ExternalDataSource {

    /**
     * Stores the total number of voxels in the data set.
     */
    protected int numVoxels;

    /**
     * The index of the next voxel to return.
     */
    protected int nextVoxel = 0;

    /**
     * The data array.
     */
    protected double[][] data;
 

    /**
     * Constructor requires the filename, the dimensions of the measurement
     * volume, the number of values in each voxel and the data type. If the
     * filename is null, the object reads from the standard input.
     * 
     * @param filename
     *            The name of the data file.
     * 
     * @param numVox
     *            The number of voxels in the measurement volume.
     * 
     * @param components
     *            The number of values in each voxel.
     * 
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     */
    public ScannerOrderDataSource(String filename, int numVox, int components, String type) {
	initFileInput(filename, false, 0);
        
	init(numVox, components, type);
    }


    /**
     * Constructor requires the filename, the dimensions of the measurement
     * volume, the number of values in each voxel and the data type. If the
     * filename is null, the object reads from the standard input.
     * 
     * @param filename
     *            The name of the data file.
     * 
     * @param numVox
     *            The number of voxels in the measurement volume.
     * 
     * @param components
     *            The number of values in each voxel.
     * 
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     *
     * @param intelByteOrder 
     *            true if the ordering of the input is little-endian.
     */
    public ScannerOrderDataSource(String filename, int numVox, int components, 
				  String type, boolean intelByteOrder) {


	initFileInput(filename, intelByteOrder, 0);

        init(numVox, components, type);
    }



    /**
     * Constructor requires the filename, the number of values in each voxel and
     * the data type. If the filename is null, the object reads from the
     * standard input.
     * 
     * @param filename
     *            The name of the data file.
     * 
     * @param numVox
     *            The number of voxels in the measurement volume.
     *
     * @param components
     *            The number of values in each voxel.
     * 
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     *
     * @param intelByteOrder 
     *            true if the byte ordering of the input stream is little-endian.
     *
     * @param offset
     *            read and discard this many bytes; used to skip headers.

     */
    public ScannerOrderDataSource(String filename, int numVox, int components, String type, 
				boolean intelByteOrder, int offset) {

	initFileInput(filename, intelByteOrder, offset);

        init(numVox, components, type);
    }


    /**
     * Constructs a data source from an array.
     *
     * @param data where <code>data[x][y]</code> is component y for voxel x.
     */
    public ScannerOrderDataSource(double[][] data) {
	this.data = data;
	numVoxels = data.length;
	numComponents = data[0].length;
    }


    /**
     * Creates the storage array and reads all the data into it.
     * 
     * @param numVox
     *            The number of voxels in the measurement volume.
     * 
     * @param components
     *            The number of values in each voxel.
     * 
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     */
    protected void init(int numVox, int components, String type) {

        numComponents = components;
        numVoxels = numVox;

        // Create the data array.
        data = new double[numVoxels][numComponents];

        try {

            // Set the data type. This can fail if the data type
            // String is not recognized.
            datatype = getDataTypeCode(type);

        }
        catch (Exception e) {
            throw new LoggedException(e);
        }

        // Read the data into the array.
        for (int i = 0; i < numComponents; i++) {
            for (int j = 0; j < numVoxels; j++)
                try {
                    if (datatype == BYTE) {
                        data[j][i] = (double) dataIn.readByte();
                    }
                    else if (datatype == UBYTE) {
                        data[j][i] = (double)dataIn.readUnsignedByte();
                    }
                    else if (datatype == SHORT) {
                        data[j][i] = (double) dataIn.readShort();
                    }
                    else if (datatype == INT) {
                        data[j][i] = (double) dataIn.readInt();
                    }
                    else if (datatype == USHORT) {
                        data[j][i] = (double) dataIn.readUnsignedShort();
                    }
                    else if (datatype == UINT) {
                        data[j][i] = (double) dataIn.readUnsignedInt();
                    }
                    else if (datatype == LONG) {
                        data[j][i] = (double) dataIn.readLong();
                    }
                    else if (datatype == FLOAT) {
                        data[j][i] = (double) dataIn.readFloat();
                    }
                    else if (datatype == DOUBLE) {
                        data[j][i] = (double) dataIn.readDouble();
                    }

                }
                catch (EOFException e) {
                    throw new LoggedException("\nRan out of data.  Check inputdatatype.\nCurrently reading input data as " + typeString(datatype) + ". \n" + e);
                }
                catch (Exception e) {
                    throw new LoggedException(e);
                }
        }
    }

    public double[] nextVoxel() throws DataSourceException {

        if (nextVoxel >= numVoxels) {
            throw new DataSourceException("No more voxels.");
        }

        double[] next = data[nextVoxel];
        nextVoxel += 1;

        return next;
    }

    public boolean more() {
        return nextVoxel < numVoxels;
    }

}

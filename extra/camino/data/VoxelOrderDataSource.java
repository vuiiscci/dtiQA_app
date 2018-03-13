package data;

import java.io.*;

import misc.*;
import tools.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Source of voxel order data from a data file or standard input stream.
 * 
 * <dt>Description:
 * 
 * <dd>Reads in the data voxel by voxel and returns arrays of data from each
 * individual voxel in turn.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class VoxelOrderDataSource extends ExternalDataSource {

    /**
     * An array that always contains the data for the next voxel.
     */
    protected double[] next;

    /**
     * Flags use to detect when the end of the input stream is reached and when
     * there is no more voxel data to return.
     */
    private boolean reachedEndOfFile;

    private boolean noMoreVoxels;

    /**
     * Default constructor required by inherited classes.
     */
    protected VoxelOrderDataSource() {
    }

    /**
     * Constructor requires the filename, the number of values in each voxel and
     * the data type. If the filename is null, the object reads from the
     * standard input.
     * 
     * @param filename
     *            The name of the data file.
     * 
     * @param components
     *            The number of values in each voxel.
     * 
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     */
    public VoxelOrderDataSource(String filename, int components, String type) {
	initFileInput(filename, false, 0);
	
        init(components, type);
    }


    /**
     * 
     * @param filename
     *            The name of the data file.
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
     */
    public VoxelOrderDataSource(String filename, int components, String type, boolean intelByteOrder) {
	initFileInput(filename, intelByteOrder, 0);

        init(components, type);
    }

    /**
     * 
     * @param filename
     *            The name of the data file.
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
    public VoxelOrderDataSource(String filename, int components, String type, boolean intelByteOrder, int offset) {
	
	initFileInput(filename, intelByteOrder, offset);
	
	init(components, type);
    }


 
    /**
     * Initialize a VoxelOrderDataSource from an input stream.
     *
     * @param input
     *            The stream, which should be ready to read data.
     *
     * @param components
     *            The number of values in each voxel.
     * 
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     *
     *
     */
    public VoxelOrderDataSource(EndianNeutralDataInputStream input, int components, String type) {
	
	dataIn = input;
	
	init(components, type);
    }
  


    /**
     * Does the work of the constructor.
     * 
     * @param components
     *            The number of values in each voxel.
     * 
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     */
    protected void init(int components, String type) {

        numComponents = components;
        next = new double[numComponents];

        try {

            // Set the data type. This can fail if the data type
            // String is not recognized.
            datatype = getDataTypeCode(type);

        }
        catch (Exception e) {
            throw new LoggedException(e);
        }

        // Initialise the flags indicating end of input file and whether
        // more voxel data can be requested.
        reachedEndOfFile = false;
        noMoreVoxels = false;

        // Read in the first voxel ready for the first data request.
        try {
            readNextVoxel();
        }
        catch (Exception e) {
            throw new LoggedException(e);
        }
    }

    public double[] nextVoxel() throws DataSourceException {

        // Copy the current voxel data into an array for return.
        double[] vox = retrieveCurrentVoxelData();

        // Read in the data from the next voxel ready for the next
        // call to this method.
        readNextVoxel();

        // The previous read may have been the last voxel, in which case
        // the current read failed indicating that this is the last
        // voxel's worth of data to return.
        if (reachedEndOfFile) {
            noMoreVoxels = true;
        }

        // Return the current voxel's data.
        return vox;
    }

    public boolean more() {
        return !noMoreVoxels;
    }

    /**
     * Puts the data for the current voxel into an array ready for return by
     * nextVoxel.
     * 
     * @return An array containing the voxel data.
     */
    protected double[] retrieveCurrentVoxelData() {
        double[] vox = new double[numComponents];
        for (int i = 0; i < numComponents; i++) {
            vox[i] = next[i];
        }
        return vox;
    }

    /**
     * Reads the next voxel's worth of data into the array next. If the read
     * fails, the method sets a flag to indicate the end of the data file. If
     * the read fails half way through, the method reports a separate error
     * indicating that the file did not contain a whole number of voxels.
     */
    private void readNextVoxel() throws DataSourceException {

        if (reachedEndOfFile) {
            noMoreVoxels = true;
            throw new DataSourceException("No more voxels in data source.");
        }

        // Need this in the Exception handling.
        int indexReached = 0;

        for (int i = 0; i < numComponents; i++)
            try {

                indexReached = i;

                if (datatype == BYTE) {
                    next[i] = (double) dataIn.readByte();
                }
                else if (datatype == UBYTE) {
                    next[i] = (double)dataIn.readUnsignedByte();
                }
                else if (datatype == SHORT) {
                    next[i] = (double) dataIn.readShort();
                }
                else if (datatype == USHORT) {
                    next[i] = (double) dataIn.readUnsignedShort();
                }
                else if (datatype == INT) {
                    next[i] = (double) dataIn.readInt();
                }
                else if (datatype == UINT) {
                    next[i] = (double) dataIn.readUnsignedInt();
                }
                else if (datatype == LONG) {
                    next[i] = (double) dataIn.readLong();
                }
                else if (datatype == FLOAT) {
                    next[i] = (double) dataIn.readFloat();
                }
                else if (datatype == DOUBLE) {
                    next[i] = (double) dataIn.readDouble();
                }

            }
            catch (IOException e) {

                reachedEndOfFile = true;

                if (indexReached > 0) {

                    // We have reached the end of the file part way
                    // through one voxel.
                    noMoreVoxels = true;
                    throw new DataSourceException(
			      "End of file reached without completing voxel.");
                }

                // Close the data stream.
                try {
                    dataIn.close();
                }
                catch (Exception e2) {
                    System.err.println("Failed to close the data input stream in VoxelOrderDataSource. Trying to continue...");
                    
                }

                return;
            }

    }

}

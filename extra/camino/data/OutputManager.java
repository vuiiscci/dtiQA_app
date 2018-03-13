package data;

import java.io.*;
import java.util.logging.Logger;

import imaging.*;
import tools.*;
import misc.LoggedException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Handles output of data from various applications.
 * 
 * <dt>Description:
 * 
 * <dd>Implements various output modes, such as output to files or
 * standard output, or stored output for return as an array.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: ModelFit.java 58 2006-10-26 11:38:18Z ucacmgh $
 *  
 */
public class OutputManager {


    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.OutputManager");


    public static final int FILEBUFFERSIZE = 1024*1024*24;


    /**
     * The stream to which to write the output.
     */
    private DataOutputStream out = new DataOutputStream(new BufferedOutputStream(System.out, FILEBUFFERSIZE));


    /**
     * The type of the output data
     */
    public static String outputDataType = "double";


    /**
     * Name of the output file.
     */
    public static String outputFile = null;


    /**
     * Write GZIP output.
     */
    public static boolean gzipOut = false;


    /**
     * Specifies whether to store the output as an array for collection
     * by some external routine, such as Matlab.
     */
    private boolean useOutputArray = false;


    /**
     * True if we output to an image. Implies useOutputArray. Image will be written once
     * the output is complete and the OutputManager is closed.
     *
     */
    private boolean useOutputImage = false;


    /**
     * This is for storing output.
     */
    private double[][][][] outputArray;


    /**
     * Size of output array.
     */
    private int outputArrayX = -1;
    private int outputArrayY = -1;
    private int outputArrayZ = -1;


    /**
     * Indices of output array.
     */
    private int nextX = 0;
    private int nextY = 0;
    private int nextZ = 0;



    /**
     * Sets up an output stream for an output file if the output file
     * has been specified, otherwise stays with the default output
     * stream.
     */
    public OutputManager() {

        // Set up the output stream.
        if (outputFile != null) { 
            
            try {

                if (!ImageHeader.getFileRoot(outputFile).equals(outputFile)) {

                    // if an extension was found and stripped off imageFileRoot, output is to an image
                    useOutputImage = true;

                    if (CL_Initializer.headerTemplate == null) {
                        CL_Initializer.initInputSpaceAndHeaderOptions();
                    }

                    int[] dataDims = CL_Initializer.headerTemplate.getDataDims();
                    
                    setOutputArray(dataDims[0], dataDims[1], dataDims[2]);
                    
                }

                
                
                
                if(gzipOut) {
                    out = new DataOutputStream
                        (new java.util.zip.GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile), 
                                                                                     FILEBUFFERSIZE)));
                }
                else {
                    out = new DataOutputStream
                        (new BufferedOutputStream(new FileOutputStream(outputFile), FILEBUFFERSIZE));
                }
            }
            catch (Exception e) {
                throw new LoggedException(e);
            }
        }
	else if(gzipOut) { 
            try {
                out = new DataOutputStream
                    (new java.util.zip.GZIPOutputStream(new BufferedOutputStream(System.out, FILEBUFFERSIZE)));
            }
            catch (Exception e) {
                throw new LoggedException(e);
            }
        }

    }


    /**
     * Outputs the data for one voxel to the selected output stream.
     * 
     * @param data
     *            The array of data values to output.
     * 
     * @param f
     *            The file output stream. If this is null, the output goes to
     *            standard out.
     * 
     * @param outputDataType
     *            The type of the data output.
     */
    public void output(double[] data) {


        // If user specifies that output should be returned as an
        // array, just call this method to store the latest voxel
        // output.  This is for the Matlab interface.
        if(useOutputArray) {
            addToOutputArray(data);
            return;
        }


        // Otherwise output the next voxel to a file or output stream
        // in the right format.

        // Would be much more efficient to do this once, when datatype is set

        int datatype = ExternalDataSource.DOUBLE;
        try {
            datatype = ExternalDataSource.getDataTypeCode(outputDataType);
        }
        catch (Exception e) {
            LoggedException.logExceptionWarning(e,Thread.currentThread().getName());  
            logger.warning("Outputting data as doubles.");

        }

        for (int i = 0; i < data.length; i++)
            try {
                if (datatype == ExternalDataSource.BYTE) {
                    out.writeByte((byte) data[i]);
                }
                if (datatype == ExternalDataSource.UBYTE) {
                    
                    int value = (int)data[i];


                    if (value < 0) {
                        value = 0;
                    }
                    else if (value > 255) {
                        value = 255;
                    }

                    // last eight bits
                    out.writeByte((byte)(value & 0xff));
                }
                else if (datatype == ExternalDataSource.SHORT) {
                    out.writeShort((short) data[i]);
                }
		else if (datatype == ExternalDataSource.USHORT) {
		    int val = (int)data[i];
	    
		    int maxUShort = (2 * Short.MAX_VALUE + 1) - 1;
		    
		    if (val > maxUShort) {
			val = maxUShort;
		    }
		    if (val < 0) {
			val = 0;
		    }
		    
		    out.writeShort((short)val);
		    
		}
                else if (datatype == ExternalDataSource.INT) {
                    out.writeInt((int) data[i]);
                }
		else if (datatype == ExternalDataSource.UINT) {
		    long val = (long)(data[i]);
		    
		    long maxUInt = 2L * (Integer.MAX_VALUE + 1L) - 1;
		    
		    if (val > maxUInt) {
			val = maxUInt;
		    }
		    if (val < 0L) {
			val = 0L;
		    }
		    
		    out.writeInt((int)(val));
		}
                else if (datatype == ExternalDataSource.LONG) {
                    out.writeLong((long) data[i]);
                }
                else if (datatype == ExternalDataSource.FLOAT) {
                    out.writeFloat((float) data[i]);
                }
                else {
                    out.writeDouble(data[i]);
                }
            }
            catch (Exception e) {
                throw new LoggedException(e);
            }
    }


    /**
     * Returns the output stream.
     */
    public DataOutputStream getOutputStream() {
        return out;
    }


    /**
     * Closes and tidies up.
     */
    public void close() {
        if (out != null) try {
            out.close();
        }
        catch(IOException e) {
            LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
        }

        if (useOutputImage) {
            ImageHeader outputHdr = CL_Initializer.headerTemplate.writeVectorImage(outputArray, ImageHeader.getFileRoot(outputFile));
        }
    }


    /**
     * Specifies that the output should be stored for collection
     * by an external process, such as Matlab.
     *
     * @param x X-size of output array
     *
     * @param y Y-size of output array
     *
     * @param z Z-size of output array
     *
     */
    public void setOutputArray(int x, int y, int z) {
        useOutputArray = true;
        outputArrayX = x;
        outputArrayY = y;
        outputArrayZ = z;

        outputArray = new double[outputArrayX][outputArrayY][outputArrayZ][0];
        nextX = 0;
        nextY = 0;
        nextZ = 0;
    }
    
    
    /**
     * Adds another voxel's worth of data to the output
     * array.
     *
     * @param data
     *            The array of data values to output.
     *
     */
    private void addToOutputArray(double[] data) {

        if(nextZ >= outputArrayZ)
            throw new LoggedException("Output array exhausted.");

        outputArray[nextX][nextY][nextZ] = data;
        nextX+=1;
        if(nextX==outputArrayX) {
            nextY+=1;
            nextX=0;
        }
        if(nextY==outputArrayY) {
            nextZ+=1;
            nextY=0;
        }
    }


    /**
     * Returns the output array.  Note that the array is in scanner
     * rather than voxel order so the last index specifies individual
     * components of one voxel.
     */
    public double[][][][] getOutputArray() {
        return outputArray;
    }

}


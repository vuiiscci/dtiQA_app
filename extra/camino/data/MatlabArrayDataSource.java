package data;

import misc.LoggedException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Base class for data sources of data from files and input streams.
 * 
 * <dt>Description:
 * 
 * <dd>Contains some basic methods for controlling data types.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: ExternalDataSource.java 176 2007-02-16 23:01:28Z ucacpco $
 */

public class MatlabArrayDataSource implements DataSource {

    // Size of array
    protected int xSize;
    protected int ySize;
    protected int zSize;
    protected int components;

    // Data array
    protected double[][][][] data;

    // Current indices
    protected int nextX;
    protected int nextY;
    protected int nextZ;


    /**
     * Constructor receives a data array from matlab process.  Note
     * that the array is in scanner rather than voxel order.
     *
     * @param matlabData Input data array.
     */
    public MatlabArrayDataSource(double[][][][] matlabData) {
        data = matlabData;
        xSize = data.length;
        ySize = data[0].length;
        zSize = data[0][0].length;
        components = data[0][0][0].length;

        nextX = 0;
        nextY = 0;
        nextZ = 0;

    }


    public double[] nextVoxel() throws DataSourceException {

        if(!more()) {
            throw new DataSourceException("No more voxels.");
        }

        double[] nextVox = data[nextX][nextY][nextZ];
        
        // Update next index
        nextX += 1;
        if(nextX >= xSize) {
            nextY += 1;
            nextX = 0;

            if(nextY >= ySize) {
                nextZ += 1;
                nextY = 0;
            }
        }

        return nextVox;
    }       


    public boolean more() {
        return (nextZ < zSize);
    }


}

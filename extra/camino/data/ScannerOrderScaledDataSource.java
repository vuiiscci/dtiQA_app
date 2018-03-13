package data;

import java.io.*;

import misc.*;
import tools.*;

/**
 *
 * Like a regular ScannerOrderDataSource, except that it scales the input data. Useful for
 * reading Analyze / NIFTI data.
 * 
 * @author Philip Cook
 * @version $Id$
 *  
 */
public class ScannerOrderScaledDataSource extends ScannerOrderDataSource {

    private double scaleSlope;

    private double scaleInter;

  
    /**
     * Constructor requires the filename, the number of values in each voxel and
     * the data type. If the filename is null, the object reads from the
     * standard input.
     * 
     * @param filename
     *            The name of the data file.
     * 
     * @param numVoxels
     *             The number of voxels in the file.            
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
     *
     */
    public ScannerOrderScaledDataSource(String filename, int numVoxels, int components, String type, boolean intelByteOrder, int offset,
				      double scaleSlope, double scaleInter) {
	
	super(filename, numVoxels, components, type, intelByteOrder, offset);

	this.scaleSlope = scaleSlope;
	this.scaleInter = scaleInter;
    }
    

    public double[] nextVoxel() throws DataSourceException {
	
	double[] voxel = super.nextVoxel();
	
	for (int i = 0; i < voxel.length; i++) {
	    voxel[i] = voxel[i] * scaleSlope + scaleInter;
	}

	return voxel;

    }


}

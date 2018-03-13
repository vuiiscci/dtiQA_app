package data;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General interface for source of multi-component voxel data, e.g.,
 * diffusion-weighted image data.
 * 
 * <dt>Description:
 * 
 * <dd>Has a single public interface method that returns a full list of the
 * measurements in the next voxel.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public interface DataSource {

    /**
     * Returns the data in the next voxel.
     * 
     * @return An array of values in the next voxel.
     */
    public double[] nextVoxel() throws DataSourceException;

    /**
     * Tests whether there are unprocessed voxels in the data source.
     * 
     * @return Boolean result of the test.
     */
    public boolean more();

}
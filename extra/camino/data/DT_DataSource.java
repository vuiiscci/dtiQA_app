package data;


/**
 *
 * Produces Camino DT data from a six-component DT input source. Much like calling nii2dt on the command line,
 * only this way users can use their actual DT images instead of going through that conversion step.
 * 
 */
public class DT_DataSource implements DataSource {

    private final boolean swap;

    private final DataSource data;


    /**
     * Construct a DT_DataSource from a tensor input source.
     *
     * @param rawData the tensor components as they are read from the file.
     *
     * @param lowerTriangular if true, re-ordering is required.
     *
     * 
     *
     */
    public DT_DataSource(DataSource rawData, boolean lowerTriangular) {
        
        swap = lowerTriangular;

        data = rawData;
        
    }
    




    /**
     * Returns the next voxel of tensor data. the exit code and lnS0 is always 0.
     * 
     * 
     * @return DT data in Camino format.
     */
    public double[] nextVoxel() throws DataSourceException {

        double[] next = new double[8];

        double[] input = data.nextVoxel();

        next[2] = input[0];
        next[3] = input[1];

        if (swap) {
            next[5] = input[2];
            next[4] = input[3];
        }
        else {
            next[4] = input[2];
            next[5] = input[3];
        }

        next[6] = input[4];
        next[7] = input[5];

        return next;
        
    }


    /**
     * Tests whether there are unprocessed voxels in the data source.
     * 
     * @return Boolean result of the test.
     */
    public boolean more() {
        return data.more();
    }






}
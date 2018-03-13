package data;

import imaging.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Source of voxel order diffusion-weighted MRI data from a data file or
 * standard input.
 * 
 * <dt>Description:
 * 
 * <dd>Doesn't do anything right now except read the data. It gets the number of components
 * in each voxel from the scheme file. This class can probably be deleted, since its original 
 * purpose was to support the now deprecated ignore measurements option.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: VoxelOrderDW_DataSource.java,v 1.5 2005/08/18 10:59:36 ucacmgh
 *          Exp $
 */
public class VoxelOrderDW_DataSource extends VoxelOrderDataSource {

    /**
     * The imaging parameters of the sequence used to acquire the data.
     */
    private DW_Scheme ip;

    /**
     * Constructor requires the filename, the number of values in each voxel and
     * the data type. If the filename is null, the input stream provides the
     * data.
     * 
     * @param filename
     *            The name of the data file.
     * 
     * @param type
     *            A string indicating the data type: either "char", "short",
     *            "int", "long", "float" or "double".
     * 
     * @param imParams
     *            The sequence defining the DW data.
     */
    public VoxelOrderDW_DataSource(String filename, String type,
            DW_Scheme imParams) {
        ip = imParams;
        int components = ip.numMeasurements();

	initFileInput(filename, false, 0);
	
        init(components, type);
    }

  
}

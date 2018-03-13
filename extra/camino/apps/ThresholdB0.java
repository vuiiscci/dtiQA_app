package apps;

import java.io.*;
import java.util.logging.Logger;

import tools.*;
import imaging.*;
import data.*;
import misc.LoggedException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Creates a background mask.
 * 
 * <dt>Description:
 * 
 * <dd> The program thresholds the average b=0 measurements to create
 * a background mask.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: ModelFit.java 58 2006-10-26 11:38:18Z ucacmgh $
 *  
 */
public class ThresholdB0 {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.ThresholdB0");

    /**
     * Output manager
     */
    private static OutputManager om;


    public static void main(String[] args) {

        // Default output format is char.
        OutputManager.outputDataType = "char";

        // Parse the command line arguments
        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
        CL_Initializer.initImagingScheme();
        CL_Initializer.initDataSynthesizer();

        om = new OutputManager();

        // Loop over the data
        int voxelNumber = 0;
        double[] outArr = new double[1];
        while (CL_Initializer.data.more())
            try {

                double[] nextVoxel = CL_Initializer.data.nextVoxel();

                // Fit or output background default.
                double backgroundB0 = CL_Initializer.imPars.geoMeanZeroMeas(nextVoxel);
                if(ModelFit.isBG(backgroundB0) || 
		   (CL_Initializer.CSFTHRESHOLD > 0.0 && backgroundB0 > CL_Initializer.CSFTHRESHOLD)) {
                    outArr[0] = 0;
                    om.output(outArr);

                }
                else {
                    outArr[0] = 1;
                    om.output(outArr);
                }
                voxelNumber++;
                logger.fine("Completed voxel: " + voxelNumber);

            }
            catch (DataSourceException e) {
                throw new LoggedException("The data file does not contain a whole number of voxels." +
                                          "Check the scheme file. Got Exception " + e);
            }

        om.close();
    }


}

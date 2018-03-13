package apps;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Date;
import java.util.logging.Logger;

import sphfunc.MaxEntProfile;
import data.*;
import tools.*;
import mesd.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Runs maximum entropy spherical deconvolution on each voxel's worth of
 * data on the input stream.
 * 
 * <dt>Description:
 * 
 * <dd>The output in each voxel is {exitcode, ln \bar{A}^\star(0), \lambda_0,
 * ..., \lambda_N}.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class MESD {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.MESD");


    /**
     * Output manager
     */
    private static OutputManager om;


    public static void main(String[] args) {

        int counter = 0;

        // Parse the command line arguments
        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
        CL_Initializer.initImagingScheme();

        CL_Initializer.initMaxEnt();
        CL_Initializer.initDataSynthesizer();

        om = new OutputManager();

        // Choose the inversion to run.
        MESD_Inversion inv = new MESD_Inversion(CL_Initializer.imPars,
                CL_Initializer.kernelParams, MaxEntProfile.getNumParams());

        // Loop over the data
        while (CL_Initializer.data.more())
            try {

                double[] nextVoxel = CL_Initializer.data.nextVoxel();


                // Fit or output background default.
                double backgroundB0 = CL_Initializer.imPars.geoMeanZeroMeas(nextVoxel);
                if (ModelFit.isBG(backgroundB0)) {
                    double[] voxOut = new double[inv.itemsPerVoxel()];

                    // Set the exitcode to -1 to indicate background.
                    voxOut[0] = -1;
                    voxOut[1] = (backgroundB0<=0)?0:Math.log(backgroundB0);

                    om.output(voxOut);
                }
                else {

                    // Fit the model and output the result.
                    double[] fittedData = inv.invert(nextVoxel);
                                       
                    om.output(fittedData);
                }

            }
            catch (Exception e) {
                throw new RuntimeException(e);
            }

        om.close();
    }

}

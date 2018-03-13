package apps;

import java.io.*;
import java.util.logging.Logger;

import tools.CL_Initializer;
import numerics.*;
import data.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Converts the spherical harmonic (SH) output of several third party programs to the format
 * expected in Camino.
 * 
 * <dt>Description:
 * 
 * <dd>Reorders and outputs SH coefficients into the format expected by Camino.
 * Camino uses the following ordering:
 * 
 *  [c00, c20, Re(c21), Im(c21), Re(c22), Im(c22), c40, Re(c41), Im(c41), ...].
 *
 * This program currently supports Tournier's Spherical Deconvolution
 * programs (both Matlab and C).
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: ImportSH_Coeffst.java 333 2007-08-01 12:02:10Z ucacdxa $
 *  
 */
public class ImportSH_Coeffs {

    /**
     * logger object
     */
    private static Logger logger = Logger.getLogger("camino.apps.ImportSH_Coeffs");

    /**
     * The maximum order of the spherical harmonic series in the input data.
     */
    private static int maxOrder = 8;

    /**
     * The type of the input data.
     */
    private static String inputDataType = "float";

    /**
     * The type of the output data.
     */
    private static String outputDataType = "double";

    /**
     * Input data file. If left null, the program reads data from the standard
     * input.
     */
    private static String inputFile = null;

    /**
     * The stream to which to write the output.
     */
    private static DataOutputStream out = new DataOutputStream(System.out);
    
    /**
     * Used by coeffSource.  Input data from Tournier's Spherical Deconvolution code.
     */
    private static final int TOURNIER=1;
    /**
     * The ordering of the input data.  The options are:
     * 
     * tournier - files from Tournier's Spherical Deconvolution code
     * 
     */
    private static int coeffSource = TOURNIER;

    /**
     * Output manager
     */
    private static OutputManager om;


    public static void main(String[] args) {

        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
	
	// parse technique specific args
	for(int i=0; i<args.length; i++) {
	    if(args[i].equals("-source")) {
		if(args[i+1].equalsIgnoreCase("tournier")){
		    coeffSource=TOURNIER;
		}
		   else{
		       logger.warning("-source value: " + args[i+1] + " not recognized.  Defaulting to \"tournier\" conversion");
		       coeffSource=TOURNIER;
		   }
		CL_Initializer.markAsParsed(i, 2);
	    }
	}

	maxOrder=CL_Initializer.maxOrder;
	inputDataType =CL_Initializer.inputDataType;
	
        om = new OutputManager();

        // Construct the data source.
        DataSource data = ExternalDataSource.getDataSource(inputFile,
                SphericalHarmonics.evenFuncsUpTo(maxOrder), inputDataType);

	if(coeffSource==TOURNIER) {
	    tournierConversion(data);
	}

        // Tidy up.
        om.close();

    }

    private static void tournierConversion(DataSource data) {
	// Loop over the data
        while (data.more()){
            try {

                // Read in the coefficients.
                double[] coeffs = data.nextVoxel();

                // The output data contains two extra values.
                double[] reordered = new double[coeffs.length + 2];

                // The first two values are an exitcode and ln A^\star(0).
                // We just put in dummy values.
                reordered[0] = 0.0;
                reordered[1] = 0.0;

                // Now reorder the coefficients.
                int i = 2;
                for (int l = 0; l <= maxOrder; l += 2) {
                    reordered[i] = coeffs[i - 2 + l];
                    for (int m = 1; m <= l; m++) {
                        reordered[i + 2 * m - 1] = coeffs[i - 2 + l + m];
                        reordered[i + 2 * m] = -coeffs[i - 2 + l - m];
                    }
                    i += 2 * l + 1;
                }

                om.output(reordered);

            }
            catch (Exception e) {
                System.err.println(e);
            }
	}
    }
}

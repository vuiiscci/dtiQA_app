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
 * <dd>Converts the output of Donald Tournier's sig2fod program to the format
 * expected in the Cameno suite.
 * 
 * <dt>Description:
 * 
 * <dd>sig2fod outputs spherical harmonic coefficients in the following order:
 * [c00, -Im(c22), -Im(c21), c20, Re(c21), Re(c22), c40, -Im(c44), ...]. The
 * Cameno tools use the following ordering: [c00, c20, Re(c21), Im(c21),
 * Re(c22), Im(c22), c40, Re(c41), Im(c41), ...]. This program reads input in
 * the former format and outputs the same data in the latter format.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class DeconvToCamenoFormat {

    /**
     * logger object
     */
    private static Logger logger = Logger.getLogger("camino.apps.DeconvToCamenoFormat");

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
     * Output manager
     */
    private static OutputManager om;


    public static void main(String[] args) {

        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);

        om = new OutputManager();

        // Construct the data source.
        VoxelOrderDataSource data = new VoxelOrderDataSource(inputFile,
                SphericalHarmonics.evenFuncsUpTo(maxOrder), inputDataType);

        // Loop over the data
        while (data.more())
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

        // Tidy up.
        om.close();

    }


}

package apps;

import java.util.logging.Logger;

import tools.*;
import misc.*;
import data.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Write binary zeros (or any other number) to the standard
 * output.  The first value can optionally be set to a different value
 * using the command line option -exitcode <val> , for example -1,
 * which is the exit code indicating background.
 * 
 * <dt>Description:
 * 
 * <dd>Used to create files full of zeros, useful for creating data
 * files for background voxels for example.
 * 
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class WriteZeros {

    private static Logger logger = Logger.getLogger("camino.apps.SphHarmFitter");

    /**
     * Output manager
     */
    private static OutputManager om;


    public static void main(String[] args) {

        if(args.length < 1) {
            System.err.println("Usage: WriteZeros <num zeros> [ options ]");
            System.exit(1);
        }

        int components = Integer.parseInt(args[0]);
 
        // Check the number of components is sensible.
        if(components<=0) {
            throw new LoggedException("Number of components must be greater than 0.");
        }

        String[] arg2 = new String[args.length - 1];
        for(int i=0; i<args.length - 1; i++) {
            arg2[i] = args[i+1];
        }

        // Set defaults
        boolean flagAsBackground = false;
        double val = 0.0;
        double exitcode = 0.0;

        // Parse the command line arguments
        CL_Initializer.CL_init(arg2);

        for (int i = 0; i < arg2.length; i++) {
            if (arg2[i].equals("-exitcode")) {
                flagAsBackground = true;
                exitcode = Double.parseDouble(arg2[i+1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i+1);
           }
            if (arg2[i].equals("-value")) {
                val = Double.parseDouble(arg2[i+1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i+1);
           }
       }

        CL_Initializer.checkParsing(arg2);

        om = new OutputManager();

        // Construct the output array.
        double[] voxOut = new double[components];
        for(int i=0; i<voxOut.length; i++) {
            voxOut[i] = val;
        }
        if(flagAsBackground) {
            voxOut[0] = exitcode;
        }

        // Output the data
        om.output(voxOut);


        // Tidy up.
        om.close();
    }

}

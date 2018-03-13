package apps;

import java.io.*;
import java.util.logging.Logger;

import tools.CL_Initializer;
import data.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Converts scanner order data to the voxel order format used by the Camino
 * suite.
 * 
 * <dt>Description:
 * 
 * <dd>Reads the whole scanner-order data set into memory and outputs the data
 * in voxel order as big-endian floating point data.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class ScannerToVoxel extends Executable{

    public ScannerToVoxel(String[] args){
        super(args);
    }
     private static Logger logger = Logger.getLogger("camino.apps.ScannerToVoxel");

    /**
     * The number of measurements in each voxel.
     */
    private static int components;
	ScannerOrderDataSource data;

    /**
     * Output manager
     */
   // private static OutputManager om;

   
     public void initDefaultVals(){
		components = 1;
	 }

	public void initOptions(String[] args){
//    public static void main(String[] args) {

        // Set default input and output data types to float.
        CL_Initializer.inputDataType = "float";
        OutputManager.outputDataType = "float";

        CL_Initializer.CL_init(args);

        // Parse the command line arguments
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-components")) {
                components = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
        }

        CL_Initializer.checkParsing(args);
	}
	
	// set up output 
    //    om = new OutputManager();

	public void initVariables(){
        // Construct the data source.
        data = new ScannerOrderDataSource(CL_Initializer.inputFile, CL_Initializer.numVoxels,
                components, CL_Initializer.inputDataType);
	}
	

	public void execute(OutputManager om){
        // Loop over the data
        while (data.more())
            try {

                // Get the data from the next voxel.
                double[] reordered = data.nextVoxel();

                // Output it.
                om.output(reordered);

            }
            catch (Exception e) {
                logger.warning(e.toString() + "  (Program will continue)");
            }

        // Tidy up.
        om.close();
    }

}

package apps;

import java.io.*;
import java.util.logging.Logger;

import misc.LoggedException;

import tools.CL_Initializer;
import data.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Converts voxel order data to scanner order.
 * 
 * <dt>Description:
 * 
 * <dd>Reads the whole voxel-order data set into memory and outputs the data in
 * scanner order as big-endian floating point data.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class VoxelToScanner extends Executable{

    public VoxelToScanner(String[] args){
        super(args);
    }
	
    private static Logger logger = Logger.getLogger("camino.apps.VoxelToScanner");

    /**
     * The number of measurements in each voxel.
     */
    private static int components;
	VoxelOrderDataSource data;
	private double[][] d;
	private int nv;

    /**
     * Output manager
     */
    //private static OutputManager om;

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
  //      om = new OutputManager();

	public void initVariables(){
        // Construct the data source.
		data = new VoxelOrderDataSource(CL_Initializer.inputFile, components,
                CL_Initializer.inputDataType); 
 /*       VoxelOrderDataSource data = new VoxelOrderDataSource(CL_Initializer.inputFile, components,
                CL_Initializer.inputDataType); */

        // Create an array to store all the data.
        d = new double[components][CL_Initializer.numVoxels];
		nv = CL_Initializer.numVoxels;
	}
		
	public void execute(OutputManager om){
        // Read all the data into the array.
        //double[] voxData = new double[CL_Initializer.numVoxels];
		for (int i = 0; i < nv; i++)
            try {

                // Get the data from the next voxel.
                //voxData[i] = data.nextVoxel();
				double[] voxData = data.nextVoxel();

                // Copy into the data array.
                for (int j = 0; j < components; j++) {
                    d[j][i] = voxData[j];
                }

            }
            catch (Exception e) {
                LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
            }

        // Now output all the data in scanner order.
        for (int i = 0; i < components; i++) {
            om.output(d[i]);
        }

        // Tidy up.
        om.close();
    }

}

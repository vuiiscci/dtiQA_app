package apps;

import java.util.logging.Logger;

import misc.*;
import tools.*;
import data.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Creates a stream of synthetic data
 * 
 * <dt>Description:
 * 
 * <dd>Outputs a specified number of voxels worth of data synthesized from the
 * specified test function with the specified imaging parameters and noise
 * level.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class SyntheticData extends Executable {

    private static Logger logger = Logger.getLogger("camino.apps.SyntheticData");
    
    
    public SyntheticData(String[] args) {
        super(args);
    }

    public void initDefaultVals() {

    }

    public void initOptions(String[] args) {

        // The output defaults to float type.
        OutputManager.outputDataType = "float";
        // The input defaults to double type.
        CL_Initializer.inputDataType = "double";
        
        // Parse the command line arguments
        CL_Initializer.CL_init(args);
        
        CL_Initializer.checkParsing(args);
        
        // Now initialize the acquisition scheme and model parameters.
        CL_Initializer.initImagingScheme();
        CL_Initializer.initDataSynthesizer();
     
    }

    public void initVariables() {

    }

    public void execute(OutputManager om) {

        // Check that a test function was specified.
        if (CL_Initializer.p == null && CL_Initializer.inputModel == null && CL_Initializer.brownianSimulation==false 
            && (CL_Initializer.bootstrap < 0) && (CL_Initializer.compartmentModel==false)) {
            throw new LoggedException("No test function or input model or input data specified.");
        }

        // Loop over the data outputting it
        while (CL_Initializer.data.more())
            try {
                om.output(CL_Initializer.data.nextVoxel());
                
            }
            catch (DataSourceException e) {
                throw new LoggedException(e);
            }

        om.close();
    }

}

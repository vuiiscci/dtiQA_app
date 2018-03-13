package apps;

import java.util.logging.Logger;

import misc.LoggedException;

import data.*;
import tools.*;
import inverters.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Classifies each voxel as background, isotropic, anisotropic Gaussian or
 * non-Gaussian.
 * 
 * <dt>Description:
 * 
 * <dd>Computes the spherical harmonic expansion of log(S(q)) in each voxel and
 * uses an F-test to determine the simplest model that gives a good fit to the
 * data.
 * 
 * Uses standard input and output streams for input and output data.
 * 
 * If no f-test thresholds are specified on the command line, the program
 * outputs all the f-statistics in each voxel for independent processing.
 * 
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class VoxelClassify extends Executable{

    private static Logger logger = Logger.getLogger("camino.apps.VoxelClassify");

    /**
     * Output manager
     */
    //private static OutputManager om;

	public VoxelClassify(String[] args){
    	super(args);
    }

	private EvenSphHarmFitter eshf;

	public void initDefaultVals(){
	}
	
	public void initOptions(String[] args){
		// Parse the command line arguments
        //CL_Initializer cl = new CL_Initializer(args);
        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
        CL_Initializer.initImagingScheme();
        CL_Initializer.initDataSynthesizer();
	}
    
	public void initVariables(){
        // If f-test thresholds are specified, we write out
        // integer classifications.
        if(CL_Initializer.f1 > 0.0)
            OutputManager.outputDataType = "int";
			
	     // Create the spherical harmonic fitter object.
		 eshf = new EvenSphHarmFitter(CL_Initializer.imPars, CL_Initializer.maxOrder);
	}
	
	public void execute(OutputManager om){	
   // public static void main(String[] args) {

        //om = new OutputManager();


        // Loop over the data
        while (CL_Initializer.data.more())
            try {

                // Fit the models and compute the f-statistics.
                double[] voxelData = CL_Initializer.data.nextVoxel();
                double[] coeffs = eshf.fit(voxelData);

                double exitCode = coeffs[0];
                double q0 = Math.exp(coeffs[1]);
                double[] ps = eshf.getF_TestProbabilities(coeffs, voxelData);

                // Write out the desired result.

                // If f-tests were specified, use them to do the model
                // selection.
                if (CL_Initializer.f1 > 0.0)
                    try {
                        int val;
                        // Apply the background and CSF thresholds first.
                        if (ModelFit.isBG(q0) || exitCode < 0.0) {
                            val = -1;
                        }
                        else if (CL_Initializer.CSFTHRESHOLD > 0.0
                                && q0 > CL_Initializer.CSFTHRESHOLD) {
                            val = 0;
                        }
                        else {
                            val = eshf.selectModel(ps,
                                    CL_Initializer.f1, CL_Initializer.f2,
                                    CL_Initializer.f3);
                        }

                        om.output(new double[] {val});
                    }
                    catch (Exception e) {
                        LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
                    }

                // Otherwise output the exit code, the average q=0 and all
                // the f-statistics.
                else {
                    // (PAC) respect the outputDataType
                    // CL_Initializer.out.writeDouble(coeffs[0]);
                    // CL_Initializer.out.writeDouble(coeffs[1]);

                    double[] coeffOut = new double[2];

                    coeffOut[0] = coeffs[0];
                    coeffOut[1] = coeffs[1];

                    om.output(coeffOut);
                    om.output(ps);
                }

            }
            catch (Exception e) {
                StackTraceElement[] stackTrace = e.getStackTrace();
                String stString = new String();

                // log the exceptions message
                logger.warning(e.toString() + "  (Program will continue)");

                // log the stack trace
                for (int i = 0; i < stackTrace.length; i++) {
                    stString += stackTrace[i] + "\n";
                }
                logger.warning(stString);

                //System.err.println("WARNING: "+e);
            }

        // Tidy up.
        om.close();
    }

}

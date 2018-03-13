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
 * <dd>Computes the trace of each diffusion tensor.
 * 
 * <dt>Description:
 * 
 * <dd>Uses standard input and output streams for input and output data.
 * 
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class TraceD extends Executable{

    private static Logger logger = Logger.getLogger("camino.apps.TraceD");

    /**
     * Output manager
     */
    //private static OutputManager om;

    public TraceD(String[] args){
    	super(args);
    }
	
	public void initDefaultVals(){
	
	}
	
	public void initOptions(String[] args){
//   public static void main(String[] args) {

        // Input type defaults to single DT. Input datatype defaults to
        // double.
        CL_Initializer.inputModel = "dt";
        CL_Initializer.inputDataType = "double";

        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);

        // This program only deals with tensor data.
        CL_Initializer.initTensorDataSource();

	}
	
	public void initVariables(){
	}
	
	public void execute(OutputManager om){	
 
//        om = new OutputManager();

        // Loop over the data
        while (CL_Initializer.data.more())
            try {

                double[] vox = CL_Initializer.data.nextVoxel();
                DT[] dts = FracAnis.getTensorList(vox, CL_Initializer.inputModel);
                double[] trs = new double[dts.length];

                // Check if the voxel is background and if so just
                // output zero.  Check the exit code in the input data
                // as well as the background threshold.
		boolean notBG = !ModelFit.isBG(Math.exp(vox[1]));

                if (vox[0]>=0 && notBG) {
                    for (int i = 0; i < dts.length; i++) {
                        trs[i] = dts[i].trace();
                    }
                }
                else {
                    for (int i = 0; i < dts.length; i++) {
                        trs[i] = 0.0;
                    }
                }

                om.output(trs);

            }
            catch (Exception e) {
                
                LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
            }

        // Tidy up.
        om.close();
    }

}

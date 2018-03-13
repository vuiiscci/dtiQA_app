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
 * <dd>Computes the eigensystem of each diffusion tensor.
 * 
 * <dt>Description:
 * 
 * <dd>Uses standard input and output streams for input and output data. The
 * output in each voxel is: {\lambda_1, e_1x, e_1y, e_1z, \lambda_2, e_2x, e_2y,
 * e_2z, \lambda_3, e_3x, e_3y, e_3z}, where \lambda_1 >= \lambda_2 >=
 * \lambda_3, for each DT.
 * 
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class DT_EigenSystem extends Executable{

    private static Logger logger = Logger.getLogger("camino.apps.DT_EigenSystem");

    /**
     * Output manager
     */
    //private static OutputManager om;

	public DT_EigenSystem(String[] args){
    	super(args);
    }
	
	public void initDefaultVals(){
	
	}
	
	public void initOptions(String[] args){
        // Input type defaults to single DT. Input datatype defaults to
        // double.
        CL_Initializer.inputModel = new String("dt");
        CL_Initializer.inputDataType = new String("double");

        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
        // This program only deals with tensor data.
        CL_Initializer.initTensorDataSource();
	}
	
	public void initVariables(){
	}
	
	public void execute(OutputManager om){	
//    public static void main(String[] args) {


	// set up output 
//	om = new OutputManager();

        // Loop over the data
        while (CL_Initializer.data.more()) {
            double[] vox = CL_Initializer.data.nextVoxel();
            DT[] dts = FracAnis.getTensorList(vox, CL_Initializer.inputModel);
            
            // Check if the voxel is background and if so just
            // output zero.  Check the exit code in the input data
            // as well as the background threshold.

			boolean notBG = !ModelFit.isBG(Math.exp(vox[1]));

            if (vox[0]>=0 && notBG) {
                for (int i = 0; i < dts.length; i++) {
                    try {
                        double[][] sEig = dts[i].sortedEigenSystem();
                        
                        // Flatten into an array
                        double[] eigFlat = new double[12];
                        for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 4; k++) {
                            eigFlat[4 * j + k] = sEig[k][j];
                        }
                        }
                        // Output.
                        om.output(eigFlat);
                    }
                    catch(Exception e) {
                        LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
                        double[] zeros = new double[12];
                        // output zeros
                        om.output(zeros); 
                    }
                    
                }
            }
            else {
                
                om.output(new double[dts.length * 12]);
            }
            
        }
        
        om.close();

    }

}

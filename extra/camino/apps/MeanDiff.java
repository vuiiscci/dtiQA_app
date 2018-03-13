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
 * <dd>Computes the mean diffusivity (trace/3) of each diffusion tensor.
 *
 * <dt>Description:
 *
 * <dd>Uses standard input and output streams for input and output data.
 *
 *
 * </dl>
 *
 * @author Danny Alexander
 * @version $Id: TraceD.java 58 2006-10-26 11:38:18Z ucacmgh $
 *
 */
public class MeanDiff extends Executable{
    
    public MeanDiff(String[] args){
        super(args);
    }
    
    private static Logger logger = Logger.getLogger("camino.apps.TraceD");
    
    /**
     * Output manager
     */
    //private static OutputManager om;
    
    
    //  public static void main(String[] args) {
    
    public void initDefaultVals(){
    }
    
    public void initOptions(String[] args){
        
        // Input type defaults to single DT. Input datatype defaults to
        // double.
        CL_Initializer.inputModel = "dt";
        CL_Initializer.inputDataType = "double";
        
        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
        
        // This program only deals with tensor data.
        CL_Initializer.initTensorDataSource();
    }
    //om = new OutputManager();
    
    public void initVariables(){
    }
    
    public void execute(OutputManager om){
        // Loop over the data
        while (CL_Initializer.data.more())
            try {
                
                double[] vox = CL_Initializer.data.nextVoxel();
                DT[] dts = FracAnis.getTensorList(vox, CL_Initializer.inputModel);
                double[] mds = new double[dts.length];
                
                // Check if the voxel is background and if so just
                // output zero.  Check the exit code in the input data
                // as well as the background threshold.
                boolean notBG = !ModelFit.isBG(Math.exp(vox[1]));
                
                if (vox[0]>=0 && notBG) {
                    for (int i = 0; i < dts.length; i++) {
                        mds[i] = dts[i].trace()/3.0;
                    }
                }
                else {
                    for (int i = 0; i < dts.length; i++) {
                        mds[i] = 0.0;
                    }
                }
                
                om.output(mds);
                
            }
            catch (Exception e) {
                
                LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
            }
        
        // Tidy up.
        om.close();
    }
    
}

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
 * <dd>Computes the fractional anisotropy of each diffusion tensor.
 *
 * <dt>Description:
 *
 * <dd>For a diffusion tensor D, the fractional anisotropy is
 *
 * \fracanis = \left(\frac32 \sum_{i=1}^3 \left(\dteval_i - \frac13
 * \trace(\difften) \right)^2 \right)^\frac12 \left(\sum_{i=1}^3 \dteval_i^2
 * \right)^{-\frac12}.
 *
 *
 *
 * </dl>
 *
 * @author Danny Alexander
 * @version $Id$
 *
 */
public class FracAnis extends Executable{
    
    private static Logger logger = Logger.getLogger("camino.apps.FracAnis");
    
    
    /**
     * Output manager
     */
    //private static OutputManager om;
    
    public FracAnis(String[] args){
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
        
        // set up output
        // om = new OutputManager();
        
        
        // Loop over the data
        while (CL_Initializer.data.more())
            try {
                
                double[] vox = CL_Initializer.data.nextVoxel();
                DT[] dts = getTensorList(vox, CL_Initializer.inputModel);
                double[] fas = new double[dts.length];
                
                // Check if the voxel is background and if so just
                // output zero.  Check the exit code in the input data
                // as well as the background threshold.
                boolean notBG = !ModelFit.isBG(Math.exp(vox[1]));
                
                if (vox[0]>=0 && notBG) {
                    for (int i = 0; i < dts.length; i++) {
                        fas[i] = dts[i].fa();
						//System.err.println(fas[i]);
                    }
                }
                else {
                    for (int i = 0; i < dts.length; i++) {
                        fas[i] = 0.0;
                    }
                }
                
                om.output(fas);
                
            }
            catch (Exception e) {
                
                LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
                
            }
        
        // Tidy up.
        om.close();
    }
    
    /**
     * Constructs a list of diffusion tensors from the input data in a single
     * voxel.
     *
     * @param vox
     *            The data.
     *
     * @param model
     *            The model.
     *
     * @return The array of diffusion tensors.
     */
    public static DT[] getTensorList(double[] vox, String model) {
        //public DT[] getTensorList(double[] vox, String model) {
        
        DT[] dts = null;
        
        if (model.equals("dt")) {
            dts = new DT[1];
            dts[0] = new DT(vox[2], vox[3], vox[4], vox[5], vox[6], vox[7]);
        }
        else if (model.equals("twotensor")) {
            dts = new DT[2];
            dts[0] = new DT(vox[4], vox[5], vox[6], vox[7], vox[8], vox[9]);
            dts[1] = new DT(vox[11], vox[12], vox[13], vox[14], vox[15], vox[16]);
        }
        else if (model.equals("threetensor")) {
            dts = new DT[3];
            dts[0] = new DT(vox[4], vox[5], vox[6], vox[7], vox[8], vox[9]);
            dts[1] = new DT(vox[11], vox[12], vox[13], vox[14], vox[15], vox[16]);
            dts[2] = new DT(vox[18], vox[19], vox[20], vox[21], vox[22], vox[23]);
        }
        else if (model.equals("multitensor")) {
            dts = new DT[CL_Initializer.maxTensorComponents];
            for(int i=0; i<dts.length; i++) {
                dts[i] = new DT(vox[4+7*i], vox[5+7*i], vox[6+7*i], vox[7+7*i], vox[8+7*i], vox[9+7*i]);
            }
        }
        else {
            throw new RuntimeException("Cannot extract tensors from model: " + model);
        }
        
        return dts;
    }
    
}

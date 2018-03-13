package apps;

import java.util.logging.Logger;

import misc.LoggedException;

import tools.*;
import data.*;
import sphfunc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Computes the anisotropy of an input stream of spherical functions.
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
public class SphFuncAnisotropy extends Executable{

    private static Logger logger = Logger.getLogger("camino.apps.SphFuncAnisotropy");


    /**
     * Output manager
     */
    //private static OutputManager om;

	public SphFuncAnisotropy(String[] args){
    	super(args);
    }

	public void initDefaultVals(){
	}
	
	public void initOptions(String[] args){
        // Parse the command line arguments
        CL_Initializer.inputDataType = "double";
        //CL_Initializer cl = new CL_Initializer(args);
        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);

        CL_Initializer.initImagingScheme();
        CL_Initializer.initMaxEnt();
        CL_Initializer.initSphFuncDataSource();	
	}
	
	public void initVariables(){
	}
	
	public void execute(OutputManager om){	
//    public static void main(String[] args) {

//        om = new OutputManager();

        // Loop over the data
        while (CL_Initializer.data.more())
            try {

                // Read in the coefficients.
                double[] coeffs = CL_Initializer.data.nextVoxel();

                // Construct a spherical function using the
                // coefficients.
                SphericalFunction sf = null;
                if (CL_Initializer.inputModel.equals("sh")) {
                    sf = new EvenSHS(coeffs, CL_Initializer.maxOrder);
                }
                else if (CL_Initializer.inputModel.equals("maxent")) {
                    sf = new MaxEntProfile(coeffs, CL_Initializer.kernelParams);
                }
                else {
                    sf = new TuchRBF_Sum(coeffs);
                }

                double[] anis = {0.0};

		boolean notBG = !ModelFit.isBG(Math.exp(coeffs[1]));

                if (coeffs[0]>=0 && notBG) {
                    anis[0] = sf.anisotropy();
                }

                om.output(anis);

            }
            catch (Exception e) {
                
                LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
            }

        // Tidy up.
        om.close();
    }

}

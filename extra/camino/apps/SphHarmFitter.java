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
 * <dd>Computes the anisotropy based on the spherical harmonic fit.
 * 
 * <dt>Description:
 * 
 * <dd>Computes the spherical harmonic expansion of log(S(q)) in each voxel and
 * outputs the coefficients of the series up to the specified order.
 * 
 * Uses standard input and output streams for input and output data.
 * 
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class SphHarmFitter extends Executable{

    private static Logger logger = Logger.getLogger("camino.apps.SphHarmFitter");

    /**
     * Output manager
     */
    //private static OutputManager om;
	private EvenSphHarmFitter eshf;
	
    public SphHarmFitter(String[] args){
    	super(args);
    }
	
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
       // Create the spherical harmonic fitter object.
        eshf = new EvenSphHarmFitter(CL_Initializer.imPars, CL_Initializer.maxOrder);
	}
	
	public void execute(OutputManager om){	
//    public static void main(String[] args) {
 //       om = new OutputManager();
 
        // Loop over the data
        while (CL_Initializer.data.more())
            try {

                double[] nextVoxel = CL_Initializer.data.nextVoxel();
                double[] nextGradAdj = {};
                if (CL_Initializer.gradAdj!= null) {
                        nextGradAdj = CL_Initializer.gradAdj.nextVoxel();
                }
                // Fit or output background default.
                double backgroundB0 = CL_Initializer.imPars.geoMeanZeroMeas(nextVoxel);
                if (ModelFit.isBG(backgroundB0)) {
                    double[] voxOut = new double[eshf.itemsPerVoxel()];

                    // Set the exitcode to -1 to indicate background.
                    voxOut[0] = -1;
                    voxOut[1] = Math.log(backgroundB0);

                    om.output(voxOut);
                }
                else {
                    // Fit and output the result.
                    double[] coeffs;
                    if (CL_Initializer.gradAdj!=null) {
                        // if gradient adjust then fit with gradadj
                        coeffs = eshf.fit(nextVoxel, nextGradAdj);

                    }
                    else{
                        // else fit normally
                        coeffs = eshf.fit(nextVoxel);
                    }
                    om.output(coeffs);
                }
            }
            catch (Exception e) {
                
                throw new LoggedException(e);
            }

        // Tidy up.
        om.close();
    }

    /**
     * Computes the spherical harmonic anisotropy.
     * 
     * @param coeffs
     *            The array returned by the array returned by
     *            EvenSphHarmFitter.fit.
     * 
     * @return The anisotropy.
     */
    //public static double getSphHarmAnisotropy(double[] coeffs) {
	public double getSphHarmAnisotropy(double[] coeffs) {

        // The first two elements of coeffs are the exit code and
        // log(A^\star(0)).
        int ind = 2;

        // Start at order zero.
        int ord = 0;

        // Initialize the anisotropy to zero.
        double anis = 0.0;

        // Keep going until we run out of coefficients. We do not need
        // to figure out the order of the series explicitly.
        while (ind < coeffs.length) {
            anis += coeffs[ind] * coeffs[ind];
            ind += 1;
            for (int i = 0; i < 2 * ord; i++) {
                anis += 2.0 * coeffs[ind] * coeffs[ind];
                ind += 1;
            }

            // Move on to the next even order.
            ord += 2;
        }

        // Finally normalize by the first moment (sqrt(4 \pi c00))
        // squared and take the square root.
        anis = Math.sqrt((anis - coeffs[2] * coeffs[2])
                / (4.0 * Math.PI * coeffs[2] * coeffs[2]));

        return anis;

    }

}

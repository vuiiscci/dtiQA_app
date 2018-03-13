package apps;

import java.util.logging.Logger;

import data.*;
import tools.*;
import inverters.*;
import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Creates a matrix for fitting the diffusion tensor via a linear
 * transform.
 * 
 * <dt>Description:
 * 
 * <dd>Constructs the matrix required for linear fitting of the diffusion
 * tensor to log measurements and outputs it.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class DT_FitMatrix extends Executable{

	public DT_FitMatrix(String[] args){
		super(args);
	}
	
	//private LinearDT_Inversion inv;
	private RealMatrix m;

    /**
     * Logger object
     */
    private static Logger logger = Logger.getLogger("camino.apps.DT_FitMatrix");

    /**
     * Output manager
     */
   // private static OutputManager om;

    public void initDefaultVals() {
		//inv = null;
		m = null;
	}

	public void initOptions(String[] args){
        // Parse the command line arguments
        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
        CL_Initializer.initImagingScheme();		
	}
	
	public void initVariables() {
       // Create the linear inversion
       LinearDT_Inversion inv = new LinearDT_Inversion(CL_Initializer.imPars);

        // Get the matrix.
        m = inv.getMatrix();
	}
	
	public void execute(OutputManager om){

        // Output the matrix
        double[] row = new double[m.columns()];
        for (int j = 0; j < m.rows(); j++) {
            for (int i = 0; i < m.columns(); i++) {
                row[i] = m.entries[j][i];
            }
            om.output(row);
        }
        om.close();
    }

}

package apps;

import java.util.logging.Logger;

import data.*;
import tools.*;
import inverters.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Performs a linear reconstruction in each voxel.
 * 
 * <dt>Description:
 * 
 * <dd>The program reads in a linear transformation matrix and uses it to
 * transform the data in each voxel.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class LinearRecon extends ModelFit {

    public LinearRecon(String[] args){
        super(args);
    }
 
    public void initVariables(){
        if (CL_Initializer.matrixFile == null) {
            logger.severe("LinearRecon requires a matrix file.");
            System.exit(1);
        }

        inv = new LinearInversion(CL_Initializer.imPars,
		CL_Initializer.matrixFile, CL_Initializer.lrNormalize, CL_Initializer.lrLog);
 	}
	

}

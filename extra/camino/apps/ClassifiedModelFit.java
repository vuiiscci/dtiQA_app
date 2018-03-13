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
 * <dd>Classified model fitting program.
 * 
 * <dt>Description:
 * 
 * <dd>The program fits models to each voxel depending on a
 * precomputed classification.
 * 
 * </dl>
 * 
 * @see inverters.ModelIndex 
 *
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class ClassifiedModelFit extends ModelFit {



    public ClassifiedModelFit(String[] args){
        super(args);
    }


    public void initVariables() {

        inv = new MultiTensorInversion(CL_Initializer.imPars, CL_Initializer.voxelClassMap, CL_Initializer.classifiedModelIndices, CL_Initializer.maxTensorComponents);

    }


    public void initOptions(String[] args) {

        super.initOptions(args);

        // Check that the essential information is provided.
        if(CL_Initializer.voxelClassMap == null) {
            logger.severe("Need a voxel classification map.");
            System.exit(1);
        }

    }


}

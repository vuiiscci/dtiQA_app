package tractography;

import misc.*;

import java.util.Random;

/**
 * Data model types: "cylsymmdt", "ballsticks".
 */
public enum BayesDataModel {
    
    CYL_SYMM_DT("cylsymmdt"),
    BALL_STICK("ballstick");
    

    BayesDataModel(String modelName) {
	name = modelName;
    }
    
    
    public String toString() {
	return name;
    }
    

    public static BayesDataModel getModel(String s) {

	for (BayesDataModel model : BayesDataModel.values()) {
	    if (s.equals(model.name)) {
		return model;
	    }
	}
	
	throw new LoggedException("Unsupported data model type " + s);
    }
    
    
    /**
     * The name of the model.
     */
    public final String name;

    
}

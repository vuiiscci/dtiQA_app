package models.compartments;

import numerics.RealMatrix;
import misc.LoggedException;
import models.ParametricModel;
import imaging.DW_Scheme;
import imaging.SimulableScheme;


/**
 * Implements the compartment interface for an isotropic restricted compartment,
 * Spherical boundary compartment with zero radius.
 * 
 * 
 * @author laura (panagio@cs.ucl.ac.uk)
 *
 */
public class Dot extends ParametricModel {



	/** constructor needs array of params. 
	 * in this case, an array of 4 parameters, the difffusivity and the angles theta and phi
	 * that determine the fibre orientation and the radius R.
	 * 
	 * @param params array
	 */
	public Dot() {

		super(CompartmentType.DOT.numParams);  

	}
	
	
	/**
	 * return a single signal from the scheme from the dot model
	 * 
	 * @return 1.0
	 */
	public double getSignal(double[] params, DW_Scheme scheme, int i){
	    
	    return 1.0;
	}



	/**
	 * generates signals from this compartment for each line
	 * in the scheme given.
	 * 
	 * @param scheme scan specifics
	 */
	public RealMatrix getSignals(double[] params, DW_Scheme scheme) {

		RealMatrix signals = new RealMatrix(scheme.numMeasurements(), 1);

		for (int i = 0; i < scheme.numMeasurements(); i++) {
			// the restricted signal
			//signals[i] = 1;
		    signals.setEntry(i, 0, 1.0);
		}
		return signals;
	}
	

}
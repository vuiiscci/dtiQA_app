package models.compartments;

import misc.LoggedException;
import models.ParametricModel;
import imaging.DW_Scheme;
import imaging.SimulableScheme;
import imaging.StejskalTannerScheme;
import numerics.ErrorFunction;
import numerics.ErrorFunctionException;
import numerics.RealMatrix;

/**
 * Implements the compartment interface for a restricted isotropic compartment,
 * sticks with distributed orientations.
 * 
 * 
 * @author laura (panagio@cs.ucl.ac.uk)
 *
 */
public class Astrosticks extends ParametricModel {

	/** gyromagnetic ratio */
	private final double GAMMA = DW_Scheme.GAMMA;

	/** constructor. needs array of params. 
	 * in this case, an array of 1 parameter the diffusivity.
	 * 
	 * @param params array
	 */
	public Astrosticks() {

		super(CompartmentType.ASTROSTICKS.numParams);
	}

	/**
	 * generates signals from this compartment for each line
	 * in the scheme given.
	 * 
	 * @param scheme scan specifics
	 */
	public RealMatrix getSignals(double[] params, DW_Scheme rawScheme) {

		RealMatrix signals = new RealMatrix(rawScheme.numMeasurements(), 1);

		StejskalTannerScheme scheme;

		try {
			scheme = (StejskalTannerScheme) rawScheme;
		} catch (ClassCastException cce) {
			throw new LoggedException(
					"scheme object passed to cylinder compartment is not a StejskalTanner sequence");
		}

		for (int i = 0; i < scheme.numMeasurements(); i++) {
			signals.setEntry(i, 0, getSignal(params, scheme, i));
		}
		return signals;
	}

	/**
	 * generates signals from this compartment for each line
	 * in the scheme given.
	 * 
	 * @param scheme scan specifics
	 */
	public double getSignal(double[] params, DW_Scheme rawScheme, int i) {

		StejskalTannerScheme scheme;

		try {
			scheme = (StejskalTannerScheme) rawScheme;
		} catch (ClassCastException cce) {
			throw new LoggedException(
					"scheme object passed to astrostick compartment is not a StejskalTanner sequence");
		}

		return getSignal(params, scheme, i);
	}

	/**
	 * evaluates the model for the i-th line of a given Stejskal-Tanner scheme 
	 * and a set of model parameters 
	 * 
	 * @param params array of model params
	 * @param scheme scan specifics (must be StejskalTanner)
	 * @param i line of scheme to evaluate
	 * 
	 * @return astrostick model at line
	 */
	private final double getSignal(double[] params,
			StejskalTannerScheme scheme, int i) {

		if (scheme.zero(i)) {
			return 1;
		}
		double signal;

		double modG = scheme.getModG(i);
		double b = scheme.getB_Value(i);

		// the perpendicular component
		double lperp = 0;

		//parallel component
		double lpar = -b / (modG*modG) * params[0];

		// the isotropic restricted signal
		try {
			signal = Math.sqrt(Math.PI) * 1
					/ (2 * modG * Math.sqrt(lperp - lpar))
					* Math.exp(modG * modG * lperp)
					* ErrorFunction.erf(modG * Math.sqrt((lperp - lpar)));
		} catch (ErrorFunctionException erfe) {
			throw new LoggedException(erfe);
		}

		return signal;
	}

}
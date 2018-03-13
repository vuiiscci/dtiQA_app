package models.compartments;

import java.util.logging.Logger;

import misc.LoggedException;
import models.ParametricModel;
import imaging.DW_Scheme;
import imaging.StejskalTannerScheme;
import numerics.ErrorFunction;
import numerics.ErrorFunctionException;
import numerics.RealMatrix;

/**
 * Implents the compartment interface for a restricted isotropic compartment,
 * cylinders with distributed orientations.
 * 
 * 
 * @author laura (panagio@cs.ucl.ac.uk)
 * 
 */
public class Astrocylinders extends ParametricModel {

	/** logging object */
	private final Logger logger = Logger.getLogger(this.getClass().getName());

	/** gyromagnetic ratio */
	private final double GAMMA = DW_Scheme.GAMMA;

	/** 60 first roots from the equation j'1(am*x)=0 */
	private final double[] am = { 1.84118307861360, 5.33144196877749,
			8.53631578218074, 11.7060038949077, 14.8635881488839,
			18.0155278304879, 21.1643671187891, 24.3113254834588,
			27.4570501848623, 30.6019229722078, 33.7461812269726,
			36.8899866873805, 40.0334439409610, 43.1766274212415,
			46.3195966792621, 49.4623908440429, 52.6050411092602,
			55.7475709551533, 58.8900018651876, 62.0323477967829,
			65.1746202084584, 68.3168306640438, 71.4589869258787,
			74.6010956133729, 77.7431620631416, 80.8851921057280,
			84.0271895462953, 87.1691575709855, 90.3110993488875,
			93.4530179063458, 96.5949155953313, 99.7367932203820,
			102.878653768715, 106.020498619541, 109.162329055405,
			112.304145672561, 115.445950418834, 118.587744574512,
			121.729527118091, 124.871300497614, 128.013065217171,
			131.154821965250, 134.296570328107, 137.438311926144,
			140.580047659913, 143.721775748727, 146.863498476739,
			150.005215971725, 153.146928691331, 156.288635801966,
			159.430338769213, 162.572038308643, 165.713732347338,
			168.855423073845, 171.997111729391, 175.138794734935,
			178.280475036977, 181.422152668422, 184.563828222242,
			187.705499575101 };

	/**
	 * constructor. needs array of params. in this case, an array of 2
	 * parameters, the diffusivity and the radius R.
	 * 
	 * @param params
	 *            array
	 */
	public Astrocylinders() {
		super(CompartmentType.ASTROCYLINDERS.numParams);
	}

	/**
	 * generates signals from this compartment for each line in the scheme
	 * given.
	 * 
	 * @param scheme
	 *            scan specifics
	 */
	public RealMatrix getSignals(double[] params, DW_Scheme rawScheme) {

		RealMatrix signals = new RealMatrix(rawScheme.numMeasurements(), 1);

		StejskalTannerScheme scheme;

		try {
			scheme = (StejskalTannerScheme) rawScheme;
		} catch (ClassCastException cce) {
			throw new LoggedException("scheme object passed to astrocylinder "
					+ "compartment is not a StejskalTanner sequence");
		}

		for (int i = 0; i < scheme.numMeasurements(); i++) {
			signals.setEntry(i, 0, getSignal(params, scheme, i));
		}
		return signals;
	}

	/**
	 * generates signals from this compartment for the given line in the scheme
	 * given.
	 * 
	 * @param params
	 *            model parameters
	 * @param scheme
	 *            scan specifics
	 * @param i
	 *            line in the scheme to evaluate
	 */
	public double getSignal(double[] params, DW_Scheme rawScheme, int i) {

		StejskalTannerScheme scheme;

		try {
			scheme = (StejskalTannerScheme) rawScheme;
		} catch (ClassCastException cce) {
			throw new LoggedException(
					"scheme object passed to cylinder compartment is not a StejskalTanner sequence");
		}

		return getSignal(params, scheme, i);
	}

	/**
	 * generates signals from this compartment for the given line in the scheme
	 * given. scheme must stejskal-tanner. this is a helper method to reduce
	 * code-duplication.
	 * 
	 * @param params
	 *            model parameters
	 * @param scheme
	 *            scan specifics
	 * @param i
	 *            line in the scheme to evaluate
	 */
	private double getSignal(double[] params, StejskalTannerScheme scheme, int i) {

		if (scheme.zero(i)) {
			return 1;
		}

		double signal;

		double R = params[1];
		double b = scheme.getB_Value(i);
		double DELTA = scheme.getDELTA(i);
		double delta = scheme.getDelta(i);
    	double modG = scheme.getModG(i);
		double[] ghat = scheme.getG_Dir(i);
		double[] G = new double[3];
		for (int j = 0; j < G.length; j++) {
			G[j] = ghat[j] * modG;
		}

		double[] am1 = new double[am.length];

		for (int i1 = 0; i1 < am.length; i1++) {
			am1[i1] = am[i1] / R;

		}

		/** calculating the sum for the perpendicular intra-cellular signal */
		double sum = 0;
		for (int i1 = 0; i1 < am1.length; i1++) {
			// d*am^2
			double dam = params[0] * am1[i1] * am1[i1];
			// -d*am^2*delta
			double e11 = -dam * delta;
			// -d*am^2*DELTA
			double e2 = -dam * DELTA;
			// -d*am^2*(DELTA-delta)
			double dif = DELTA - delta;
			double e3 = -dam * dif;
			// -d*am^2*(DELTA+delta)
			double plus = DELTA + delta;
			double e4 = -dam * plus;
			// numerator of the fraction
			double nom = 2 * dam * delta - 2 + (2 * Math.exp(e11))
					+ (2 * Math.exp(e2)) - Math.exp(e3) - Math.exp(e4);

			// denominator
			double denom = dam * dam * am1[i1] * am1[i1]
					* (R * R * am1[i1] * am1[i1] - 1);

			// the sum of the fraction
			sum += (nom / denom);
		}
		// the perpendicular component
		double lperp = (-2 * GAMMA * GAMMA * sum);

		// parallel component
		double lpar = -b * 1/(modG*modG) * params[0];

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
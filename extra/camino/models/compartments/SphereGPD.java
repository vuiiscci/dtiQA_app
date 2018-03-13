package models.compartments;

import java.util.logging.Logger;

import numerics.ErrorFunction;
import numerics.ErrorFunctionException;
import numerics.RealMatrix;
import misc.LoggedException;
import models.ParametricModel;
import imaging.DW_Scheme;
import imaging.SimulableScheme;
import imaging.StejskalTannerScheme;

/**
 * Implements the compartment interface for an isotropic restricted compartment,
 * Spherical boundary compartment using the GPD approximation.
 * 
 * 
 * @author laura (panagio@cs.ucl.ac.uk)
 *
 */
public class SphereGPD extends ParametricModel {

	/** gyromagnetic ratio */
	private final double GAMMA = DW_Scheme.GAMMA;

	  /**
     * Logging object
     */
    protected static Logger logger = Logger.getLogger("camino.models.compartments.CylinderGPD");



    /** limits required precision in sum. Increase order of magnitude to speed up. */
    private double REQPREC = 1E-7;

    // Speed up for HARDI protocols to avoid repeating calculations
    private double lastDelta = 0;
    private double lastDELTA = 0;
    private double lastDiff = 0;
    private double lastR = 0;
    private double lastSum = 0;

    /** 60 first roots from the equation (am*x)j3/2'(am*x)- 1/2 J3/2(am*x)=0 */
	double[] am = { 2.08157597781810, 5.94036999057271, 9.20584014293667,
			12.4044450219020, 15.5792364103872, 18.7426455847748,
			21.8996964794928, 25.0528252809930, 28.2033610039524,
			31.3520917265645, 34.4995149213670, 37.6459603230864,
			40.7916552312719, 43.9367614714198, 47.0813974121542,
			50.2256516491831, 53.3695918204908, 56.5132704621986,
			59.6567290035279, 62.8000005565198, 65.9431119046553,
			69.0860849466452, 72.2289377620154, 75.3716854092873,
			78.5143405319308, 81.6569138240367, 84.7994143922025,
			87.9418500396598, 91.0842274914688, 94.2265525745684,
			97.3688303629010, 100.511065295271, 103.653261271734,
			106.795421732944, 109.937549725876, 113.079647958579,
			116.221718846033, 116.221718846033, 119.363764548757,
			122.505787005472, 125.647787960854, 128.789768989223,
			131.931731514843, 135.073676829384, 138.215606107009,
			141.357520417437, 144.499420737305, 147.641307960079,
			150.783182904724, 153.925046323312, 157.066898907715,
			166.492397790874, 169.634212946261, 172.776020008465,
			175.917819411203, 179.059611557741, 182.201396823524,
			185.343175558534, 188.484948089409, 191.626714721361 };

    
    /** constructor. needs array of params. 
     * in this case, an array of 2 parameters, the diffusivity and the radius R.
     * 
     * @param params array
     */
    public SphereGPD() {
        super(CompartmentType.SPHEREGPD.numParams);
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
        
        try{
            scheme=(StejskalTannerScheme)rawScheme;
        }
        catch(ClassCastException cce){
            throw new LoggedException("scheme object passed to sphere compartment is not a StejskalTanner sequence");
        }
        
        for (int i = 0; i < scheme.numMeasurements(); i++) {
            signals.setEntry(i, 0, getSignal(params, scheme, i));
        }
        
        return signals;
    }

    
    /**
     * generates signals from this compartment for one particular line
     * in the scheme given.
     * 
     * @param scheme scan specifics
     */
    public double getSignal(double[] params, DW_Scheme rawScheme, int i) {

        double signal;

        StejskalTannerScheme scheme;
        
        try{
            scheme=(StejskalTannerScheme)rawScheme;
        }
        catch(ClassCastException cce){
            throw new LoggedException("scheme object passed to sphere compartment is not a StejskalTanner sequence");
        }
        
        return getSignal(params, scheme, i);
    }

    
    private double getSignal(double[] params, StejskalTannerScheme scheme, int i){

        if (scheme.zero(i))
        {
        	return 1;
        }
    	
        double DELTA = scheme.getDELTA(i);
        double delta = scheme.getDelta(i);
        
        double modG = scheme.getModG(i);
        double[] ghat = scheme.getG_Dir(i);
        double[] G = new double[3];
        for (int j = 0; j < G.length; j++) {
            G[j] = ghat[j] * modG;
        }

        return getSignal(params, delta, DELTA, G);
    }

        
    /**
     * generates signals from this compartment for one particular line
     * in the scheme given.
     * 
     * @param params The model parameters [diff,  R]
     *
     * @param delta Pulse width
     *
     * @param DELTA Pulse separation
     *
     * @param G gradient vector
     */
    public double getSignal(double[] params, double delta, double DELTA, double[] G) {

        double diff = params[0];
        double R = params[1];

        double t= DELTA-delta/3;
        
	double modG = 0;
	for(int i=0; i<G.length; i++)
	    modG += G[i]*G[i];
		modG = Math.sqrt(modG);
	
        double[] am1 = new double[am.length];

        for (int i1 = 0; i1 < am.length; i1++) {
            am1[i1] = am[i1] / R;

        }

        
        double sum = computeGPD_Sum(am1, delta, DELTA, diff, R);

        
        // the restricted signal
        double signal = Math.exp(-2 * GAMMA * GAMMA * modG * modG * sum);
        return signal;
        
    }
    
    /**
     * Function calculating the signal from inside the sphere.
     * @param am1 roots of the bessel function.
     * @param delta gradient duration 
     * @param DELTA diffusion time
     * @param diff parallel diffusivity  in  sphere
     * @param R radius of cylinder
     * @return
     */
    protected double computeGPD_Sum(double[] am1, double delta, double DELTA, double diff, double R) {

	// Check we haven't just computed the same thing.
	if(lastDelta == delta && lastDELTA == DELTA && lastDiff == diff && lastR == R) {
	    return lastSum;
	}

        /** calculating the sum for the perpendicular intra-cellular signal */
        double sum = 0;
	int i1;
        for (i1 = 0; i1 < am1.length; i1++) {
            // d*am^2
            double dam = diff * am1[i1] * am1[i1];
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
                    * (R * R * am1[i1] * am1[i1] - 2);

            // the sum of the fraction
	    double term = nom/denom;
            sum += term;
	    //System.err.println(R+ " " + i1 + " "  + sum + " " + Math.exp(-2 * GAMMA * GAMMA * modG * modG * Math.sin(omega) * Math.sin(omega) * sum));
	    if(term<REQPREC*sum)
		break;
        }
	//System.err.println(i1);
	if(i1==am1.length)
	    logger.warning("Last term in SphereGPD sum greater than " + REQPREC + " times sum.  May lack precision. R=" + R + ". diff=" + diff + ". delta=" + delta + ". DELTA=" + DELTA + ".");

	lastDelta = delta;
	lastDELTA = DELTA;
	lastDiff = diff;
	lastR = R;
	lastSum = sum;

	return sum;
    }

    
    /*    
    public RealMatrix getJacobian(final double[] modParams, final DW_Scheme scheme) {
    	RealMatrix jac = new RealMatrix(scheme.numMeasurements(), modParams.length);
    	
    	//diff
    	RealMatrix diffJ = super.getJacobian(modParams, scheme, 0);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 0, diffJ.entry(r, 0));
        } 
        
        //theta
        double theta = modParams[1];
        
        if (theta < 1e-6)
        {
        	modParams[1] = 1e-6;
        }
        
    	RealMatrix thetaJ = super.getJacobian(modParams, scheme, 1);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 1, thetaJ.entry(r, 1));
        } 
        
        //phi
        
        RealMatrix phiJ = super.getJacobian(modParams, scheme, 2);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 2, theta > 1e-6 ? phiJ.entry(r, 2) : 1);
        } 
        
        //R
        
        RealMatrix rJ = super.getJacobian(modParams, scheme, 3);
        
        for(int r=0; r < scheme.numMeasurements();r++){
        	jac.setEntry(r, 3,  rJ.entry(r, 3));
        } 
    	
    	return jac;
    	
    }
    */
    
}
package models.compartments;

import numerics.RealMatrix;
import misc.LoggedException;
import java.util.logging.*;
import models.ParametricModel;
import imaging.DW_Scheme;
import imaging.StejskalTannerScheme;

public class CylinderGPD extends ParametricModel {

    /**
     * Logging object
     */
    protected static Logger logger = Logger.getLogger("camino.models.compartments.CylinderGPD");

    /** gyromagnetic ratio */
    private final double GAMMA = DW_Scheme.GAMMA;

    /** limits required precision in sum. Increase order of magnitude to speed up. */
    private double REQPREC = 1E-7;

    // Speed up for HARDI protocols to avoid repeating calculations
    private double lastDelta = 0;
    private double lastDELTA = 0;
    private double lastDiff = 0;
    private double lastR = 0;
    private double lastSum = 0;

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

    
    
    /** constructor. needs array of params. 
     * in this case, an array of 4 parameters, the diffusivity, the angles theta and phi
     * that determine the fibre orientation and the radius R.
     * 
     * @param params array
     */
    public CylinderGPD() {
        super(CompartmentType.CYLINDERGPD.numParams);
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
            throw new LoggedException("scheme object passed to cylinder compartment is not a StejskalTanner sequence");
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
            throw new LoggedException("scheme object passed to cylinder compartment is not a StejskalTanner sequence");
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
     * @param params The model parameters [diff, theta, phi, R]
     *
     * @param delta Pulse width
     *
     * @param DELTA Pulse separation
     *
     * @param G gradient vector
     */
    public double getSignal(double[] params, double delta, double DELTA, double[] G) {

        double diff = params[0];
        double theta = params[1];
        double phi = params[2];
        double R = params[3];

        double t= DELTA-delta/3;
        
	double modG = 0;
	for(int i=0; i<G.length; i++)
	    modG += G[i]*G[i];
	modG = Math.sqrt(modG);
	
        double[] am1 = new double[am.length];

        for (int i1 = 0; i1 < am.length; i1++) {
            am1[i1] = am[i1] / R;

        }

        double sinT = Math.sin(theta);
        double cosT = Math.cos(theta);

        double sinP = Math.sin(phi);
        double cosP = Math.cos(phi);
        
      

        // fibre orientation
        double[] n = new double[3];
        n[0] = cosP * sinT;
        n[1] = sinP * sinT;
        n[2] = cosT;

        double modn = 0.0;
        for (int i1 = 0; i1 < n.length; i1++) {
            modn += n[i1] * n[i1];
        }
        modn = Math.sqrt(modn);

        /** omega is the angle between n and G */

        double dotprodGn = 0.0;
        double unitGn = 0.0;
        dotprodGn = (n[0] * G[0]) + (n[1] * G[1]) + (n[2] * G[2]);
        double omega = 0.0;

        unitGn = (modG == 0.0) ? 0.0 : dotprodGn / (modG * modn);

        omega = Math.acos(unitGn);

	double sum = computeGPD_Sum(am1, delta, DELTA, diff, R);

        double sRperp = Math.exp(-2 * GAMMA * GAMMA * modG * modG
                * Math.sin(omega) * Math.sin(omega) * sum);

        // the restricted parallel signal
        double sRpar = Math.exp(-t
                * (GAMMA * delta * modG * Math.cos(omega) * (GAMMA * delta
                        * modG * Math.cos(omega))) * diff);

        // the restricted signal
        double signal = sRperp * sRpar;
        return signal;
        
    }
    
    /**
     * Function calculating the signal from inside the cylinder.
     * @param am1 roots of the bessel function.
     * @param delta gradient duration 
     * @param DELTA diffusion time
     * @param diff parallel diffusivity  in  cylinder
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
                    * (R * R * am1[i1] * am1[i1] - 1);

            // the sum of the fraction
	    double term = nom/denom;
            sum += term;
	    //System.err.println(R+ " " + i1 + " "  + sum + " " + Math.exp(-2 * GAMMA * GAMMA * modG * modG * Math.sin(omega) * Math.sin(omega) * sum));
	    if(term<REQPREC*sum)
		break;
        }
	//System.err.println(i1);
	if(i1==am1.length)
	    logger.warning("Last term in CylinderGPD sum greater than " + REQPREC + " times sum.  May lack precision. R=" + R + ". diff=" + diff + ". delta=" + delta + ". DELTA=" + DELTA + ".");

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
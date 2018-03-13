package imaging;

import misc.*;

import java.util.logging.Logger;


/**
 * 
 * Superclass for Stejskal-Tanner imaging schemes. This superclass class assumes 
 * no knowledge of the shape of the gradient pulses. It defines the wavenumber q based
 * on the average gradient strength |g| over the duration of the pulse delta. The b-value
 * is not computed by this class because different pulse configurations produce different
 * effective diffusion times.
 *
 * @author Philip Cook
 *
 * @version $Id$
 *  
 */
public abstract class StejskalTannerScheme extends DW_Scheme implements SimulableScheme {




    /**
     * logging object
     */
    private static final Logger logger = Logger.getLogger("camino.imaging.StejskalTannerScheme");
    
   
    /**
     * The average gradient strength of the diffusion-weighting pulses, such that q = gamma * delta * modG.
     */
    protected double[] modG;

    /**
     * The pulse separations.
     */
    private final double[] bigDel;

    /**
     * The pulse widths.
     */
    private final double[] smallDel;

    /**
     * The echo times.
     */
    private final double[] TE;

  
    /**
     * Construct a scheme from the imaging parameters.
     *
     * @param initG_Dir gradient directions.
     * @param initModG the average gradient strengths.
     * @param initBigDel gradient separations.
     * @param initSmallDel gradient durations.
     * @param initTE echo times.
     *
     */
    protected StejskalTannerScheme(double[][] initG_Dir, double[] initModG, double[] initBigDel, 
				   double[] initSmallDel, double[] initTE) {
        
	super(initG_Dir);

	// Create the pulse-sequence parameter arrays.
	modG = new double[numMeas];
        bigDel = new double[numMeas];
        smallDel = new double[numMeas];
        TE = new double[numMeas];
	
	System.arraycopy(initModG, 0, modG, 0, numMeas);
	System.arraycopy(initBigDel, 0, bigDel, 0, numMeas);
	System.arraycopy(initSmallDel, 0, smallDel, 0, numMeas);
	System.arraycopy(initTE, 0, TE, 0, numMeas);
	
   
    }

    
    /**
     * Copy constructor used by flipping / gradient direction swapping methods. 
     *
     * @param initG_Dir gradient directions for the new scheme. The number of directions
     * should equal the number of directions in the source scheme.
     *
     * @param source all imaging parameters <b>except</b> the gradient directions are taken from 
     * this object.
     *
     */
    protected StejskalTannerScheme(double[][] initG_Dir, StejskalTannerScheme source) {
	this(initG_Dir, source.modG, source.bigDel, source.smallDel, source.TE);

	if (numMeas != source.numMeas) {
	    throw new LoggedException("Error in scheme copy construction, number of measurements was " + 
				      numMeas + ", expected " + source.numMeas); 
	}
    }
    
    /**
     * Copy constructor used by copy scheme methods.
     *
     *
     * @param source all imaging parameters are taken from
     * this object.
     *
     */
    protected StejskalTannerScheme(StejskalTannerScheme source) {
        this(source.gDir, source.modG, source.bigDel, source.smallDel, source.TE);
        
        if (numMeas != source.numMeas) {
            throw new LoggedException("Error in scheme copy construction, number of measurements was " +
                                      numMeas + ", expected " + source.numMeas);
        }
    }
      
    /**
     * Calculates the mean |q| in the set containing all measurements with |q| > 0.
     *
     * @return the mean |q| of the non-zero measurements.
     *
     */
    public final double getMeanNonZeroModQ() {
        double sumNZMQ = 0.0;
        for(int i=0; i<numMeas; i++) {
            if(!zero(i)) sumNZMQ += getModQ(i);
        }

        return sumNZMQ / (numMeas - numZeroMeas);
    }


    /**
     *
     */
    public double[] getModG() {
        return modG;
    }


    /**
     * Gets all non-zero q normalized by <code>getMeanNonZeroModQ()</code>.
     *
     * @return a list of normalized non-zero q vectors.
     */
    public final double[][] getNormNonZeroQs() {

        double[][] nzq = new double[numMeas - numZeroMeas][3];
        int nextInd = 0;

	double normMod = getMeanNonZeroModQ();

	for(int i=0; i<numMeas; i++) {
            if(!zero(i)) {

		double modQ = getModQ(i);

		double[] gDir = getG_Dir(i);

                for(int j=0; j<3; j++) {
                    nzq[nextInd][j] = gDir[j] * modQ / normMod;
                }
                nextInd += 1;
            }
        }

        return nzq;
    }
 

    /**
     * Gets all non-zero q.
     *
     * @return all q vectors that are not zero.
     */
    public final double[][] getNonZeroQs() {
        double[][] nzq = new double[numMeas - numZeroMeas][3];
        int nextInd = 0;
        for(int i=0; i<numMeas; i++) {
            if(!zero(i)) {
		
		double modQ = getModQ(i);
		
		double[] gDir = getG_Dir(i);
                
		for(int j=0; j<3; j++) {
                    nzq[nextInd][j] = modQ*gDir[j];
                }
                nextInd += 1;
            }
        }

        return nzq;
    }


    /**
     * Gets all non-zero gradient vectors. 
     *
     * @return the un-normalized, non-zero diffusion gradient vectors.
     */
    public final double[][] getNonZeroGs() {
        double[][] nzg = new double[numMeas - numZeroMeas][3];
        int nextInd = 0;
        for(int i=0; i<numMeas; i++) {
            if(!zero(i)) {
		
		double[] gDir = getG_Dir(i);
                
		for(int j=0; j<3; j++) {
                    nzg[nextInd][j] = modG[i]*gDir[j];
                }
                nextInd += 1;
            }
        }

        return nzg;
    }


    /**
     * Gets DELTA (gradient pulse separation) for all non-zero measurements. 
     *
     * @return an array of DELTA for measurements with q > 0.
     */
    public final double[] getNonZeroBigDeltas() {
	return getNonZeroParam(bigDel);
    }
 
   
    /**
     * Gets delta (gradient pulse duration) for all non-zero measurements. 
     *
     * @return an array of delta for measurements with q > 0.
     */
    public final double[] getNonZeroSmallDeltas() {
    	return getNonZeroParam(smallDel);
    }
    

    /**
     * Gets gradient strengths for all non-zero measurements.
     * 
     * @return an array of modG for measurements with q > 0.
     *
     * @see #getModG(int)
     */
    public final double[] getNonZeroModGs() {
	return getNonZeroParam(modG);
    }


    /**
     * Gets TE for all non-zero measurements.
     * 
     * @return an array of TE for measurements with q > 0.
     *
     */
    public final double[] getNonZeroTEs() {
	return getNonZeroParam(TE);
    }

    /**
     * Gets the gradient pulse duration for measurement i.
     *
     * @return the gradient pulse duration for measurement i.
     *
     */
    public final double getDelta(int i) {
        return smallDel[i];
    }


    /**
     * Gets the gradient pulse separation for measurement i.
     *
     * @return the gradient pulse separation for measurement i.
     *
     */
    public final double getDELTA(int i) {
        return bigDel[i];
    }


    /**
     * Gets the average gradient strength for measurement i. 
     * The average is defined over the time delta.
     *
     * @return the average gradient strength for measurement i.
     *
     */
    public final double getModG(int i) {
        return modG[i];
    }


    /**
     * Gets the gradient vector g for measurement i, which is the gradient direction scaled by 
     * the average gradient magnitude.
     *
     *
     * @return the diffusion-weighting gradient vector for measurement i.
     */
    public final double[] getG(int i) {
	double[] g = getG_Dir(i);

	for (int j = 0; j < 3; j++) {
	    g[j] *= modG[i];
	}

	return g;
    }



    /**
     * Gets the radial wavenumber q for measurement i.
     *
     * @return the radial wavenumber q.
     */
    public final double[] getQ(int i) {
	double[] q = getG_Dir(i);

	double modQ = getModQ(i);
	
	for (int j = 0; j < 3; j++) {
	    q[j] *= modQ;
	}

	return q;
    }


    /**
     * Gets the echo time for measurement i. 
     *
     * @return the echo time for measurement i.
     *
     */
    public final double getTE(int i) {
        return TE[i];
    }


    /**
     * Gets the wavenumber of measurement i.
     *
     * @return |q| for measurement i.
     */
    public final double getModQ(int i) {
	return GAMMA * smallDel[i] * modG[i];
    }

 
    /**
     * returns the gradient impulse between this and most recent call
     * 
     * @param i line in scheme file
     * @param t current time
     * @param tLast time of most recent previous call
     */
    public abstract double[] getGradImpulse(int i, double t, double tLast);

    
    /**
     * returns the duration necessary for simulation dynamics to cover the
     * whole specified sequence
     * 
     * @return max TE in the scheme in seconds
     */
    public abstract double getDuration();



}

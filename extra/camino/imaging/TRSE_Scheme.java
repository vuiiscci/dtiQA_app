package imaging;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.logging.Logger;

import misc.LoggedException;


/**
 * Scheme file for twice-refocused pulse sequences, which use four diffusion-weighting
 * pulses, which begin at time t_delN and last for delN seconds. 
 * <p>
 * Note that the time of the start of the third pulse is NOT specified in the scheme file, 
 * as it is assumed to start at the end of the second pulse. Also, the duration of the 
 * fourth pulse is not specified either but calculated so as to enforce 
 * del1 + del2 = del3 + del4. If del3 is greater than del1 + del2 an exception will be 
 * thrown.
 * <p> 
 * This implementation also assumes that the gradient strength in each
 * pulse is the same, and in the same direction.
 * <p>
 * For more information see Reese et al, 
 * http://www.ncbi.nlm.nih.gov/pubmed/12509835
 *
 * @author Matt Hall, Philip Cook
 * @version $Id$
 *
 *
 */
public abstract class TRSE_Scheme extends DW_Scheme implements SimulableScheme {

    /** logging object */
    private static final Logger logger= Logger.getLogger("camino.imaging.TRSE_Scheme");
    
    /** mod G in each dir */
    protected double[] modG;

    /** big delta in each dir */
    private final double[] bigDel;
    
    /** duration of first grad block in each dir */
    private final double[] del1;
    
    /** time of first grad block in each dir */
    private final double[] t_del1;
    
    /** duration of second grad block in each dir */
    private final double[] del2;

    /** time of second grad block in each dir */
    private final double[] t_del2;
    
    /** duration of third grad block in each dir */
    private final double[] del3;

    /** time of third grad block in each dir */
    private final double[] t_del3;
        
    /** duration of fourth grad block in each dir */
    private final double[] del4;

    /** time of fourth grad block in each dir */
    private final double[] t_del4;
    
    /** echo time in each dir */
    private final double[] TE;

    
      
		
    /**
     * @param gDir gradient directions.
     *
     *
     */    
    protected TRSE_Scheme(double[][] gDir, double[] modG, double[] del1, double[] t_del1, 
		       double[] del2, double[] t_del2, double[] del3, double[] t_del4, double[] te) {
	    
	super(gDir);
	   
	// should make a def copy here
 
	this.modG = new double[numMeas];
        this.bigDel = new double[numMeas];
	this.del1 = new double[numMeas];
	this.t_del1 = new double[numMeas];
	this.del2 = new double[numMeas];
	this.t_del2 = new double[numMeas];
	this.del3 = new double[numMeas];
	this.t_del3 = new double[numMeas];
	this.del4 = new double[numMeas];
	this.t_del4 = new double[numMeas];
	this.TE = new double[numMeas];

	System.arraycopy(modG, 0, this.modG, 0, numMeas);
	System.arraycopy(del1, 0, this.del1, 0, numMeas);
	System.arraycopy(t_del1, 0, this.t_del1, 0, numMeas);
	System.arraycopy(del2, 0, this.del2, 0, numMeas);
	System.arraycopy(t_del2, 0, this.t_del2, 0, numMeas);
	System.arraycopy(del3, 0, this.del3, 0, numMeas);
	System.arraycopy(t_del4, 0, this.t_del4, 0, numMeas);

	System.arraycopy(te, 0, this.TE, 0, numMeas);
	    
	for (int i = 0; i < numMeas; i++) {
	    t_del3[i] = t_del2[i] + del2[i];

	    bigDel[i] = t_del3[i] - t_del1[i];

            if(del3[i]>(del1[i]+del2[i])){
            	throw new LoggedException(" third gradient pulse in direction "
					  +i+" of schemefile is longer than sum of first two.");
            }

	    del4[i] = (del1[i]+del2[i])-del3[i];
	}
	    
    }
	


    /**
     * Copy constructor needed by flip / gradient swap methods.
     */
    protected TRSE_Scheme(double[][] gDir, TRSE_Scheme source) {
	    
	this(gDir, source.modG, source.del1, source.t_del1, source.del2, 
	     source.t_del2, source.del3, source.t_del4, source.TE);
  
    }
	
        /**
     * Copy constructor used by copy scheme methods.
     *
     *
     * @param source all imaging parameters are taken from
     * this object.
     *
     */
    protected TRSE_Scheme(TRSE_Scheme source) {
        this(source.gDir, source.modG, source.del1, source.t_del1, source.del2, 
         source.t_del2, source.del3, source.t_del4, source.TE);
    }

    /**
     *
     * @return all non-zero gradient vectors. 
     *
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
     * Gets DELTA (separation of the onset of the first and third pulse) for all non-zero measurements. 
     *
     * @return an array of DELTA for measurements with b > 0.
     */
    public final double[] getNonZeroBigDeltas() {
	return getNonZeroParam(bigDel);
    }
 
   
    /**
     * Gets delta1 for all non-zero measurements. 
     *
     * @return an array of t_delta1 for measurements with b > 0.
     */
    public final double[] getNonZeroDelta1s() {
    	return getNonZeroParam(del1);
    }

    /**
     * Gets t_delta1 for all non-zero measurements. 
     *
     * @return an array of t_delta1 for measurements with b > 0.
     */
    public final double[] getNonZeroT_Delta1s() {
    	return getNonZeroParam(t_del1);
    }


    /**
     * Gets delta2 for all non-zero measurements. 
     *
     * @return an array of t_delta2 for measurements with b > 0.
     */
    public final double[] getNonZeroDelta2s() {
    	return getNonZeroParam(del2);
    }

    /**
     * Gets t_delta2 for all non-zero measurements. 
     *
     * @return an array of t_delta2 for measurements with b > 0.
     */
    public final double[] getNonZeroT_Delta2s() {
    	return getNonZeroParam(t_del2);
    }

    /**
     * Gets delta3 for all non-zero measurements. 
     *
     * @return an array of t_delta3 for measurements with b > 0.
     */
    public final double[] getNonZeroDelta3s() {
    	return getNonZeroParam(del3);
    }

    /**
     * Gets t_delta3 for all non-zero measurements. 
     *
     * @return an array of t_delta3 for measurements with b > 0.
     */
    public final double[] getNonZeroT_Delta3s() {
    	return getNonZeroParam(t_del1);
    }

    /**
     * Gets delta4 for all non-zero measurements. 
     *
     * @return an array of t_delta4 for measurements with b > 0.
     */
    public final double[] getNonZeroDelta4s() {
    	return getNonZeroParam(del4);
    }

    /**
     * Gets t_delta4 for all non-zero measurements. 
     *
     * @return an array of t_delta4 for measurements with b > 0.
     */
    public final double[] getNonZeroT_Delta4s() {
    	return getNonZeroParam(t_del4);
    }
    
    /**
     * Gets TE for all non-zero measurements. 
     *
     * @return an array of TE for measurements with b > 0.
     */
    public final double[] getNonZeroTEs() {
    	return getNonZeroParam(TE);
    }
     
   

    /**
     * Gets all gradient strengths for non-zero measurements.
     * 
     * @return an array of modG for measurements with q > 0.
     *
     * @see #getModG(int)
     */
    public final double[] getNonZeroModGs() {
	return getNonZeroParam(modG);
    }



    /**
     * @return the gradient vector g for measurement i, which is the gradient direction scaled by 
     * the gradient magnitude.
     *
     */
    public final double[] getG(int i) {
	double[] g = getG_Dir(i);

	for (int j = 0; j < 3; j++) {
	    g[j] *= modG[i];
	}

	return g;
    }



    /** 
     *
     * @return big delta for measurement i.
     *
     */
    public final double getDELTA(int i) {
	return bigDel[i];
    }

    	
    /**
     *
     * @return gradient strength for measurement i.
     *
     */
    public final double getModG(int i) {
	return modG[i];
    }
   
    /**
     * @return time of start of first grad block for measurement i.
     */
    public final double getT_Del1(int i){
	return t_del1[i];
    }

    /**
     * @return time of start of second grad block for measurement i.
     */
    public final double getT_Del2(int i){
	return t_del2[i];
    }

    /**
     * @return time of start of third grad block for measurement i.
     */
    public final double getT_Del3(int i){
	return t_del3[i];
    }

    /**
     * @return time of start of fourth grad block for measurement i.
     */
    public final double getT_Del4(int i){
	return t_del4[i];
    }

    /**
     * @return the duration of the first grad block for measurement i.
     */
    public final double getDel1(int i){
	return del1[i];
    }

    /**
     * @return the duration of the second grad block for measurement i.
     */
    public final double getDel2(int i){
	return del2[i];
    }

    /**
     * @return the duration of the third grad block for measurement i.
     */
    public final double getDel3(int i){
	return del3[i];
    }

    /**
     * @return the duration of the fourth grad block for measurement i.
     */
    public final double getDel4(int i){
	return del4[i];
    }


    /**
     * @return the echo time for measurement i.
     */
    public final double getTE(int i) {
        return TE[i];
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

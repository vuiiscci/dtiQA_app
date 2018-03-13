package imaging;

import misc.*;
import java.util.logging.Logger;

/**
 * Super class for double Pulsed Field Gradients scheme files.
 * @author gmorgan
 * @version $Id$
 *
 */

public abstract class DPFG_Scheme extends DW_Scheme implements SimulableScheme {
	
	/**
	 * logging object
	 */
    private static final Logger logger = Logger.getLogger("camino.imaging.DPFG_Scheme");
    
    /**
     * The average gradient strength of the first set of diffusion-weighting pulses, such that q = gamma * delta * modG.
     */
    private final double[] modG1;

    /**
     * The pulse separations of the first PGSE block.
     */
    private final double[] bigDel1;

    /**
     * The pulse widths of the first PGSE block.
     */
    private final double[] smallDel1;

    /**
     * The mixing time between the two PGSE blocks
     */
    private final double[] TM;
    
    /**
     * The gradient direction of the second PGSE block
     */
    private final double[][] gDir2;
    
    /**
     * The average gradient strength of the second set of diffusion-weighting pulses, such that q = gamma * delta * modG.
     */
    private final double[] modG2;

    /**
     * The pulse separations of the second PGSE block.
     */
    private final double[] bigDel2;

    /**
     * The pulse widths of the second PGSE block.
     */
    private final double[] smallDel2;

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
    protected DPFG_Scheme(double[][] initG_Dir1, double[] initModG1, double[] initBigDel1, 
    		double[] initSmallDel1, double[] initTM, double[][] initG_Dir2, double[] initModG2, double[] initBigDel2,
    		double[] initSmallDel2, double[] initTE) {

    	super(initG_Dir1);

    	// Create the pulse-sequence parameter arrays.
    	this.modG1 = new double[numMeas];
    	this.bigDel1 = new double[numMeas];
    	this.smallDel1 = new double[numMeas];
    	this.TM = new double[numMeas];
    	this.gDir2 = new double[numMeas][3];
    	this.modG2 = new double[numMeas];
    	this.bigDel2 = new double[numMeas];
    	this.smallDel2 = new double[numMeas];
    	this.TE = new double[numMeas];

    	System.arraycopy(initModG1, 0, modG1, 0, numMeas);
    	System.arraycopy(initBigDel1, 0, bigDel1, 0, numMeas);
    	System.arraycopy(initSmallDel1, 0, smallDel1, 0, numMeas);
    	System.arraycopy(initTM, 0, TM, 0, numMeas);
    	System.arraycopy(initG_Dir2, 0, gDir2, 0, numMeas);
    	System.arraycopy(initModG2, 0, modG2, 0, numMeas);
    	System.arraycopy(initBigDel2, 0, bigDel2, 0, numMeas);
    	System.arraycopy(initSmallDel2, 0, smallDel2, 0, numMeas);
    	System.arraycopy(initTE, 0, TE, 0, numMeas);

    }
    
    /**
     * Copy constructor needed by flip / gradient swap methods.
     */
    protected DPFG_Scheme(double[][] gDir, double[][] gDir2, DPFG_Scheme source) {
	    
	     this(gDir, source.modG1, source.bigDel1, source.smallDel1, source.TM, source.gDir2, 
	     source.modG2, source.bigDel2, source.smallDel2, source.TE);
  
    }

    /**
     * Copy constructor needed to copy scheme
     */
    protected DPFG_Scheme(DPFG_Scheme source) {
      
       this(source.gDir, source.modG1, source.bigDel1, source.smallDel1, source.TM, source.gDir2, 
       source.modG2, source.bigDel2, source.smallDel2, source.TE);
  
    }

    /**
    *
    * @return all non-zero gradient vectors. 
    *
    */
   public final double[][] getNonZeroG1s() {
       double[][] nzg = new double[numMeas - numZeroMeas][3];
       int nextInd = 0;
       for(int i=0; i<numMeas; i++) {
           if(!zero(i)) {

		double[] gDir1 = getG_Dir(i);
               
		for(int j=0; j<3; j++) {
                   nzg[nextInd][j] = modG1[i]*gDir1[j];
               }
               nextInd += 1;
           }
       }

       return nzg;
   }


   /**
    * Gets bigDel1 (separation of the onset of the first and second pulse) for all non-zero measurements. 
    *
    * @return an array of DELTA for measurements with b > 0.
    */
   public final double[] getNonZeroBigDelta1s() {
	return getNonZeroParam(bigDel1);
   }

  
   /**
    * Gets smallDel1 for all non-zero measurements. 
    *
    * @return an array of smallDel1 for measurements with b > 0.
    */
   public final double[] getNonZeroDelta1s() {
   	return getNonZeroParam(smallDel1);
   }

   /**
    * Gets TM for all non-zero measurements. 
    *
    * @return an array of TM for measurements with b > 0.
    */
   public final double[] getNonZeroTMs() {
   	return getNonZeroParam(TM);
   }

   /**
   *
   * @return all non-zero gradient vectors for the second PGSE block. 
   *
   */
  public final double[][] getNonZeroG2s() {
      double[][] nzg = new double[numMeas - numZeroMeas][3];
      int nextInd = 0;
      for(int i=0; i<numMeas; i++) {
          if(!zero(i)) {

		double[] gDir2 = getG2(i);
              
		for(int j=0; j<3; j++) {
                  nzg[nextInd][j] = modG2[i]*gDir2[j];
              }
              nextInd += 1;
          }
      }

      return nzg;
  }

   /**
    * Gets bigDelta2 for all non-zero measurements. 
    *
    * @return an array of bigDelta2 for measurements with b > 0.
    */
  public final double[] getNonZeroBigDelta2s() {
	  return getNonZeroParam(bigDel2);
  }

   /**
    * Gets smallDelta2 for all non-zero measurements. 
    *
    * @return an array of smallDelta2 for measurements with b > 0.
    */
   public final double[] getNonZeroSmallDelta2s() {
   	return getNonZeroParam(smallDel2);
   }

   /**
    * Gets TE for all non-zero measurements. 
    *
    * @return an array of TE for measurements with b > 0.
    */
   public final double[] getNonZeroTE() {
   	return getNonZeroParam(TE);
   }

   /**
    * Gets all gradient strengths for non-zero measurements.
    * 
    * @return an array of modG1 for measurements with q > 0.
    *
    * @see #getModG(int)
    */
   public final double[] getNonZeroModG1s() {
	return getNonZeroParam(modG1);
   }
   
   /**
    * Gets all gradient strengths for non-zero measurements.
    * 
    * @return an array of modG2 for measurements with q > 0.
    *
    * @see #getModG(int)
    */
   public final double[] getNonZeroModG2s() {
	return getNonZeroParam(modG2);
   }

   /**
    * @return the gradient vector g for the first PGSE block for measurement i, which is the gradient direction scaled by 
    * the gradient magnitude.
    *
    */
   public final double[] getG1(int i) {
	double[] g = getG_Dir(i);

	for (int j = 0; j < 3; j++) {
	    g[j] *= modG1[i];
	}

	return g;
   }
   	
   /**
    *
    * @return gradient strength for the first PGSE block for measurement i.
    *
    */
   public final double getModG1(int i) {
	return modG1[i];
   }
  
   /**
    * @return bigDelta for first PGSE block for measurement i.
    */
   public final double getBigDel1(int i){
	return bigDel1[i];
   }

   /**
    * @return smallDelta for first PGSE block for measurement i.
    */
   public final double getSmallDel1(int i){
	return smallDel1[i];
   }
   
   /** 
   *
   * @return mixing time between the PGSE blocks for measurement i.
   *
   */
  public final double getTM(int i) {
	return TM[i];
  }
   
   /**
    * @return the gradient vector g for the second PGSE block for measurement i, which is the gradient direction scaled by 
    * the gradient magnitude.
    *
    */
   public final double[] getG2(int i) {
	   double[] g2=new double[3];

       g2[0]=gDir2[i][0];
       g2[1]=gDir2[i][1];
       g2[2]=gDir2[i][2];

       return g2;
   }
   	
   /**
    *
    * @return gradient strength for the second PGSE block for measurement i.
    *
    */
   public final double getModG2(int i) {
	return modG2[i];
   }
  
   /**
    * @return bigDelta for second PGSE block for measurement i.
    */
   public final double getBigDel2(int i){
	return bigDel2[i];
   }

   /**
    * @return smallDelta for first PGSE block for measurement i.
    */
   public final double getSmallDel2(int i){
	return smallDel2[i];
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

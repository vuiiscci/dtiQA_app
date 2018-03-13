package imaging;

import java.util.Scanner;
import java.util.Vector;
import java.util.logging.Logger;

import simulation.DiffusionSimulation;

import misc.LoggedException;

/**
 *  scheme object for general gradient waveforms. these are defined
 *  step-wise and can potentially take on arbitrary pulse shapes, 
 *  including pulses that vary in direction as well as intensity.
 *  
 *  scheme file format is:
 *  VERSION: GENERAL_WAVEFORM
 *  K dt G1x G1y G1z G2x G2y G2z ... GKx GKy GKz
 *  
 *  where K is the number of gradient-slices specified on the line and
 *  dt is the time increment associated with each gradient slice. The rest 
 *  of the line consistes of 3K entries describing the gradient strength
 *  and orientation in each slice.
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class GradientWaveform_Scheme extends DW_Scheme implements SimulableScheme {

    /** logging object */
    private final Logger logger= Logger.getLogger(this.getClass().getName()); 
    
    /** version string for scheme */
    protected static final String VERSION = "GRADIENT_WAVEFORM";
    
    /** waveforms arrays */
    private final double[][][] gradWave;

    /** time increment array */
    private final double[] dt;
    
    /** array of zero flags (overriding superclass) */
    protected final boolean[] zero;
            
    /** number of zero measurements */
    private final int numZeros;
    
    /** 
     * constructor. takes waveforms and time increments.
     * 
     * 
     * @param gradwave array of complete (i.e. both halves) gradient waveforms
     * @param dt array of time increments. one for each line in scheme
     */
    private GradientWaveform_Scheme(double[][][] gradWave, double[] dt){
        
        super(new double[0][0]);
        
        this.dt=dt;
        this.gradWave=gradWave;
        this.zero= new boolean[gradWave.length];
        
        int numZeros=0;
        
        for(int i=0; i<gradWave.length; i++){
            this.zero[i]=true;
            
            for(int j=0; j<gradWave[i].length; j++){
                for(int k=0; k<3; k++){
                    if(gradWave[i][j][k]!=0.0){
                        this.zero[i]=false;
                        break;
                    }
                }
                if(!this.zero[i]){
                    break;
                }
            }
            if(zero[i]){
                numZeros++;
            }
        }
        
        this.numZeros=numZeros;
    }

    /**
     * Copy constructor
     *
     * @param GradientWaveform_Scheme that we want to copy.
     *
     */
    public GradientWaveform_Scheme(GradientWaveform_Scheme source) {
        this(source.gradWave, source.dt);
    }

    // /**
    // * Resests the values in scheme to be the same as in source.
    // */
    // public void resetScheme(HCPScheme source){
    //     //Reset scheme 
    // }

    // /**
    // *
    // */
    // public void modifyScheme(double[] gradAdj){
    //     // modifies the scheme with gradient adjustment.
    // }
    
    // /**
    //  * @return the b-value array.
    //  */
    // public double[] getB_Values() {
    //    return null;
    // }

    // /**
    //  * @return ModG
    //  */
    // public double[] getModG() {
    //    return null;
    // }

    /**
     * overrides the number of measurements method
     * 
     * @return number of measurements in the gradient wave file
     */
    public int numMeasurements(){
        
        return gradWave.length;
    }
    
    /**
     * overrides the number of zeros measurements method to 
     * use local count.
     */
    public int numZeroMeasurements(){
        
        return this.numZeros;
    }
    
    /**
     * flips the current entry in the 
     * @param i
     */
    private GradientWaveform_Scheme flippedScheme(int flipInd){
        
        double[][][] flippedGradWave= new double[gradWave.length][][];
        
        for(int i=0; i<flippedGradWave.length; i++){
            flippedGradWave[i]= new double[gradWave[i].length][3];
            
            for(int j=0; j<flippedGradWave[i].length; j++){
                for(int k=0; k<3; k++){
                    if(k==flipInd){
                        flippedGradWave[i][j][k]=-gradWave[i][j][k];
                    }
                    else{
                        flippedGradWave[i][j][k]= gradWave[i][j][k];
                    }
                }
            }
        }
        
        return new GradientWaveform_Scheme(flippedGradWave, dt);
    }
    
    
    /**
     * returns a new copy of the current scheme
     * with gradients flipped in the x-direction
     */    
    public DW_Scheme flipX() {
        return flippedScheme(0);
    }

    /**
     * returns a new copy of the current scheme
     * with gradients flipped in the y-direction
     */    
    public DW_Scheme flipY() {
        return flippedScheme(1);
    }

    /**
     * returns a new copy of the current scheme
     * with gradients flipped in the z-direction
     */    
    public DW_Scheme flipZ() {
        return flippedScheme(2);
    }

    // /**
    //  * Copies the scheme.
    //  */
    // public HCPScheme copyScheme(){
    //     return new GradientWaveform_Scheme(this);
    // }

    /** 
     * calculates the b-factor of a general waveform using
     * Price's framework (Price, 1997). This defines b-value 
     * as 
     *      b= \gamma \int_0^t F^2(t')dt'
     * where
     *      F(t')= \int_0^t' \mathbf{G}(t'')dt''
     * 
     * i.e. the integral of the square of the integral over
     * the waveform in 3D. Here we use the same time increment
     * for both integrals.
     * 
     * @param the index of the direction we want the b-value in
     * 
     * @return the b-value, as the result of both integrals 
     *         in SI units
     * 
     * @see imaging.DW_Scheme#getB_Value(int)
     */
    public double getB_Value(int i) {
        
        // accumulate b value along wave
        double b= 0.0;
        
        for(int t=0; t<gradWave[i].length; t++){
            
            // get gradient integral
            double[] F = getGradientIntegral(i, t);
            
            // square the F vector
            double Fsquared=0.0;
            for(int k=0; k<3; k++){
                Fsquared+=F[k]*F[k];
            }
            
            // accumulate F^2
            b+= Fsquared;
            
        }
        
        // multiply by dt
        b*=dt[i];
        
        // multiply by gyromag ratio
        b*= DW_Scheme.GAMMA*DW_Scheme.GAMMA;
        
        // ... and we're done.
        return b;
    }

    
    /**
     * calculates the integral over a specified gradient
     * waveform between zero and a specified index (inclusive)
     * 
     * @param i index of waveform
     * @param t index to stop at (includes this indexes contribution)
     * 
     * @return vector of integrated gradients
     */
    private double[] getGradientIntegral(int i, int t){
        
        final double[] F= new double[3];
        
        
        for(int k=0; k<3; k++){
         // reset summation value
            F[k]=0.0;
        
            // sum over histogram values
            for(int j=0; j<=t; j++){
                F[k]+=gradWave[i][j][k];
            }

            // multiply by dt
            F[k]*=dt[i];
        }
        
        return F;
        
    }
    
    /**
     * generates a new scheme from a specified subset of the acquisitions
     * in the current scheme.
     * 
     * @param indeices the subset desired
     * 
     * @return new GradientWaveform_Scheme from subset
     */
    public DW_Scheme getSubsetScheme(int[] indices) {
        
        // space for subset of waveforms
        double[][][] subsetGradWave= new double[indices.length][][];
        
        // space for subset of dts
        double[] subsetDt= new double[indices.length];
        
        // counter for subsets
        int subsetIndex=0;
        
        for(int i=0; i<indices.length; i++){
        
            // are we in the subset?
            if(i==indices[i]){
                
                // allocate space for waveform
                subsetGradWave[subsetIndex]= new double[gradWave[i].length][3];
                
                // copy waveform
                for(int j=0; j<subsetGradWave[subsetIndex].length; j++){
                    for(int k=0; k<3; k++){
                        subsetGradWave[subsetIndex][j][k]=gradWave[i][j][k];
                    }
                }
                
                // copy dt
                subsetDt[subsetIndex]= dt[i];
                
                // increment subset index
                subsetIndex++;
            }
        
        }
        
        // construct and return subset scheme 
        return new GradientWaveform_Scheme(subsetGradWave, subsetDt);
    }

    /**
     * constructs a new set of waveforms in an order specified by
     * the given array of indices. This permuted scheme uses fresh
     * instances of all the arrays, not the originals, so changes
     * made will not propagate back to the original scheme.
     * 
     * note that only lines specified in the order array will be 
     * copied into the permuted scheme, so if you leave any out 
     * they will not be present in the output scheme.
     * 
     * @param order set of indices defining the new order
     * 
     * @return new scheme with permuted acquisition order
     * 
     */
    public DW_Scheme gradOrder(int[] order) {
        
        double[][][] reorderedGradWave= new double[order.length][][];
        double[] reorderedDt= new double[order.length];
        
        
        for(int i=0; i<order.length; i++){
            int index= order[i];
            
            reorderedGradWave[i]= new double[gradWave[index].length][3];
            for(int j=0; j<reorderedGradWave[i].length; j++){
                for(int k=0; k<3; k++){
                    reorderedGradWave[i][j][k]= gradWave[index][j][k];
                }
                
                reorderedDt[i]= dt[index];
            }
        }
         
        
        return new GradientWaveform_Scheme(reorderedGradWave, reorderedDt);
    }


    /**
     * returns the duration necessary from the dynamics.
     * this is the length of the longest acquisition
     * waveform in the scheme, so we calculate the duration
     * of each wave and return it.
     * 
     * @return required duration of dynamics in seconds
     */
    public double getDuration() {
        
        double duration= 0.0;
        
        for(int i=0; i<gradWave.length; i++){
            
            double duration_i= gradWave[i].length*dt[i];
            
            if(duration< duration_i){
                duration= duration_i;
            }
            
        }
        
        return duration;
    }

    /**
     * calculates the net area under the gradient pulse. In this case
     * this is effectively an intergal over the duration of the simulation
     * timestep. Since in general this will not match the increment 
     * associated with steps in the waveform, care must be taken to
     * include partial steps at the beginning and end and to catch the case 
     * where simulation timesteps are completely contained within a single
     * waveform update.
     * 
     * @param i acquisition index
     * @param t current time (end of simulation timestep)
     * @param last previous time called (beginning of simulation timestep)
     *  
     *  
     * @see imaging.SimulableScheme#getGradStrength(int, double, double)
     */
    public double[] getGradImpulse(int i, double t, double last) {
        // return array
        final double[] Gdt=new double[DiffusionSimulation.D];
       
        // index of first waveform block intersected by timestep
        int iStart= (int)Math.floor(last/dt[i]);
        
        // index of last waveform block intersected by timestep
        int iFinish= (int)Math.floor(t/dt[i]);
        
        if(iStart>=gradWave[i].length){
        	for(int j=0; j<Gdt.length; j++){
        		Gdt[j]=0.0;
        	}
        	return Gdt;
        }
        
        // is the timestep completely contained within a single block?
        if(iStart==iFinish){
        
            double len= t-last;
            for(int j=0; j< Gdt.length; j++){
                
                Gdt[j]=gradWave[i][iStart][j]*len;
                
            }
            
            return Gdt;
        }
        
        // reset return array values
        for(int j=0; j<Gdt.length; j++){
            Gdt[j]=0.0;
        }
        
        for(int j=0; j<Gdt.length; j++){
            // add first partial block            
            Gdt[j]= ((iStart+1)*dt[i] - last)*gradWave[i][iStart][j];
            
            
            // add end partial block
            Gdt[j]+= (t - (iFinish*dt[i]))*gradWave[i][iFinish][j];
        }
        
        // if these are not in adjacent waveform blocks, add the ones in the middle
        if(iFinish>iStart+1){
            
        	// include blocks between first and last (2nd-(n-1)th)
            for(int block=iStart+1; block<iFinish; block++){
                for(int j=0; j<Gdt.length; j++){   
                    Gdt[j]+= gradWave[i][block][j]*dt[i];
                }
            }
        }
        
        
        return Gdt;
    }

    
    /**
     * reads the lines from the scheme file and constructs the 
     * scheme object containing the waveforms specified.
     * 
     * @param lines raw data from scheme file
     * 
     * @return new scheme object
     */
    protected static GradientWaveform_Scheme readScheme(Vector<String> lines){
        
        int numMeas= lines.size();
        
        if(numMeas==0){
            throw new LoggedException("no data in scheme file.");
        }
        
        // initialise waveforms array
        double[][][] gradWave= new double[numMeas][][];
        
        // initialise time increment array
        double[] dt= new double[numMeas];
        
        for(int i=0; i<numMeas; i++){
            
            // get next line
            String line= lines.get(i);
            
            // new line scanner object
            Scanner lineScanner = new Scanner(line);
            
            // number of slices
            int K= (int)lineScanner.nextDouble();
            
            // time increment
            dt[i] = lineScanner.nextDouble();
            
            gradWave[i]= new double[K][3];
            
            for(int j=0; j<K; j++){
                for(int k=0; k<3; k++){
                    gradWave[i][j][k]= lineScanner.nextDouble();
                }
            }
            
        }
        
        return new GradientWaveform_Scheme(gradWave, dt);
        
    }
    
    
    /**
     * testing entrypoint. reads a GEN scheme and calculates b-values
     */
    public static void main(String[] args){
        
        System.err.println("Testing generalised gradient scheme files");
        
        System.err.print("Reading scheme...");
        DW_Scheme scheme= DW_Scheme.readScheme("gen_test.scheme");
        System.err.println("done.");
        
        System.err.println("Checking b value calulation");
        
        int numAqs= scheme.numMeasurements();
        int numZeros= scheme.numZeroMeasurements();
        
        System.err.println("parent scheme expects "+numAqs+" measurements, of which "+numZeros+" are zero");
        System.err.println("b-values are as follows:");
        for(int i=0; i<numAqs; i++){
            System.err.println(i+") "+scheme.getB_Value(i));
        }
    }
}

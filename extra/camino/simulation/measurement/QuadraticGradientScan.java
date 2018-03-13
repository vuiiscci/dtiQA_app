package simulation.measurement;

import java.io.BufferedWriter;

//import com.sun.istack.internal.logging.Logger;
import java.util.logging.Logger;

import imaging.DW_Scheme;
import imaging.RectQuadraticGradSteTanScheme;
import simulation.DiffusionSimulation;
import simulation.dynamics.Walker;
import simulation.geometry.substrates.Substrate;
import tools.CL_Initializer;

/**
 * alternate version of the AgnosticScan that can cope with 
 * quadratic gradients. Same as superclass but with overridden
 * phase shift and signal calulation.
 * 
 * @author matt (matt.hall@ucl.ac.uk)
 *
 */
public class QuadraticGradientScan extends AgnosticScan {

	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** coords of centre of substrate */
	private final double[] C= new double[DiffusionSimulation.D];
	
	/** array of initial phases */
	private final double[][] dphi0;
	
	/** local copy of scheme - enforces type */
	private final RectQuadraticGradSteTanScheme scheme;
	
    /** gyromagnetic ratio */
    private final double GAMMA= DW_Scheme.GAMMA;
    
    /** flag indicating if we've initialised the initial phase offsets */
    private boolean dphi0Set= false;
    
    /** array of walkers in the sim */
    private final Walker[] walker;
    
    /** substrate */
    private final Substrate substrate;
    
    
	public QuadraticGradientScan(RectQuadraticGradSteTanScheme scheme, Walker[] walker,
			Substrate substrate) {
		super(scheme, walker, substrate);
		
		this.scheme=scheme;
		
		this.walker=walker;
		
		this.substrate= substrate;
		
		// set the substrate size
		double[] subsize= substrate.getSubstrateSize();
		
		for(int i=0; i<DiffusionSimulation.D; i++){
			C[i]=subsize[i]/2;
		}
		
		// initialise dphi0 array
		this.dphi0= new double[scheme.numMeasurements()][walker.length];
	}

	
	
	
    /**
     * calculates a phase shift for the given spin at the given time in
     * the given direction
     * 
     * @param walker the spin
     * @param t the current time
     * @param dir direction index (line in scheme file)
     * @param last time of last call
     * 
     * @return phase increment for spin
     */
    public double getPhaseShift(Walker walker, double t, int dir, double last) {    	
    	
        // gradient strength times duration in gradient pulse
        double Gdt = scheme.getGradientWeight(dir, t, last);
        double[] Gdir= scheme.getG_Dir(dir);
                
        double gradDotPos = 0.0;
        //double gradDotR0= 0.0;
        
                
        for(int i=0; i<Gdir.length; i++){
            gradDotPos+= Gdir[i]*(walker.r[i]-C[i]);
          //gradDotR0+= Gdir*(r0[i]-C[i]);
        }
        
        // mod operation here is rare, but PI is correct, not 2PI (see mapToCircle() comment)
        //return (GAMMA*gradDotPos*gradDotPos*Gdt)%(2*Math.PI);
        return (GAMMA*gradDotPos*gradDotPos*Gdt)%(2*Math.PI);
    }


    /**
     * returns the phase shift using current and initial position rather
     * than a Walker object. This is used in conjuction with trajfiles
     * and the Scan command.
     * 
     * @param r current walker position
     * @param r0 initial walker position
     * @param t current time
     * @param dir gradient direction index
     * @param last time of last call
     * 
     * @return the phase shift for the walker during the specified update
     */
    public double getPhaseShift(double[] r, double[] r0, double t, int dir, double last) {
        
        // gradient strength times duration in gradient pulse
        double Gdt = scheme.getGradientWeight(dir, t, last);
        double[] Gdir= scheme.getG_Dir(dir);
        
        double gradDotPos = 0.0;
        //double gradDotR0= 0.0;
        
        for(int i=0; i<Gdir.length; i++){
            gradDotPos+= Gdir[i]*(r[i]-C[i]);
            //gradDotR0+= Gdir*(r0[i]-C[i]);
        }
        
        // mod operation here is rare, but PI is correct, not 2PI (see mapToCircle() comment)
        return (GAMMA*gradDotPos*gradDotPos*Gdt)%(2*Math.PI);
    }

    
    
    public double[] getSignals() {
        
    	int numMeas= scheme.numMeasurements();
    	
        double[] signal=new double[numMeas];
        
        double S=0.0;
        
        /*BufferedWriter phaseWriter= null;
        
        try{
        	phaseWriter= new BufferedWriter(new FileWriter("phases.csv")); 
        }
        catch(IOException ioe){
        	throw new LoggedException(ioe);
        }*/
        
        
        logger.info("generating "+numMeas+" signals");
        for(int i=0; i<numMeas; i++){
            
            // reset the signal sums
            S=0.0;
            
            for(int j=0; j<walker.length; j++){
                
                // if this walker isn't in the voxel, skip it
                if(!substrate.voxelContains(walker[j].r)){
                    continue;
                }
                
                double phi=walker[j].getPhaseShift(i);
                
                /*try{
                	phaseWriter.write(phi+",");
                }
                catch(IOException ioe){
                	throw new LoggedException(ioe);
                }*/
                
                // add to sum of signals
                if(scheme.isMagnitude(i)){
                	S+=Math.cos(phi);
                }
                else{
                	S+=Math.sin(phi);
                }
               
            }

            
            /*double noiseTerm=0.0;
            
            double snr= CL_Initializer.SNR;
            double noiseDev= 0.0;
            
            if(snr>0.0){
                noiseDev=((double)walker.length)/snr;
            }
            
            // add noise to net signal
            if(snr>0.0){
                noiseTerm=twister.nextGaussian()*noiseDev;
            }
            Sreal+=noiseTerm;

            if(snr>0.0){
                noiseTerm=twister.nextGaussian()*noiseDev;              
            }
            
            // take modulus of real and imaginary parts to provide measured signal
            if(snr>0.0){
            	signal[i]=Math.sqrt(Sreal*Sreal + noiseTerm*noiseTerm);
            }
            else{
            	signal[i]=Sreal;
            }*/
            
            //double signalInt= SintReal;
            //double signalExt= SextReal;

            //double b= ((DW_Scheme)scheme).getB_Value(i);
            //double Sfree= Math.exp(-b*CL_Initializer.DIFF_CONST);
            
            signal[i]=S;
            
            logger.info("signal ["+i+"] = "+signal[i]);
            
            /*try{
            	phaseWriter.write("\n");
            }
            catch(IOException ioe){
            	throw new LoggedException(ioe);
            }*/

        }
        
        /*try{
        	phaseWriter.flush();
        	phaseWriter.close();
        }
        catch(IOException ioe){
        	throw new LoggedException(ioe);
        }*/
        
        
        return signal;
    }

    
}

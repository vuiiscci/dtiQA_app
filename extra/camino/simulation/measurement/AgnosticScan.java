package simulation.measurement;

import imaging.DW_Scheme;
import imaging.RectGradSteTanScheme;
import imaging.SimulableScheme;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.logging.Logger;

import apps.Executable;

import data.DataSource;
import data.OutputManager;

import misc.LoggedException;
import numerics.MTRandom;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.Walker;
import simulation.geometry.substrates.Substrate;
import simulation.measurement.ScanFactory.ScanType;
import tools.CL_Initializer;

/**
 * agnostic scan is a passive scan module that is 
 * completely agnostic towards the actual pulse-sequence 
 * defined in the scheme file. it simply accumulates
 * phase shifts for all spins according to the getGradStrength()
 * method in SimulableScheme and then reads out the signals
 * by summing over contributions from all spins at the end of
 * the scan. 
 * 
 * It can work with any SimulableScheme object, including those
 * with non-linear gradients.
 *
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class AgnosticScan implements SyntheticScan {

    /** logging object */
    private static final Logger logger= Logger.getLogger("simulation.measurement.AgnosticScan");
    
    /** dimensionality of space */
    private final int D=DiffusionSimulation.D;
    
    /** scheme to use */
    private final SimulableScheme scheme;
    
    /** gyromagnetic ratio */
    private final double GAMMA= DW_Scheme.GAMMA;
    
    /** number of measurements */
    private final int numMeas;
    
    /** array of walkers updated by simulation main loop */
    private final Walker[] walker;
    
    /** the substrate that walkers are contained in */
    private final Substrate substrate;
    
    /** random number generator */
    private final MTRandom twister=new MTRandom((1736401757<<32)|(CL_Initializer.seed));

    /** trajfile reader */
    private DataInputStream trajReader=null;
    
    /** trajfile header - duration */
    private double traj_dur;
    
    /** trajfile header - num walkers */
    private int N_walkers;
    
    /** trajfile header - number of timesteps */
    private double traj_tmax;
    
    /** trajfile header - time increment */
    private double traj_dt;
    
    /** a little class that hold a line from a trajfile */
    private final class TrajFileLine{
    	
    	private double t;
    	private int i;
    	private double[] r= new double[D];
    	
    }
    
    
    /**
     * constructor. needs scheme.
     * 
     * @param scheme SimulableScheme to ask for gradient strengths
     * @param walker array of spins
     */
    public AgnosticScan(SimulableScheme scheme, Walker[] walker, Substrate substrate){
        
        this.scheme=scheme;
        
        numMeas= ((DW_Scheme)scheme).numMeasurements();
               
        this.walker=walker;
        
        this.substrate= substrate;
        
    }
    
    /**
     * constructor for trajfile data synthesis. needs scheme and trajfile name
     * 
     * @param scheme simulable scheme to use for synthesis
     * @param trajfile trajectories to use
     * 
     */
    public AgnosticScan(SimulableScheme scheme, String trajfile){
    	
    	this.scheme=scheme;
    	
    	numMeas= ((DW_Scheme)scheme).numMeasurements();
    	
    	this.walker=null;
    	
    	this.substrate=null;
    	
    	// open the trajfile
    	FileInputStream fStream;
    	try{
    		fStream= new FileInputStream(trajfile);
    	}
    	catch(IOException ioe){
    		throw new LoggedException(ioe);
    	}
    	
    	this.trajReader= new DataInputStream(new BufferedInputStream(fStream, 1048576));
    	
    	// read the trajfile header
    	try {
    		
            traj_dur = trajReader.readDouble();
            N_walkers= (int)(trajReader.readDouble());
            traj_tmax= (int)(trajReader.readDouble());
            traj_dt= traj_dur/traj_tmax;
            
        } catch (IOException ioe) {
            throw new LoggedException(ioe);
        }
    	
    }
    
    
    public double[] getCompartmentalSignals(boolean intra) {
        // TODO Auto-generated method stub
        double[] signal=new double[numMeas];
        
        double Sreal;        

        int numIn=0;
        int numExt=0;
        
        double noiseTerm;
                
                
        // now calculate the remainder of the signals
        logger.info("generating "+numMeas+" signals");
        for(int i=0; i<numMeas; i++){
            
            // reset the signal sums
            Sreal=0.0;

            numIn=0;
            numExt=0;
            
            
            double min=Double.MAX_VALUE;
            double max=-Double.MAX_VALUE;
                        
            for(int j=0; j<walker.length; j++){
                
                double phi=walker[j].getPhaseShift(i);
                
                if(phi<min){
                    min=phi;
                }
                
                if(phi>max){
                    max=phi;
                }
                
                // add to sum of signals 
                Sreal+=Math.cos(phi);
               
                // check if the signals come from intra
                // or extra cellular compartments and
                // update the appropriate signal
                if((substrate.intracellular(walker[j]))==intra){
                    Sreal+=Math.cos(phi);
                    
                    numIn++;
                }
            }

            
            noiseTerm=0.0;
            
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
            signal[i]=Math.sqrt(Sreal*Sreal + noiseTerm*noiseTerm);
            
        }
        
        return signal;
    }
    
    

    

    /** 
     * number of measurements in the scan. same as length of gDirs array
     * 
     * @return number of measurements
     */
    public int getNumMeasurements() {
        return numMeas;
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
        double[] Gdt = scheme.getGradImpulse(dir, t, last);
        
        double gradDotPos = 0.0;
                
        for(int i=0; i<Gdt.length; i++){
            gradDotPos+= Gdt[i]*(walker.r[i]-walker.r0[i]);
        }
        
        // mod operation here is rare, but PI is correct, not 2PI (see mapToCircle() comment)
        return (GAMMA*gradDotPos)%(2*Math.PI);
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
        double[] Gdt = scheme.getGradImpulse(dir, t, last);
        
        double gradDotPos = 0.0;
                
        for(int i=0; i<Gdt.length; i++){
            gradDotPos+= Gdt[i]*(r[i]-r0[i]);
        }
        
        // mod operation here is rare, but PI is correct, not 2PI (see mapToCircle() comment)
        return (GAMMA*gradDotPos)%(2*Math.PI);
    }

    /** 
     * maps a phase into the interval [-2pi, 2pi). This was originally a helper method in {@link Walker}
     * and is required to avoid horrendous rounding errors in accumulating phase shifts. Moved it here
     * so that getSignalsFromTrajectories() and Walker.update() use the same piece of code.
     * 
     * @param rawPhase number on unbounded interval
     * 
     * @return rawPhase mod 2pi
     */
    public static final double mapToCircle(double rawPhase){
    	
    	/** 
    	 * OLD COMMENT: mod pi, not mod 2pi - mods of negative numbers are negative
    	 * (at least in Java they are - this behaviour is ambiguous in
    	 * the IEEE floating point standard and is not conserved across 
    	 * languages - beware!), this operation maps all phases into 
    	 * the range [-pi,pi], which is what we want.
    	 * 
    	 * EDIT 26/3/13: Changed the range of the Modulus operation and this was causing problems when
    	 * mapping numbers greater than PI but less than 2*PI. They were getting mapped to between
    	 * zero and PI rather than -PI and -0, which resulted in phases with the wrong sign.
    	 * 
    	 */
    	return rawPhase%(2*Math.PI);
    	
    }
    

    
    /**
     * returns an agnostic identifier
     */
    public ScanType getScanType() {
        return ScanType.AGNOSTIC;
    }

    
    public double[] getSignals() {
        
        double[] signal=new double[numMeas];
        
        double Sreal;
        double SintReal;
        double SextReal;
        
        int numIn;
        int numExt;
        
        BufferedWriter phaseWriter= null;
        
        /*try{
        	phaseWriter= new BufferedWriter(new FileWriter("phases.csv")); 
        }
        catch(IOException ioe){
        	throw new LoggedException(ioe);
        }*/
        
        
        logger.info("generating "+numMeas+" signals");
        for(int i=0; i<numMeas; i++){
            
            // reset the signal sums
            Sreal=0.0;

            SintReal=0.0;
            
            SextReal=0.0;
            
            numIn=0;
            numExt=0;
            
            for(int j=0; j<walker.length; j++){
                
                // if this walker isn't in the voxel, skip it
                if(!substrate.voxelContains(walker[j].r)){
                    continue;
                }
                
                double phi=walker[j].getPhaseShift(i);
                double M= Math.exp(walker[j].getLogMagnetisation(i));
                
                /*try{
                	phaseWriter.write(phi+",");
                }
                catch(IOException ioe){
                	throw new LoggedException(ioe);
                }*/
                
                // add to sum of signals 
                Sreal+=M*Math.cos(phi);
               
                // check if the signals come from intra
                // or extra cellular compartments and
                // update the appropriate signal
                if(substrate.intracellular(walker[j])){
                    SintReal+=M*Math.cos(phi);
                    
                    numIn++;
                }
                else{
                    SextReal+=M*Math.cos(phi);
                    
                    numExt++;
                }
            }

            
            double noiseTerm=0.0;
            
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
            }
            
            double signalInt= SintReal;
            double signalExt= SextReal;

            //double b= (GAMMA*delta[i]*G[i])*(GAMMA*delta[i]*G[i])*(DELTA[i]-delta[i]/3);
            double b= ((DW_Scheme)scheme).getB_Value(i);
            double Sfree= Math.exp(-b*CL_Initializer.DIFF_CONST);
            
            logger.info("signal = "+    signal[i]+"  intra= "+signalInt+
                    " extra= "+signalExt+" ("+numIn+" in, "+numExt+" out)"+" free = "+Sfree);
            
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

    /**
     * Calculates an array of signals from a trajfile and scheme.
     * This is just a matter of accumulating a phase shift for each walker
     * by iterating along its trajectory.
     * 
     * No distinction is made between weighted and unweighted acquisitions.
     * In the case of zero diffusion weighting all phase shifts will be zero 
     * and the calculated signal equal to the number of walkers.
     * 
     * TODO: data synthesis from trajectories has no substrate and so can't do
     * magnetisation vector manipulation. This isn't implemented properly elsewhere
     * anyway, but this will be a shortcoming if we do improve the scan physics
     * 
     * TODO: with current trajfile format there is no information about inner&outer 
     * voxel dimensions. Need to add this to traj file header.
     * 
     * Returns one signal for each line in the scheme file
     * 
     * @return array of signals
     */
    public double[] getSignalsFromTrajectories(){
    	
    	// return array for signals
    	double[] signal= new double[numMeas];
    	
    	// array of initial walker positions (N_w x D)
    	double[][] r0 = new double[N_walkers][D];
    	
    	// array of phase shifts - one per walker per signal
    	double[][] dphi= new double[numMeas][N_walkers];
    	
    	// storage for next line in trajfile
    	TrajFileLine tfl= new TrajFileLine();
    	
    	double schemeDuration= scheme.getDuration();
    	int scheme_tmax= (int)(Math.floor(schemeDuration/traj_dt));
    	
    	double tLast=0.0;
    	
    	BufferedWriter phaseWriter;
    	
        try{
        	phaseWriter= new BufferedWriter(new FileWriter("phases_traj.csv")); 
        }
        catch(IOException ioe){
        	throw new LoggedException(ioe);
        }
        
    	
    	// initial read initial positions and do first update
    	for(int i=0; i<N_walkers; i++){
    		readTrajFileLine(tfl);
    		
    		if(i!=tfl.i){
    			throw new LoggedException("traj file reader and contents are out of sync. aborting.");
    		}
    		
    		for(int j=0; j<D; j++){
    			r0[i][j]=tfl.r[j];
    		}
    		
    		// set initial value of phase shifts
    		for(int j=0; j<numMeas; j++){
    			dphi[j][i]=getPhaseShift(tfl.r, r0[i], tfl.t, j, 0.0);						// accumulate raw phase shift
    			dphi[j][i]=AgnosticScan.mapToCircle(dphi[j][i]);							// avoid rounding issues due to winding
    		}
    		
    		tLast=tfl.t;
    	}
    	
    	// accumulate phase sifts for all walkers across their trajectories
    	//for(int t=1; t<traj_tmax; t++){
    	for(int t=1; t<scheme_tmax; t++){
    	    		
    		for(int i=0; i<N_walkers; i++){
    			
    			readTrajFileLine(tfl);
        		
        		if(i!=tfl.i){
        			throw new LoggedException("traj file reader and contents are out of sync. aborting.");
        		}
        		
        		// accumulate phase shifts
        		for(int j=0; j<numMeas; j++){
        			dphi[j][i]+=getPhaseShift(tfl.r, r0[i], tfl.t, j, tLast);
        			dphi[j][i]=AgnosticScan.mapToCircle(dphi[j][i]);
        		}
    			
    		}
    		
    		// now we've read all the walkers in this timestep, update tLast for next time 
    		tLast=tfl.t;
    	}
    	
    	// finally, calculate signals from walker phase shifts
    	for(int i=0; i<numMeas; i++){
    		
    		// reset signal
    		signal[i]=0.0;

    		for(int j=0; j<N_walkers; j++){
    			
                try{
                	phaseWriter.write(dphi[i][j]+",");
                }
                catch(IOException ioe){
                	throw new LoggedException(ioe);
                }
    			
    			// assume signal is symmetric
    			signal[i]+= Math.cos(dphi[i][j]);
    		}
    		
    		
    		//signal[i]= Math.sqrt(signal[i]*signal[i]);
    		
    		double b= ((DW_Scheme)scheme).getB_Value(i);
            double Sfree= Math.exp(-b*CL_Initializer.DIFF_CONST);
            
            logger.info("signal = "+    signal[i]+", free = "+Sfree);	

            try{
            	phaseWriter.write("\n");
            }
            catch(IOException ioe){
            	throw new LoggedException(ioe);
            }

    	
    	}

        try{
        	phaseWriter.flush();
        	phaseWriter.close();
        }
        catch(IOException ioe){
        	throw new LoggedException(ioe);
        }

    	
    	// return the signals
    	return signal;
    	
    }
    
    /** 
     * reads one record from the trajfile into a TrajFileLine object
     * 
     * @param tfl space to store the variables
     */
    private final void readTrajFileLine(TrajFileLine tfl){
    	
    	try{
    		
    		tfl.t= trajReader.readDouble();
    		tfl.i= trajReader.readInt();
    		for(int i=0; i<D; i++){
    			tfl.r[i]= trajReader.readDouble();
    		}
    		
    	}
    	catch(IOException ioe){
    		throw new LoggedException(ioe);
    	}
    	
    }
    
    
    public void update(int t) {
 
    }

    
}

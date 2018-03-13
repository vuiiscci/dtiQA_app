package simulation.measurement;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Logger;

import data.OutputManager;

import misc.LoggedException;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.Walker;
import tools.CL_Initializer;

/**
 * generates mean-squared displacements
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class MSdispStatsModule extends StatisticsModule {

	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());

	/** dimensionality of space */
	private static final int D= DiffusionSimulation.D;



	/**
	 * constructor, takes array of walkers
	 */
	public MSdispStatsModule(Walker[] walker){

		super(walker);

		Ns= 2*D;

		stats= new double[Ns];
	}

	/**
	 * overrides abstract runtime stats method. this is called at every update
	 * of the simulation and here is used to generate the mean-squared displacements
	 * of spins in each of the three cardinal directions at the current time
	 * 
	 * @param t time in sencods
	 */
	public double[] getRuntimeStats(double t){

		/* NB -- this zeros Ns values, which happens
		 * currently to be the same as D but in future this
		 * may change. don't forget!
		 */
		for(int j=0; j<Ns; j++){
			stats[j]=0.0;
		}
		// accumulate mean-squared displacements
		for(int i=0; i<walker.length; i++){
			for(int j=0; j<D; j++){
				stats[j]+= (walker[i].r[j]- walker[i].r0[j])*(walker[i].r[j]- walker[i].r0[j]);
			}
			// Calculating the kurtosis here as well
			for(int j=D; j<2*D; j++){
				stats[j]+= (walker[i].r[j-D]- walker[i].r0[j-D])*(walker[i].r[j-D]- walker[i].r0[j-D])*(walker[i].r[j-D]- walker[i].r0[j-D])*(walker[i].r[j-D]- walker[i].r0[j-D]);
			}
		}

		// divide by number of walkers to give mean squared displacement
		for(int j=0; j<D; j++){
			stats[j]/=(walker.length);
		}

		// Now dividing by number of walkers and MSD to get kurtosis 
		for (int j = D; j < 2*D; j++) {
			stats[j]/=(walker.length*stats[j-D]*stats[j-D]);
			stats[j]-=3;
		}

		// Now divide the MSD by 2*time to get the diffusion constant
		//for(int j=0; j<D; j++){
		//	stats[j]/=(2*t);
		//}

		return stats;

	}


	/**
	 * overrides the post simulation abstract method. in this instance 
	 * we do nothing and the method returns null.
	 * 
	 * @return null
	 */
	public double[] getPostSimulationStats(){

		return null;
	}

	
	/**
	 * parse a stats file into human readable form
	 * 
	 */
	public static void main(String[] args){
	    
	    String delimiter=" ";
	    
	    CL_Initializer.CL_init(args);
	    
	    if(SimulationParams.sim_statsfile==null){
	        throw new LoggedException("no stats file specified");
	    }
	    
	    
	    double[] line = new double[7];
	    
	    try{
	        DataInputStream inFile= new DataInputStream(new FileInputStream(SimulationParams.sim_statsfile));
	        
	        String outfilename= new String(SimulationParams.sim_statsfile.substring(0, SimulationParams.sim_statsfile.indexOf('.'))+".dat");
	        
	        FileWriter outfile= new FileWriter(outfilename);
	        
	        while(true){
	            // read input data
	            for(int i=0; i<line.length; i++){
	                line[i]=inFile.readDouble();
	            }
	            
	            // output line in human readable format
                outfile.write(line[0]+delimiter);
                System.out.print(line[0]+delimiter);
                
                for(int i=1; i<7; i++){
                    outfile.write(line[i]+delimiter);
                    System.out.print(line[i]+delimiter);    
                }
	            
	            outfile.write("\n");
	            System.out.println();
	        }
	        
	    }
	    catch(EOFException eof){
	        
	    }
	    catch(IOException ioe){
	        throw new LoggedException(ioe);
	    }
	    
	    System.err.println("done.");
	    
	}
	
}

package apps;

import imaging.DW_Scheme;

import java.util.logging.Logger;

//import simulation.SimulationParams;
import tools.CL_Initializer;
import data.DataSource;
import data.VoxelOrderDataSource;

public class GatherStats {

	// logging object
	private static final Logger logger= Logger.getLogger("apps.GatherStats");
	
	// (cheating) x axis values
	private static double[] xaxis= new double[]{	
			                             0.6800,
			                             0.6921,
			                             0.7026,
			                             0.7148,
			                             0.7249,
			                             0.7313,
			                             0.7388,
			                             0.7436,
			                             0.7485,
			                             0.7526
										};
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
				
        // Parse the command line arguments
        CL_Initializer.CL_init(args);
        
                
        CL_Initializer.checkParsing(args);
        int numIncs=1;

        /** 
         * assume that if brownian sim flag is set, that means the
         * increments variable has been set from the commandline
         * so we need to update the default value
         */
        if(CL_Initializer.brownianSimulation==true){
        	CL_Initializer.brownianSimulation=false;
//         	numIncs=SimulationParams.sim_inflamm_increments;
                throw new misc.LoggedException("Simulation is broken :-(");
        }
        
        // Now initialize the acquisition scheme and model parameters.
        CL_Initializer.initImagingScheme();
     
        
        // imaging scheme
        DW_Scheme scheme=CL_Initializer.imPars;
        
        // data source for noiseless signals
        DataSource data= new VoxelOrderDataSource(CL_Initializer.inputFile, 
        		       scheme.numMeasurements(), CL_Initializer.inputDataType);
        
        
        int numMeas=scheme.numMeasurements();
        
        
        
        
        double[][] mean= new double[numMeas][numIncs];
        double[][] meanSq= new double[numMeas][numIncs];
        
        for(int i=0; i<numMeas; i++){
        	for(int j=0; j<numIncs; j++){
        		mean[i][j]=0.0;
        		meanSq[i][j]=0.0;
        	}
        }
        
        int count=0;
        
        // loop over noiseless input voxels
        while(data.more()){
    		for(int j=0; j<numIncs; j++){
        		// get next voxel's worth of data
        		double[] voxel= data.nextVoxel();
        	
        		double S0=voxel[0];
        	
        		// add meassurements to mean & mean square
        		for(int i=1; i<numMeas; i++){
        			mean[i][j]+= voxel[i]/S0;
        			meanSq[i][j]+= (voxel[i]/S0)*(voxel[i]/S0);
        		}
        	}
        	
        	count++;
        }
        
        double[][] stdDev= new double[numMeas][numIncs];
        double meanStdDev=0.0;

	double minStdDev=Double.MAX_VALUE;
	double maxStdDev=0.0;
        
        // divide the means through by the number of measurements
        for(int i=1; i<numMeas; i++){
        	for(int j=0; j<numIncs; j++){
	        	mean[i][j]/=count;
	        	meanSq[i][j]/=count;
	        	
	        	stdDev[i][j]=Math.sqrt(Math.abs(mean[i][j]*mean[i][j]-meanSq[i][j]));
	        	meanStdDev+=stdDev[i][j];
	        	if(stdDev[i][j]<minStdDev){
	        		minStdDev=stdDev[i][j];
	        	}
	        	if(stdDev[i][j]>maxStdDev){
	        		maxStdDev=stdDev[i][j];
	        	}
        	}
        }
        meanStdDev/=numMeas;
        
        logger.info("V_i, mean and std devs of each signal:");
        for(int i=1; i<numMeas; i++){
        	for(int j=0; j<numIncs; j++){
        		//logger.info(""+stdDev[i]);
        		System.err.println(xaxis[j]+","+mean[i][j]+","+stdDev[i][j]);
        	}
        	System.err.println();
        }

        logger.info("mean std dev is "+meanStdDev);
        logger.info("min  std dev is "+minStdDev);
        logger.info("max  std dev is "+maxStdDev);

	}
}

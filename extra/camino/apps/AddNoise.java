package apps;

import java.io.IOException;
import java.util.Random;
import java.util.logging.Logger;

import numerics.MTRandom;

import data.*;
import imaging.DW_Scheme;
import tools.CL_Initializer;

/**
 * app for adding rician noise to noiseless synthetic noiseless 
 * measurements given SNR and scheme file
 * 
 * @author matt
 *
 */
public class AddNoise extends Executable{

	// logging object
	private static final Logger logger= Logger.getLogger("apps.AddNoise");
	
	private int voxels;

	private boolean gaussianNoise;
	
	private Random twister;
	
//	private DataSource data;
	
	
	public AddNoise(String[] args){
		super(args);
	}
	
	public void initDefaultVals(){
		
	}
	
	public void initOptions(String[] args){
	 
		// Set the default input and output type to float
	    CL_Initializer.inputDataType = "float";
	    OutputManager.outputDataType = "float";
	       
	    // Parse the command line arguments
	    CL_Initializer.CL_init(args);
	    CL_Initializer.checkParsing(args);

	    // data source for noiseless signals
	    CL_Initializer.data = new VoxelOrderDataSource(CL_Initializer.inputFile, 1, CL_Initializer.inputDataType);
	}

	
	public void initVariables() {
		
	    // random number generator
	    // Random twister = new MTRandom(CL_Initializer.seed);
		twister = new MTRandom(CL_Initializer.seed);

	    
		
	    // Check noise level has been set
	    if(CL_Initializer.sigma<=0.0){
        	String errMess= new String("Must set the noise level.  Include -sigma <sigma>.");
        	
        	logger.severe(errMess);
        	throw new RuntimeException(errMess);
	    }
        
        // Determine the noise type
		gaussianNoise = false;
	    if(CL_Initializer.noiseType.toLowerCase().equals("gaussian"))
	    	gaussianNoise = true;

	    // number of noise realisations per input voxel (at least 1)
		voxels= Math.max(CL_Initializer.numVoxels, 1);
	}
	    
	 
	public void execute(OutputManager om){
		// loop over noiseless input voxels
	    while(CL_Initializer.data.more()){
        	
        	// get next voxel's worth of data
        	double[] noiselessVoxel= CL_Initializer.data.nextVoxel();
       	
        	// add noise to voxel
        	for(int i=0; i<voxels; i++){
        		double[] noisyVoxel;
        		if(gaussianNoise) {
        			noisyVoxel=DataSynthesizer.addGaussianNoise(noiselessVoxel, CL_Initializer.sigma, twister);
        		}
        		else {
        			noisyVoxel=DataSynthesizer.addNoise(noiselessVoxel, CL_Initializer.sigma, twister);
        		}
        		        		
        		// output noisy voxel
        		om.output(noisyVoxel);
        		
		}	
	    }
	    om.close();
	}
}

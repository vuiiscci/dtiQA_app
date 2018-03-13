package apps;


import java.util.logging.Logger;

import data.*;

import tools.CL_Initializer;

/**
 * module that takes several voxels of data and generates the mean
 * in each direction.
 * 
 * note that this is different to "mean diffusivity", which is 
 * directionally averaged. this command returns the mean of the 
 * unnormalised signal in each direction spearately. a single voxel
 * of data is returned.
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class VoxelMean {

	/** logging object */
	private static final Logger logger= Logger.getLogger("apps.VoxelMean");
	
	/**
	 * reads in voxels one at a time, calculates means, outputs
	 * mean voxel. 
	 * 
	 * @param args command line args
	 */
	public static void main(String[] args) {
		
		/** set defaults */
	    CL_Initializer.inputDataType = "float";
	    OutputManager.outputDataType = "float";
		
	    /** parse commandline */
	    CL_Initializer.CL_init(args);
	    CL_Initializer.checkParsing(args);
	    
	    /** initialise scan parameters */
	    CL_Initializer.initImagingScheme();
	    //CL_Initializer.initDataSynthesizer();
	    
	    /** initialise output data object */
	    OutputManager om= new OutputManager(); 
	    
	    /** get total number of measurements in each voxel */
	    int numElements= CL_Initializer.numElements;
	    
	    /** initialise space for mean values */
	    double[] meanVoxel= new double[numElements];
	    
	    /** reset mean values to zero */
	    for(int i=0; i<meanVoxel.length; i++){
        	meanVoxel[i]=0.0;
	    }
	    
	    /** counter for number of input voxels */
	    int numVoxels=0;
	    
	    /** initialise input data object */
	    DataSource data = ExternalDataSource.getDataSource(CL_Initializer.inputFile, numElements, CL_Initializer.inputDataType);
	    
	    
	    /** read input and calculate totals */
	    while(data.more()){
        	double[] voxel= data.nextVoxel();
        	
        	for(int i=0; i<voxel.length; i++){
		    meanVoxel[i]+=voxel[i];
        	}
        	numVoxels++;
	    }
	    
	    
	    /** divide by number of voxels */
	    if(numVoxels>0){
        	for(int i=0; i<meanVoxel.length; i++){
		    meanVoxel[i]/=numVoxels;
        	}
	    }
	    

	    /** output mean values */
	    om.output(meanVoxel);
	    om.close();
	}

}

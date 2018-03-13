package apps;

import tools.CL_Initializer;
import numerics.*;
import data.BallStick;
import data.OutputManager;

/**
 * Outputs a ball and stick model fanning structure
 * 
 * @author Shahrum
 *
 */

public class FanningGrid extends Executable{

	public FanningGrid(String[] args){
        super(args);
    }
	
	private double diffusivity;
	private double volfrac;
	private int[] gridsize;
	private double centreDist;

     public void initDefaultVals(){
		diffusivity = 0;
		double volfrac = 0;
		int[] gridsize = null;
		double centreDist = 0;

	 }

	public void initOptions(String[] args){
//	public static void main(String[] args) {

		/* 
		 * We want as inputs:
		 *  gridsize: dimensions of our grid (units voxel length) 
		 *  gridcentre: vector showing where circle centre is relative to centre of grid
		 *  diffusivity: diffusivity, assumed to be the same in all voxels in the grid
		 *  volfrac: volume fraction of the fibre population, assumed to be the same in all voxels in the grid
		 *   
		 *   */
    
		// initializing command line parameters
		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
	}
	
	public void initVariables(){
		diffusivity = CL_Initializer.DIFF_CONST; //is this wise?
		volfrac = CL_Initializer.volFrac;
		gridsize = CL_Initializer.dataDims;
		centreDist = CL_Initializer.centreDist;
	}
		
	//	OutputManager om = new OutputManager();
	public void execute(OutputManager om){
	    for (int k = 0; k < gridsize[2]; k++) {
	    	for (int j = 0; j < gridsize[1]; j++) {
	    		for (int i = 0; i < gridsize[0]; i++) {
	    			double[] voxelOrientation = getOrientation(i,j,k,gridsize,centreDist);
	    			BallStick b = new BallStick(diffusivity, volfrac, voxelOrientation);
	    			
	    			double[] informationb = b.getInformationArray();
	    			
	    			double[] fullinformationb = new double[informationb.length+2];
	    			fullinformationb[0]=0;
	    			fullinformationb[1]=0;
	    			for(int ii=2; ii<fullinformationb.length; ii++){
	    				fullinformationb[ii] = informationb[ii-2];
					//				  	System.err.println(fullinformationb[ii]);
	    			}	    			
	    			om.output(fullinformationb);	
	    		}
	    	}
	    }
	    om.close();
	}
	
	private static double[] getOrientation(int i, int j, int k, int[] gridsize, double centreDist) {
	    Vector3D p = new Vector3D(-gridsize[0]/2.0+i+0.5, -gridsize[1]/2.0+j+0.5, centreDist-gridsize[2]/2.0+k+0.5);
	    
	    // (PAC) need to check that modulus is not zero (happens with 1x1x1 grid)
	    double mod = p.mod();

	    if (mod > 0.0) {
		double[] r = new double[3];
		r[0] = p.x / mod;
		r[1] = p.y / mod;
		r[2] = p.z / mod;
		return r;
	    }
	    else {
		return new double[] {1.0, 0.0, 0.0};
	    }
	}
	
}

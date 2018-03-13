package simulation.geometry.substrates;

import imaging.RectGradSteTanScheme;
import imaging.SimulableScheme;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.logging.Logger;

import misc.LoggedException;
import numerics.MTRandom;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.geometry.substrates.SubstrateFactory.SubstrateType;
import tools.CL_Initializer;

public class PercolationSubstrate extends SquashyInflammationSubstrate {

	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** the max number of attempts to place a new cylinder
	 * before giving up
	 */
	private static final int MAX_TRIES=10000;
	
	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;

	/** size of spatial sorting grid */
	private final int gridSize;
	
	/** size of grid elements */
	private double[] gridLen;
	
	/** constructor. calls super constructor and re-does cylinder positions.
	 * 
	 * cylinders are uniform randomly posiitoned in space. positions are checked
	 * such that cylinders are guarenteed not to be overlapping at initial
	 * radius rMin. This means that there is a possibility that it will
	 * not be possible to place the specified number of cylinders in the substrate 
	 * size given. In this case the constructor fits as many as possible given
	 * a max of 10000 tries per cylinder and issues a warning.
	 * 
	 * @param rmin min cylinder radius
	 * @param rmax max cylinder radius
	 * @param numCylinders number cylinders
	 * @param numVoxels number of voxels to synthesise
	 * @param numIncrements number of inflammation increments
	 * @param substrateSize size of substrate
	 */
	public PercolationSubstrate(SimulationParams simParams){
		super(simParams);
		
		double rmin =SimulationParams.sim_cyl_min_r; 
		double rmax =SimulationParams.sim_cyl_max_r;
        int numCylinders =SimulationParams.sim_cyl_dist_size;
        int numVoxels =CL_Initializer.numVoxels;
        int numIncrements =SimulationParams.sim_inflamm_increments;
        double substrateSize =SimulationParams.sim_L;
        int gridSize =SimulationParams.sim_spatial_grid_size;
		
		
		
		this.gridSize=gridSize;
		this.gridLen=new double[D];
		
		for(int i=0; i<D; i++){
			gridLen[i]=substrateSize/gridSize;
		}
		
		double[][] P=new double[numCylinders][D];
		int placed=0; // count the number of cylinders placed so far
		
		// random number generator
		MTRandom rng=new MTRandom(CL_Initializer.seed+98736);
		
		double z=super.getPeakCoord();
		
		for(int i=0; i<numCylinders; i++){
			
			int numTries=0;
			boolean overlapping=true;
			double[] centre= new double[D];
			
			centre[D-1]=z;
			while(overlapping){
				for(int j=0; j<D-1; j++){
					centre[j]=rng.nextDouble()*substrateSize;
				}
				
				numTries++;
				overlapping=false;
				for(int j=0; j<placed; j++){
					if(super.twoDdist(P[j], centre)<=2.0*rmin){
						overlapping=true;
						break;
					}
				}
				if(numTries>=MAX_TRIES){
					break;
				}
				
			}
			if(numTries>=MAX_TRIES){
				break;
			}
			
			P[i]=centre;
			placed++;
		}
		
		if(placed!=numCylinders){
			logger.warning("could only place "+placed+" of "+numCylinders+" cylinders of radius "+rmin+" in substrate of size "+substrateSize);

			double[][] Pnew= new double[placed][];
			
			for(int i=0; i<placed; i++){
				Pnew[i]=P[i];
			}
		
			P=Pnew;
		}
		
		/** replace the cylinder position array in superclass */
		super.P=P;

		/** handle the cloning for cylinders near the substrate edges */
		super.initClones();
		
		/** initialise the spatial sorting grid */
		initSpatialSorting();
		
	}
	
	/**
	 * returns the grid element from spatial coordinates
	 * 
	 * @return spatial sorting grid index
	 */
	private final int getSpatialGridElement(double[] subsCoords){
		
		int[] raw= new int[D-1];
		
		for(int i=0; i<D-1; i++){
			raw[i]=(int)Math.floor(subsCoords[i]/gridLen[i]);
		}
		
		/** WARNING! this only works for D=3 */
		return raw[0]+raw[1]*raw[1];
	}
	
	/**
	 * initialise the maps of grid element to list of intersecting cylinders
	 * 
	 */
	private void initSpatialSorting(){
		
		
	}
	
	
	/**
	 * test construction of percolation substrate
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		
		double rmin=1E-6;
		double rmax=5E-6; 
		int numCylinders=100;
		int numVoxels=10;
		int numIncrements=10; 
		double substrateSize=2.3E-5;
		int gridSize=10;
		
        SimulableScheme scheme;
        try{
            URI uri= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= uri.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }
		
		SimulationParams simParams= new SimulationParams(1, 1000, 1.0, 1, 
				SubstrateFactory.SubstrateType.CYL_1_FIXED, StepType.FIXEDLENGTH,
				1.0, scheme);
		
		SimulationParams.sim_cyl_min_r=rmin;
		SimulationParams.sim_cyl_max_r=rmax;
		SimulationParams.sim_cyl_dist_size=numCylinders;
		CL_Initializer.numVoxels=numIncrements;
		SimulationParams.sim_L=substrateSize;
		SimulationParams.sim_spatial_grid_size= gridSize;
		
		PercolationSubstrate percSubs= new PercolationSubstrate(simParams);
		
		
		
	}

}

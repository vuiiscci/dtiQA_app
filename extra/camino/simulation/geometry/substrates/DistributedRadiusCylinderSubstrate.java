package simulation.geometry.substrates;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.logging.Logger;

import misc.LoggedException;
import numerics.GammaRandom;
import numerics.MTRandom;
import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepGeneratorFactory;
import simulation.geometry.elements.Cylinder;
import simulation.geometry.elements.CylinderFactory;
import simulation.geometry.elements.CylinderFactory.CylType;
import simulation.geometry.elements.NestedCylinder;
import simulation.geometry.elements.SubstrateObject;
import tools.CL_Initializer;

/** 
 *  implements a class of substrate whereby the radii of cylinders
 * 	are drawn from a specified distribution. 
 * 
 *  Currently cylinder radii are drawn from a uniform distribution
 *  and are placed in a square packing with fixed separation R.
 *  
 *  TODO: In the future this class will also implement non-fixed 
 *  cylinder separations and radii drawn from arbitrary distributions
 *  
 * 
 * 
 * 
 * @author matt m.hall@cs.ucl.ac.uk
 *
 */
public class DistributedRadiusCylinderSubstrate extends CylinderSubstrate {

	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** dimensionality of space */
	private static final int D=DiffusionSimulation.D;
	
	/** size of substrate */
	private final double[] l;
	
	/** set of cylinder centres */
	private double[][] P;
	
	/** set of radii of cylinders */
	private double[] radius;
	
	/** border for cloning cylinders */
	//private static final double border= DiffusionSimulation.border;
	
	/** place to store clones */
	private final ArrayList<double[]> clones=new ArrayList<double[]>();
	
	/** place to store the dynamic spac opt map */
	private ArrayList<Cylinder>[] dynamicVoxMap;

	/** list of all cylinders, including clones */
	private final ArrayList<Cylinder> allCyls;


	/** 
	 * constructor.
	 * 
	 * @param subsSize size of substrate
	 * @param N number of cylinders on the substrate
	 * @param k shape param for gamma distributed cylinder radii
	 * @param beta scale param for gamma distributed cylinder radii
	 * @param simParams simulation parameters object
	 * 
	 */
	public DistributedRadiusCylinderSubstrate(SimulationParams simParams){

		super(new double[]{SimulationParams.sim_L, SimulationParams.sim_L, SimulationParams.sim_L}, simParams, true);
				
		double k=CL_Initializer.gamma_k;
		double beta=CL_Initializer.gamma_beta;
		double subsSize= SimulationParams.sim_L;

		// set substrate size
		this.l=new double[]{subsSize, subsSize, subsSize};
		
		int N=SimulationParams.sim_cyl_dist_size;
		
		cyl_p= new double[N];
		
		double p= simParams.getP();
		for(int i=0; i<cyl_p.length; i++){
			cyl_p[i]=p;
		}
		
		allCyls= new ArrayList<Cylinder>();
		
		Cylinder[] cylinder= arrangeCylinders(k, beta, subsSize, N, simParams);

		if(SimulationParams.cylinderType==CylType.NESTED){
			// switch from the basic cylinders used in positioning to nested cyls
			ArrayList<Cylinder> allCyls= new ArrayList<Cylinder>();
			
			// add all cyls to a new global set
			for(int i=0; i<cylinder.length; i++){
				// construct the nested cylinder
				NestedCylinder nestedCyl= (NestedCylinder)CylinderFactory.getCylinder(cylinder[i].getPosition(), cylinder[i].getRadius(), p, CylType.NESTED);
				
				// split it up into constituent cylinders and add them to a global set
				allCyls.addAll(nestedCyl.allCylinders());

			}

			// copy new global set to cylinder array
			cylinder= new Cylinder[allCyls.size()];
			for(int i=0; i<cylinder.length; i++){
				cylinder[i]= allCyls.get(i);
			}
			
			// calculate mean diameter -- this is aslightly funny distribution but this is a decent stab at it
			double meanSize=0.0;
			for(int i=0; i<allCyls.size(); i++){
				
				meanSize += allCyls.get(i).getRadius();
				
			}
			
			meanSize*=(2.0/allCyls.size());
			
			// remake dynamic map
			int[] n_new = new int[]{(int)(L[0]/meanSize), (int)(L[1]/meanSize), 1};
			
			initDynamicSpacOpt(n_new);
			
			for(int i=0; i<cylinder.length; i++){
				int[] cells= getIntersectingCells(cylinder[i].getPosition(), cylinder[i].getRadius());
				addSingleCylToMap(p, cells, cylinder[i]);
			}
			
		}
		
		
		// set the cylinder array in the superclass
		setCylinders(cylinder);
		
		// finalise the spatial optimisation map
		initSpacOptFromDynamicMap();

		
		// dynamic spatial optimisation performed during cylinder arrangement so this part has moved.
		/*logger.info("initialising spatial optimisation");
		int[] n= getSpatialOptGridSize(k,beta,L);
		
		initialiseSpatialOptimisation(n);
		logger.info("done.");
		*/
		
		//drawCrossSection();
	}
	
	
	/** 
	 * samples cylinder radii from the specified gamma distribution and attempts to arrange them in teh speciied 
	 * region such that they don't overlap. Also handles spatial optimisation and cloning.
	 * 
	 * @param k shape param of gamma distribution
	 * @param beta scale param of gamma distribution
	 * @param subsSize size of region to pack into
	 * @param N number of cylinders to try to pack
	 * @param simParams simulation paramteres object
	 * 
	 * @return array of cylinder objects arranged in postions
	 */
	protected final Cylinder[] arrangeCylinders(double k, double beta, double subsSize, int N, SimulationParams simParams){
		
		MTRandom rng= new MTRandom(CL_Initializer.seed);
		GammaRandom grng= new GammaRandom(CL_Initializer.seed+8273l, k, beta);
				
		
		if(N==0){
			logger.warning("number of cylinders is zero!");
		}
		
		// initialise radius array
		radius= new double[N];

		// initialise radii
		for(int i=0; i<N; i++){
			radius[i]= grng.nextGamma();
		}
		

		// check if cross sectional areas are compatible
		double cylArea=0.0;
		for(int i=0; i<radius.length; i++){
			cylArea+=Math.PI*radius[i]*radius[i];
		}
		
		double sqArea=L[0]*L[1];
		if(cylArea>=sqArea){
			logger.warning("total area of specified gamma-distributed cylinders ("+cylArea+") is greater than cross section of substrate ("+sqArea+")");
		}

		logger.info("initialising dynamic spatial optimisation");
		int[] n= getSpatialOptGridSize(k, beta, L);
		initDynamicSpacOpt(n);
		
		// sort radii into ascending order
		Arrays.sort(radius);
		
		// reverse it into descending order
		double[] tempRadius= new double[radius.length];
		for(int i=0; i<radius.length; i++){
			tempRadius[radius.length-i-1]=radius[i];
		}
		radius=tempRadius;
		
		// array of cylinder positions
		P= new double[N][D];

		int interval= N/20;
		
		logger.info("arranging cylinders");
		
		int maxInd=N;
		// place cylinders on substrate
		for(int i=0; i<radius.length; i++){
			boolean overlapping= true;
			int count=0;
			
			if(i%interval==0){
				logger.info("cylinder "+i+" ("+(100*i/N)+"% complete)");
			}
			
			int[] cells= null;
			
			ArrayList<double[]> newClones=null;
			
			while(overlapping && (count<10000)){
				P[i][0]= L[0]*rng.nextDouble();
				P[i][1]= L[1]*rng.nextDouble();
				
				overlapping=false;
				
				// dynamic spac-opt checking
				cells= getIntersectingCells(P[i], radius[i]);
				overlapping= checkIntersections(P[i], radius[i], cells, 
						StepGeneratorFactory.getStepGenerator(simParams).getWalkerRadius());
				
				// check if the new cylinder needs cloning, 
				// and check the clones for overlaps
				if(!overlapping){
					int[] cloneCells;
					
					newClones=checkForCloning(P[i], radius[i]);
	
					for(int j=0; j<newClones.size(); j++){
						cloneCells=getIntersectingCells(newClones.get(j), radius[i]);
						overlapping= checkIntersections(newClones.get(j), radius[i], cloneCells, 
								StepGeneratorFactory.getStepGenerator(simParams).getWalkerRadius());
						if(overlapping){
							break;
						}
					}
				}
	
				count++;
			}
						
			if(count==10000){
				logger.warning("could only place "+(i+1)+" of "+ N +" cylinders on substrate");
				maxInd=i;
				
				double[][] Pnew = new double[maxInd][];
				
				for(int j=0; j<maxInd; j++){
				    Pnew[j]=P[j];
				}
				
				P=Pnew;
				
				break;
			}

			// if we're down here and not overlapping then it should be ok to add the clones
			if(!overlapping){
				addToDynamicMap(P[i], radius[i], cyl_p[i], cells, newClones, CylType.BASIC);
			}
		}
		
		logger.info("dynamic cylinder placment finished");
		
		// calculate intracellular vol frac
		cylArea=0.0;
		for(int i=0; i<maxInd; i++){
			cylArea+=Math.PI*radius[i]*radius[i];
		}
		
		double V_I= cylArea/sqArea;
		
		logger.info("intracellular volume fraction "+V_I);
		
		
		
		// remember how many we've managed to place
		int Nbefore=maxInd;
		
		logger.info("constructing runtime spatial optimisation arrays");
		// set the cylinders arrays here and in superclass
		this.cylinder= new Cylinder[allCyls.size()];
		for(int i=0; i<allCyls.size(); i++){
			// bloody stupid java generics can't do this automatically via toArray(). useless.
			cylinder[i]=allCyls.get(i);
		}
		
		return cylinder;
		
		/*try{
			drawCrossSection();
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}*/
		
	}
	
	
	/*private int[] getIntersectingCells(double[] pos, double r) {
		int[] maxInd= new int[D];
		int[] minInd= new int[D];
		
		for(int i=0; i<D; i++){
			
			double max= pos[i]+r;
			double min= pos[i]-r;
			
			double mincell= Math.floor(min/super.s[i]);
			double maxcell= Math.floor(max/super.s[i]);
			
			minInd[i]= (int)Math.max(mincell, 0);
			maxInd[i]= (int)Math.min(maxcell, super.n[i]-1);
		}
		
		int numCells=1;
		
		for(int i=0; i<D; i++){
			numCells*= (maxInd[i]-minInd[i]+1);
		}
		
		int[] cells= new int[numCells];
		int count=0;
		for(int i=minInd[0]; i<=maxInd[0]; i++){
			for(int j=minInd[1]; j<=maxInd[1]; j++){
				for(int k=minInd[2]; k<=maxInd[2]; k++){
					cells[count]=super.getSubVoxelIndex(i, j, k);
					count++;
				}
			}
		}

		return cells;

	}*/


	/** returns the size of the substrate (N*R) 
	 */
	public double[] getSubstrateSize() {
		return L;
	}


	/** returns the coord for the delta peak (centre of substrate)
	 */
	public double getPeakCoord() {
		return L[0]/2;
	}

	
	/**
	 * check if a cylinder needs to be cloned
	 */
	private ArrayList<double[]> checkForCloning(double[] P, double radius){
		
		ArrayList<double[]> added= new ArrayList<double[]>();
		
		for(int j=0; j<2; j++){
			// if cylinder lower side is close enough...
			if(P[j]-radius<=border[j]){
				double[] clone= new double[D+2];
				
				for(int k=0; k<D; k++){
					if(k==j){
						clone[k]=P[k]+(double)(L[k]);
					}
					else{
						clone[k]=P[k];
					}
				}
				// clone radius
				clone[D]=radius;
				
				//clones.add(clone);
				added.add(clone);
			}

			// if cylinder upper side is close enough...
			if(L[0]-(P[j]+radius)<=border[j]){
				double[] clone= new double[D+2];
				
				for(int k=0; k<D; k++){
					if(k==j){
						clone[k]=P[k]-L[k];
					}
					else{
						clone[k]=P[k];
					}

				}
				clone[D]=radius;
				
				added.add(clone);
			}
		}
		
		// check for corners
		// top-right
		if((L[0]-(P[0]+radius)<=border[0])&&(L[1]-(P[1]+radius)<=border[1])){
			double[] clone= new double[]{P[0]-(double)(L[0]), 
										 P[1]-(double)(L[1]), 
										 P[2], 
					 					 radius};
			
			added.add(clone);
		}
		// bottom-left
		if((P[0]-radius<=border[0])&&(P[1]-radius<=border[1])){
			double[] clone= new double[]{P[0]+(double)(L[0]), 
										 P[1]+(double)(L[1]), 
										 P[2], 
					 					 radius};
			
			added.add(clone);
		}
		// top-left
		if((L[0]-(P[0]+radius)<=border[0])&&(P[1]-radius<=border[1])){
			double[] clone= new double[]{P[0]-(double)(L[0]), 
					 					 P[1]+(double)(L[1]), 
					 					 P[2], 
					 					 radius};

			added.add(clone);
		}
		// bottom-right
		if((P[0]-(radius)<=border[0])&&(L[1]-(P[1]+radius)<=border[1])){
			double[] clone= new double[]{P[0]+(double)L[0], 
					 					 P[1]-(double)L[1], 
					 					 P[2], 
					 					 radius};

			added.add(clone);
		}
		
		return added;
	}

	
	/**
	 * check if a cylinder needs to be cloned
	 */
	/*private ArrayList<double[]> checkForCloning(double[] P, double radius){
		
		ArrayList<double[]> added= new ArrayList<double[]>();
		
		for(int j=0; j<2; j++){
			// if cylinder lower side is close enough...
			if(P[j]-radius<=border){
				double[] clone= new double[D+1];
				
				for(int k=0; k<D; k++){
					if(k==j){
						clone[k]=P[k]+(double)(L[k]);
					}
					else{
						clone[k]=P[k];
					}
				}
				clone[D]=radius;
				
				//clones.add(clone);
				added.add(clone);
			}

			// if cylinder upper side is close enough...
			if(L[0]-(P[j]+radius)<=border){
				double[] clone= new double[D+1];
				
				for(int k=0; k<D; k++){
					if(k==j){
						clone[k]=P[k]-(double)(L[k]);
					}
					else{
						clone[k]=P[k];
					}

				}
				clone[D]=radius;
				
				//clones.add(clone);
				added.add(clone);
			}
		}
		
		// check for corners
		// top-right
		if((L[0]-(P[0]+radius)<=border)&&(L[1]-(P[1]+radius)<=border)){
			double[] clone= new double[]{P[0]-(double)(L[0]), 
										 P[1]-(double)(L[1]), 
										 P[2], 
					 					 radius};
			//clones.add(clone);
			added.add(clone);
		}
		// bottom-left
		if((P[0]-radius<=border)&&(P[1]-radius<=border)){
			double[] clone= new double[]{P[0]+(double)(L[0]), 
										 P[1]+(double)(L[1]), 
										 P[2], 
					 					 radius};
			//clones.add(clone);
			added.add(clone);
		}
		// top-left
		if((L[0]-(P[0]+radius)<=border)&&(P[1]-radius<=border)){
			double[] clone= new double[]{P[0]-(double)(L[0]), 
					 					 P[1]+(double)(L[1]), 
					 					 P[2], 
					 					 radius};
			//clones.add(clone);
			added.add(clone);
		}
		// bottom-right
		if((P[0]-(radius)<=border)&&(L[1]-(P[1]+radius)<=border)){
			double[] clone= new double[]{P[0]+(double)L[0], 
					 					 P[1]-(double)L[1], 
					 					 P[2], 
					 					 radius};
			//clones.add(clone);
			added.add(clone);
		}
		
		return added;
	}*/
	
	
	/** 
	 * handle the clones around the cylinder edges
	 */
	protected void initClones(){
		
		// add the clones to the array of centres
		double[][] Pfull=new double[P.length+clones.size()][];
		
		for(int i=0; i<P.length; i++){
			Pfull[i]=new double[D];
			for(int j=0; j<D; j++){
				Pfull[i][j]=P[i][j];
			}
		}
		
		for(int i=0; i<clones.size(); i++){
			double[] clone=clones.get(i);
			
			Pfull[i+P.length]=clone;
		}
		
		// copy the new array (including clones) to the
		// cylinder centres array
		P=Pfull;
		
		// initialise the cylinders array
		cylinder= new Cylinder[P.length];
	}

	/**
	 * estimates the optimal size of the spatial optimisation grid
	 * given the size of the substrate and the distribution of the cylinder
	 * radii - from experiments, the grid should be approximately the size of
	 * the mean axon radius
	 * 
	 */
	private final int[] getSpatialOptGridSize(double kappa, double beta, double[] subsSize){
		
		//calculate mean radius size
		double meanRad = kappa*beta;
		// x,y dimensions of the grid
		
		int[] n = new int[]{(int)(Math.ceil(subsSize[0]/meanRad)), (int)(Math.ceil(subsSize[1]/meanRad)), 1};
		
		return n;
	}
	
	/**
	 * writes the cross sections of each cylinder to a file
	 *
	 */
	public final void drawCrossSection(){
		
		int Nincs=20;
		
		BufferedWriter out;
		
		try{
			out= new BufferedWriter(new FileWriter("gamma_cyls.csv"));
		
			
			for(int i=0; i<cylinder.length; i++){
				
				double r= cylinder[i].getRadius();
				for(int t=0; t<Nincs; t++){
					double theta1= (((double)t)/((double)Nincs))*2.0*Math.PI;
					double theta2= (((double)(t+1))/((double)Nincs))*2.0*Math.PI;
					
					double x1= P[i][0]+r*Math.cos(theta1);
					double y1= P[i][1]+r*Math.sin(theta1);
					
					double x2= P[i][0]+r*Math.cos(theta2);
					double y2= P[i][1]+r*Math.sin(theta2);
					
					out.write(x1+","+y1+"\n");
					out.write(x2+","+y2+"\n");
				}
			}
			
			out.flush();
			out.close();
			
		}
		catch(IOException ioe){
			throw new RuntimeException(ioe);
		}
	}
		
	
	
	
	/**
	 * initialises an array of empty ArrayLists with a guess at
	 * their lengths assuming a roughly equal number of cylinders 
	 * in each grid cell.
	 * 
	 * @param n array of spac opt grid dimensions
	 */
	private final void initDynamicSpacOpt(int[] n){
		
		int numCells=1;
		
		for(int i=0; i<D; i++){
			super.n[i]=n[i];
			s[i]= SimulationParams.sim_L/n[i];
			numCells*=n[i];
		}
		
        voxToObjects= new SubstrateObject[numCells][];
        candidateSubVox= new int[voxToObjects.length];
		
		this.dynamicVoxMap= new ArrayList[numCells];
		
		for(int i=0; i<numCells; i++){
			// initialise 
			dynamicVoxMap[i]= new ArrayList<Cylinder>((int)(SimulationParams.sim_cyl_dist_size/numCells));
		}
	}
	
	
	private final int[] getIntersectingCells(double[] pos, double r){
		
		int[] maxInd= new int[D];
		int[] minInd= new int[D];
		
		for(int i=0; i<D; i++){
			
			double max= pos[i]+r;
			double min= pos[i]-r;
			
			double mincell= Math.floor(min/s[i]);
			double maxcell= Math.floor(max/s[i]);
			
			minInd[i]= (int)Math.max(mincell, 0);
			maxInd[i]= (int)Math.min(maxcell, n[i]-1);
		}
		
		int numCells=1;
		
		for(int i=0; i<D; i++){
			if(maxInd[i]>minInd[i]){
				numCells*= (maxInd[i]-minInd[i]+1);
			}
		}
		
		int[] cells= new int[numCells];
		int count=0;
		for(int i=minInd[0]; i<=maxInd[0]; i++){
			for(int j=minInd[1]; j<=maxInd[1]; j++){
				for(int k=minInd[2]; k<=maxInd[2]; k++){
					int cellInd= super.getSubVoxelIndex(i, j, k);
					//if((cellInd>=0)&&(cellInd<numCells)){
						cells[count]= cellInd;
						count++;
					//}
				}
			}
		}

		return cells;
	}
	
	
	private final boolean checkIntersections(double[] pos, double r, int[] cells, double walkerRadius){
		
		/* check cylinder proximity to boundaries - if the edge of the cylinder is within
		 * two walker radii of the the permeable voxel boundary this will cause a leak
		 * in the cylinder when walkers transition the boundary.
		 * 
		 * this isn't solved by cloning, as a cylinder that's on the far side of the boundary
		 * won't be cloned (nor should it be -- too many clones is inefficient and capturing these
		 * would require enlarging the subvoxel grid).
		 */
		for(int i=0; i<2; i++){
			
			// check lower boundary
			if((pos[i]-r>0.0) && (pos[i]-r<2.0*walkerRadius)){
				return true;
			}
			
			// check upper boundary
			if((pos[i]+r<border[i]) && (pos[i]+r-border[i]<2.0*walkerRadius)){
				return true;
			}
			
		}
		
		// check intersections with cylinders already placed on the substrate
		for(int i=0; i<cells.length; i++){
			
			int ind= cells[i];
			
			ArrayList<Cylinder> cyls= dynamicVoxMap[ind];
			
			Iterator<Cylinder> cylsIt= cyls.iterator();
			while(cylsIt.hasNext()){
				Cylinder cyl= cylsIt.next();
				
				double dist= cyl.getDistanceFrom(pos);
				//if(dist<=cyl.getRadius()+r){
				// cylinders need to be non-overlapping and a walker diameter apart 
				if(dist<=cyl.getRadius()+r+2.0*walkerRadius){
					return true;
				}
			}
		}
		
		return false;
		
	}
	
	
	
	private final void initSpacOptFromDynamicMap(){
		
		super.voxToObjects= new SubstrateObject[dynamicVoxMap.length][];
		
		for(int i=0; i<voxToObjects.length; i++){
			if(dynamicVoxMap[i].size()>0){
				// transfer objects out of the ArrayList and into a static array
				voxToObjects[i]= new SubstrateObject[dynamicVoxMap[i].size()];
				for(int j=0; j< dynamicVoxMap[i].size(); j++){
					voxToObjects[i][j]= dynamicVoxMap[i].get(j);
				}
			}
			else{
				voxToObjects[i]=null;
			}
		}
		
		checked= new SubstrateObject[allCyls.size()];
	}
	

	private final void addSingleCylToMap(double p_perc, int[] cells, Cylinder cylinder){
		
		for(int i=0; i<cells.length; i++){
			dynamicVoxMap[cells[i]].add(cylinder);
		}
		
		allCyls.add(cylinder);
		
		ArrayList<double[]> clones= checkForCloning(cylinder.getPosition(), cylinder.getRadius());
		
		
		if(clones!=null){
			for(int i=0; i<clones.size(); i++){
				int[] cloneCells= getIntersectingCells(clones.get(i), clones.get(i)[D]);
				
				//SquashyCylinder clone= new SquashyCylinder(clones.get(i), clones.get(i)[D], p_perc);
				Cylinder clone= CylinderFactory.getCylinder(clones.get(i), clones.get(i)[D],  p_perc, SimulationParams.cylinderType);
				
				
				for(int j=0; j<cloneCells.length; j++){
					dynamicVoxMap[cloneCells[j]].add(clone);
				}
				
				allCyls.add(clone);
			}
		}
		
	}
	
	
	private final void addSingleCylToMap(double[] pos, double r, double p_perc, int[] cells, Cylinder cylinder){
		
		for(int i=0; i<cells.length; i++){
			dynamicVoxMap[cells[i]].add(cylinder);
		}
		
		allCyls.add(cylinder);
		
		ArrayList<double[]> clones= checkForCloning(pos, r);
		
		
		if(clones!=null){
			for(int i=0; i<clones.size(); i++){
				int[] cloneCells= getIntersectingCells(clones.get(i), clones.get(i)[D]);
				
				//SquashyCylinder clone= new SquashyCylinder(clones.get(i), clones.get(i)[D], p_perc);
				Cylinder clone= CylinderFactory.getCylinder(clones.get(i), r,  p_perc, SimulationParams.cylinderType);
				
				
				for(int j=0; j<cloneCells.length; j++){
					dynamicVoxMap[cloneCells[j]].add(clone);
				}
				
				allCyls.add(clone);
			}
		}
		
	}
	
	
	private final void addToDynamicMap(double[] pos, double r, double p_perc, int[] cells, ArrayList<double[]> clones, CylType cylType){
				
		Cylinder cylinder= CylinderFactory.getCylinder(pos, r, p_perc, SimulationParams.cylinderType);

		addSingleCylToMap(pos, r, p_perc, cells, cylinder);
			
	}

	
	/**
	 * override the positionOK method from superclass to utilise spatial optimisation
	 * 
	 * @param r0 walker position
	 * @para R walker radius
	 */
	public boolean positionOk(double[] r0, double R){
    	
    	
    	double minDist=Double.MAX_VALUE;
    	
    	int[] cell= getIntersectingCells(r0, R);
    	
    	for(int i=0; i<cell.length; i++){
    		
    		SubstrateObject[] cyls= voxToObjects[cell[i]];
    		
    		if(cyls!=null){
	    		for(int j=0; j<cyls.length; j++){
	    			double distToObj= cyls[j].getDistanceFrom(r0);
				
		    		if(distToObj<minDist){
		    			minDist=distToObj;
		    		}
	    		}
    		}
    	}
    	  	
    	return (minDist>R);
    	
    }
		
	
	
	
}

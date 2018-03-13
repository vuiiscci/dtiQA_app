package simulation.geometry.substrates;

import imaging.SimulableScheme;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.logging.Logger;

import org.omg.CORBA.portable.OutputStream;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepAmender;
import simulation.dynamics.StepGenerator;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.dynamics.Walker;
import simulation.geometry.elements.Cylinder;
import simulation.geometry.elements.SquashyCylinder;
import simulation.geometry.elements.SquashyCylinder.Chord;
import simulation.geometry.elements.SubstrateObject;
import simulation.geometry.substrates.SubstrateFactory.SubstrateType;
import simulation.measurement.SyntheticScan;
import tools.CL_Initializer;

import misc.LoggedException;
import numerics.GammaRandom;
import numerics.MTRandom;

public class SquashyInflammationSubstrate extends CylinderSubstrate {

	/** for debugging */
	public boolean chuck=false;
	
	/** logging object */
	private final Logger logger=Logger.getLogger(this.getClass().getName());
	
	/** random number generator */
	private final MTRandom rng= new MTRandom(CL_Initializer.seed);
	
	/** dimensionality of space */
	private final int D=DiffusionSimulation.D;
	
	/** size of substrate */
	private final double[] L;
	
	/** array of cylinders */
	public SquashyCylinder[] cylinder;
	
	/** min cylinder radius */
	private double[] radius;
		
	/** radius increment at each step */
	private double[] rinc;
	
	/** number of steps in incrementing  cylinder size */
	private final int numIncrements;
	
	/** number of steps between re-sizing cylinders */
	private int inflammationIncrementModulus;
	
	/** counter for the number of times the substrate initialiser function is called */
	//private int calls=0;
	
	/** index of current increment */
	private int n=0;
	
	/** array of cylinder positions */
	public double[][] P;
	
	/** counter for miscellainious purposes */
	private int counter=0;
	
	/** records which cylinder an intersection was with */
	private int[] cylCrossed;
	
	/** list of shortest distances to barrier crossings */
	private double[] shortestDist;
	
	/** list of membrane permeabilities */
	private double[] cyl_p;
	
	/** list of d[0] values for barrier intersections */
	private double[] intOriginDist;
	
	/** list of intersection normals for barrier intersections */
	private double[][] intNormal;
	
	/** list of initially in or out flags for intersections */
	private boolean[] intInside;
	
	/** which cylinder was the most recent detected intersection with? */
	private int lastCrossed=-1;	
	
	/** membrane permiability */
	private final double p;
	
	/** expansion coeff */
	private final double K=1E4;
	
	/** gamma random number generator */
	private final GammaRandom grng;
	
	/** place to store clones */
	private final ArrayList<double[]> clones=new ArrayList<double[]>();

	/** place to store the dynamic spac opt map */
	private ArrayList<Cylinder>[] dynamicVoxMap;

	/** place to store all cylinders as we create them (needs to be dynamic because of clones) */
	private final ArrayList<SquashyCylinder> allCyls;
	
	/** 
	 * number of placed cylinders before cloning
	 */
	private final int Nbefore;
	
	/** 
	 * constructor. needs min and max radius, the number of 
	 * cylinders, and the number of increments 
	 * 
	 */
	public SquashyInflammationSubstrate(SimulationParams simParams) {

		super(new double[]{SimulationParams.sim_L, 
		                   SimulationParams.sim_L, 
		                   SimulationParams.sim_L}, simParams, true);
		
		double k=CL_Initializer.gamma_k;
		double beta=CL_Initializer.gamma_beta;
		int numCylinders= SimulationParams.sim_cyl_dist_size;
		int numVoxels=CL_Initializer.numVoxels;
		int numIncrements= SimulationParams.sim_inflamm_increments;
		double substrateSize= SimulationParams.sim_L;
		
		
		grng= new GammaRandom(CL_Initializer.seed+8273l, k, beta);
		
		this.numIncrements= numIncrements;
		this.L=new double[]{substrateSize, substrateSize, substrateSize};
		
		if(numVoxels<numIncrements){
			logger.warning("number of voxels less than number if inflammation increments.\n"+
					"inflammation increment modulus =1 will be used.");
			this.inflammationIncrementModulus=1;
		}
		else{
			this.inflammationIncrementModulus=numVoxels/numIncrements;
		}
		
		P= new double[numCylinders][D];
		
		cylCrossed= new int[numCylinders];
		cyl_p= new double[numCylinders];
		shortestDist= new double[numCylinders];
		intOriginDist= new double[numCylinders];
		intNormal= new double[numCylinders][D];
		intInside= new boolean[numCylinders];
		
		allCyls= new ArrayList<SquashyCylinder>(numCylinders);

		this.p= simParams.getP();
		for(int i=0; i<cyl_p.length; i++){
			cyl_p[i]=p;
		}

		logger.info("arranging cylinders");
		
		// initialise radius array
		radius= new double[numCylinders];

		// initialise radii
		for(int i=0; i<numCylinders; i++){
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
		// set the increments
		rinc= new double[radius.length];
		
		int interval= numCylinders/20;
		
		for(int i=0; i<radius.length; i++){
			double r= radius[i];
			
			double a= Math.PI*r*r;
			
			double inc= K*a;
			
			rinc[i]=inc;
		}
		// array of cylinder positions
		P= new double[numCylinders][D];
		
		logger.info("arranging cylinders");
		
		int maxInd=numCylinders;
		// place cylinders on substrate
		for(int i=0; i<radius.length; i++){
			boolean overlapping= true;
			int count=0;
			
			if(i%interval==0){
				logger.info("cylinder "+i+" ("+(100*i/numCylinders)+"% complete)");
			}
			
			int[] cells= null;
			
			ArrayList<double[]> newClones=null;
			
			while(overlapping && (count<10000)){
				P[i][0]= L[0]*rng.nextDouble();
				P[i][1]= L[1]*rng.nextDouble();
				
				overlapping=false;
				
				// dynamic spac-opt checking
				cells= getInersectingCells(P[i], radius[i]);
				overlapping= checkIntersections(P[i], radius[i], cells, 
						StepGeneratorFactory.getStepGenerator(simParams).getWalkerRadius());
				
				// check if the new cylinder needs cloning, 
				// and check the clones for overlaps
				if(!overlapping){
					int[] cloneCells;
					
					newClones=checkForCloning(P[i], radius[i], i);
	
					for(int j=0; j<newClones.size(); j++){
						cloneCells=getInersectingCells(newClones.get(j), radius[i]);
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
				logger.warning("could only place "+(i+1)+" of "+ numCylinders+" cylinders on substrate");
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
				addToDynamicMap(P[i], radius[i], cyl_p[i], cells, newClones);
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
		Nbefore=maxInd;
		
		logger.info("constructing runtime spatial optimisation arrays");
		// set the cylinders arrays here and in superclass
		this.cylinder= new SquashyCylinder[allCyls.size()];
		for(int i=0; i<allCyls.size(); i++){
			// bloody stupid java generics can't do this automatically via toArray(). useless.
			cylinder[i]=(SquashyCylinder)allCyls.get(i);
		}
		setCylinders(cylinder);
		
		initSpacOptFromDynamicMap();
		
		/*try{
			drawCrossSection();
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}*/
		
		// spit out all the cylinder centres and radii
		/*try{
			FileWriter cylWriter= new FileWriter("allcyls.csv");
			
			for(int i=0; i<cylinder.length; i++){
				
				double[] pos= cylinder[i].getPosition();
				double R= cylinder[i].getRadius();
				
				cylWriter.write(pos[0]+","+pos[1]+","+pos[2]+","+R+"\n");
			}
			
			cylWriter.flush();
			cylWriter.close();
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}*/
		
        logger.info("done.");
		
	}

	/**
	 * Constructs a substrate from a csv file that defines the position and radius of
	 * all cylinders. NB- NO ADDITIONAL CLONING OR CHECKING IS PERFORMED, so you're 
	 * on your own as far as placement is concerned. CSV file should be formatted in 
	 * four columns (x,y,z,R), which is the same as the debgging format in the regular
	 * constructor.
	 * 
	 * @param simParams simulation params object
	 * @param csvfilename name of the csv file to read
	 */
	public SquashyInflammationSubstrate(SimulationParams simParams, String csvfilename){
		
		// sets the substrate size
		super(new double[]{SimulationParams.sim_L, 
                SimulationParams.sim_L, 
                SimulationParams.sim_L}, simParams, true);

		logger.info("setting substrate parameters");
		
		p=SimulationParams.sim_p;
		grng=null;
		allCyls=null;
		numIncrements= simParams.sim_inflamm_increments;
		int numVoxels= CL_Initializer.numVoxels;
		inflammationIncrementModulus= inflammationIncrementModulus=numVoxels/numIncrements;
		
		
		L= new double[]{
				SimulationParams.sim_L, 
                SimulationParams.sim_L, 
                SimulationParams.sim_L};
		
		// somewhere dynamic to store the cylinders we read from the file
		ArrayList<SquashyCylinder> cyls= new ArrayList<SquashyCylinder>(SimulationParams.sim_num_cylinders);
		
		logger.info("reading cylidners file");
		
		// file reader
		Scanner csvReader= null;
		
		// cylinder counter
		int cylCount=0;
		
		// open the file
		try{
			csvReader= new Scanner(new FileInputStream(csvfilename));
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}
		
		// read the file & instantiate cylinders
		try{
			while(csvReader.hasNextLine()){			// end of file?
				String line= csvReader.nextLine();			// read the whole line
		
				StringTokenizer strTok= new StringTokenizer(line, ",");	// tokenise
				
				// construct position array
				double[] P= new double[D];
				for(int i=0; i<D; i++){
					P[i]= Double.parseDouble(strTok.nextToken());
				}
				
				// get radius
				double R= Double.parseDouble(strTok.nextToken());
				
				// add new cylinder to array
				cyls.add(new SquashyCylinder(P, R, SimulationParams.sim_p));
				
				// increment the counter
				cylCount++;
			}
		}
		catch(NoSuchElementException nsee){
			throw new LoggedException(nsee);
		}
		catch(IllegalStateException ise){
			throw new LoggedException(ise);
		}
		finally{
			// close the file reader
			csvReader.close();
		}
		
		logger.info("done. read "+cyls.size()+" cylinders");
		
		logger.info("constructing static array");
		this.cylinder= new SquashyCylinder[cyls.size()];
		
		for(int i=0; i<cyls.size(); i++){
			
			cylinder[i]=cyls.get(i);
		}
		
		Nbefore= cylinder.length;
		
		setCylinders(cylinder);
		
		logger.info("initialising spatial optimisation");
		double k= CL_Initializer.gamma_k;
		double beta= CL_Initializer.gamma_beta;
		
		int[] n= getSpatialOptGridSize(k, beta, L);
		
		initialiseSpatialOptimisation(n);
		logger.info("done");
	}
	
	
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
	}
	
	/**
	 * calculates the ditance between two given points
	 * projected into the plane perp to cylinder axis
	 * 
	 * @param p1 posn of first point
	 * @param p2 posn of second point
	 * 
	 * @return euclidean distance
	 */
	protected final double twoDdist(double[] p1, double[] p2){
		
		double sqDist=0.0;
		for(int i=0; i<2; i++){
			sqDist+=(p1[i]-p2[i])*(p1[i]-p2[i]);
		}
		
		return Math.sqrt(sqDist);
		
	}
		
	/** checks if a walker's step will take it across a membrane or not.
	 * this just involves checking every cylinder in turn. The cylinder
	 * crosses() method should fill in the blanks where necessary.
	 * 
	 * @param walker the walker to check
	 * @param stepVector the step being made
	 * @param normal space to store the normal if a barrier is crossed
	 * @param d space to store barrier distance from origin dotted with normal
	 * @param skipCurrent flag indicating we're sitting on a barrier and should 
	 *        ignore the closest one.
	 * @param originLength the original length of the step vector;
	 * 
	 * @return true or false
	 */
	/*public boolean crossesMembrane(Walker walker, double[] offset, double[] stepVector,
			double[] normal, double[] d, boolean skipCurrent, double origLength, boolean[] in, double[] p) {
		
		double len= 0.0;
		double[] walkerPos= new double[D];
		double[] intDist= new double[1];
		
		boolean[] intIn= new boolean[1];
		double[] intP= new double[1];
		
		getSubstrateCoords(walker.r, offset, walkerPos);
		
		for(int i=0; i<stepVector.length; i++){
			len += stepVector[i]*stepVector[i];
		}
		len=Math.sqrt(len);
		
		if(len/origLength<=1E-14){
			return false;
		}
			
		int numIntersections=0;
		boolean skip=skipCurrent;
		
		for(int i=0; i<cylinder.length; i++){
			
			if(skipCurrent){
				if(i==lastCrossed){
					intIn[0]=in[0];
					skip=true;
				}
				else{
					intIn[0]=cylinder[i].inside(walkerPos);
					skip= false;
				}				
			}
			
			if(cylinder[i].crosses(walkerPos, stepVector, normal, d, skip, origLength, intDist, intIn, intP)){
				
				shortestDist[numIntersections]=intDist[0];
				intOriginDist[numIntersections]=d[0];
				for(int j=0; j<D; j++){
					intNormal[numIntersections][j]=normal[j];
				}
				intInside[numIntersections]=intIn[0];
				cylCrossed[numIntersections]=i;
				cyl_p[numIntersections]=intP[0];
				
				numIntersections++;
				
			}
			
		}
		

		if(numIntersections==1){
			// if there's one interesection everything is  
			// as it should be so just copy the distance
			// and normal over and return true
			d[0]=intOriginDist[0];
			for(int i=0; i<D; i++){
				normal[i]=intNormal[0][i];
			}
			in[0]=intInside[0];
			
			p[0]=cyl_p[0];
			
			lastCrossed=cylCrossed[0];
			
			return true;
		}
		
		if(numIntersections>1){
			
			// if there are more than intersections
			// we must pick the one that would happen
			// first, and set up the distance and 
			// normal accordingly
			double closestDist=Double.MAX_VALUE;
			int closest=-1;
			
			// find closest interesection to the walker
			for(int i=0; i<numIntersections; i++){
				if(shortestDist[i]<closestDist){
					closestDist=shortestDist[i];
					closest=i;
					lastCrossed=cylCrossed[i];
				}
			}
			
			if(closest==-1){
				return false;
			}
			
			// set the distance and normal
			d[0]=intOriginDist[closest];
			for(int i=0; i<D; i++){
				normal[i]=intNormal[closest][i];
			}
			in[0]=intInside[closest];
			p[0]=cyl_p[closest];
			
			return true;
		}
		
		// in this case there are no intersections
		return false;
	}*/

	/**
	 * @return linear size of substrate
	 */
	public double[] getSubstrateSize() {
		// TODO Auto-generated method stub
		return L;
	}

	/** returns the centra of a square substrate
	 * 
	 * @return substrate size L/2.0
	 */
	public double getPeakCoord() {
		return L[0]/2;
	}

	/**
	 * initialise squashy cylinders of given radius at
	 * each centre specified by the P array
	 * 
	 */
	public void init() {
		
	    if(numIncrements!=1){
    		if((DiffusionSimulation.calls%inflammationIncrementModulus == 0)){
    		    
    			for(int i=0; i<cylinder.length; i++){
    				
    				// reinitialise cylinder either with new or same radius
    				if((cylinder[i]==null)||(cylinder[i].isExpanding())){
    					cylinder[i]=new SquashyCylinder(P[i], radius[i]+n*rinc[i], p);
    				}
    				else{
    					cylinder[i]= new SquashyCylinder(P[i], cylinder[i].getRadius(), p);
    					
    					// cylinders are expanding by default. if static in the previous
    					// iteration, must be static in this.
    					cylinder[i].stopExpanding();
    				}
    			}
    		}

			// reassemble clonal relationships
			for(int i=Nbefore; i<P.length; i++){
				int cloneOf= (int)P[i][D+1];
				
				cylinder[cloneOf].addToMyClones(cylinder[i]);
			}
			
			// propagate clone relationships
			for(int i=0; i<Nbefore; i++){
				cylinder[i].propagateCloneRelationships();
				
				// this might seem a little odd! we have to
				// propagate out the non-expanding flag to
				// the clones. when we set it before they
				// weren't affected by the call (but it needed 
				// to be made otherwise we'd have no record of
				// whether the cylinder is expanding or not)
				if(!cylinder[i].isExpanding()){
					cylinder[i].stopExpanding();
				}
			}
			
			
			
			for(int i=0; i<cylinder.length; i++){				
				// add intersections
				for(int j=0; j<i; j++){
					if(i==j){
						continue;
					}
									
					if(cylinder[i].abutts(cylinder[j])){
						if(!cylinder[i].hasIntersectionWith(cylinder[j])){
							cylinder[i].addIntersectionWith(cylinder[j]);
						}
						if(!cylinder[j].hasIntersectionWith(cylinder[i])){
							cylinder[j].addIntersectionWith(cylinder[i]);
						}
					}
				}

				for(int j=0; j<i; j++){
					Iterator<Chord> chIt=((SquashyCylinder)cylinder[j]).chords.iterator();
					
					while(chIt.hasNext()){
						Chord chord= (Chord)chIt.next();

						if(chord.tmin>=chord.tmax){
							System.err.println("added cylinder "+i+", cyl "+j+" has zero-length chord");
							System.err.println("tmin= "+chord.tmin+", tmax= "+chord.tmax);
						}
						
						
						if(Double.isNaN(chord.tmin)){
							System.err.println("added cylinder "+i+", cylinder "+j+" tmin is NaN");
						}
						if(Double.isNaN(chord.tmax)){
							System.err.println("added cylinder "+i+", cylinder "+j+" tmax is NaN");
						}
					}
				}
			}
						
			n++;
			logger.info("iteration "+DiffusionSimulation.calls+"  vol frac= "+intraCellularVolFrac());				

		}
		
		if(SimulationParams.sim_drawCrossSection){
		    try{
		        drawCrossSection();
		    }
		    catch(IOException ioe){
		        throw new LoggedException(ioe);
		    }
		}
		
		// read out the cylinder radii and positions
		if(substrateInfo){
    		logger.info("cylinder positions and radii:");
    		
    		String data= new String();
    		
    		for(int i=0; i<cylinder.length; i++){
    		    double[] pos= cylinder[i].getPosition();
    		    double radius= cylinder[i].getRadius();
    		    
    		    data+="cylinder "+i+" pos = ("+pos[0]+", "+pos[1]+", "+pos[2]+")\tradius = "+radius+"\n";
    		}
    		
    		logger.info("cylinder positions and radii:\n"+data);
    		
    		substrateInfo= false;
		}		
		
		lastCrossed=-1;
		
		DiffusionSimulation.calls++;
	}

	/** 
	 * tests if a position is intracellular or not.
	 * 
	 * @param walker walker whose position to test
	 * 
	 * @return true if inside a cylinder, otherwise false
	 */
	public boolean intracellular(Walker walker){
		
		final double[] substrateCoords= new double[D];
		final double[] noStep= new double[]{0.0, 0.0, 0.0};
		getSubstrateCoords(walker.r, new double[]{0.0,0.0,0.0}, substrateCoords);
		
		assembleSubVoxelList(substrateCoords,noStep);
		
		while(moreCandidates()){
			SquashyCylinder cyl = (SquashyCylinder) nextCandidate();
			
			if(cyl==null){
				continue;
			}
			
			if(cyl.inside(substrateCoords)){
				return true;
			}
		}
		return false;
	}

	
	
	/** numerically calculates the intracellular volume fraction of 
	 *  the substrate by checking a large number of points across the
	 *  cross section of the cylinders (xy-plane)
	 *  
	 *  @return intracellular volume fraction
	 */
	private double intraCellularVolFrac(){
		
		int pts= 100;
		
		double inc=(double)L[0]/(double)pts;
		
		int intraCellularPts=0;
		
		double[] point= new double[D];
		point[2]=0.0;
		
		for(int i=0; i<pts; i++){
			point[0]=i*inc;
			for(int j=0; j<pts; j++){
				point[1]=j*inc;
				
				for(int k=0; k<cylinder.length; k++){
					
					
					if(cylinder[k].inside(point)){
						intraCellularPts++;
						break;
					}
				}
			}
		}
		
		return ((double)intraCellularPts)/((double)(pts*pts));
	}
	
	
	/**
	 * check if a cylinder needs to be cloned
	 */
	private ArrayList<double[]> checkForCloning(double[] P, double radius, int i){
		
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
				
				// remember who you are a clone of
				clone[D+1]=i;
				
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
				clone[D+1]=i;
				
				added.add(clone);
			}
		}
		
		// check for corners
		// top-right
		if((L[0]-(P[0]+radius)<=border[0])&&(L[1]-(P[1]+radius)<=border[1])){
			double[] clone= new double[]{P[0]-(double)(L[0]), 
										 P[1]-(double)(L[1]), 
										 P[2], 
					 					 radius, 
					 					 i};
			
			added.add(clone);
		}
		// bottom-left
		if((P[0]-radius<=border[0])&&(P[1]-radius<=border[1])){
			double[] clone= new double[]{P[0]+(double)(L[0]), 
										 P[1]+(double)(L[1]), 
										 P[2], 
					 					 radius, 
					 					 i};
			
			added.add(clone);
		}
		// top-left
		if((L[0]-(P[0]+radius)<=border[0])&&(P[1]-radius<=border[1])){
			double[] clone= new double[]{P[0]-(double)(L[0]), 
					 					 P[1]+(double)(L[1]), 
					 					 P[2], 
					 					 radius,
					 					 i};

			added.add(clone);
		}
		// bottom-right
		if((P[0]-(radius)<=border[0])&&(L[1]-(P[1]+radius)<=border[1])){
			double[] clone= new double[]{P[0]+(double)L[0], 
					 					 P[1]-(double)L[1], 
					 					 P[2], 
					 					 radius,
					 					 i};

			added.add(clone);
		}
		
		return added;
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
	 * writes the cross sections of each cylinder to a file.
	 * this is done by checking the intersections along the 
	 * edges of a fine grid defined over the substrate.
	 * 
	 *
	 */
	public final void drawCrossSection() throws IOException{
		
		logger.info("generating cylinder cross section image.");
		
		int pts=SimulationParams.sim_crossSectionImageSize;

		// distance between grid points
		double dx = (this.getSubstrateSize())[0]/pts;
		
		// step vectors along each grid line
		double[][] s= new double[][]{{ dx, 0.0, 0.0},
				                     {-dx, 0.0, 0.0},
				                     {0.0,  dx, 0.0},
				                     {0.0, -dx, 0.0}};
		
		// file output
		FileOutputStream fstream = null;
    	try{
            //String fname= new String("crossSec_"+DiffusionSimulation.calls+".gray");
    		String fname= new String("debug_crossSec_"+cylinder.length+".gray");//SimulationParams.crossSectionFig;
    		fstream= new FileOutputStream(fname);
    	}
    	catch(IOException e){
    		throw new LoggedException(e);
    	}
    	
    	DataOutputStream out= new DataOutputStream(new BufferedOutputStream(fstream, 1024));
		
    	SimulationParams simParams= new SimulationParams(1, 1, 0.0, SimulationParams.UNIFORM, SubstrateType.CYL_1_INFLAM, StepType.FIXEDLENGTH, super.L[0], 1.0);
    	StepGenerator stepGen= StepGeneratorFactory.getStepGenerator(simParams);
    	
		double[] pos= new double[D];
		logger.info("instantiating non-simulation walker for cross section. ignore the following warning.");
		Walker walker= new Walker(pos, stepGen);
		logger.info("non-simulation walker instantiated.");
		
		walker.testReplacePositionVector(pos);
		
		double[] offset = new double[D];
		double[] normal = new double[D];
		double[] d = new double[1];
		boolean[] in = new boolean[1];
		double[] pFake= new double[]{1.0};
		
		for(int i=0; i<pts; i++){
			
			// set x xoord
			pos[0]=dx*(i+0.5);
			
			System.err.print("\ri= "+i+"     ");
			for(int j=0; j<pts; j++){
				
				// set y coord
				pos[1]=dx*(j+0.5);
				
				//check intersection
				boolean crosses= false;
				
				for(int k=0; k<s.length; k++){
					
					try {
						crosses=crossesMembrane(walker, offset, s[k], normal, d, false, dx, in, pFake, false, null);
					} catch (StepRejectedException e) {
						e.printStackTrace();
					}
					
					if(crosses){
						break;
					}
				}
				
				if(crosses){
					out.writeByte(0);
				}
				else{
					if(intracellular(walker)){
						out.writeByte(-1);
					}
					else{
						out.write(-128);
					}
				}
			}
		}
		
		out.flush();
		out.close();
		
		System.err.println();
		
		logger.info("cylinder cross section image generated");
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
			super.s[i]= SimulationParams.sim_L/n[i];
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
	
	
	private final int[] getInersectingCells(double[] pos, double r){
		
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
	
	
	private final void addToDynamicMap(double[] pos, double r, double p_perc, int[] cells, ArrayList<double[]> clones){
		
		SquashyCylinder cylinder= new SquashyCylinder(pos, r, p_perc);
		
		for(int i=0; i<cells.length; i++){
			dynamicVoxMap[cells[i]].add(cylinder);
		}
		
		allCyls.add(cylinder);
		
		if(clones!=null){
			for(int i=0; i<clones.size(); i++){
				int[] cloneCells= getInersectingCells(clones.get(i), clones.get(i)[D]);
				
				SquashyCylinder clone= new SquashyCylinder(clones.get(i), clones.get(i)[D], p_perc);
				
				for(int j=0; j<cloneCells.length; j++){
					dynamicVoxMap[cloneCells[j]].add(clone);
				}
				
				cylinder.addToMyClones(clone);
				allCyls.add(clone);
			}
		}
		cylinder.propagateCloneRelationships();
	}
	
	
	
	/**
	 * test function that assembles a known buggy case from a csv file substrate plus
	 * specified information about a walker and step. This is useful if we have a bug that
	 * occurs after a long period of time.  
	 */
	public static void testLongTimeBug(String[] args){
		
		// parse command line
		CL_Initializer.CL_init(args);
		
		// initialise scheme
		CL_Initializer.initImagingScheme();
		
		// get scheme object
		SimulableScheme scheme= (SimulableScheme)CL_Initializer.imPars;
		
		// construct simulation params object
        SimulationParams simParams = new SimulationParams(
                SimulationParams.sim_N_walkers,
                SimulationParams.sim_tmax, SimulationParams.sim_p,
                SimulationParams.sim_initial,
                SimulationParams.sim_geomType,
                SimulationParams.sim_stepType,
                SimulationParams.sim_voxelSize,
                scheme);
		
		// get step params array
        double[] stepParams = StepGeneratorFactory.getStepParamsArray(SimulationParams.sim_stepType, simParams);
        
        // add step params to simparams
        simParams.setStepParams(stepParams);
        
        // construct substrate
        Substrate substrate= new SquashyInflammationSubstrate(simParams, "allCyls.csv");
        
        // construct simulation
		DiffusionSimulation diffSim= new DiffusionSimulation(simParams, scheme, substrate);
		
		// get the ste generator out of the siuation object
		StepGenerator stepGen= diffSim.getStepGenerator();
		
		// get the scan object out of the simulation object
		SyntheticScan scan= diffSim.getScan(); 

		// set walker position (from debug output
		double[] r0= new double[]{2.5686876843502166E-4,2.1144633425261604E-4,2.4779049425419504E-4};
		
		// construct a test walker
		Walker walker= new Walker(r0, stepGen, substrate, scan, null);

		// construct the problematic step (from debug output)
		double[] step= new double[]{-5.899642786304271E-6,-2.129506892156036E-5,2.051132015790311E-5};
		
		// have the substrate amend the step
		boolean stepOk=false;
		int count=0;
		
		while(!stepOk){
			stepOk= substrate.amend(walker, step, 0, 0, false, null);
			if(++count>10){
				throw new LoggedException("rejected 10 steps - can't find a step that isn't rejected!");
			}
			
			step=new double[]{-5.899642786304271E-6,-2.129506892156036E-5,2.051132015790311E-5};
		}
		
	}
	
	/**
	 * override the positionOK method from superclass to utilise spatial optimisation
	 * 
	 * @param r0 walker position
	 * @para R walker radius
	 */
	public boolean positionOk(double[] r0, double R){
    	
    	
    	double minDist=Double.MAX_VALUE;
    	
    	int[] cell= getInersectingCells(r0, R);
    	
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
	
	/**
	 * entrypoint for testing
	 */
	public static void main(String[] args){
		
		testLongTimeBug(args);
	}
}

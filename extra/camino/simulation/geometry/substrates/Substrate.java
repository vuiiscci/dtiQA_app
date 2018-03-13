/** 
 * The substrates package contains classes that represent
 * complete substrates. A substrate defines a spatial region,
 * handles boundary conditions, barrier intersention and step
 * amendment. Substrates extend the abstract class Substrate.
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation.geometry.substrates;


import imaging.RectGradSteTanScheme;
import imaging.SimulableScheme;

import java.io.FileWriter;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Logger;

import misc.LoggedException;
import numerics.MTRandom;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepAmender;
import simulation.dynamics.StepAmenderFactory;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.StepGenerator;
import simulation.dynamics.FixedLengthStepGenerator;
import simulation.dynamics.Walker;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.geometry.elements.SubstrateObject;
import simulation.geometry.elements.Triangle;
import simulation.geometry.substrates.SubstrateFactory.SubstrateType;
import tools.CL_Initializer;


/**
 *  Abstract root class of all substrates.
 * 
 * 
 * 
 *
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public abstract class Substrate{

    /**
     * logging object
     */
    Logger logger = Logger.getLogger(this.getClass().getName());
    
    /**
     * the dimensionality of the substrate
     */
    private final int D=DiffusionSimulation.D;
    
    /**
     * membrane transition probability
     */
    private final double[] p= new double[1];
        
    /**
     * substrate dimensions
     */
    protected final double[] L;
    
    /**
     * space to store normal
     */
    private final double[] normal;
        
    /**
     * space to store amended step
     */
    private final double[] amended;
    
    /**
     * space to store unamended step
     */
    private final double[] unamended;

    /**
     * space to store step to barrier
     */
    private final double[] toBarrier;
    
    /**
     * point of intersection
     */
    private final double[] intPoint;
    
    /**
     * space for amended step
     */
    private final double[] newStep;

    /**
     * space to store membrane distance
     */
    private final double[] d;

    /**
     * space to store intracellular flag
     */
    private final boolean[] in;
    
    /**
     * space to store position of walker on substrate
     */
    protected final double[] subsCoords= new double[D];
    
    /**
     * space to store barrier transition steps
     */
    private final double[] transition= new double[D];
    
    /**
     * rechecking offset
     */
    private final double[] offset= new double[D];
    
    /**
     * random number generator
     */
    private final MTRandom twister;
    
    /**
     * step amender
     */
    private final StepAmender amender;

    /** 
     * fraction of substrate size to distance of test boundary from actual 
     */
    public final double a= 1000000000;
    
    /** 
     * boundary bottom-left 
     */
    private final double[] xmin= new double[D];
    
    /** 
     * boundary top-right 
     */
    private final double[] xmax= new double[D];
    
    /**
     * boundary widths in each direction
     */
    public double[] border= new double[D];
    
    /**
     * flag to read-out cylinder info
     */
    protected boolean substrateInfo= SimulationParams.substrateInfo;
    
    /**
     * list of objects on the substrate
     */
    protected SubstrateObject[] subsObj;
    
    /**
     * number of subvoxels in each direction for spatial optimisation
     */
    protected final int[] n=new int[D];
    
    /**
     * flag indicating if spatial opt grid has been set or not
     */
    public boolean spatialOptInitialised= false;
    
    /**
     * dimensions of subvoxels in spatial optimisation grid
     */
    protected final double[] s= new double[D];

    /** initial pos */
    private final double[] px= new double[D];
    
    /** final pos */
    private final double[] qx= new double[D];
    
    /** gradient vector */
    private final double[] du= new double[D];
    
    /** increment test vector */
    private final double[] e= new double[D];
    
    /** increment factor vector */
    private final double[] r= new double[D];
        
    /** initial subvoxel */
    private final int[] cp= new int[D];
    
    /** unreflected subvoxel coords */
    private final int[] ctrue= new int[D];
    
    /** final subvoxel */
    private final int[] cq= new int[D];
    
    /** flags indicating reflection symmetry true= flipped */
    private final boolean[] flipped= new boolean[D];  
    

    
    
    /**
     * map of subvoxels to arrays of substrate objects for 
     * spatial optimisation
     */
    public SubstrateObject[][] voxToObjects=null;
    
    
    /**
     * space to store latest list of candidates for intersection testing
     */
    protected int[] candidateSubVox;
    
    /**
     * if we're iterating through a list of candidates, this is the line we're on
     */
    private int currentSubVox=-1;
    
    /**
     * number of candidates in current list
     */
    private int subVoxListLength=0;
    
    /**
     * object index in current voxel
     */
    private int objIndex=-1;

    
    /**
     * to avoid repeated checks of objects 
     */
    protected SubstrateObject[] checked;
    
    /**
     * counter for length of list of checked objects
     */
    private int checkedLength=-1;
    
    /**
     * flag to say whether an intersection is with the cell boundary or not
     */
    public boolean intersectsBoundary=false;
    
  
    /**
     * space to store a subvox index array
     */
    private final int[] subvox= new int[D];
    
    /** 
     * constructor with p specified
     * 
     * @param simParams simulation parameters
     * @param substrateDims substrate dimensions
     */
    public Substrate(SimulationParams simParams, double[] substrateDims){

        this.L=substrateDims;
        
        this.normal=new double[D];
        this.amended= new double[D];
        this.unamended= new double[D];
        this.d= new double[1];
        this.in= new boolean[1];
        this.toBarrier= new double[D];
        this.intPoint= new double[D];
        this.newStep= new double[D];
        this.p[0]=simParams.getP(); 
        
        this.twister=new MTRandom(CL_Initializer.seed);
        
        this.amender=StepAmenderFactory.getStepAmender(SimulationParams.sim_amender_type, this);
        
        // set intersection check coords
        initBoundaryIntersectionArrays();
    }
        
    /** 
     * constructor for testing. only dimensions specified
     * 
     * @param simParams simulation parameters
     * @param substrateDims substrate dimensions
     */
    public Substrate(double[] substrateDims){

    	logger.warning("Substrate constructed without SimulationParams obj. No step amender will be constructed.");
    	
        this.L=substrateDims;
        
        this.normal=new double[D];
        this.amended= new double[D];
        this.unamended= new double[D];
        this.d= new double[1];
        this.in= new boolean[1];
        this.toBarrier= new double[D];
        this.intPoint= new double[D];
        this.newStep= new double[D];
        
        this.twister= new MTRandom(CL_Initializer.seed);
        
        //this.amender=StepAmenderFactory.getStepAmender(SimulationParams.sim_amender_type, this);
        this.amender= null;
        
        // set intersection check coords
        initBoundaryIntersectionArrays();
    }
    /** 
     * sets the values of the min and max values to use when
     * checking if a spin walks off the substrate. these are slightly 
     * wider than the substrate boundaries themselves, which is so that
     * numerical nonsense with mapping into the substrate.
     * 
     */
    public void initBoundaryIntersectionArrays(){
        for(int i=0; i<D; i++){
            xmin[i]=-L[i]/a;
            xmax[i]= L[i]*(1+1/a);
            
            border[i]=L[i]/a;
        }

    }
    
    
    public final boolean amend(Walker walker, double[] step, double t, int n, boolean report, FileWriter debugWriter){
   	
        boolean crosses=false;
    	
        boolean isAmended= false;
        
        boolean stepMade=false;
        
        double origLength=0.0;
    	
        FileWriter cylWriter=null;
        
        if(report){
        	try{
        		cylWriter= new FileWriter("cylsTested.csv");
        	}
        	catch(IOException ioe){
        		throw new LoggedException(ioe);
        	}
        }
        
        // reset the new step
        for(int i=0; i<D; i++){
        	newStep[i]=0.0;
        }
        
        // check if step will cross a barrier
        try{
        	crosses=crossesMembrane(walker, newStep, step, normal, d, stepMade, origLength, in, p, report, cylWriter);
        }
        catch(StepRejectedException sre){
        	// reject if step is pathological
        	return false;
        }
        
        int count=0;
    	
        // if it crosses, we have to amend it and check for more crossings
        while(crosses){                    
        	
        	getSubstrateCoords(walker.r, newStep, subsCoords);
        	
        	count++;
        	
        	if(count==1000){
        	    logger.warning("t= "+t+" walker= "+n+" amendment loop count= "+count);
        	}
        	
        	if(count==10000000){
        		throw new LoggedException("substrate amendment has exceeded max iterations allowed t="+t+" i="+n);
        	}
        	
            // amend the step -- this gives the step the barrier 
            // and the amended and unamended steps away from it
            amender.amendStep(walker, subsCoords, newStep, step, normal, d, origLength, toBarrier, amended, unamended, L, t, n);
                        
            // set the amended flag
            isAmended= true;
            
            // adjust step lengths to account for finite walker size
            double R= walker.R;
            
            
            double toBarrierStepLength=0.0;
            double unamendedStepLength=0.0;
            double amendedStepLength=0.0;
            
            for(int i=0; i<D; i++){
            	toBarrierStepLength+=toBarrier[i]*toBarrier[i];
            	unamendedStepLength+=unamended[i]*unamended[i];
            	amendedStepLength+=amended[i]*amended[i];
            }
            
            toBarrierStepLength=Math.sqrt(toBarrierStepLength);
            unamendedStepLength=Math.sqrt(unamendedStepLength);
            amendedStepLength=Math.sqrt(amendedStepLength);

            if(toBarrierStepLength<R){
                return false;
            }
            
            double lengthRatToBarrier= 1 - R/toBarrierStepLength;
            double lengthRatAmended= 1 + R/amendedStepLength;
            double lengthRatUnamnded= 1 - R/unamendedStepLength;
            
            //logger.info("Original to barrier: "+toBarrier[0]+", "+toBarrier[1]+", "+toBarrier[2]);
            
            
            for(int i=0; i<D; i++){
            	
            	// shorten step to barrier
            	toBarrier[i]*=lengthRatToBarrier;
            	
            	// shorten unamended step
            	unamended[i]*=lengthRatUnamnded;
            	
            	// lengthen amended step
            	amended[i]*=lengthRatAmended;
            	
            	// construct barrier transition step
            	//transition[i]=2*R*toBarrier[i]/toBarrierStepLength;
            	
            }
            
            double newToBarrierStepLength=0.0;
            double newAmendedStepLength=0.0;
            
            for(int i=0; i<D; i++){            
		        //calculate new, shorter newToBarrierStepLength;
		        newToBarrierStepLength+=toBarrier[i]*toBarrier[i];
		        newAmendedStepLength+=amended[i]*amended[i];
            }
            
            newToBarrierStepLength=Math.sqrt(newToBarrierStepLength);
            newAmendedStepLength=Math.sqrt(newAmendedStepLength);
            
            if(newAmendedStepLength<R){
            	logger.warning("amended step length is less than walker radius!");	
            }
            
            for(int i=0; i<D; i++){
                //construct barrier transition step
                transition[i]=2*R*toBarrier[i]/newToBarrierStepLength;
                
                //change amended step to take into account the new transition step
                //amended[i]=transition[i]+amended[i];
            }
            
            //logger.info("To barrier: "+toBarrier[0]+", "+toBarrier[1]+", "+toBarrier[2]);
            //logger.info("amended: "+amended[0]+", "+amended[1]+", "+amended[2]);
            //logger.info("unamended: "+unamended[0]+", "+unamended[1]+", "+unamended[2]);
            //logger.info("transition: "+transition[0]+", "+transition[1]+", "+transition[2]);
            
            // always step to the barrier
            //walker.makeStep(toBarrier);                    
            for(int j=0; j<D; j++){
            	newStep[j]+=toBarrier[j];
            }
            
            if(report){
            	double[] pos= new double[D];
            	
            	getSubstrateCoords(walker.r, newStep, pos);
            	
            	try{
            		debugWriter.write(pos[0]+","+pos[1]+","+pos[2]+"\n");
            	}
            	catch(IOException ioe){
            		throw new LoggedException();
            	}
            }
            
        	// pick a number
            double prob=twister.nextDouble();
            //logger.info("Rand number: "+prob);
            
            
            if(prob>=p[0]){
            	// we're reflecting and new step is the 
            	// amended step reflected away from the barrier
            	for(int j=0; j<D; j++){
            		step[j]=amended[j];
            	}
            }
            else{
            	// we're not reflecting, and the new step 
            	// is the unamended step through the barrier
            	// first transition the membrane
            	//walker.makeStep(transition);
            	
            	
            	// set the next step to be queried as the namended step
            	for(int j=0; j<D; j++){
            		newStep[j]+=transition[j];
            		step[j]=unamended[j];
            	}
            }

            
            //logger.info("offset: "+newStep[0]+", "+newStep[1]+", "+newStep[2]);
            //logger.info("step: "+step[0]+", "+step[1]+", "+step[2]);
            
            stepMade= true;
            
            // in either case, must check if the new step also crosses a barrier 
            try{
            	crosses=crossesMembrane(walker, newStep, step, normal, d, stepMade, origLength, in, p, report, cylWriter);
            }
            catch(StepRejectedException sre){
            	// if we're in here, the step is pathological, so reject it and start again.
            	return false;
            }
        }
        
        // if step is amended
        if(isAmended){
        	for(int i=0; i<D; i++){
        		// add the final part of the step to the newStep
        		newStep[i]+=step[i];
        		
        		// replace the step value for returning
        		step[i]=newStep[i];
        	}
        }

        if(report){
        	try{
        		cylWriter.flush();
        		cylWriter.close();
        		debugWriter.write(step[0]+","+step[1]+","+step[2]+"\n");
        	}
        	catch(IOException ioe){
        		throw new LoggedException(ioe);
        	}
        }
        
        
        return true;

    }
    
    
    /**
     * checks intersection with the boundaries of the substrate. The
     * actual barriers that are checked against are a very tiny amount 
     * wider than those of the actual substrate. This is to ensure that 
     * there's no issue with floating point accuracy when the walker is 
     * mapped back into substrate coords the next time we go round the 
     * barrier check loop.
     * 
     * 
     * @param walkerPos the position of the walker in substrate coords
     * @param offset current offset from starting point
     * @param stepVector proposed step vector
     * @param normal space to store boundary normal
     * @param d distance to boundary
     * @param skipCurrent should we skip the current boundary
     * @param origLength original step length
     * @param in initially intracellular?
     * @param p permeability of boundary (always one here)
     * 
     * @return true is boundary crossed
     */
    protected boolean checkBoundaryIntersection(double[] walkerPos, double[] offset, double[] stepVector,
            double[] normal, double[] d, boolean skipCurrent, double origLength, double[] intDist,  
            boolean[] in, double[] p){
        
       
        //getSubstrateCoords(walkerPos, offset, subsCoords);
        
        double[] dist= new double[D];
        int[] crossInd= new int[D];
        boolean[] minMax= new boolean[D]; // true for min, false for max
        
        int crossCount=0;
        
        
        for(int i=0; i<D; i++){
            double startPoint=walkerPos[i];
            double endPoint=walkerPos[i]+stepVector[i];
            if(endPoint<xmin[i]){
                // we're crossing on the lower end
                dist[crossCount]=(xmin[i]-startPoint)/stepVector[i];
                crossInd[crossCount]=i;
                minMax[crossCount]=true;
                
                crossCount++;
            }
            else if(endPoint>xmax[i]){
                // we're crossing at the upper end
                
                dist[crossCount]= (xmax[i]-startPoint)/stepVector[i];
                crossInd[crossCount]= i;
                minMax[crossCount]=false;
                
                crossCount++;
            }
        }
        
        // if no intersection return false
        if(crossCount==0){
            return false;
        }
        
        // otherwise find closest intersection
        int closest=0;
        double minDist=dist[0];
        int minInd=crossInd[0];
        
        for(int i=1; i<crossCount; i++){
            if(dist[i]<minDist){
                minDist=dist[i];
                minInd=crossInd[i];
                closest=i;
            }
        }
        
        // assemble geometric data 
        double stepLen=0.0;
        for(int i=0; i<D; i++){
            stepLen+=stepVector[i]*stepVector[i];
        }
        stepLen=Math.sqrt(stepLen);
        
        if(minMax[closest]){
            d[0]= -L[minInd]/a;
        }
        else{
            d[0]=L[minInd]*(1+1/a);
        }

        intDist[0]=minDist;
        
        // always totally permeable
        p[0]=1.0;
        
        // normal parallel to appropriate axis, sign not important
        for(int i=0; i<D; i++){
            if(i==minInd){
                normal[i]=1.0;
            }
            else{
                normal[i]=0.0;
            }
        }
        
        return true;
    }
    
    /** 
     * checks if a step will cross a barrier or not
     * 
     * @param walker the walker
     * @param step the step vector
     * @param normal container for the normal
     * @param skipCurrent flag saying whether 
     *                     to ignore the current 
     *                     normal and distance 
     *                     or not
     * 
     * @return the geometry's answer
     * @throws StepRejectedException 
     */
    public abstract boolean crossesMembrane(Walker walker, double[] offset, double[] step, double[] normal, double[] d, boolean skipCurrent, 
    		double origLength, boolean[] in, double[] p, boolean report, FileWriter debugWriter) throws StepRejectedException;
    
    /**
     * @return the size of the substrate
     */
    public abstract double[] getSubstrateSize();
    
    
    /** 
     * @return substrate peak coord
     */
    public abstract double getPeakCoord();

        
    public final int[] getN() {
        return n;
    }

    /**
     *  initialiser for substrate. usually does nothing, but it's
     *  there if we need it for whatever reason.
     *  
     */
    public abstract void init();
    
    
    /**
     * checks if a walker is in intracellular (true) or
     * extracellular (false) space
     * 
     * @param walker the walker to check
     * 
     * @return true if intracellular, false otherwise
     */
    public abstract boolean intracellular(Walker walker);
    
    /**
     * returns the diffusivity at a given location on the substrate
     * 
     * @param location in world coords
     * @return diffusivity at that location
     */
    public double getDiffusivityAt(double[] walkerPos){
    	
    	return CL_Initializer.DIFF_CONST;
    }
    
    
    /**
     * returns the change in magnetisation at a given location and time
     * 
     * @param walker the walker
     * @param t current time (end of timestep)
     * @param tLast the last time we were here (beginning of timestep)
     * 
     * @return zero by default. override to implement relaxation effects on substrate.
     */
    public double getLogMagnetisationChange(Walker walker, double t, double tLast) {
        
        return 0.0;
    }
    
    
    /** 
     * maps physical walker location into substrate cell.
     * this basic version assumes the cell is square.
     * 
     * @param walker the walker whose position we are mapping
     * 
     */
    public void getSubstrateCoords(double[] walkerPos, double[] offset, double[] newPos){
    	
    	for(int i=0; i<D; i++){
    		newPos[i]=walkerPos[i]+offset[i];
    		double windingNum=Math.floor(newPos[i]/L[i]);
    		newPos[i]=newPos[i]-windingNum*L[i];
    	}

    }
    
    public final SubstrateObject[] getSubsObj() {
        return subsObj;
    }

    /**
     * maps the given step into the substrate coordinates.
     * usually, this will not need to do anything at all, since
     * the periodic constraints on position don't apply to step
     * vectors.
     * 
     * this is here to be overridden in any cases where unusual
     * coordinate systems are used, such as those that employ a
     * rotation of different coordinate system to the step generator.
     * 
     * this implementation just returns the original step.
     * 
     * @param walker the walk
     * @param offset offset from current walker position
     * @param step the step
     * 
     * @return the new step
     * 
     */
    public double[] mapStepIntoSubstrate(Walker walker, double[] offset, double[] step){
       
        // GNDN (for all you trekkies...)
        return step;
        
    }
    
    /**
     * unmaps a step that has been transformed using mapStepIntoSubstrate().
     * Usually this does nothing.
     * 
     *  @param subsCoords the position of the walker on a substrate
     *  @param offset offset from current position
     *  @param step the vector to transform
     *  
     *  @return the un-transformed vector
     */
    public double[] unmapStepFromSubstrate(double[] subsCoords, double[] offset, double[] step){
        
        // default case does nothing
        return step;
    }
    
    
    
    /** wrapper for index method that takes array
     * 
     * @param c array of sub vox coords of length D
     * 
     * @return linear index
     */
    private int getSubVoxelIndex(int[] c){
        return getSubVoxelIndex(c[0], c[1], c[2]);
    }
    
    /**
     * calculates linear index from three subvoxel coordinates.
     * 
     * no checking is performed for out of bounds coords.
     * 
     * @param i x-coord
     * @param j y-coord
     * @param k z-coord
     * 
     * @return i+n[0]*j+(n[0]*n[1])*k
     */
    public int getSubVoxelIndex(int i, int j, int k){
    	
        return i + n[0]*j + (n[0]*n[1])*k;
        
        
    }
    
    /**
     * calculates the spatial coordinates of the lower corner
     * of the subvoxel with given integer indices. as this 
     * is only used during initialisation, anew array is 
     * instantiated each time.
     * 
     * @param i x index of subvox
     * @param j y index of subvox
     * @param k z index of subvox
     * 
     * @return spatial coord of lowest corner
     */
    public double[] getBottomLeft(int i, int j, int k){
        
        double[] bottomLeft= new double[D];
        
        bottomLeft[0]= i*s[0];
        bottomLeft[1]= j*s[1];
        bottomLeft[2]= k*s[2];
        
        return bottomLeft;
    }
    
    /**
     * calculates the spatial coordinates of the upper right
     * corner of the subvoxel given. as this is used is only 
     * used during initialisation, a new array is instantiated 
     * each call.
     * 
     * @param i x index of subvox
     * @param j y index of subvox
     * @param k z index of subvox
     * 
     * @return spatial coord of upper corner
     */
    public double[] getUpperRight(int i, int j, int k){
        
        return getBottomLeft(i+1, j+1, k+1);
    }
    
    protected int initSpatialOptArrays(int[] n){
    	
    	// initalise subvoxel counter
    	int numSubVoxels= 1;
        
        for(int i=0; i<D; i++){
            if(n[i]==0){
                throw new LoggedException("number of spatial optimisation voxels in direction "+i+" is zero");
            }
            
            // count subvoxels
            numSubVoxels*=n[i];
            
            // initialise subvoxel dimension
            s[i]=L[i]/n[i];
            
            // set global voxel number array
            this.n[i]=n[i];
        }
        
        voxToObjects= new SubstrateObject[numSubVoxels][];
        candidateSubVox= new int[voxToObjects.length];
        
        checked= new SubstrateObject[subsObj.length];
        
        return numSubVoxels;
    }
    
    
    /**
     * initialises the map of substrate objects to subvoxels
     * in the spatial optimisation grid.
     * 
     * This is a general purpose version that checks every object 
     * against every box.
     * 
     * @param n number of subvoxels in each direction
     */
    public void initialiseSpatialOptimisation(int[] n){
        
        if(spatialOptInitialised){
            logger.warning("spatial optimisation grid being initialised with "+n[0]+","+n[1]+","+n[2]);
            logger.warning("currently n= "+n[0]+","+n[1]+","+n[2]);
        }
        
        spatialOptInitialised=true;
       
        logger.info("initialising spatial optimisation...");
        
        initSpatialOptArrays(n);
        
        for(int i=0; i<n[0]; i++){
            for(int j=0; j<n[1]; j++){
                for(int k=0; k<n[2]; k++){
                    int index= getSubVoxelIndex(i, j, k);
                    
                    double[] bottomLeft= getBottomLeft(i, j, k);
                    double[] topRight= getUpperRight(i, j, k);
                    
                    ArrayList<SubstrateObject> intersectingObjects= new ArrayList<SubstrateObject>();
                    
                    for(int o=0; o<subsObj.length; o++){
                        SubstrateObject so = subsObj[o];
                        if(so.intersectsCubicRegion(bottomLeft, topRight)){
                            intersectingObjects.add(so);                            
                        }
                    }
                    
                    if(intersectingObjects.size()>0){
                        voxToObjects[index]= new SubstrateObject[intersectingObjects.size()];
                        for(int o=0; o<intersectingObjects.size(); o++){
                            voxToObjects[index][o]= intersectingObjects.get(o);
                        }
                    }
                    else{
                        voxToObjects[index]=null;
                    }           
                }
            }
        }
        logger.info("spatial optimisation initialised.");
    }
    
    
    /** 
     * returns the subvoxel for a given location on the substrate.
     * 
     * @param subscoords location to check
     * @param subvox space to store return values
     */
    protected final void getSubVoxel(double[] subscoords, int[] subvox){
    	
    	for(int i=0; i<D; i++){
    		subvox[i]=(int)Math.floor(subsCoords[i]/s[i]);
    	}
    	
    }
    
    
    /**
     * having found the index of a subvoxel to add to the 
     * candidates list, this method does the adding. New
     * entries are checked to make sure that the particular
     * subvoxel contains objects. if no objects, they aren't 
     * added to the list.
     * 
     * @param index the linear cell index to add
     */
    private final void report(int[] c){
        
        // ignore cells that are off the edges of the substrate
        for(int i=0; i<c.length; i++){
            if(c[i]>=n[i]){
                return;
            }
            if(c[i]<0){
                return;
            }
        }
        
        int index= getSubVoxelIndex(c);
        
        if(voxToObjects[index]!=null){
            candidateSubVox[subVoxListLength]=index;
            subVoxListLength++;
        }
        
    }
    
    
    /**
     * checks if a walker is in the current voxel. This is designed
     * to be used by the synthetic scan so that walkers outside of
     * a central voxel can be eliminated from from data generation.
     * 
     * This is the default version and just returns true all the time.
     * Override to get different behaviour.
     * 
     * @param pos the position to check
     * 
     * @return always true here
     */
    public boolean voxelContains(double[] pos){
        
        return true;
    }
    
    /**
     * calculates a list of indicies of subvoxels that are
     * intersected by a given step made by a walker. This
     * uses an algorithm from "Computer graphics and virtual
     * environments" by Slater, Stted and Chrysathou (pg 370)
     * which is based on Bresenham's 2D line-drawing algorithm.
     * 
     * Since the version in the book assumes that all 
     * components of the gradient are positive and non-zero, so 
     * this implementation reflects the coordinates of the 
     * step through the plane perpendicular to any negative 
     * comps of the step gradient and then re-reflects the 
     * subvoxel indices accordingly in the results.
     * 
     * 
     * 
     * @param subsCoords coordinates of the walker in substrate
     * @param offset net previous steps taken so far
     * @param step step vector to test
     * 
     * @return list of int indices to intersected subvoxels 
     */
    protected final void assembleSubVoxelList(double[] subsCoords, double[] step){
        
        subVoxListLength=0;
        checkedLength=0;
        
        // initialise vectors
        for(int j=0; j<D; j++){
            
            // initialise starting pos and gradients, reflecting if negative
            if(step[j]<0){
                flipped[j]=true;
                px[j]=L[j]-subsCoords[j];
                qx[j]= px[j]-step[j];
                du[j]=-step[j];
            }
            else{
                flipped[j]=false;
                px[j]=subsCoords[j];
                qx[j]= px[j]+step[j];
                du[j]=step[j];
            }

            // initialise starting cell
            cp[j]= (int)(px[j]/s[j]);
            
            // initialise end cell
            cq[j]= (int)(qx[j]/s[j]);
            
            // initialise unflipped cell
            ctrue[j]=(int)(subsCoords[j]/s[j]);
            
            if(flipped[j]){
            	/* if we're flipped, it's possible that rounding error can introduce
            	 * a discrepancy between the flipped and unflipped box coordinates.
            	 * In this case we trust the unflipped coordinate, ctrue. 
            	 */
            	if(cp[j]!=(n[j]-1)-ctrue[j]){
            		cp[j]=(n[j]-1)-ctrue[j];
            	}
        	}
            
            // initialise update test vector
            e[j]= ((cp[j]+1)*s[j]-px[j])/du[j];
            
            // need to take care of the -Infinities that can occur if du[j]=-0.0
            if(Double.isInfinite(e[j])||Double.isNaN(e[j])){
            	e[j]=Double.MAX_VALUE;
            }
            
            r[j]=s[j]/du[j];
            // need to take care of the -Infinities that can occur if du[j]=-0.0
            if(Double.isInfinite(r[j])||Double.isNaN(r[j])){
            	e[j]=Double.MAX_VALUE;
            }
            
        }
        
        // report ctrue;
        report(ctrue);
        
        // main loop
        while(!((cp[0]==cq[0])&&(cp[1]==cq[1])&&(cp[2]==cq[2]))){

            // find smallest e[j]
            double smallest =e[0];
            int smallestJ =0;
            for(int j=1; j<D; j++){
                if(e[j]<smallest){
                    smallest=e[j];
                    smallestJ= j;
                }
            }
            
            // update subvoxel coord
            cp[smallestJ]++;
            
            // the true (unflipped) coord
            if(flipped[smallestJ]){
                // if flipped, we're going downwards
                ctrue[smallestJ]--;
            }
            else{
                // otherwise we're just doing the same as cp
                ctrue[smallestJ]++;
            }

            e[smallestJ]+=r[smallestJ];

            
            // move
            report(ctrue);
        }
        
        // set the list length
        currentSubVox=0;
        objIndex=0;
        
    }
    
    
    public void initCandidates(Walker walker, double[] offset, double[] step){
        
        // map position to substrate
        getSubstrateCoords(walker.r, offset, subsCoords);
        
        // get intersecting subvoxels
        assembleSubVoxelList(subsCoords, step);

    }
    
    /**
     * are there more candidates to check?
     * 
     * @return
     */
    public final boolean moreCandidates(){
        
        
        if(currentSubVox<subVoxListLength){
            int index= candidateSubVox[currentSubVox];
            if(objIndex<voxToObjects[index].length){
                return true;
            }
        }
        
        return false;
    }
    
    /** 
     * get the next substrate object in the candidate list
     * 
     * @return next object
     */
    public SubstrateObject nextCandidate(){
        
        // get line to read along
        int index= candidateSubVox[currentSubVox];
        
        // space to store return value
        SubstrateObject lastObj= voxToObjects[index][objIndex];
        SubstrateObject nextObj=null;
        boolean alreadyChecked=true;
        
        while(alreadyChecked){
            
            
            nextObj=voxToObjects[index][objIndex++];
        
            
            
            if(objIndex>=voxToObjects[index].length){
                currentSubVox++;
                objIndex=0;
            }
            
            // check if we've already done this one
            alreadyChecked=false;
            for(int i=0; i<checkedLength; i++){
                if(checked[i]==nextObj){
                	nextObj=null;
                    break;
                }
            }
            
            // this is here to catch the case where the last n objects 
            // in the list have all been checked before. 
            if(currentSubVox==subVoxListLength){
                //nextObj=lastObj;
                break;
            }

            index= candidateSubVox[currentSubVox];
            
            // check if we've already done this one
            alreadyChecked=false;
            for(int i=0; i<checkedLength; i++){
                if(checked[i]==nextObj){
                    alreadyChecked=true;
                }
            }
        }
        
        if(nextObj!=null){
        	checked[checkedLength++]=nextObj;
        }
        
        return nextObj;
    }
    
    
    /**
     * test if a proposed location for an initial walker position is OK or not
     * 
     * @param r0 the postion to check
     * @param R the radius of a walker
     * 
     * @return true if the position is further from a barrier than R
     */
    public boolean positionOk(double[] r0, double R){
    	
    	getSubVoxel(r0, subvox);
    	
    	int imin=Math.max(subvox[0]-1,0);
    	int imax=subvox[0]+1;
    	
    	int jmin=Math.max(subvox[1]-1,0);
    	int jmax=subvox[1]+1;
    	
    	int kmin=Math.max(subvox[2]-1,0);
    	int kmax=subvox[2]+1;
    	
    	double minDist=Double.MAX_VALUE;
    	for(int i=imin; i<=imax; i++){
    		for(int j=jmin; j<=jmax; j++){
    			for(int k=kmin; k<=kmax; k++){
    				
    				int index= getSubVoxelIndex(i, j, k);
    				
    				if(voxToObjects==null){
    					continue;
    				}
    				
    				SubstrateObject[] subsObj= voxToObjects[index];
    				
    				if(subsObj==null){
    					continue;
    				}
    				
    				for(int o=0; o<subsObj.length; o++){
    					double distToObj= subsObj[o].getDistanceFrom(r0);
    					
    					if(distToObj<minDist){
    						minDist=distToObj;
    					}
    				}
    			}
    		}
    	}
    	
    	return (minDist>R);
    	
    }
    
    
    /**
     * test access to the amender. this is used in various test cases
     * across the various substrate classes
     * 
     * @param walker walker making the step
     * @param step step being made
     * @param normal barrier normal
     * @param d distance to barrier
     * @param origLength original length opf step
     * @param toBarrier step to the barrier
     * @param amended amended portion of step
     * @param unamended unamended portion of step
     */
    public void testAmendment(Walker walker, double[] step, 
    		double[] normal, double[] d, double origLength, 
    		double[] toBarrier, double[] amended, double[] unamended){
    	
    	double[] offset= new double[]{0.0, 0.0, 0.0};
    	
    	getSubstrateCoords(walker.r, offset, subsCoords);
    	
    	// pass everything directly to the amender
    	amender.amendStep(walker, subsCoords, offset, step, normal, d, origLength, toBarrier, 
    					amended, unamended, L, 0.0, 1);
    	
    }
    
    
    /** 
     * test of step amending code.
     * 
     * creates a walker and a step with known geometry and tests what comes out
     * to see that the geometry is amending steps in the correct way.
     * 
     * walker at coords (1.5, 1.5), plane membrane at p=1.0 from origin with
     * normal (0.0, 1.0) and a step of length 2.0, with vector (1.0, -1.0).
     * 
     * should mean that the step is amended to (1.0, 0.0) and the final coords
     * of the walker are (2.5, 1.5).
     */
    public static void testStepAmendment() {

        // the next four lines are dummies that allow us to 
        // instantiate a substrate. The details of the Geometry
        // object and the simulation parameters are all ignored
        // in this test code.
        double[] stepParams = new double[] {0.1};

        SimulableScheme scheme;
        try{
            URI uri= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= uri.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }
      
        SimulationParams simParams= new SimulationParams(1, 1000, 0.0, 
                SimulationParams.SPIKE, SubstrateFactory.SubstrateType.CYL_1_INFLAM,
                StepType.FIXEDLENGTH, 1.5, scheme);
 
        Substrate substrate = SubstrateFactory.getSubstrate(SubstrateFactory.SubstrateType.CYL_1_INFLAM, simParams);
        
        
        double[] r0= new double[] {1.5, 1.5, 0.0};      // (initial) walker position
        double[] step= new double[] {1.0, -1.0, 0.0};   // step vector
        
        double[] normal= new double[] {0.0, 1.0, 0.0};  // plane normal
        double[] d=new double[] {1.0};             // distance of plane from origin
        
        Walker walker = new Walker(r0);
        
        // test the new step size
/*        assertTrue(Math.abs(step[0]-1.0)<=1E-12);
        assertTrue(Math.abs(step[1])<=1E-12);
        assertTrue(Math.abs(step[2])<=1E-12);
        
        // test the walker's ability to make the step
        assertTrue(Math.abs(walker.r[0]-2.5)<=1E-12);
        assertTrue(Math.abs(walker.r[1]-1.5)<=1E-12);
        assertTrue(Math.abs(walker.r[2])<1E-12);
*/        
            
    }
    
    
    protected static final void testBoundaryIntersection(){
        
        SimulationParams simParams= new SimulationParams(1, 1000, 0.0, SimulationParams.UNIFORM, SubstrateFactory.SubstrateType.CYL_1_FIXED, StepType.FIXEDLENGTH, 1E-5, 0.1);
        
        SimulationParams.sim_R=4E-6;
        SimulationParams.sim_r=1E-6;
        SimulationParams.sim_cyl_pack= ParallelCylinderSubstrate.SQUARE;
        
        CylinderSubstrate cylSubs = new ParallelCylinderSubstrate(simParams);
        
        
        // construct walker and step for lower x boundary
        double[] r0=new double[]{2E-6, 2E-6, 2E-6};
        Walker walker = new Walker(r0);
        
        double[] step= new double[]{-3E-6, 0.0, 0.0};
        
        cylSubs.amend(walker, step, 0.0, 0, false, null);
        
        // construct walker and step for upper x boundary
        r0= new double[]{6E-6, 6E-6, 6E-6};
        walker= new Walker(r0);
        
        step = new double[]{3E-6, 0.0, 0.0};
        
        cylSubs.amend(walker, step, 0.0, 0, false, null);
        
    }
    
    
    private static final void testIntersectionOrdering(){
        
        SimulationParams simParams= new SimulationParams(1, 1000, 0.0, SimulationParams.UNIFORM, SubstrateFactory.SubstrateType.CYL_1_FIXED, StepType.FIXEDLENGTH, 1E-5, 0.1);

        SimulationParams.sim_R=4E-6;
        SimulationParams.sim_r=1E-6;
        SimulationParams.sim_cyl_pack= ParallelCylinderSubstrate.SQUARE;
        

        CylinderSubstrate cylSubs = new ParallelCylinderSubstrate(simParams);
        
        // test multiple intersections. cylinder then boundary
        double[] r0= new double[]{4E-6, 6E-6, 4E-6};
        Walker walker= new Walker(r0);
        
        double[] step= new double[]{0.0, 2.5E-6, 0.0};
        
        cylSubs.amend(walker, step, 0.0, 0, false, null);
        
        
        // boundary then cylinder
        r0= new double[]{0.5E-6, 4E-6, 4E-6};
        walker = new Walker(r0);
        
        step= new double[]{-3E-6, 0.0, 0.0};
        
        cylSubs.amend(walker, step, 0.0, 0, false, null);
        
    }
    
    
    
    private static final void testVeryLongStep(){
        
        SimulationParams simParams= new SimulationParams(1, 1000, 0.0, SimulationParams.UNIFORM, SubstrateFactory.SubstrateType.CYL_1_FIXED, StepType.FIXEDLENGTH, 1E-5, 0.1);

        SimulationParams.sim_R=4E-6; //separation
        SimulationParams.sim_r=1E-6; //radius
        SimulationParams.sim_cyl_pack= ParallelCylinderSubstrate.SQUARE;
        
        double[] stepLength= new double[]{20e-6};
        simParams.setStepParams(stepLength);

        CylinderSubstrate cylSubs = new ParallelCylinderSubstrate(simParams);
        
        FixedLengthStepGenerator stepGen = new FixedLengthStepGenerator(simParams);
        
        double walkerSize=1e-12;
        //double stepLength= 3e-6;
        
        
        // step that wraps the substrate more than once (downwards)
        double[] r0= new double[]{5E-6+(walkerSize*stepLength[0]), 2E-6, 2E-6};
        //double[] r0= new double[]{6E-6, 2E-6, 2E-6};
        Walker walker= new Walker(r0, stepGen);
        
        double[] step= new double[]{0.0, -1*stepLength[0], 0.0};
        //double[] step= new double[]{0.0, -20E-6, 0.0};
        
        cylSubs.amend(walker, step, 0.0, 0, false, null);
        
        
        // now do same upwards
        step= new double[]{0.0, stepLength[0], 0.0};
        //step= new double[]{0.0, 20E-6, 0.0};
        
        cylSubs.amend(walker, step, 0.0, 0, false, null);
        
    }

    
    /**
     * tests regular grid spatial optimisation
     */
    private static void testRegularGridSpatialOpt(){
        
        // construct simulation params
        SimulationParams simParams=new SimulationParams(1, 1000, 0.0, SimulationParams.UNIFORM, SubstrateFactory.SubstrateType.CYL_1_FIXED, StepGeneratorFactory.StepType.FIXEDLENGTH, 5E-5, 0.1);
        
        double r=1E-6;
        double R=3E-6;
        
        SimulationParams.sim_r=r;
        SimulationParams.sim_R=R;
        SimulationParams.sim_cyl_pack= ParallelCylinderSubstrate.SQUARE;
        
        // construct substrate
        Substrate substrate= SubstrateFactory.getSubstrate(SimulationParams.sim_geomType, simParams);
        
        // set number of subVoxels
        int[] n= new int[]{3, 3, 1};
        
        substrate.initialiseSpatialOptimisation(n);
        
        // check contents of spatial optimisation map
        SubstrateObject[][] optMap= substrate.voxToObjects;
        
        System.err.println("optMap has "+optMap.length+" entries:");
        for(int i=0; i<optMap.length; i++){
            System.err.print(i+") ");
            if(optMap[i]!=null){
                System.err.print("["+optMap[i].length+"] ");
                for(int j=0; j<optMap[i].length; j++){
                    System.err.print(optMap[i][j]+" ");
                }
                System.err.println();
            }
            else{
                System.err.println("(empty)");
            }
        }
        
        // now construct position and step
        double[] walkerPos= new double[]{R/2, R/2, 0.0};
        double[] step= new double[]{1.25*R, 0.0, 0.0};
        double[] offset= new double[]{0.0, 0.0, 0.0};
        
        Walker walker= new Walker(walkerPos);
        
        // check candidates list calculated
        substrate.initCandidates(walker, offset, step);
        
        int count=0;
        
        System.err.println("candidates are:");
        while(substrate.moreCandidates()){
            SubstrateObject sobj= substrate.nextCandidate();
            System.err.println(count+") "+ sobj);
            count++;
        }
        
    }
    
    public abstract Collection<Triangle> getTriangles();
    	    	

    
    public static void main(String[] args){
        
        //testBoundaryIntersection();
        
        //testIntersectionOrdering();
        
        testVeryLongStep();
        
        //testRegularGridSpatialOpt();
    }
}

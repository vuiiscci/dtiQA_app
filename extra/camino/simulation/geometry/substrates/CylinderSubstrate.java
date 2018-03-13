package simulation.geometry.substrates;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Logger;

import misc.LoggedException;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.Walker;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.dynamics.exceptions.TooDamnCloseException;
import simulation.geometry.elements.Cylinder;
import simulation.geometry.elements.SquashyCylinder;
import simulation.geometry.elements.SubstrateObject;
import simulation.geometry.elements.Triangle;

public abstract class CylinderSubstrate extends Substrate{

	
	// logging object
	Logger logger = Logger.getLogger(this.getClass().getName());

	/** size of substrate */
	protected double[] L;
	
	// the cylinders 
	public Cylinder cylinder[];
	
	// dimensionality of space
	protected final int D= DiffusionSimulation.D;

	/** which cylinder was crossed */
	protected int[] cylCrossed;
	
	/** object ref to cylinder crossed */
	protected SubstrateObject[] cylCrossedRef;
	
	/** membrane permeability of crossed cylinder */
	protected double[] cyl_p;
	
	/** list of shortest distances to barrier crossings */
	protected double[] shortestDist;
	
	/** list of d[0] values for barrier intersections */
	protected double[] intOriginDist;
	
	/** list of intersection normals for barrier intersections */
	protected double[][] intNormal;
	
	/** were we initially inside or outside the cylinder crossed? */
	protected boolean[] intInside;

	/** was the intersection with a boundary or an object? */
	protected boolean[] isBoundary;
	
	/** which cylinder was the most recent detected intersection with? */
	protected int lastCrossed=-1;	

	/** for spatially optimised cylinders we can't use index to identify 
	 * cylinder, so use an object ref instead
	 */ 
	protected SubstrateObject toSkip=null;
	
	/** use spatial optimisation or not? */
	private final boolean useSpaceOpt;

	
	
	public CylinderSubstrate(double[] L, SimulationParams simParams, boolean useSpaceOpt){
		super(simParams, L);

		this.L=L;
		
		this.useSpaceOpt= useSpaceOpt;
	}
	
	
	

	/** 
	 * checks if a walker's step will take it across a membrane or not.
	 * this just involves checking every cylinder in turn. The cylinder
	 * crosses() method should fill in the blanks where necessary.
	 * 
	 * @param walker the walker to check
	 * @param offset offset from current walker position for beginning of step
	 * @param stepVector the step being made
	 * @param normal space to store the normal if a barrier is crossed
	 * @param d space to store barrier distance from origin dotted with normal
	 * @param skipCurrent flag indicating we're sitting on a barrier and should 
	 *        ignore the closest one.
	 * @param originLength the original length of the step vector;
	 * 
	 * @return true or false
	 */
	public boolean crossesMembrane(Walker walker, double[] offset, double[] stepVector,
			double[] normal, double[] d, boolean skipCurrent, double origLength, 
			boolean[] in, double[] p, boolean report, FileWriter debugWriter) throws StepRejectedException {
		double len= 0.0;
		double[] walkerPos= new double[D];
		double[] intDist= new double[1];
		
		boolean[] intIn= new boolean[1];
		double[] intP= new double[1];
		
		for(int i=0; i<stepVector.length; i++){
			len += stepVector[i]*stepVector[i];
		}
		getSubstrateCoords(walker.r, offset, walkerPos);
		len=Math.sqrt(len);
		
		if(len/origLength<=1E-14){
			return false;
		}
		
		int numIntersections=0;
        boolean skip=skipCurrent;
    
		
		if(useSpaceOpt){
		    
		    initCandidates(walker, offset, stepVector);
	        
	        int count =0;
	        // check interactions with all triangles
	        //while(triIt.hasNext()){
	        while(moreCandidates()){
	            // get the next triangle
	            Cylinder cyl=(Cylinder)nextCandidate();
	            
	            if(cyl==null){
	            	continue;
	            }
	            
	            // check skipping
	            if(skipCurrent){
	                //if(i==lastCrossed){
	                if(cyl==toSkip){
                        intIn[0]=in[0];
                        skip=true;
                    }
                    else{
                        intIn[0]=cyl.inside(walkerPos);
                        skip= false;
                    }	                
	            }

    			boolean crosses=false;
    			try{
    				crosses=cyl.crosses(walkerPos, stepVector, normal, d, skip, origLength, intDist, intIn, intP, walker.R);
    			}
	            catch(TooDamnCloseException tdce){
	            	throw new StepRejectedException(tdce.getMessage());
	            }
	            
	            // check intersection 
	           if(crosses){
                    shortestDist[numIntersections]=intDist[0];
                    intOriginDist[numIntersections]=d[0];
                    for(int j=0; j<D; j++){
                        intNormal[numIntersections][j]=normal[j];
                    }
                    intInside[numIntersections]=intIn[0];
                    cylCrossed[numIntersections]=count;
                    cylCrossedRef[numIntersections]=cyl;
    
                    cyl_p[numIntersections]=intP[0];
                    
                    
                    isBoundary[numIntersections]=false;
                    
                    /*try{
                        BufferedWriter out= new BufferedWriter(new FileWriter("trickyCyl.csv"));
                        
                        ((SquashyCylinder)cyl).drawCrossSection(out);
                        
                        out.flush();
                        out.close();
                    }
                    catch(IOException ioe){
                        throw new LoggedException(ioe);
                    }*/
                    
                    if(report){
                    	cyl.toFile(debugWriter);
                    }
                    
                    numIntersections++;
                    
	            }
	            count++;

	        }
		}
		else{
		    
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
    			
    			
    			boolean crosses=false;
    			try{
    				crosses=cylinder[i].crosses(walkerPos, stepVector, normal, d, skip, origLength, intDist, intIn, intP, walker.R);
    			}
    			catch(TooDamnCloseException tdce){
    				throw new StepRejectedException(tdce.getMessage());
    			}
    			
    			if(crosses){
    
    				shortestDist[numIntersections]=intDist[0];
    				intOriginDist[numIntersections]=d[0];
    				for(int j=0; j<D; j++){
    					intNormal[numIntersections][j]=normal[j];
    				}
    				intInside[numIntersections]=intIn[0];
    				cylCrossed[numIntersections]=i;
    
    				cyl_p[numIntersections]=intP[0];
    				
    				
    				isBoundary[numIntersections]=false;
    				
    				numIntersections++;
    				
    			}
    		}			
		}
		
		// check if we've intersected the boundaries of the substrate
		if(checkBoundaryIntersection(walkerPos, offset, stepVector, normal, d, skipCurrent, origLength, intDist, intIn, intP)){
            shortestDist[numIntersections]=intDist[0];
		    intOriginDist[numIntersections]=d[0];
		    for(int j=0; j<D; j++){
		        intNormal[numIntersections][j]=normal[j];
		    }
		    intInside[numIntersections]=intIn[0];
		    
		    // in this case we haven't crossed a cylinder, so must check them all next time
		    cylCrossed[numIntersections]=-1; 
		    cylCrossedRef[numIntersections]=null;

		    cyl_p[numIntersections]=intP[0];
		    
		    
		    isBoundary[numIntersections]= true;
		    
		    numIntersections++;
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
			
			lastCrossed=cylCrossed[0];
			toSkip=cylCrossedRef[0];
			
			p[0]=cyl_p[0];
			
			
			intersectsBoundary= isBoundary[0];
			
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
					toSkip=cylCrossedRef[i];
				}
			}
			
			// find the second closest
			double secondClosestDist= Double.MAX_VALUE;
			for(int i=0; i<numIntersections; i++){
				if(i==closest){
					continue;
				}
				
				if(shortestDist[i]<secondClosestDist){
					secondClosestDist= shortestDist[i];
				}
			}
			
			// check that there aren't two barrier within walker radius
			if((secondClosestDist-closestDist)*len<=walker.R){
				throw new StepRejectedException("Two barriers within walker radius. rejecting step");
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
			
			
			intersectsBoundary= isBoundary[closest];
			
			return true;
		}
		
		// in this case there are no intersections
		return false;
	}


	public abstract double[] getSubstrateSize();

	public abstract double getPeakCoord();

	public void init(){
		
	}

	/**
	 * change the substrate dimensions
	 * 
	 * @param L the new dimensions
	 */
	public void setSubstrateDims(double[] L){
	    this.L=L;
	}
	
	public boolean intracellular(Walker walker){
		
		double[] substrateCoords= new double[D];
		final double[] noStep= new double[]{0.0, 0.0, 0.0};
		getSubstrateCoords(walker.r, new double[]{0.0,0.0,0.0}, substrateCoords);
		
		//
		if(useSpaceOpt){
			assembleSubVoxelList(substrateCoords,noStep);
			
			while(moreCandidates()){
				Cylinder cyl = (Cylinder) nextCandidate();
				
				if(cyl==null){
					continue;
				}
				
				if(cyl.inside(substrateCoords)){
					return true;
				}
			}
			
			return false;
		}
		else{
			for(int i=0; i<cylinder.length; i++){
				if(cylinder[i].inside(substrateCoords)){
					return true;
				}
			}
			
			return false;
		}
	}

	
	/**
	 * 
	 * 
	 * @return array of cyli
	 * nders from substarte
	 */
	public final Cylinder[] getCylinders(){
		if(cylinder==null){
			throw new LoggedException("attaempt to retreive uninitialised array of cylinders");
		}
		
		return cylinder;
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
    	
    	
    	double minDist=Double.MAX_VALUE;
    	
    				
    	for(int o=0; o<cylinder.length; o++){
    		double distToObj= cylinder[o].getDistanceFrom(r0);
    					
    		if(distToObj<minDist){
    			minDist=distToObj;
    		}
    	}
    	
    	return (minDist>R);
    	
    }
	
	
	final void setCylinders(Cylinder[] cylinder) {
		this.cylinder = cylinder;
		super.subsObj= cylinder;
		
		cylCrossed= new int[cylinder.length+2];
		cylCrossedRef= new SubstrateObject[cylinder.length+2];
		cyl_p= new double[cylinder.length+2];
		shortestDist= new double[cylinder.length+2];
		intOriginDist= new double[cylinder.length+2];
		intNormal= new double[cylinder.length+2][D];
		intInside= new boolean[cylinder.length+2];
		isBoundary= new boolean[cylinder.length+2];
		
		
		
		/*
		 * debug code. outputs cylinder cross sections to csv file.
		 */
		/*Writer cylWriter;
		int N=20;
		
		try{
			cylWriter= new FileWriter(cylinder.length+"_cyls_crossing.csv");
			
			for(int i=0; i<cylinder.length; i++){
				double[] cylP= cylinder[i].getPosition();
				double cylR= cylinder[i].getRadius();
				
				for(int j=0; j<N; j++){
					double t1= (((double)j)/(double)N)*2.0*Math.PI;
					double t2= (((double)(j+1))/(double)N)*2.0*Math.PI;
					
					double x1= cylP[1]+cylR*Math.cos(t1);
					double y1= cylP[2]+cylR*Math.sin(t1);
					
					double x2= cylP[1]+cylR*Math.cos(t2);
					double y2= cylP[2]+cylR*Math.sin(t2);
					
					cylWriter.write(x1+","+y1+"\n"+x2+","+y2+"\n");
				}
			}
			cylWriter.flush();
			cylWriter.close();
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}*/
	}
	
	private static final void testSpatialOptInitialisation(int[] n){
	    
	    double L= 1E-5;
	    
	    SimulationParams simParams= new SimulationParams(1, 1, 0.0, 
	            SimulationParams.UNIFORM, SubstrateFactory.SubstrateType.CYL_1_INFLAM, 
	            StepGeneratorFactory.StepType.FIXEDLENGTH, L, 0.1);
	    
	    SimulationParams.sim_L=L;
	    
	    // construct substrate
	    CylinderSubstrate cylSubs= new SquashyInflammationSubstrate(simParams);

	    // construct cylinders
	    int sqrtN=10;

	    double radius= 9E-7;
	    double spacing= L/sqrtN;
	    
	    
	    Cylinder[] cyl= new Cylinder[sqrtN*sqrtN];
	    
	    for(int i=0; i<sqrtN; i++){
	        for(int j=0; j<sqrtN; j++){
	            int index =i+sqrtN*j;
	            
	            double[] P= new double[]{(i+0.5)*spacing, (j+0.5)*spacing, 0.5*L};
	            
	            cyl[index]= new SquashyCylinder(P, radius, 0.0);
	        }
	    }
	    // set cylinders on substrate
	    cylSubs.setCylinders(cyl);
	    
	    // initialise spatial optimisation
	    cylSubs.initialiseSpatialOptimisation(n);
	    
	    // check sub-voxel mappings
	    
	    
	}
	
	public Collection<Triangle> getTriangles(){
		
		Collection<Triangle> tris= new ArrayList<Triangle>(subsObj.length);
		
		for(int i=0; i<subsObj.length; i++){
			
			tris.addAll(subsObj[i].getTriangles());
			
		}
		
		return tris;
	}
	
	
	public static void main(String[] args){
	    
	    int n[] = new int[]{10, 10, 1};
	    
	    testSpatialOptInitialisation(n);
	    
	}
	
}

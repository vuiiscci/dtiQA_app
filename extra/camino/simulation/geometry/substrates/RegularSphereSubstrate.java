package simulation.geometry.substrates;

import java.io.FileWriter;
import java.util.Collection;
import java.util.logging.Logger;

import misc.LoggedException;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.Walker;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.geometry.elements.Cylinder;
import simulation.geometry.elements.Sphere;
import simulation.geometry.elements.SubstrateObject;
import simulation.geometry.elements.Triangle;

/**
 * Substrate with regularly square-packed spheres of constant radius.
 * 
 * Spheres are primitive cubic-packed (radius r, sep R >= 2r) and non-overlapping.
 * 
 * @author matt hall (matt.hall@ucl.ac.uk)
 *
 */
public class RegularSphereSubstrate extends Substrate {

	/** logging object */
	private final Logger logger = Logger.getLogger(this.getClass().getName());
	
	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;
	
	/** sphere radius */
	private final double r;
	
	/** sphere sep */
	private final double R;
	
	/** substrate size array */
	private final double[] substrateSize;
	
	/** permeability */
	private final double p;
	
	/** the sphere itself */
	private final Sphere sphere;
	
	
	/** storage for 3D coords */
	private final double[] coords= new double[D];
	
	/** null vector */
	private final double[] nullVector= new double[]{0.0,0.0,0.0};
	

	
	/** which cylinder was crossed */
	protected final int[] cylCrossed;
	
	/** object ref to cylinder crossed */
	protected final SubstrateObject[] cylCrossedRef;
	
	/** membrane permeability of crossed cylinder */
	protected final double[] cyl_p;
	
	/** list of shortest distances to barrier crossings */
	protected final double[] shortestDist;
	
	/** list of d[0] values for barrier intersections */
	protected final double[] intOriginDist;
	
	/** list of intersection normals for barrier intersections */
	protected final double[][] intNormal;
	
	/** were we initially inside or outside the cylinder crossed? */
	protected final boolean[] intInside;

	/** was the intersection with a boundary or an object? */
	protected final boolean[] isBoundary;
	
	/** which cylinder was the most recent detected intersection with? */
	protected int lastCrossed=-1;	

	/** for spatially optimised cylinders we can't use index to identify 
	 * cylinder, so use an object ref instead
	 */ 
	protected SubstrateObject toSkip=null;
	

	
	
	/**
	 * Instantiate substrate with simulation parameters object
	 * 
	 * @param simParams
	 */
	public RegularSphereSubstrate(SimulationParams simParams){
		super(simParams, new double[]{SimulationParams.sim_R,SimulationParams.sim_R,SimulationParams.sim_R});
		
		this.r=SimulationParams.sim_r;
		
		this.R=SimulationParams.sim_R;
		
		if(R<2.0*r){
			logger.warning("Sphere separation "+R+" is less than sphere diameter +"
						+(2*r)+" overlapping object may cause unexpected behaviour");
		}
		
		this.p=SimulationParams.sim_p;
		
		double[] P= new double[]{R/2, R/2, R/2};

		substrateSize= new double[]{R, R, R};
		
		this.sphere= new Sphere(r, P, p);

		// that right across the sphere plus a boundary, plus one more for luck
		int maxInts=4;
		
		cylCrossed= new int[maxInts];
		
		cylCrossedRef= new SubstrateObject[maxInts];
		
		cyl_p= new double[maxInts];
		
		shortestDist= new double[maxInts];
		
		intOriginDist= new double[maxInts];
		
		intNormal= new double[maxInts][D];
		
		intInside= new boolean[maxInts];

		isBoundary= new boolean[maxInts];

	}
	

	
	@Override
	public boolean crossesMembrane(Walker walker, double[] offset,
			double[] step, double[] normal, double[] d, boolean skipCurrent,
			double origLength, boolean[] in, double[] p, boolean report, FileWriter debugWriter)
			throws StepRejectedException {
		
		/*double[] walkerPos= new double[D];
		double[] intDist= new double[1];
		
		getSubstrateCoords(walker.r, offset, walkerPos);
		
		return sphere.crosses(walkerPos, step, normal, d, skipCurrent, origLength, intDist, in, p);
		*/
		
		double len= 0.0;
		double[] walkerPos= new double[D];
		double[] intDist= new double[1];
		
		boolean[] intIn= new boolean[1];
		double[] intP= new double[1];
		
		for(int i=0; i<step.length; i++){
			len += step[i]*step[i];
		}
		getSubstrateCoords(walker.r, offset, walkerPos);
		len=Math.sqrt(len);
		
		if(len/origLength<=1E-14){
			return false;
		}
		
		int numIntersections=0;
        boolean skip=skipCurrent;
    
		
		    
		/*if(skipCurrent){
			if(i==lastCrossed){
				intIn[0]=in[0];
				skip=true;
			}
			else{
				intIn[0]=cylinder[i].inside(walkerPos);
				skip= false;
			}				
		}*/
		
		if(sphere.crosses(walkerPos, step, normal, d, skip, origLength, intDist, intIn, intP, walker.R)){

			shortestDist[numIntersections]=intDist[0];
			intOriginDist[numIntersections]=d[0];
			for(int j=0; j<D; j++){
				intNormal[numIntersections][j]=normal[j];
			}
			intInside[numIntersections]=intIn[0];
			cylCrossed[numIntersections]=0;

			cyl_p[numIntersections]=intP[0];
			
			
			isBoundary[numIntersections]=false;
			
			numIntersections++;
			
		}
		
			
		// check if we've intersected the boundaries of the substrate
		if(checkBoundaryIntersection(walkerPos, offset, step, normal, d, skipCurrent, origLength, intDist, intIn, intP)){
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


	public double[] getSubstrateSize() {
		return substrateSize;
	}


	public double getPeakCoord() {
		return R/2;
	}

	
	public void init() {

	}



	
	public boolean intracellular(Walker walker) {
		getSubstrateCoords(walker.r,  nullVector, coords);
		
		return sphere.inside(coords);
	}


	public Collection<Triangle> getTriangles() {
			return sphere.getTriangles();
	}


	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

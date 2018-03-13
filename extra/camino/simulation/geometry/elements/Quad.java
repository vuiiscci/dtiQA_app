package simulation.geometry.elements;

import java.util.logging.Logger;

import simulation.DiffusionSimulation;

/**
 * a quadrilateral defined by two vectors and 
 * parameter min and max in each. Can have an 
 * arbitrary oreientation in 3D space but is
 * an essentiall 2D object.
 * 
 * These are useful as bounding box walls.
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class Quad {
	
	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());

	/** dimensionality of space (assumed to be 3) */
	private final int D= DiffusionSimulation.D;
	
	/** vectors defining the plane of the quad */
	private final double[][] v= new double[2][D];

	/** min params in each vector */
	private final double[] min= new double[2];
	
	/** max params in each vector */
	private final double[] max= new double[2];
	
	/** normal (from cross product) */
	private final double[] n;
	
	/** a point in the quad */
	private final double[] V0= new double[D];
	
	
	/** 
	 * constructor. take two vectors, mins and maxes
	 * 
	 * @param v 2xD array storing two plane vectors
	 * @param min array of 2 min params
	 * @param max array of 2 max params
	 * @param V0 a point in the plane (specifies distance from origin)
	 */
	protected Quad(double[][] v, double[] min, double[] max, double[] V0){
		
		for(int i=0; i<2; i++){
			for(int j=0; j<D; j++){
				this.v[i][j]=v[i][j];
			}
			this.min[i]=min[i];
			this.max[i]=max[i];
		}

		for(int i=0; i<D; i++){
			this.V0[i]=V0[i];
		}
		
		n=Triangle.normCrossProd(v[0], v[1]);
		
		
	}
	
	
	/** 
	 * check plane-step intersection
	 * 
	 * @param walkerPos position of walker (substrate coords)
	 * @param step step vector
	 */
	protected final boolean crossedBy(double[] walkerPos, double[] step){
		
		
    	// start by checking intersection with plane containing triangle
    	double nVmP=0.0;       // normal dot (V0 - walkerpos) 
    	double nStep=0.0;	   // normal dot step
    	
    	
    	// calculate lengyel's 4D dot product of vector & plane
    	// this equals dp of normal and vector plus plane constant
    	for(int i=0; i<D; i++){
    		nVmP+=n[i]*(V0[i]-walkerPos[i]);
    		nStep+=n[i]*step[i];
    	}

    	
    	// if n dot step is zero (within tolerence) there's no plane intersection
    	if(Math.abs(nStep)<=1E-14){
    		return false;
    	}

    	
    	double tInt=nVmP/nStep;

    	// if tInt is out of range, there's no interaction
    	if(tInt<0.0){
    		return false;
    	}
    	    	
    	//if(tInt>=stepLen){
    	if(tInt>=1.0){
    		return false;
    	}
    	
    	// if we cross the plane, project intersection coords into plane and
    	// check ranges
    	double[] projectedPos= new double[D];
		double[] projectedStep= new double[D];
    	
		for(int i=0; i<2; i++){
			double coord=0.0;
			double disp=0.0;
			for(int j=0; j<D; j++){
				coord+=walkerPos[j]*v[i][j];
				disp+=step[j]*v[i][j];
			}
			
			projectedPos[i]=coord;
			projectedStep[i]=disp;
		}
    	/*
		// third coord comes from projection onto plane normal
		projectedPos[2]=0.0;
		
		
		// get intersection coords in box frame 
		for(int j=0; j<D; j++){
			projectedPos[2]+=walkerPos[j]*n[j];
			projectedStep[2]+=step[j]*n[j];
		}
		*/
		// construct intersection coords in plane coords (only need first two coords)
		double[] intCoords= new double[2];
		
		for(int i=0; i<2; i++){
			intCoords[i]=projectedPos[i]+tInt*projectedStep[i];			                                            
		}
		
		
		// check if coordinates are in plane range
		boolean inQuad=true;
		for(int i=0; i<2; i++){
			if((intCoords[i]<min[i])||(intCoords[i]>max[i])){
				inQuad=false;
				break;
			}
		}
		
		
    	return inQuad;
	}
		
	
	private static final void testQuadIntersection(){
		
		double[][] v= new double[][]{{1.0, 0.0, 0.0},
									 {0.0, 1.0, 0.0}};
		double[] vmin= new double[]{0.0, 0.0};
		double[] vmax= new double[]{1.0, 1.0};
		
		double[] V0= new double[]{0.0, 0.0, 0.5};
		
		// construct quad in xy plane
		Quad quad= new Quad(v, vmin, vmax, V0);
		
		// construct step that intersects the quad
		double[] pos= new double[]{0.5, 0.5, 0.4};
		double[] step= new double[]{0.0, 0.0, 0.2};
		
		// check crossing
		boolean crosses=quad.crossedBy(pos, step);
		System.err.println("intersecting step. crosses is "+crosses+" (true)");
		
		// construct non-intersecting step
		step[2]=0.05;
		
		// check crossing again
		crosses=quad.crossedBy(pos, step);
		System.err.println("non-intersecting step. crosses is "+crosses+" (true)");
		
	}
	
	/**
	 * entrypoint. test quad intersection
	 * 
	 * @param args ignored.
	 */
	public static void main(String[] args){
		
		testQuadIntersection();
	}
	
}

package simulation.geometry.elements;

import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Logger;

import numerics.Rotations;

import simulation.DiffusionSimulation;
import simulation.geometry.substrates.CellularLattice;


/**
 * Sphere object for substrates. This is a thin-walled sphere
 * analogous to BasicCylinder, defined by a position and a
 * radius.
 * 
 * @author Matt Hall (matt.hall@ucl.ac.uk)
 *
 */
public class Sphere implements SubstrateObject {
	
	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** numerical tolerance */
	private static final double TINYNUM= CellularLattice.TINYNUM;
	
	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;
	
	/** radius of sphere */
	private final double R;
	
	/** position of sphere centre in 3D space */
	private final double[] P;
	
	/** diffusivity inside the sphere */
	private final double diffInside= DiffusionSimulation.DIFF_CONST;
	
	/** permeability of sphere */
	private final double p;

	/** index of initial position of walker before reflection loop */
	private double[] initialPos= new double[D];
	
	/** records whether the initial walker position is inside or outside the cylinder */
	private boolean initiallyIn;

	/** array of principle axes of the sphere - use identity */
	private final double[][] v;
	
	/** array of min values for bounding box */
	private final double[] vmin;
	
	/** array of max values for bounding box */
	private final double[] vmax;
	
	/** bounding box for sphere */
	private final BoundingBox bbox;
	
	/** number of sub-divisions for triangulation */
	private final int M=30;
	
	
	/** 
	 * constructor. takes radius, position and permeability.
	 *  
	 * @param R sphere radius
	 * @param P sphere position
	 * @param p permeability
	 */
	public Sphere(double R, double[] P, double p){
		
		this.R=R;
		
		this.P= new double[D];
		for(int i=0; i<D; i++){
			this.P[i]=P[i];
		}
		
		this.p= p;
		
		double left= P[0]-R;
		double right= P[0]+R;
		
		double bottom= P[1]-R;
		double top= P[1]-R;
		
		double back= P[2]-R;
		double front= P[2]-R;
		
		// principle axes
		v= new double[][]{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
		
		// vector of minimum values
		vmin= new double[]{left, bottom, back};
		vmax= new double[]{right, top, front};
		
		bbox= new BoundingBox(v, vmin, vmax);
	}
	
	
	/**
	 * Checks is a given step intersects the sphere. If it does, it returns permeability, 
	 * normal and origin distance information for step amendment.
	 * 
	 * @see BasicCylinder.crosses() for more information
	 * 
	 * @param walkerPos position of walker
	 * @param step the step to make
	 * @param normal space to store the normal
	 * @param d space to store distance to origin
	 * @param skipCurrent flag to say if we're skipping the intersection test
	 * @param origLength the length fot eh original step
	 * @param intDist space to store arc-length to intersection
	 * @param flag indicating if we're inside the sphere or not
	 * @param p space to store pereambility
	 * 
	 * @return true if crossing, otherwise false
	 */
	public boolean crosses(double[] walkerPos, double[] step, double[] normal,
			double[] d, boolean skipCurrent, double origLength,
			double[] intDist, boolean[] in, double[] p, double walkerRad) {
		
		
				
		// start and end points of current step projected into plane in axis coords
		double[] w1= new double[]{walkerPos[0]-P[0], walkerPos[1]-P[1], walkerPos[2]-P[2]};

		double[] w2= new double[]{walkerPos[0]-P[0]+step[0], 
								  walkerPos[1]-P[1]+step[1], 
								  walkerPos[2]-P[2]+step[2]};
		
		
		double[] newPos= new double[D];
		
		// distance of start and end points from sphere centre
		double d1sq=(w1[0]*w1[0] + w1[1]*w1[1] + w1[2]*w1[2]);
		double d2sq=(w2[0]*w2[0] + w2[1]*w2[1] + w2[2]*w2[2]);
		
		double Rsq= R*R;
		
		// if skipCurrent is false, this is an inital call
		// before the reflection loop, so we need to update 
		// the initial position of walker and the initially
		// in flag for future skip checks
		if(!skipCurrent){
			initiallyIn= (d1sq<=R*R);
			in[0]=initiallyIn;
		}

				
		// if start and end points and both inside the cylinder, 
		// we're definitely not crossing the circumference
		if((d1sq<=R*R)&&(d2sq<=R*R)){
			return false;
		}
		
		
		// otherwise there potentially is an intersection
		// get the intersection parameters
		double[] root=getIntersectionRoots(walkerPos, step);
		
		// no real roots = no intersections
		if(root==null){
			return false;
		}

		// both roots less than zero = no intersections
		if((root[0]<0.0) && (root[1]<0.0)){
			return false;
		}

		// if first root is out of range but second is in, swap them
		if(!((root[0]>=0.0)&&(root[0]<=1.0))){
			if((root[1]>=0.0)&&(root[1]<=1.0)){
				double temp=root[0];
				
				root[0]=root[1];
				root[1]=temp;				
			}
		}
		
		// if both roots in range, make sure lowest one
		// is first so that correct intersection is used
		if((root[0]>=0.0)&&(root[0]<=1.0)){
			if((root[1]>=0.0)&&(root[1]<=1.0)){
				if(root[0]>root[1]){
					double temp=root[0];
					
					root[0]=root[1];
					root[1]=temp;
				}
			}
		}
				

		// calculate step length
		double stepLen=0.0;
		for(int i=0; i<D; i++){
			stepLen+=step[i]*step[i];
		}
		stepLen=Math.sqrt(stepLen);
		
		for(int i=0; i<D; i++){
			newPos[i]=walkerPos[i]+step[i];
		}
		

		
		// construct intersection point and normal
		for(int i=0; i<root.length; i++){
			if((root[i]>0.0)&&(root[i]<=1.0)){
				// set normal and distance from origin in normal direction
				double[] intPoint= new double[]{walkerPos[0]+root[i]*step[0]-P[0], 
												walkerPos[1]+root[i]*step[1]-P[1], 
												walkerPos[2]+root[i]*step[2]-P[2]}; 

				
				//double r= Math.sqrt(intPoint[0]*intPoint[0] + intPoint[1]*intPoint[1] + intPoint[2]*intPoint[2]);
				//double phi= Math.atan2(intPoint[1], intPoint[0]);
				//double theta= Math.acos(intPoint[2]/R);
				
				double newNormal[] = new double[normal.length];
				
				double newD;
				
				// set normal using radial unit vector
				//newNormal[0]=Math.sin(theta)*Math.cos(phi);
				//newNormal[1]=Math.sin(theta)*Math.sin(phi);
				//newNormal[2]=Math.cos(theta);
				for(int j=0; j<D; j++){
					newNormal[j]= intPoint[j]/R;
				}
				
				
				// check if we need to skip this intersection
				if(skipCurrent){
	
					boolean skipIt=checkSkipping(newPos);
					
					// if we're skipping, jump to the next root in the i loop
					if(skipIt){
						continue;
					}
				}

				// distance from origin in normal dir is rotation-invariant
                newD=0.0;
                for(int j=0; j<D; j++){
                    newD+=(intPoint[j]+P[j])*newNormal[j];
                }

								
				// if we've got here then we need to return 
				// the distance and normal to the amending
				// routines
				d[0]=newD;
				intDist[0]=0.0;
				
								
				for(int j=0; j<normal.length; j++){
					normal[j]=newNormal[j];
				}
				intDist[0]=root[i];
				
				
				// if skipCurrent is false, this is an inital call
				// before the reflection loop, so we need to update 
				// the initial position of walker and the initially
				// in flag for future skip checks
				if(!skipCurrent){
					for(int j=0; j<D; j++){
						initialPos[j]=w1[j];
					}
					if(inside(walkerPos)){
						initiallyIn=true;
					}
					else{
						initiallyIn=false;
					}
				}
				p[0]=this.p;
				return true;
			}
		}
		
		return false;

	}

	/** 
	 * @returns the diffusivity inside the sphere
	 */
	public double getDiffusivityAt(double[] subsCoords) {
		
		return diffInside;
	}

	
	/**
	 * @return bounding box for the sphere, we have a 
	 * free choice about alignment, so it's axis aligned
	 */
	public BoundingBox getBoundingBox() {
		
		return bbox;
	}

	/**
	 * checks if a given step intersects the bounding box (unused)
	 */
	public boolean boundingBoxIntersects(double[] subsCoords, double[] step) {
		
		return bbox.intersectedBy(subsCoords, step);
	}

	
	
	/** 
	 * checks if a position is inside the sphere. Uses distance to 
	 * sphere centre squared less equal radius squared.
	 * 
	 * @param coords on the substrate
	 * 
	 * @return true if inside
	 * 
	 */
	public boolean inside(double[] subsCoords) {
		
		double posSq= ((subsCoords[0]-P[0])*(subsCoords[0]-P[0]) +
					   (subsCoords[1]-P[1])*(subsCoords[1]-P[1]) +
					   (subsCoords[2]-P[2])*(subsCoords[2]-P[2]));
		
		return (posSq<=R*R);
	}

	
	
	/**
	 * checks if the sphere intersects a given cubic region. uses a clamping
	 * procedure to get the distance-squared between the axis-aligned cubic region
	 * and the centre of the sphere and the cubic region. If this is less-equal than
	 * the radius squared then there is an intersection
	 * 
	 * @param bottomLeft lower corner of cube
	 * @param topRight upper corner of cube
	 * 
	 * @return true if intersection
	 * 
	 */
	public boolean intersectsCubicRegion(double[] bottomLeft, double[] topRight) {
		
		double distSq=0.0;
		
		double clampP=0.0;
		
		for(int i=0; i<D; i++){
			/*if(P[i]<bottomLeft[i]){
				clampP=bottomLeft[i];
			}
			else if(P[i]>topRight[i]){
				clampP=topRight[i];
			}
			else{
				clampP=P[i];
			}*/
			
			clampP= BasicCylinder.clamp(P[i], bottomLeft[i], topRight[i]);
			
			distSq+=(clampP-P[i])*(clampP-P[i]);
		}
		
		return distSq<=R*R;
		
	}

	
	/**
	 * returns the distance from the centre of the sphere
	 * 
	 * @param point
	 * 
	 * @return distance (not squared)
	 * 
	 */
	public double getDistanceFrom(double[] r0) {
		double distSq= 0.0;
		
		for(int i=0; i<D; i++){
			distSq+= (r0[i]-P[i])*(r0[i]-P[i]);
		}
		
		return Math.sqrt(distSq);
	}

	
	/**
	 * tesselates the sphere into a set of triangles for visualisation.
	 * Each polar angle is subdivided into a number of intervals given 
	 * by the constant M defined above.
	 * 
	 * @return a Collection of Triangle objects
	 */
	public ArrayList<Triangle> getTriangles() {
		
		ArrayList<Triangle> tris= new ArrayList<Triangle>();
		
		double[] v1;
		double[] v2;
		double[] v3;
		
		// construct the cap at the north pole
		v1= new double[]{0.0, 0.0, P[2]+R};
		
		double theta= Math.PI/M;
		
		for(int i=0; i<M; i++){
			
			double phi= 2.0*i*Math.PI/M;
			double phip1= 2.0*(i+1)*Math.PI/M;
			
			v2= new double[]{P[0]+R*Math.cos(phi)*Math.sin(theta),
					 P[1]+R*Math.sin(phi)*Math.sin(theta),
					 P[2]+R*Math.cos(theta)};
	
			v3= new double[]{P[0]+R*Math.cos(phip1)*Math.sin(theta),
					 P[1]+R*Math.sin(phip1)*Math.sin(theta),
					 P[2]+R*Math.cos(theta)};
	
			
			tris.add(new Triangle(v1, v3, v2, p));
		}
		
		
		// construct main part of sphere
		for(int i=1; i<M-1; i++){
			
			theta=i*Math.PI/M;
			double thetap1= (i+1)*Math.PI/M;
			
			for(int j=0; j<M; j++){
				double phi= 2.0*j*Math.PI/M;
				double phip1= 2.0*(j+1)*Math.PI/M;
				
				// first triangle
				v1= new double[]{ 	P[0]+R*Math.cos(phi)*Math.sin(theta),
									P[1]+R*Math.sin(phi)*Math.sin(theta),
									P[2]+R*Math.cos(theta)};
				v2= new double[]{ 	P[0]+R*Math.cos(phi)*Math.sin(thetap1),
									P[1]+R*Math.sin(phi)*Math.sin(thetap1),
									P[2]+R*Math.cos(thetap1)};
				v3= new double[]{ 	P[0]+R*Math.cos(phip1)*Math.sin(theta),
									P[1]+R*Math.sin(phip1)*Math.sin(theta),
									P[2]+R*Math.cos(theta)};
				
				tris.add(new Triangle(v1, v2, v3, p));

				// second triangle
				v1= new double[]{ 	P[0]+R*Math.cos(phip1)*Math.sin(theta),
									P[1]+R*Math.sin(phip1)*Math.sin(theta),
									P[2]+R*Math.cos(theta)};
				v2= new double[]{ 	P[0]+R*Math.cos(phi)*Math.sin(thetap1),
									P[1]+R*Math.sin(phi)*Math.sin(thetap1),
									P[2]+R*Math.cos(thetap1)};
				v3= new double[]{ 	P[0]+R*Math.cos(phip1)*Math.sin(thetap1),
									P[1]+R*Math.sin(phip1)*Math.sin(thetap1),
									P[2]+R*Math.cos(thetap1)};
	
				tris.add(new Triangle(v1, v2, v3, p));
			}
			
		}

		
		// construct the cap at the south pole
		v1=new double[]{0.0, 0.0, P[2]-R};
		
		theta= (M-1)*Math.PI/M;
		
		for(int i=0; i<M; i++){
			
			double phi= 2.0*i*Math.PI/M;
			double phip1= 2.0*(i+1)*Math.PI/M;
			
			v2= new double[]{P[0]+R*Math.cos(phi)*Math.sin(theta),
					 P[1]+R*Math.sin(phi)*Math.sin(theta),
					 P[2]+R*Math.cos(theta)};
	
			v3= new double[]{P[0]+R*Math.cos(phip1)*Math.sin(theta),
					 P[1]+R*Math.sin(phip1)*Math.sin(theta),
					 P[2]+R*Math.cos(theta)};
	
			
			tris.add(new Triangle(v1, v3, v2, p));
		}

		
		return tris;
	}

	/**
	 * finds the arc-lengths of the points of intersection on the 
	 * step with the sphere
	 * 
	 * @param pos initial position
	 * @param s step vector
	 * 
	 * @return array containing two roots in ascending order
	 */
	private final double[] getIntersectionRoots(double[] pos, double[] step){

		
		double[] x0= pos;
		
		final double[] x1= new double[]{pos[0]+step[0], 
										pos[1]+step[1], 
										pos[2]+step[2]};
		

		double A= (x0[0]-P[0])*(x0[0]-P[0]) + (x0[1]-P[1])*(x0[1]-P[1]) + (x0[2]-P[2])*(x0[2]-P[2]) - R*R;
		double C= dot(step, step);
		double B= A + C - (x1[0]-P[0])*(x1[0]-P[0]) + (x1[1]-P[1])*(x1[1]-P[1]) + (x1[2]-P[2])*(x1[2]-P[2]);
		
		
		
		double det= B*B - 4*A*C;
		
		if(det<0){
			return null;
		}
		
		
		
		double[] t= new double[2];
		
		t[0]= (-B + Math.sqrt(det))/(2*A);
		t[1]= (-B - Math.sqrt(det))/(2*A);
		
		return t;
		
		
		/*final double[] l= new double[]{step[0], step[1], step[2]};
		
		
		//Normalize l into a unit vector
	    double norml = Math.sqrt(l[0]*l[0] + l[1]*l[1] + l[2]*l[2]);
	    l[0] /= norml;
	    l[1] /= norml;
	    l[2] /= norml;

	    //We're going to treat the line as though it went through (0, 0, 0)
	    //So we modify the sphere centre accordingly
	    double[] c = new double[3];
	    c[0] = P[0] - pos[0];
	    c[1] = P[1] - pos[1];
	    c[2] = P[2] - pos[1];

	    //The part under the radical
	    double inside = dot(l, c)*dot(l, c) - dot(c, c) + R*R;

	    //inside < 0 means no solution (no intersection)
	    if(inside < 0)
	        return null;

		double[] t= new double[2];
		

	    //Get solutions
	    t[0] = dot(l, c) + Math.sqrt(inside);
	    t[1] = dot(l, c) - Math.sqrt(inside);
		
	    // normalise by step length so that root are on the interval [0,1]
		for(int i=0; i<2; i++){
			t[i]/=norml;
			
			// round to zero if extremely small
			if(Math.abs(t[i])<TINYNUM){
				t[i]=0.0;
			}
		}
		
		return t;*/

	}

	private final double dot(double[] u, double[] v){
		
		double dp=0.0;
		
		for(int i=0; i<D; i++){
			dp+=(u[i]*v[i]);
		}
		
		return dp;
	}
	
	/**
	 * check if we're skipping a particular intersection (prevents
	 * rounding errors)
	 * 
	 * @param newPos position after step
	 * @return true or false
	 */
	private final boolean checkSkipping(double[] newPos){
		
		if(initiallyIn==inside(newPos)){
			return true;
		}
		else{
			return false;
		}
	}

	
	private static final void testSpheres(){
		
		int D= DiffusionSimulation.D;
		
		double[] P= new double[]{1.0, 1.0, 1.0};
		
		double R= 0.5;
		
		Sphere sphere= new Sphere(R, P, 0.0);
		
		double[] normal= new double[D];
		double[] d= new double[D];
		double[] intDist= new double[2];
		double[] p= new double[1];
		boolean[] in= new boolean[1];
		
		
		double[] posInside= new double[]{0.75, 1.0, 1.0};
		double[] posOutside= new double[]{1.9, 1.0, 1.0};
		
		double[] step1= new double[]{0.5, 0.0, 0.0};
		double[] step2= new double[]{-0.5, 0.0, 0.0};
		
		boolean inside;
		
		inside=sphere.inside(posInside);
		
		System.out.println("pos1 inside is :"+inside);
		
		inside= sphere.inside(posOutside);
		
		System.out.println("pos2 inside is :"+inside);
		
		boolean crosses=sphere.crosses(posInside,  step2, normal, d, false, 0.5, intDist, in, p, 0.005);
		
		System.out.println("pos1, step out, ("+(posInside[0]+step2[0])+","+(posInside[1]+step2[1])+","+(posInside[2]+step2[2])+") crosses= "+crosses);
		
		crosses= sphere.crosses(posInside, step1, normal, posOutside, false, 0.5, intDist, in, p, 0.005);
		
		System.out.println("pos1, step in, crosses= "+crosses);
		
	}
	
	
	public final double[] getPosition(){
		
		return new double[]{P[0], P[1], P[2]};
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		testSpheres();

	}

}

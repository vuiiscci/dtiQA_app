package simulation.geometry.elements;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.exceptions.TooDamnCloseException;
import simulation.geometry.substrates.CellularLattice;
import tools.CL_Initializer;

import misc.LoggedException;
import numerics.Point3D;
import numerics.RealMatrix;
import numerics.Rotations;
import numerics.Vector3D;

/** implements a cylindrical geometry. Cylinder is defined 
 *  like a ray: a vector V defines its axis, another P its position
 *  in space and a scalar radius r provides its thickness.
 *  
 *  This class also provides methods for getting intersections
 *  and normals for reflecting steps.
 * 
 * @author matt
 *
 */
public class BasicCylinder implements Cylinder, SubstrateObject{

	/** numerical tolerence */
	private static final double TINYNUM=CellularLattice.TINYNUM;
	
	/** number of subdivisions around circumferences when constructing triangle set */
	protected static int C=60;
	
	/** position vector */
	public final Point3D P;
	
	/** position vector in rotated cylinder coords */
	private final Point3D Prot;
	
	/** pos as point */
	public final double[] Parr;
	
	   /** pos as point */
    public final double[] ProtArr;

	/** orientation vector */
	private final Vector3D V;
	
	/** orientation vector as array */
	private final double[] Varr;
	
	/** radius */
	private final double r;

	/** 
	 * rotation matrix that rotates the cylinders 
	 * orientation onto the coordinate z-axis
	 */
	private final RealMatrix rotMat;
	
	/** rotation matrix that inverts the above rotation */
	private final RealMatrix invRotMat;
	
	/** permeability of membrane -- moved here from substrate */
	private final double p;
	
	/** 
	 * diffusivity inside.
	 * 
	 * TODO: tie into constructor
	 */
	private final double d=CL_Initializer.DIFF_CONST;
	
	/** dimensionality of space */
	private final int D=DiffusionSimulation.D;

	/** index of initial positon of walker before reflection loop */
	private double[] initialPos= new double[D];
	
	/** records whether the initial walker position is inside or outside the cylinder */
	private boolean initiallyIn;

	/** number of steps when outputting cross section to file */
	private int N=20;
	
	/** bounding box for cylinder */
	private final BoundingBox bbox;

	/** 
	 * constructor without specifying direction. constructs a 
	 * cylinder parallel to the z-axis at the position specified
	 * 
	 * @param P position of cylinder
	 * @param r radius of cylinder
	 * @param p permeability of cylinder
	 */
	public BasicCylinder(double[] P,  double r, double p){
		
		// set cylinder position and point... 
		this.P= new Point3D(P);

		// ... and as array...
		this.Parr= new double[D];
		for(int i=0; i<D; i++){
			Parr[i]=P[i];
		}
		
		// here no rotation of position is neccessary
		Prot=this.P;
		ProtArr= this.Parr;
		
		// set radius
		this.r=r;
		
		// set permeability
		this.p=p;

		// set default orientation parallel to z-axis
		Varr= new double[]{0.0, 0.0, 1.0};
		
		// set orientation
		this.V= new Vector3D(Varr);
		
		
		// construct bounding box
		this.bbox = null;
				
		// rotation matrices not used in this case
		rotMat=null;
		invRotMat=null;
	}

	
	
	/** constructor that specifies orientation.
	 * 
	 * @param P position of cylinder
	 * @param V orientation of cylinder axis
	 * @param r radius of cylinder
	 * @param p permeability of cylinder
	 */
	public BasicCylinder(double[] P, double[] V, double r, double p){
		
		// set cylinder position and point... 
		this.P= new Point3D(P);

		// ... and as array...
		this.Parr= new double[D];
		for(int i=0; i<P.length; i++){
			Parr[i]=P[i];
		}
		
		
		// set radius
		this.r=r;
		
		// set permeability
		this.p=p;
		
		// normalise V vector
		double modV=0.0;
		for(int i=0; i<D; i++){
			modV+=V[i]*V[i];
		}
		
		modV=Math.sqrt(modV);
		
		for(int i=0; i<D; i++){
			V[i]=V[i]/modV;
		}
		
		// set orientation
		this.V= new Vector3D(V);
		this.Varr=V;
		
		// construct bounding box
		//this.bbox = new BoundingBox(v, vmin, vmax);
		this.bbox=null;
		
		// angle between z-axis and cylinder axis
		double theta= Math.acos(V[2]);
		
		
		// set rotations for cylinders that are not parallel to z-axis
		Vector3D zAxis= new Vector3D(0.0, 0.0, 1.0);
		Vector3D axel= this.V.cross(zAxis);
		axel=axel.normalized();
		
		rotMat=Rotations.getRotMat(axel, -theta);
		invRotMat=Rotations.getRotMat(axel, theta);
		
		
		// finally, set the rotated cylinder position
		// ... and rotated into cylinder coords
		Prot=new Point3D(Rotations.transformPoint(invRotMat, Parr));
		ProtArr= new double[]{Prot.x, Prot.y, Prot.z};
	}
	
	/**
	 * checks if step will intersect the cylinder. if so, provides normal 
	 * and distances accordingly.
	 * 
	 * @param rawPos position of walker
	 * @param rawStep step vector
	 * @param normal space to store intersection normal
	 * @param d space to store distance from origin in normal direction
	 * @param origLength the length of the original step
	 * @param intDist space to store arc-length to intersection point
	 * @param in space to store if walker was initially in or out of cylinder
	 * @param p space to store permeability for cylinder if intersection
	 * 
	 * @return true if crossing, with geometry calculated other false and 
	 * geometry arrays untouched
	 */
	public boolean crosses(double[] rawPos, double[] rawStep,
			double[] normal, double[] d, boolean skipCurrent, double origLength, 
			double[] intDist, boolean[] in, double[] p, double walkerRad)
			throws TooDamnCloseException {
		
		// position and radius of the cylinder
		double[] Praw= getPosition();
		double[] P= Praw;
 
		double r= getRadius();
	
		
		double[] walkerPos=rawPos;
		double[] step=rawStep;
		
		if(invRotMat!=null){
			/* 
			 * if cylinder is not parallel to z-axis,
			 * rotate walker, step and position of cylinder
		     * so that they are.
		     * we do all the geometry there and then rotate 
		     * the intersection normal back at the end.
		     */
			walkerPos=Rotations.transformPoint(invRotMat, rawPos);
			step=Rotations.transformPoint(invRotMat, rawStep);
			P=Rotations.transformPoint(invRotMat, Praw);
		    
		}
		
		if(rotMat!=null){
		    double[] finalPos = new double[D];
		    for(int i=0; i<D; i++){
		        finalPos[i]=walkerPos[i]+step[i];
		    }
		    
		    double[] backRotated=Rotations.transformPoint(rotMat, finalPos);
		    
		    boolean inside= inside(backRotated);
		    boolean flumdah=true;
		}
		
		
		// start and end points of current step projected into plane in axis coords
		double[] w1= new double[]{walkerPos[0]-P[0], walkerPos[1]-P[1], 0.0};

		double[] w2= new double[]{walkerPos[0]-P[0]+step[0], 
								  walkerPos[1]-P[1]+step[1], 
								  0.0};
		
		
		double[] newPos= new double[D];
		
		// distance of start and end points from cylinder axis
		double d1=Math.sqrt(w1[0]*w1[0] + w1[1]*w1[1]);
		double d2=Math.sqrt(w2[0]*w2[0] + w2[1]*w2[1]);
		
		
		// if skipCurrent is false, this is an inital call
		// before the reflection loop, so we need to update 
		// the initial position of walker and the initially
		// in flag for future skip checks
		if(!skipCurrent){
			initiallyIn=inside(rawPos);
			in[0]=initiallyIn;
		}

		// check for close approach
		if(Math.abs(d2-r)<walkerRad){
			throw new TooDamnCloseException("step would take walker too close to cylinder surface. rejecting.");
		}
		
		
		// if start and end points and both inside the cylinder, 
		// we're definitely not crossing the circumference
		if((d1<=r)&&(d2<=r)){
			// we're inside the circle
			
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
		for(int i=0; i<2; i++){
			stepLen+=step[i]*step[i];
		}
		stepLen=Math.sqrt(stepLen);
		
		for(int i=0; i<D; i++){
			newPos[i]=rawPos[i]+rawStep[i];
		}
				

		// if we're not in a chord interval, check intersection
		// with circumference.
		for(int i=0; i<root.length; i++){
			if((root[i]>0.0)&&(root[i]<=1.0)){
				// set normal and distance from origin in normal direction
				double[] intPoint= new double[]{walkerPos[0]+root[i]*step[0]-P[0], 
												walkerPos[1]+root[i]*step[1]-P[1], 
												walkerPos[2]+root[i]*step[2]-P[2]}; 

			    
				double theta=Math.atan2(intPoint[1], intPoint[0]);

				double newNormal[] = new double[normal.length];
				
				double newD;
				
				// set normal using radial unit vector
				newNormal[0]=Math.cos(theta);
				newNormal[1]=Math.sin(theta);
				newNormal[2]=0.0;
				
				
				
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
                for(int j=0; j<2; j++){
                    newD+=(intPoint[j]+P[j])*newNormal[j];
                }

				
				
				// in case of arbitrary orientation, need to
                // back-rotate the normal into the voxel frame
                // before returning it
                if(rotMat!=null){
                    double[] tempNormal= Rotations.transformPoint(rotMat, newNormal);
                    newNormal=tempNormal;
                }

				
				

								
				// if we've got here then we need to return 
				// the distance and normal to the amending
				// routines
				d[0]=newD;
				intDist[0]=0.0;
				
								
				for(int j=0; j<normal.length; j++){
					normal[j]=newNormal[j];
					//intDist[0]+=(root[i]*step[i])*(root[i]*step[i]);
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
					if(inside(rawPos)){
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
	 * checks if the current intersection should be ignored.
	 * The safest way to do this is to use the fact that the
	 * cylinder is circular and to check the initially inside
	 * flag, which is true if the walker starts inside the 
	 * cylinder, and compare it to whether the end of the step 
	 * is inside or outside.
	 * 
	 * Given that skipping is always true when this method is
	 * called, and the fact that the cylinder is circular it means
	 * that the current intersection should be ignored if the
	 * flags are equal (start and finish inside or start and 
	 * finish outside) because in that case we're sitting on the
	 * circumference when this method is called.
	 * 
	 * simple routine, very complex reasoning!
	 * 
	 *  @param newPos the new position after the step is made
	 */
	private final boolean checkSkipping(double[] newPos){
		
		if(initiallyIn==inside(newPos)){
			return true;
		}
		else{
			return false;
		}
	}

	
	
	/** 
	 * calculates the distance d of a given point from the central
	 * axis of the cylinder using the formula
	 * 
	 * d = sqrt( (Q-P)^2 - [(Q-P).V]^2)
	 * 
	 * @param Q position vector or point
	 * @return the shortest distance from the point Q to the central 
	 *         axis of the cylinder.
	 */
	public double getDistanceFrom(double[] Q){
		
		// vector Q-P
		Vector3D QminusP= new Vector3D(new Point3D(Q), P);
		
		// length of vector Q-P squared
		double QminusPsquared=QminusP.modSquared();
		
		// projection of Q-P onto orientation vector V
		double QminusPdotV=QminusP.dot(V);
		
		// distance is square root of difference between the two above
		double dist=Math.sqrt(QminusPsquared-QminusPdotV*QminusPdotV);
		
		return dist;
	}
	
	/**
	 * decides whether a point Q is inside the cylinder or not by comparing
	 * the shortest disance of the point to the central axis to the radius of 
	 * the cylinder
	 * 
	 * @param Q the point to test
	 * @return true if dist from axis <= radius, otherwise false
	 */
	public boolean inside(double[] Q){
		
		double dist =getDistanceFrom(Q);
		
		if(dist<=r){
			return true;
		}
		else{
			return false;
		}
	}
	
	/**
	 * solves the quadratic ray-cylinder intersection equation to give
	 * parametric lengths along the step. both roots are returned unless 
	 * they are complex, in which case returns null.
	 * 
	 * TODO: assumes cylinder is parallel to z-axis
	 * 
	 * @param pos walker pos in cell coords
	 * @param s normalised step vector
	 * 
	 * @return both roots, possibly the same. null if complex
	 */
	protected double[] getIntersectionRoots(double[] pos, double[] s){
		
		// space for roots
		double[] t= new double[2];
		
		// position coords in cylinder's coord frame
		double[] Q=new double[D];
		
		//double[] ProtArr= {Prot.x, Prot.y, Prot.z};
		
		// transform cell coords into cylinder's coord frame
		// transate into coord frame
		for(int i=0; i<D; i++){  
			Q[i]=pos[i]-ProtArr[i];
		}
		// rotate so the z-axes coincide
		
		
		// construct coefficients
		double A=(s[0]*s[0] + s[1]*s[1]);
		double B=2.0*(Q[0]*s[0] + Q[1]*s[1]);
		double C=Q[0]*Q[0] + Q[1]*Q[1] - r*r;
		
		// catch complex case 
		if(B*B<4.0*A*C){
			return null;
		}
		else{
			// otherwise construct solutions
			t[0]= (-B + Math.sqrt(B*B - 4*A*C))/(2.0*A);
			t[1]= (-B - Math.sqrt(B*B - 4*A*C))/(2.0*A);

			for(int i=0; i<2; i++){
				if(Math.abs(t[i])<TINYNUM){
					t[i]=0.0;
				}
			}
			
			
			return t;
		}
		
	}
	
	/**
	 * calculates the surface normal at the intersection of
	 * a step and the cylinder.
	 * 
	 * TODO: assumes cylinder is parallel to the z-axis
	 * 
	 * @param pos the walker pos in cell coords
	 * @param t length param along step
	 * @param s normalised step vector
	 * @param normal space for surface normal
	 */
	protected void getIntersectionNormal(double[] pos, double t, double[] s, double[] normal){
		
		double[] Q= new double[D];
		double[] Praw= {P.x, P.y, P.z};
		
		Vector3D sRot=new Vector3D(s);
		
		// rotate vector into new coord from
		sRot=Rotations.rotateVector(sRot, rotMat);
		
		for(int i=0; i<D; i++){
			Q[i] = (pos[i] - Praw[i]) + t*s[i];
		}
				
		for(int i=0; i<D-1; i++){
			normal[i] = Q[i]/r;
		}
		
		normal[D-1]=0.0;
	}
	
	/**
	 * @return diffusivity inside the cylinder
	 */
	public double getDiffusivityAt(double[] subsCoords){
		if(inside(subsCoords)){
			return d;
		}
		
		return CL_Initializer.DIFF_CONST;
	}

	/**
	 * @return oriented bounding box for cylinder
	 */
	public BoundingBox getBoundingBox(){
		return bbox;
	}

	/**
	 * checks if bounding box is intersected by a step
	 * 
	 * @param subsCoords position mapped onto substrate
	 * @param step step vector
	 * 
	 * @return true of false
	 */
	public boolean boundingBoxIntersects(double[] subsCoords, double[] step) {
		
		return bbox.intersectedBy(subsCoords, step);
	}

	
	/**
	 * clamps variable to a given range. used in intersection-checking.
	 * 
	 * @param x variable to clamp
	 * @param A range minimum
	 * @param B range maximum (A<B unchecked)
	 * 
	 * @return a number with minumum A and maximum B 
	 */
	protected static final double clamp(double x, double A, double B){
		
		/* 
		 * Couldn't resist doing it this way. nested question
		 * operators:
		 * 	  is x>A? 
		 * 		if no, return A,
		 * 		if yes, is x<B?
		 * 			if no, return B,
		 * 			else return x
		 * 
		 * which is a speedy clamp with no function calls
		 */
		return (x>A)? ((x<B)? x:B) : A;
	}
	
	/**
	 * checks if cylinder intersects a specified axis-aligned cubic region.
	 * 
	 * THIS ASSUMES THAT CYLINDERS ARE ALIGNED WITH THE LOCAL Z AXIS!!
	 * 
	 * In an axis-aligned bounding box, the closest point in the box to
	 * the cylinder axis is the axis position clamped to the range of the
	 * bounding box.
	 * 
	 * If this closest point is closer to the cylinder position than the
	 * the surface of the cylinder, the cylinder and box intersect.
	 * 
	 * Check uses squares of distances to avoid the square root.
	 * 
	 * Nice, simple, fast. No memory allocation, no arrays.
	 * 
	 * @param bottomLeft coords of bottom-left corner of cubic region
	 * @param topRight coords of top-right corner of cubic region
	 * 
	 * @return true or false
	 */
	public final boolean intersectsCubicRegion(double[] bottomLeft, double[] topRight){
		
		// find closest point in cube to axis of cylinder (in cross-section)
		double cx=clamp(P.x, bottomLeft[0], topRight[0]);
		double cy=clamp(P.y, bottomLeft[1], topRight[1]);
		
		// distance of closest point to cylinder axis
		double dx = cx-P.x;
		double dy = cy-P.y;
		
		// return values is same as distance squared is less equal to the cylinder radius squared
		return (dx*dx+dy*dy) < r*r;
		
	}
	
	/**
	 * returns the position vector for the cylinder.
	 * If the cylinder is not aligned with the z-axis
	 * the position will be in the cylinder's rotated
	 * coordinate frame.
	 * 
	 * @return P vector as array of doubles
	 */
	public double[] getPosition(){
		return Parr;
	}
	
	/** returns the axis vector of the cylinder
	 * 
	 * @return V vector (new array)
	 */
	public double[] getAxis(){
		return new double[]{V.x, V.y, V.z};
	}
	
	/**
	 * returns the radius of the cylinder
	 * 
	 * @return cylinder radius in meters
	 */
	public double getRadius(){
		return r;
	}
	
	/**
	 * returns the set of triangles for this cylinder. 
	 * Number of subdivisions around the circumference are
	 * defined as constant C static field of this class.
	 * 
	 * TODO: assumes triangle parallel to z-axis
	 * 
	 * @return list of triangles for this cylinder
	 */
	public ArrayList<Triangle> getTriangles(){
		
		ArrayList<Triangle> triangles= new ArrayList<Triangle>(2*C);
		
		double L= SimulationParams.sim_L;
		
		double[] c_i= new double[D];
		double[] c_ip1= new double[D];
		
		double[] d_i= new double[D];
		double[] d_ip1= new double[D];
		
		for(int i=0; i<C; i++){
			
			// angles
			double theta_i= 2.0*i*Math.PI/C;
			double theta_ip1= 2.0*(i+1)*Math.PI/C;
			
			// cosines
			double sinTh= Math.sin(theta_i);
			double cosTh= Math.cos(theta_i);
			
			double sinThp1= Math.sin(theta_ip1);
			double cosThp1= Math.cos(theta_ip1);
			
			// set coords of vertices on both circumferences
			c_i[0]= P.x + r*cosTh;
			c_i[1]= P.y + r*sinTh;
			c_i[2]= -L;
			
			c_ip1[0]= P.x + r*cosThp1;
			c_ip1[1]= P.y + r*sinThp1;
			c_ip1[2]= -L;
			
			d_i[0]= P.x + r*cosTh;
			d_i[1]= P.y + r*sinTh;
			d_i[2]= L;

			d_ip1[0]= P.x + r*cosThp1;
			d_ip1[1]= P.y + r*sinThp1;
			d_ip1[2]= L;

			// construct triangles
			Triangle T1= new Triangle(c_i, c_ip1, d_i, p);
			Triangle T2= new Triangle(d_i, c_ip1, d_ip1, p);
			
			triangles.add(T1);
			triangles.add(T2);
		}
		
		return triangles;
	}

	
	
	public void toFile(FileWriter writer) {
		
		try{
			
			writer.write("\n\n");
			
			for(int i=0; i<=N; i++){
				double theta= 2.0*Math.PI*((double)i)/N;
				
				double x= P.x + r*Math.cos(theta);
				double y= P.y + r*Math.sin(theta);
				
				writer.write(x+","+y+"\n");
			}
		}
		catch(IOException ioe){
			throw new LoggedException();
		}
		
	}
	
	
	/**
	 * test rotation int arbitrarily orientated
	 * cylinder.
	 * 
	 */
	private static void testPerpOrientation(){
		
		double[] P= new double[]{0.0, 2E-6, 2E-6};
		double[] V= new double[]{1.0, 0.0, 0.0};
		double r= 1E-6;
		double p= 0.0;
		
		// cylinder parallel to x-axis
		BasicCylinder xCyl= new BasicCylinder(P, V, r, p);
		
		// test position and step. intersect with step in z-dir
		double[] rawPos= new double[]{0.0,
									  2E-6, 
									  5E-7};
		double[] rawStep= new double[]{0.0,
									   0.0,
									   1E-6};
		
		double[] rawFinalPos= new double[rawPos.length];
		for(int i=0; i<rawFinalPos.length; i++){
			rawFinalPos[i]=rawPos[i]+rawStep[i];
		}
		
		//check the start and end points to see if they are inside the cylinder
		System.err.println("start point inside is "+xCyl.inside(rawPos));
		System.err.println("end point inside is "+xCyl.inside(rawFinalPos));
		
		
		// return spaces
		double[] normal= new double[rawPos.length];
		double[] d= new double[1];
		double[] pSpace= new double[1];
		double[] intDist= new double[1];
		boolean[] in= new boolean[1];
		
		// call the crossing code
		boolean crosses= false;
		
		try{
			xCyl.crosses(rawPos, rawStep, normal, d, false, 0.5, intDist, in, pSpace, 1E-8);
			crosses= xCyl.crosses(rawPos, rawStep, normal, d, true, 0.5, intDist, in, pSpace, 1E-8);
		}
		catch(TooDamnCloseException tdce){
			throw new LoggedException(tdce);
		}
		
		// check results
		System.err.println("crosses is "+crosses+" (expected true)");
		System.err.println("ind dist is "+intDist[0]+" (expected 0.1)");
		System.err.println("normal is ("+normal[0]+","+normal[1]+","+normal[2]+") (expected (0.0,0.0,1.0)");
	}
	
	
	public static final void testArbOrientation(){
	    
	    // construct a system parallel to the z axis with a known intersection
	    double[] P= new double[]{2.1E-6, 2.1E-6, 2.1E-6};
	    
	    double r = 1E-6;
	    
	    double[] walkerPos= new double[]{3.2E-6, 2.1E-6, 2.1E-6};
	    double[] step= new double[]{-1E-6, 0.0, 0.0};
        int D= DiffusionSimulation.D;
        
        double[] normal= new double[D];
        double[] d= new double[1];
        double[] intDist= new double[1];
        boolean[] in= new boolean[1];
        double[] p= new double[1];
        
        System.err.println("Unrotated:");
        System.err.println("walkerPos= ("+walkerPos[0]+", "+walkerPos[1]+", "+walkerPos[2]+")");
        System.err.println("step= ("+step[0]+", "+step[1]+", "+step[2]+")");
                
	    
	    BasicCylinder bc= new BasicCylinder(P, r, 0.0);
	    
	    boolean crosses=false;
	    
	    try{
	    	crosses= bc.crosses(walkerPos, step, normal, d, false, 1E-6, intDist, in, p, 1E-8);
	    }
		catch(TooDamnCloseException tdce){
			throw new LoggedException(tdce);
		}

	    System.err.println("unrotated check:");
	    System.err.println("unrotated cylinder crossing= "+crosses);
	    System.err.println("normal= ("+normal[0]+", "+normal[1]+", "+normal[2]+")");
	    System.err.println("d= "+d[0]);
	    System.err.println("intDist= "+intDist[0]);
	    
	    
	    // now define a cylinder with a different orientation
	    double[] V= new double[]{0.707106, 0.0, 0.707106}; 
	    //double[] Prot= Rotations.transformPoint(rotMat, P);
	    
	    
	    BasicCylinder cyl= new BasicCylinder(P, V, r, 0.0);
	    
	    // use the cylinder's rotation matrix to rotate the walker and step 
	    RealMatrix rotMat= cyl.rotMat;
	    RealMatrix invRotMat= cyl.invRotMat;

	    // transform position
	    double[] Prot= Rotations.transformPoint(rotMat, P);
	    
	    cyl= new BasicCylinder(Prot, V, r, 0.0);
	    
	    rotMat= cyl.rotMat;
        invRotMat= cyl.invRotMat;

	    
	    
	       // check rotation
        double[] zAxis= new double[]{0.0, 0.0, 1.0};
        double[] vCheck= Rotations.transformPoint(rotMat, zAxis);
        
        System.err.println("rotation check: ");
        System.err.println("V     = ("+V[0]+", "+V[1]+", "+V[2]+")");
        System.err.println("Vcheck= ("+vCheck[0]+", "+vCheck[1]+", "+vCheck[2]+")");

        double[] axCheck= Rotations.transformPoint(invRotMat, V);
        
        System.err.println("inverse rotation check: ");
        System.err.println("z-axis: ("+zAxis[0]+", "+zAxis[2]+", "+zAxis[2]+")");
        System.err.println("check:  ("+axCheck[0]+", "+axCheck[1]+", "+axCheck[2]+")");
	    
	    
	    walkerPos=Rotations.transformPoint(rotMat, walkerPos);
	    step=Rotations.transformPoint(rotMat, step);

	    System.err.println("rotated:");
	    System.err.println("walkerPos= ("+walkerPos[0]+", "+walkerPos[1]+", "+walkerPos[2]+")");
	    System.err.println("step= ("+step[0]+", "+step[1]+", "+step[2]+")");
	        

	    
	    try{
	    	crosses= cyl.crosses(walkerPos, step, normal, d, false, 1E-6, intDist, in, p, 1E-8);
	    }
		catch(TooDamnCloseException tdce){
			throw new LoggedException(tdce);
		}

	    System.err.println("rotated cylinder crossing");
	    System.err.println("crosses= "+crosses);
	    System.err.println("normal= ("+normal[0]+", "+normal[1]+", "+normal[2]+")");
        System.err.println("d= "+d[0]);
        System.err.println("intDist= "+intDist[0]);
        
	    
	}
	
	
	public static final void testAnnoyingCase(){
        
        // construct a system parallel to the z axis with a known intersection
        double[] P= new double[]{4.2E-6, 2.1E-6, 4.2E-6};
        double[] V= new double[]{-0.7071080798594737, 0.0, 0.7071054825112364}; 

        
        double r = 1E-6;
        
        double[] walkerPos= new double[]{3.079652947112408E-6, 2.4201214556896393E-6, 3.98082367145738E-6};
        double[] step= new double[]{-4.472501558330728E-7, 4.920544958220048E-7, 6.561171170215304E-7};
        
        int D= DiffusionSimulation.D;
        
        double[] normal= new double[D];
        double[] d= new double[1];
        double[] intDist= new double[1];
        boolean[] in= new boolean[1];
        double[] p= new double[1];
        
        System.err.println("Unrotated:");
        System.err.println("walkerPos= ("+walkerPos[0]+", "+walkerPos[1]+", "+walkerPos[2]+")");
        System.err.println("step= ("+step[0]+", "+step[1]+", "+step[2]+")");
                
        
        BasicCylinder bc9= new BasicCylinder(P, V, r, 0.0);
        
        // new construct the other one
        double[] P5= new double[]{2.1E-6, 2.1E-6, 2.1E-6};
        double[] V5= new double[]{-0.7071080798594737, 0.0, 0.7071054825112364}; 

        
        BasicCylinder bc5= new BasicCylinder(P5, V5, r, 0.0);
        
        
        boolean crosses= false;
        
        try{
        	crosses=bc9.crosses(walkerPos, step, normal, d, false, 1E-6, intDist, in, p, 1E-8);
        }
		catch(TooDamnCloseException tdce){
			throw new LoggedException(tdce);
		}


        System.err.println("cylinder 9 crossing= "+crosses);
        System.err.println("normal= ("+normal[0]+", "+normal[1]+", "+normal[2]+")");
        System.err.println("d= "+d[0]);
        System.err.println("intDist= "+intDist[0]);

        
        walkerPos[2]-=4.2E-6;
        
        try{
        	crosses= bc5.crosses(walkerPos, step, normal, d, false, 1E-6, intDist, in, p, 1E-8);
        }
		catch(TooDamnCloseException tdce){
			throw new LoggedException(tdce);
		}

        System.err.println("cylinder 5 crossing= "+crosses);
        System.err.println("normal= ("+normal[0]+", "+normal[1]+", "+normal[2]+")");
        System.err.println("d= "+d[0]);
        System.err.println("intDist= "+intDist[0]);
        
        
        // now define a cylinder with a different orientation
        //double[] Prot= Rotations.transformPoint(rotMat, P);
        
        
        BasicCylinder cyl= new BasicCylinder(P, V, r, 0.0);
        
        // use the cylinder's rotation matrix to rotate the walker and step 
        RealMatrix rotMat= cyl.rotMat;
        RealMatrix invRotMat= cyl.invRotMat;

        // transform position
        double[] Prot= Rotations.transformPoint(rotMat, P);
        
        cyl= new BasicCylinder(Prot, V, r, 0.0);
        
        rotMat= cyl.rotMat;
        invRotMat= cyl.invRotMat;

        
        
           // check rotation
        double[] zAxis= new double[]{0.0, 0.0, 1.0};
        double[] vCheck= Rotations.transformPoint(rotMat, zAxis);
        
        System.err.println("rotation check: ");
        System.err.println("V     = ("+V[0]+", "+V[1]+", "+V[2]+")");
        System.err.println("Vcheck= ("+vCheck[0]+", "+vCheck[1]+", "+vCheck[2]+")");

        double[] axCheck= Rotations.transformPoint(invRotMat, V);
        
        System.err.println("inverse rotation check: ");
        System.err.println("z-axis: ("+zAxis[0]+", "+zAxis[2]+", "+zAxis[2]+")");
        System.err.println("check:  ("+axCheck[0]+", "+axCheck[1]+", "+axCheck[2]+")");
        
        
        walkerPos=Rotations.transformPoint(rotMat, walkerPos);
        step=Rotations.transformPoint(rotMat, step);

        System.err.println("rotated:");
        System.err.println("walkerPos= ("+walkerPos[0]+", "+walkerPos[1]+", "+walkerPos[2]+")");
        System.err.println("step= ("+step[0]+", "+step[1]+", "+step[2]+")");
            

        
        try{
        	crosses= cyl.crosses(walkerPos, step, normal, d, false, 1E-6, intDist, in, p, 1E-8);
        }
		catch(TooDamnCloseException tdce){
			throw new LoggedException(tdce);
		}

        System.err.println("rotated cylinder crossing");
        System.err.println("crosses= "+crosses);
        System.err.println("normal= ("+normal[0]+", "+normal[1]+", "+normal[2]+")");
        System.err.println("d= "+d[0]);
        System.err.println("intDist= "+intDist[0]);
        
        
    }
	
	
	
	/** 
	 * entrypoint. calls test code.
	 * 
	 * @param args
	 */
	public static void main(String[] args){
		
		testAnnoyingCase();
		
	}
	
}

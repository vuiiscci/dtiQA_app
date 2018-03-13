/**
 * 
 */
package simulation.geometry.elements;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;

import simulation.DiffusionSimulation;
import tools.CL_Initializer;

import misc.LoggedException;
import numerics.Point3D;
import numerics.Vector3D;


/**
 * Cylinder implementation with polygonal cross-section.
 * Number of "facets" (i.e. number of sides in the polygonal
 * cross-section) is specified in the constructor. The facets
 * approximate the circular cross section of the smooth
 * cylindrical surface but suffer from fewer of the complications
 * that come with a smooth surface. essentially they behave like
 * the chords of the squashy cylinder arranged in a regular
 * pattern around the circumference of the cross-section.
 * 
 * The more facets, the better the approximation.
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class FacetCylinder implements Cylinder, SubstrateObject{

	private final int D= DiffusionSimulation.D;
	
	private final double d=CL_Initializer.DIFF_CONST;
	
	private final BoundingBox bbox=null;
	
	/** number of facets around cylinder */
	private final int N;
	
	/** Array of Chords for facets */
	private final Line[] facet;
	
	/** Alignment vector */
	private final Vector3D V;
	
	/** position vector */
	protected final Vector3D P;
	
	/** position as point */
	private final Point3D Ppt;
	
	/** position vector as array */
	private final double[] Parr;
	
	/** set of reflected midpoints for inside() function */
	private final double[][] Q;
	
	/** membrane permeability */
	private final double p;
	
	/** radius of cylinder */
	private final double r;

	/** last facet crossed */
	private int lastCrossed=-1;

	/** default orientation */
	private static final Vector3D Vdef= new Vector3D(0.0, 0.0, 1.0);
	
		
	/** constructor.  
	 * @param V alignment
	 * @param P position
	 * @param r radius
	 * @param N number of facets 
	 */
	public FacetCylinder(Vector3D V, Vector3D P, double r, int N, double p){
		
		this.V=V;
		
		this.P=P;
		
		this.Parr=new double[]{P.x, P.y, P.z};
		
		this.Ppt=new Point3D(Parr);
		
		this.r=r;
		
		this.N=N;
		
		this.p=p;
		
		this.facet= new Line[N];

				
		for(int i=0; i<N; i++){
			
			double thetaMin = 2.0*Math.PI*i/(double)N;
			double thetaMax = 2.0*Math.PI*(i+1)/(double)N;
			
			double lowerX= P.x + r*Math.cos(thetaMin);
			double lowerY= P.y + r*Math.sin(thetaMin);
			
			double upperX= P.x + r*Math.cos(thetaMax);
			double upperY= P.y + r*Math.sin(thetaMax);
			
			facet[i]= new Line(lowerX, lowerY, upperX, upperY);
			
		}
		
		
		
		this.Q= new double[N][D-1];
	
		/** the test for whether a point is inside the cylinder
		 * or not requires a set of points that are the cylinder 
		 * centre reflected over each of the facets. This is 
		 * initialised here
		 */
		for(int i=0; i<N; i++){
			// midpoint of facet
			double[] M=facet[i].midpoint();
			
			// vector of centre to midpoint
			double[] PM= new double[]{M[0]-P.x, M[1]-P.y};
			
			// reflected point is centre + twice midpoint vector
			Q[i][0]=P.x+2.0*PM[0];
			Q[i][1]=P.y+2.0*PM[1];
		}
		
	}
	
	/** constructor with default alignment
     * 
     * @param P position
     * @param r radius
     * @param N number of facets
     * @param p permeability
     */
    public FacetCylinder(Vector3D P, double r, int N, double p){
        this(Vdef, P, r, N, p);
    }
    

	
	
	
	
	/**
	 *  checks cylinder crossing by identifying candidate facet
	 *  and checking line-line intersection.
	 *  
	 * @param walkerPos position of walker
	 * @param step step vector
	 * @param normal space to store normal to vector
	 * @param d space to store distance to membrane
	 * @param skipCurrent are we sitting on the membrane right now?
	 * @param origLength original length of step
	 * @param intDist space for arclength to interaction point
	 * @param in inside or outside?
	 * @param space to store membrane permeability
	 * 
	 * @return true if crossing with other bits filled in
	 *                otherwise false with other bits untouched
	 */
	public boolean crosses(double[] walkerPos, double[] step,
			double[] normal, double[] d, boolean skipCurrent, double origLength, double[] intDist, 
			boolean[] in, double[] p, double walkerRad){
		
		
		double[] dists= new double[N];
		double[] intD= new double[N];
		double[][] intNormals= new double[N][D];
		int[] facetNum= new int[N];
		int numInts=0;
		
		
		
		for(int i=0; i<N; i++){
			if(skipCurrent && (i==lastCrossed)){
				continue;
			}
			
			boolean crossing=facet[i].crossedBy(walkerPos, step, normal, d, skipCurrent, origLength, intDist);
			if(crossing){
				//lastCrossed=i;
				//return true;
				
				dists[numInts]=intDist[0];
				facetNum[numInts]=i;
				intD[numInts]=d[0];
				for(int j=0; j<D; j++){
					intNormals[numInts][j]=normal[j];
				}
				numInts++;
			}
		}
		
		// check the number of interactions
		if(numInts==1){
			// if only one, then set the last crossed flag and exit
			lastCrossed=facetNum[0];
			return true;
		}
		
		if(numInts>1){
			// if more than one interaction, pick the closest one
			// then set the appropriate parameters 
			double minIntDist=dists[0];
			int minIntDistIndex=0;
			
			for(int i=1; i<numInts; i++){
				if(dists[i]<minIntDist){
					minIntDist=dists[i];
					minIntDistIndex=i;
				}
			}
			
			lastCrossed=minIntDistIndex;
			d[0]=intD[minIntDistIndex];
			for(int i=0; i<D; i++){
				normal[i]=intNormals[minIntDistIndex][i];
			}

			p[0]=this.p;
			
			return true;
		}
		
		// otherwise there are no interactions, so return false
		return false;
	}
	
	/** calculates the distance d of a given point from the central
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
		Vector3D QminusP= new Vector3D(new Point3D(Q), Ppt);
		
		// length of vector Q-P squared
		double QminusPsquared=QminusP.modSquared();
		
		// projection of Q-P onto orientation vector V
		double QminusPdotV=QminusP.dot(V);
		
		// distance is square root of difference between the two above
		return Math.sqrt(QminusPsquared-QminusPdotV*QminusPdotV);
	}
	
	/**
	 * write out cross section to file
	 * 
	 * @param writer file to writer to
	 */
	
	
	
	
	/** 
	 * checks if a position is inside or outside the 
	 * cylinder. this uses a rather spiffy algorithm i
	 * found on everything2.com which works as follows:
	 * 
	 * 1. calculate a set of points {q} representing the centre 
	 * of the cylinder reflected about each facet of the
	 * cylinder.
	 * 
	 * 2. calculate distance d of point from cylinder centre
	 * 
	 * 3. calculate distance of point from each q
	 * 
	 * 4. if d < q for all q, we're inside the cylinder
	 * 
	 * This routine uses squared-distances for numerical
	 * optimisation.
	 * 
	 * @param coord in cylinder frame
	 * 
	 * @return true or false, in or out
	 */
	public final boolean inside(double[] pos){
		
		double dp= (pos[0]-P.x)*(pos[0]-P.x) + (pos[1]-P.y)*(pos[1]-P.y);
		
		for(int i=0; i<Q.length; i++){
			double dq= (pos[0]-Q[i][0])*(pos[0]-Q[i][0]) + (pos[1]-Q[i][1])*(pos[1]-Q[i][1]);
			
			if(dq<dp){
				return false;
			}
		}
		
		return true;
	}
	
	/**
	 * @return free diffusivity
	 */
	public double getDiffusivityAt(double[] r){
		return d;
	}
	
	/**
	 * @return position of cylinder
	 */
	public double[] getPosition(){
		return Parr;
	}
	
	/**
	 * @return cylinder radius
	 */
	public double getRadius(){
		return r;
	}
	
	
	/**
	 * TODO: not implement. bbox always null.
	 * @return bounding box for cylinder
	 */
	public BoundingBox getBoundingBox(){
		return bbox;
	}
	
	/**
	 * TODO: always returns true
	 * 
	 * checks if a step will intersect the bounding box
	 * 
	 * @param subsCoords walker pos on substrate
	 * @param step step vector
	 * 
	 * @return true or false if box intersected by step
	 */
	public boolean boundingBoxIntersects(double[] subsCoords, double[] step){
		return true;
	}
	
	public final boolean intersectsCubicRegion(double[] bottomLeft, double[] topRight){
		
		return true;
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
		
		ArrayList<Triangle> triangles= new ArrayList<Triangle>(2*N);
		
		double L= 4*r;
		
		double[] c_i= new double[D];
		double[] c_ip1= new double[D];
		
		double[] d_i= new double[D];
		double[] d_ip1= new double[D];
		
		for(int i=0; i<N; i++){
			
			// angles
			double theta_i= 2.0*i*Math.PI/N;
			double theta_ip1= 2.0*(i+1)*Math.PI/N;
			
			// cosines
			double sinTh= Math.sin(theta_i);
			double cosTh= Math.cos(theta_i);
			
			double sinThp1= Math.sin(theta_ip1);
			double cosThp1= Math.cos(theta_ip1);
			
			// set coords of vertices on both circumferences
			c_i[0]= P.x + r*cosTh;
			c_i[1]= P.y + r*sinTh;
			c_i[2]= 0.0;
			
			c_ip1[0]= P.x + r*cosThp1;
			c_ip1[1]= P.y + r*sinThp1;
			c_ip1[2]= 0.0;
			
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

	
	/**
	 * write out cross section to file
	 * 
	 */
	public final void toFile(FileWriter writer){
		
		try{
			for(int i=0; i<=N; i++){
				
				double thetaMin = 2.0*Math.PI*i/(double)N;
				
				double lowerX= P.x + r*Math.cos(thetaMin);
				double lowerY= P.y + r*Math.sin(thetaMin);
				
				writer.write(lowerX+" "+lowerY+"\n");
				
			}
		}
		catch(IOException ioe){
			throw new LoggedException();
		}
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Vector3D V= new Vector3D(new double[]{0.0, 0.0, 1.0});
		Vector3D P= new Vector3D(new double[]{0.5, 0.3, 0.0});
		
		double r=1.0;
		
		int N=10;
		
		FacetCylinder cyl= new FacetCylinder(V, P, r, N, 0.0);
		
		double[] pos= new double[]{0.0, 1.0, 0.0};
		
		boolean inside= cyl.inside(pos);
		System.err.println("pos inside is "+inside);
		
	}

}

package simulation.geometry.elements;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.logging.Logger;

import misc.LoggedException;

import simulation.DiffusionSimulation;
import simulation.dynamics.exceptions.TooDamnCloseException;
import tools.CL_Initializer;

/**
 * Thick-walled cylinder. This two concentric
 * cylinders with a different diffusivities
 * in the central region and in the annular
 * cross sectional region between them
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class MyelinatedCylinder implements Cylinder, SubstrateObject {

	/** logging object */
	private final Logger logger = Logger.getLogger(this.getClass().getName());
	
	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;
	
	/** inner and outer membrane cylinders */
	private final BasicCylinder[] cylinder= new BasicCylinder[2]; 
	
	/** diffusivities */
	private final double[] d;
	
	/** position in space */
	private final double[] P;
	
	/** orientation vector */
	private final double[] V;
	
	/** bounding box for object */
	private final BoundingBox bbox;
	
	/** which cylinder was the most recent detected intersection with? */
	private int lastCrossed=-1;
	
	/** which cylinder was crossed */
	private int[] cylCrossed= new int[cylinder.length];
	
	/** list of shortest distances to barrier crossings */
	private double[] shortestDist= new double[cylinder.length];
	
	/** list of d[0] values for barrier intersections */
	private double[] intOriginDist= new double[cylinder.length];
	
	/** list of intersection normals for barrier intersections */
	private double[][] intNormal= new double[cylinder.length][D];
	
	/** were we initially inside or outside the cylinder crossed? */
	private boolean[] intInside= new boolean[cylinder.length];
	
	/** space to store permeabilities */
	private double[] cyl_p= new double[cylinder.length];
	
	public MyelinatedCylinder(double[] P, double r1, double r2, double d1, double d2, double p){
	    this(new double[]{0.0, 0.0, 1.0}, P, r1, r2, d1, d2, p);
	}
	
	public MyelinatedCylinder(double[] V, double[] P, double r1, double r2, double d1, double d2, double p){
		
		/** set position */
		this.P=P;
		
		/** set orientation */
		this.V=V;
		
		/** set diffusivities */
		this.d=new double[]{d1, d2};
		
		/** inner cylinder */
		cylinder[0]= new BasicCylinder(P, r1, p);

		/** outer cylinder */
		cylinder[1]= new BasicCylinder(P, r2, p);
		
		/** calculate coords for bounding box */
		double[][] v= new double[][]{{1,0,0},{0,1,0}};
		double[] vmin= new double[]{P[0]-r2, P[1]-r2};
		double[] vmax= new double[]{P[0]+r2, P[1]+r2}; 
		
		this.bbox = new BoundingBox(v, vmin, vmax);
	}
	
	
	public boolean crosses(double[] walkerPos, double[] step, double[] normal,
			double[] d, boolean skipCurrent, double origLength,
			double[] intDist, boolean[] in, double[] p, double walkerRad) {

		double len=0.0;		
		boolean[] intIn= new boolean[1];
		double[] intP= new double[1];

		
		for(int i=0; i<step.length; i++){
			len += step[i]*step[i];
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
			
			boolean crosses=false;
			try{
				crosses=cylinder[i].crosses(walkerPos, step, normal, d, skip, origLength, intDist, intIn, intP, walkerRad);
			}
			catch(TooDamnCloseException tdce){
				throw new LoggedException(tdce);
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
			
			lastCrossed=cylCrossed[0];
			
			p[0]= cyl_p[0];
			
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

	}

	/**
	 * @return outer cylinder radius
	 */
	public final double getRadius(){
		return cylinder[1].getRadius();
	}
	
	/**
	 * @return position of inner cylinder (same as outer)
	 */
	public final double[] getPosition(){
		return cylinder[0].getPosition();
	}
	
	/**
	 * @return distance from axis of inner cylinder (same as outer)
	 */
	public final double getDistanceFrom(double[] Q){
		return cylinder[0].getDistanceFrom(Q);
	}
	
	/**
	 * checks if a postion is inside the cylinder.
	 * this is the same as the check against the
	 * outer cylinder, to the check is passed there.
	 * 
	 * @return true or false
	 */
	public final boolean inside(double[] pos){
		return cylinder[1].inside(pos);
	}
	
	/**
	 * returns the oriented bounding box for this object 
	 * 
	 * @return pre-calculated bounding box
	 */
	public BoundingBox getBoundingBox() {
	
		return bbox;
	}

	/**
	 * @return diffusivity at given coordinates
	 */
	public double getDiffusivityAt(double[] subsCoords) {
		if(cylinder[1].inside(subsCoords)){
			if(cylinder[0].inside(subsCoords)){
				// inside both cylinders
				return d[0];
			}
			else{
				// inside annual region between cylinders
				return d[1];
			}
		}
		
		// outside both cylinders -- shouldn't have been called
		return CL_Initializer.DIFF_CONST;
	}

	/**
	 * checks if bounding box is intersected by a step
	 * 
	 * @param subsCoords position mapped onto substrate
	 * @param step step vector
	 * 
	 * @return true of false
	 */
	public final boolean boundingBoxIntersects(double[] subsCoords, double[] step) {
		// TODO Auto-generated method stub
		return bbox.intersectedBy(subsCoords, step);
	}

	
	/** 
	 * checks if cylinder intersects axis-aligned cubic region
	 * 
	 * @param bottomLeft lowest corner of region
	 * @param topRight highest corner of region
	 * 
	 * @return true or false
	 */
	public final boolean intersectsCubicRegion(double[] bottomLeft, double[] topRight){
		
		return cylinder[0].intersectsCubicRegion(bottomLeft, topRight)||
			   cylinder[1].intersectsCubicRegion(bottomLeft, topRight);
	}
	
	/**
	 * spits out cylinder cross section
	 */
	public void toFile(FileWriter writer){
		cylinder[0].toFile(writer);
		cylinder[1].toFile(writer);
	}
	

	/**
	 * return triangles for visualisation/meshisation
	 * 
	 * @return triangles for inner & outer cylinders
	 */
	public ArrayList<Triangle> getTriangles(){
		
		ArrayList<Triangle> triangles= new ArrayList<Triangle>(4*BasicCylinder.C);
		
		ArrayList<Triangle> cylTris= cylinder[0].getTriangles();
		
		triangles.addAll(cylTris);
		
		cylTris= cylinder[1].getTriangles();
		
		triangles.addAll(cylTris);
		
		return triangles;
		
	}

}

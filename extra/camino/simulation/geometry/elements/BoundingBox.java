package simulation.geometry.elements;

import java.util.logging.Logger;

import Jama.Matrix;

import numerics.RealMatrix;
import numerics.Rotations;

import misc.LoggedException;

import simulation.DiffusionSimulation;

/**
 * An optimally oriented bounding box that
 * minimally surrounds a Substrate object
 * 
 * a bounding box is defined by three vectors
 * defining a coordinate system centred on
 * the substrate origin, and max and min values
 * of each of them  for points on the object.
 * 
 * it can also be defined in cross section by
 * only supplying 2 vectors to the constructor
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class BoundingBox {
	
	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;
	
	/** 
	 * set of quads making up the box. a 2d box has 1,
	 * a 3d box has 6 (because it's cuboidal)
	 */
	private final Quad[] face;
	
	
	/**
	 * principle axes of box
	 */
	private final double[][] v;
	
	/**
	 * parameter mins along principle axes
	 */
	private final double[] min;
	
	/**
	 * parameter maxes along principle axes
	 */
	private final double[] max;
	
	
	/**
	 * constuct 3D bounding box from three (perpendicular) vectors
	 * and the max and min distance along each of them defining the
	 * corners of the cuboid
	 * 
	 * @param v three principle axis vectors
	 * @param vmin min param in each direction
	 * @param vmax max param in each direction
	 */
	protected BoundingBox(double[][] v, double[] vmin, double[] vmax){

		if(v.length!=3){
			throw new LoggedException("cuboidal bounding box requires 3 vectors but was handed "+v.length);
		}

		this.v=new double[D][D];
		this.min=new double[D];
		this.max= new double[D];
		
		for(int i=0; i<D; i++){
			for(int j=0; j<D; j++){
				this.v[i][j]= v[i][j];
			}
			this.min[i]=vmin[i];
			this.max[i]=vmax[i];
		}
		
		this.face=new Quad[6];
		
		// construct faces
		double[][] vFace= new double[2][D];
		double[] minFace= new double[2];
		double[] maxFace= new double[2];
		double[] V0= new double[D];
		
		// front face
		vFace[0]=v[0];
		vFace[1]=v[1];
		
		minFace[0]=vmin[0];
		minFace[1]=vmin[1];
		
		maxFace[0]=vmax[0];
		maxFace[1]=vmax[1];
		
		for(int i=0; i<D; i++){
			for(int j=0; j<3; j++){
				V0[i]=vmin[j]*v[j][i];
			}
		}
		
		face[0]= new Quad(vFace, minFace, maxFace, V0);		
		
		// back face
		for(int i=0; i<D; i++){
			for(int j=0; j<3; j++){
				V0[i]=vmax[j]*v[j][i];
			}
		}
		
		face[1]= new Quad(vFace, minFace, maxFace, V0);
		
		// right face
		vFace[0]=v[0];
		vFace[1]=v[2];
		
		minFace[0]=vmin[0];
		minFace[1]=vmin[2];
		
		maxFace[0]=vmax[0];
		maxFace[1]=vmax[2];
		
		face[2]= new Quad(vFace, minFace, maxFace, V0);
		
		// left face
		for(int i=0; i<D; i++){
			for(int j=0; j<3; j++){
				V0[i]=vmin[j]*v[j][i];
			}
		}

		face[3]= new Quad(vFace, minFace, maxFace, V0);
		
		// bottom face
		vFace[0]=v[1];
		vFace[1]=v[2];
		
		minFace[0]=vmin[1];
		minFace[1]=vmin[2];
		
		maxFace[0]=vmax[1];
		maxFace[1]=vmax[2];
		
		face[4]= new Quad(vFace, minFace, maxFace, V0);
		
		// top face
		for(int i=0; i<D; i++){
			for(int j=0; j<3; j++){
				V0[i]=vmax[j]*v[j][i];
			}
		}
			
		face[5]= new Quad(vFace, minFace, maxFace, V0);
		
	}

	
	protected final boolean intersectedBy(double[] walkerPos, double[] step){
		
		for(int i=0; i<face.length; i++){
			if(face[i].crossedBy(walkerPos, step)){
				return true;
			}
		}
		
		// check if step is inside box 
		if(face.length>1){
			double[] newPos= new double[D];
			for(int i=0; i<D; i++){
				newPos[i]=walkerPos[i]+step[i];
			}
			
			boolean walkerIn=contains(walkerPos);
			boolean newPosIn=contains(newPos);
			
			return walkerIn||newPosIn;
		}
		
		return false;
		
	}
	
	
	
	/**
	 * checks if a given point is inside the box or not by
	 * comparing parameter ranges
	 * 
	 * @param pos point in box coords
	 * 
	 * @return true or false
	 */
	private final boolean contains(double[] pos){

		double[] projectedPos= new double[D];
		
		for(int i=0; i<3; i++){
			double coord=0.0;
			for(int j=0; j<D; j++){
				coord+=pos[j]*v[i][j];
			}
			
			projectedPos[i]=coord;
		}
		
		boolean in= true;
		
		for(int i=0; i<D; i++){
			if((projectedPos[i]<min[i])||(projectedPos[i]>max[i])){
				in=false;
				break;
			}
		}
		
		return in;
	}
	
	
	
	/**
	 * constructs and returns the optimal bounding box for an array of
	 * triangles. The bounding box is constructed using principle 
	 * component analyis of all triangle vertices and so will be optimally
	 * orientated to be as small as possible whilst still containing all points.
	 * 
	 * PCA involves constructing the covarience matrix for the points in 
	 * the object and diagonalising. The eigenvectors are the principle
	 * axes of the box, so by projecting into this frame and finding the
	 * max and min in each component we define the optimal box.
	 * 
	 * @param triangle list of all triangles in the object
	 * 
	 * @return optimally orientated minimal bounding box
	 */
	public static BoundingBox getBoundingBoxFromTriangles(Triangle[] triangle){
		
		// find centre of mass
		int D=DiffusionSimulation.D;
		
		double[] m = new double[D];
		
		for(int i=0; i<D; i++){
			m[0]=0.0;
		}
		
		int count=0;
		
		for(int i=0; i<triangle.length; i++){
			for(int j=0; j<3; j++){
				double[] v= triangle[i].getVertex(j);
				
				for(int k=0; k<D; k++){
					m[k]+=v[k];
				}
				count++;
			}
		}
		
		for(int i=0; i<D; i++){
			m[i]/=count;
		}
		
		// construct covarience matrix
		double[][] s= new double[D][D];
		
		// construct matrix elements
		for(int i=0; i<D; i++){
			for(int j=0; j<D; j++){
				for(int k=0; k<triangle.length; k++){
					for(int l=0; l<3; l++){
						double[] v= triangle[k].getVertex(l);
						s[i][j]+=(v[i]-m[i])*(v[j]-m[j]);
					}
				}
			}
		}
		// using a diffusion tensor is a bit odd, 
		// but it's got the neccessary machinery
		Jama.Matrix S=new Jama.Matrix(s);
		
		// eigen-decomposition of S
		Jama.EigenvalueDecomposition e = new Jama.EigenvalueDecomposition(S);
            
		// extract eigenvectors
		double[][] v= e.getV().getArrayCopy();
				
		// get max and min in each direction
		double[] min= new double[D];
		double[] max= new double[D];
		for(int i=0; i<D; i++){
			min[i]=Double.MAX_VALUE;
			max[i]=-Double.MAX_VALUE;
		}
		
		for(int i=0; i<triangle.length; i++){  // each triangle
			for(int j=0; j<3; j++){                 // each vertex
				double[] vertex= triangle[i].getVertex(j);
					
				for(int k=0; k<3; k++){                    // each eVec
					double coord=0.0;
					for(int l=0; l<D; l++){                    // project vertex onto eVec
						coord+=vertex[l]*v[k][l];
					}
					
					// check against max and min
					if(min[k]>coord){
						min[k]=coord;
					}
					
					if(max[k]<coord){
						max[k]=coord;
					}
				}
			}
		}
				
		return new BoundingBox(v, min, max);
	}
	
	/**
	 * return the vertices of the corners of the cube
	 * 
	 * @return new array of vertex coords in substrate coordinates
	 */
	protected double[][] getVertices(){
		
		double[] min= new double[D];
		double[] max= new double[D];
		
		// get min and max coords in substrate coords
		for(int i=0; i<D; i++){
			for(int j=0; j<D; j++){
				min[i]+=this.min[i]*v[j][i];
				max[i]+=this.max[i]*v[j][i];
			}
		}
		
		// corners of cube are permutations of mins and maxes
		return new double[][]{{min[0], min[1], min[2]},
				              {max[0], min[1], min[2]},
				              {min[0], max[1], min[2]},
				              {max[0], max[1], min[2]},
				              {min[0], min[1], max[2]},
				              {max[0], min[1], max[2]},
				              {min[0], max[1], max[2]},
				              {max[0], max[1], max[2]}};
		
	}
	
	/** 
	 * test construction of simple cubic bounding box and
	 * intersection with various faces.
	 * 
	 */
	private static void testCubicBox(){
		
		double[][] v= new double[][]{{1.0, 0.0, 0.0},
									 {0.0, 1.0, 0.0},
									 {0.0, 0.0, 1.0}};

		double[] vmin= new double[]{0.0, 0.0, 0.0};
		double[] vmax= new double[]{1.0, 1.0, 1.0};
		
		//construct box with axes parallel to substrate frame and unit volume
		BoundingBox bbox= new BoundingBox(v, vmin, vmax);
		
		// construct walker at centre of box
		double[] pos= new double[]{0.5, 0.5, 0.5};
		
		// construct step that will end inside box
		double[] step= new double[]{0.1, 0.2, 0.3};
		
		boolean crosses= bbox.intersectedBy(pos, step);
		
		System.err.println("wholey contained step intersection is "+crosses+" (true)");
		
		// construct step that ends up outside box
		step[0]=0.6;
		
		crosses=bbox.intersectedBy(pos, step);

		System.err.println("step starting inside intersection is "+crosses+" (true)");

		// construct step that starts outside and finishes inside
		step[0]=-0.4;
		pos[0]=1.2;
		
		crosses=bbox.intersectedBy(pos, step);

		System.err.println("step starting outside intersection is "+crosses+" (true)");
		
		// construct non-intersecting step
		step[0]=0.0;
		pos[0]=1.2;
		
		crosses=bbox.intersectedBy(pos, step);

		System.err.println("non-intersecting step intersection is "+crosses+" (false)");

		// construct intersecting step that starts and finishes outside box
		pos=new double[]{0.5, 1.25, 0.5};
		step= new double[]{0.0, -1.0, 1.0};
		
		crosses=bbox.intersectedBy(pos, step);
		
		System.err.println("intersecting step that starts an finishes outside box gives "+crosses+" (true)");
	}
	
	
	private static void testRotatedBox(){
		
		double[][] vRaw= new double[][]{{1.0, 0.0, 0.0},
				                        {0.0, 1.0, 0.0},
				                        {0.0, 0.0, 1.0}};

		double[] vmin= new double[]{0.0, 0.0, 0.0};
		double[] vmax= new double[]{1.0, 1.0, 1.0};

		// rotation angle
		double theta = Math.PI/3;
		
		// rotation axis
		double[] axel = new double[]{0.0, 1/Math.sqrt(2), 1/Math.sqrt(2)};
		
		// Rotation matrix
		RealMatrix rotMat= Rotations.getRotMat(axel, theta);
		
		// space for rotated box vectors
		double[][] v= new double[3][];
		
		// transform the box axes
		for(int i=0; i<3; i++){
			v[i]= Rotations.transformPoint(rotMat, vRaw[i]);
			//v[i]=vRaw[i];
		}
		
		//construct box with axes parallel to substrate frame and unit volume
		BoundingBox bbox= new BoundingBox(v, vmin, vmax);
		
		// use same tests but rotated, so should get the same results as in cubic case
		//construct walker at centre of box
		double[] pos= new double[]{0.5, 0.5, 0.5};
		pos=Rotations.transformPoint(rotMat, pos);
		
		// construct step that will end inside box
		double[] step= new double[]{0.1, 0.2, 0.3};
		step=Rotations.transformPoint(rotMat, step);
		
		boolean crosses= bbox.intersectedBy(pos, step);
		
		System.err.println("wholey contained step intersection is "+crosses+" (true)");
		
		// construct step that ends up outside box
		step=new double[]{0.6, 0.2, 0.3};
		step=Rotations.transformPoint(rotMat, step);
		
		crosses=bbox.intersectedBy(pos, step);

		System.err.println("step starting inside intersection is "+crosses+" (true)");

		// construct step that starts outside and finishes inside
		step=new double[]{-0.4, 0.2, 0.3};
		step=Rotations.transformPoint(rotMat, step);

		pos= new double[]{1.2, 0.5, 0.5};
		pos=Rotations.transformPoint(rotMat, pos);
		
		crosses=bbox.intersectedBy(pos, step);

		System.err.println("step starting outside intersection is "+crosses+" (true)");
		
		// construct non-intersecting step
		step=new double[]{0.0, 0.2, 0.3};
		step=Rotations.transformPoint(rotMat, step);

		pos= new double[]{1.2, 0.5, 0.5};
		pos=Rotations.transformPoint(rotMat, pos);
		
		crosses=bbox.intersectedBy(pos, step);

		System.err.println("non-intersecting step intersection is "+crosses+" (false)");

		// construct intersecting step that starts and finishes outside box
		pos=new double[]{0.5, 1.25, 0.5};
		pos=Rotations.transformPoint(rotMat, pos);
		
		step= new double[]{0.0, -1.0, 1.0};
		step=Rotations.transformPoint(rotMat, step);
		
		crosses=bbox.intersectedBy(pos, step);
		
		System.err.println("intersecting step that starts an finishes outside box gives "+crosses+" (true)");
	
	}
	
	
	
	public static void main(String[] args){		
		
		testCubicBox();
		
		testRotatedBox();
	}
		
		
}

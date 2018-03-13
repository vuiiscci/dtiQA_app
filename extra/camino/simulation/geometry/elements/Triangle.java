package simulation.geometry.elements;

import java.util.ArrayList;
import java.util.logging.Logger;

import simulation.DiffusionSimulation;
import simulation.dynamics.exceptions.TooDamnCloseException;
import tools.CL_Initializer;
import misc.LoggedException;

/**
 *  
 */
public class Triangle implements SubstrateObject {
	
	/** dimensionality of space (should be 3) */
	private static final int D= DiffusionSimulation.D;
	
	/** logging object */
    private Logger logger= Logger.getLogger(this.getClass().getName());
    
    /** surface normal */
    private double[] normal;
    
    /** array of 3 vertices */
    private final double[][] vertex;

    /** vertices displaced so that box centre is origin (for box overlap check)*/
    private final double[][] w= new double[3][D];
    
    /** box half lengths (for box overlap check) */
    private final double[] boxhalfsize= new double[3];
    
    /** centre of box */
    private final double[] boxcenter= new double[D]; 
    
    /** v1-v0 */
    private final double[] u= new double[D];

    /** v2-v0 */
    private final double[] v= new double[D];
    
    /** dot product of normal and first vertex */
    private double nP;
    
    /** dot product of u with self */
    private double uu;
    
    /** dot product of v with self */
    private double vv;
    
    /** dot prod of u with v */
    private double uv;
    
    /** permeability */
    private final double p;
    
    /**
     * define a triangle from three vertices. normal calulated internally
     *
     * @param vert1 the first vertex of the triangle
     * @param vert2 the second vertex of the triangle
     * @param vert3 the third vertex of the triangle
     */
    public Triangle(double [] vert1, double [] vert2, double [] vert3, double p){
    	vertex = new double [3][D];

       	for(int j=0; j<3; j++){
       		vertex[0][j]=vert1[j];
       		vertex[1][j]=vert2[j];
       		vertex[2][j]=vert3[j];
       	}

       	this.p=p;
       	
       	initVectors();
    }
    
    /**
     * reset vertex i with the given vector
     * 
     * @param i which vertex do we set?
     * @param vert new vertex
     */
    public void setVertex(int i, double[] v){
        
        for(int j=0; j<D; j++){
            vertex[i][j]=v[j];
        }
    }
    
    /**
     * replaces all vertices with the values given
     * 
     * @param vert 3xD array of new vertices
     */
    public void setVertices(double[][] v){
        
        for(int i=0; i<3; i++){
            setVertex(i, v[i]);
        }
    }
    
    
    
    /**
     * recalculate the normal and the u and v vectors used for
     * triangle intersection 
     * 
     */
    public void initVectors(){
       	for(int i=0; i<D; i++){
       		u[i]= vertex[1][i]-vertex[0][i];
       		v[i]= vertex[0][i]-vertex[2][i];
       	}
       	
       	normal=Triangle.normCrossProd(u, v);
       	
       	nP=0.0;
       	uu=0.0;
       	vv=0.0;
       	uv=0.0;
       	for(int i=0; i<D; i++){
       		v[i]*=-1.0;
       		nP+=normal[i]*vertex[0][i];
       		uu+=u[i]*u[i];
       		vv+=v[i]*v[i];
       		uv+=u[i]*v[i];
       	}
       	

    }
    
    /**
     * returns the normalised cross product of two vectors
     * 
     * @param u a vector
     * @param v another vector
     * 
     * @return (u x v)/|u x v|
     */
    protected static final double[] normCrossProd(double[] u, double[] v){
    	
    	double[] w = new double[u.length];
    	
    	w[0]= (u[1]*v[2] - u[2]*v[1]);
    	w[1]=-(u[0]*v[2] - u[2]*v[0]);
    	w[2]= (u[0]*v[1] - u[1]*v[0]);
    	
    	double modw=0.0;
    	for(int i=0; i<w.length; i++){
    		modw+=w[i]*w[i];
    	}
    	
    	modw=Math.sqrt(modw);
    	for(int i=0; i<w.length; i++){
    		w[i]/=modw;
    	}
    	
    	return w;
    }
    
    /*
     * returns the specified vertex.
     * @param i the index of the vertex in the range 0-2
     * @return the vertex (as an array of doubles)
     */
    public double [] getVertex(int i){
        if(i>=0 || i<3){
        	return vertex[i];
        }
        else{
        	throw new LoggedException("Index to verts array out of range 0-2.  Value given: " + i);
        }
    }

    /** returns the normal of the triangle
     * @return the normal
     */
    public double[] getNormal()
    {
    	return normal;
    }

    /**
     * checks intersection of a given step with the triangle.
     * This method is from  
     * <a href='http://www.geometryalgorithms.com/Archive/algorithm_0105/algorithm_0105.htm'>
     * GeometryAlgorithms.com
     * </a>
     * 
     * it's very speedy as it only requires evaluating 5 dot products 
     * after the intersection with the plane is found.
     * 
     * it works by transforming the intersection point into plane-parameteric
     * cooordinates, with the origin at the first vertex of the 
     * triangle P0, and the basis vectors being u=P1-P0 and v=P2-P0.
     * Any point in the plane may be represented as P=P0+su+tv where
     * s and t are scalar parameters.
     * 
     * solving for sI and tI representing the point of intersection then
     * allows a straightforward test for whether the intersection is in
     * the triangle (sI>=0, tI>=0, sI+tI<=1). very neat!
     * 
     * expressions for sI and tI involve only dot products and are very 
     * easy to implement. see above link for details.
     * 
     * 
     * @param walkerPos start point of walker
     * @param step step vector
     * @param normal space to store traingle normal
     * @param d distance of interrsection from plane
     * @param origLength original length of step
     * @param intDist arcLength to interaction
     * 
     * @return true if intersection is inside triangle
     */
    public boolean crosses(double[] walkerPos, double[] step, double[] normal, double[] d, boolean ignore, 
    		double origLength, double[] intDist, boolean[] in, double[] p, double walkerRad) throws TooDamnCloseException {
    	
    	// start by checking intersection with plane containing triangle
    	double nVmP=0.0;       // normal dot (V0 - walkerpos) 
    	double nStep=0.0;	   // normal dot step
    	
    	double nDotPosPlusStep=0.0;
    	
    	// calculate lengyel's 4D dot product of vector & plane
    	// this equals dp of normal and vector plus plane constant
    	for(int i=0; i<D; i++){
    		nVmP+=this.normal[i]*(vertex[0][i]-walkerPos[i]);
    		nStep+=this.normal[i]*step[i];
    		
    		nDotPosPlusStep+=this.normal[i]*(vertex[0][i]-(walkerPos[i]+step[i]));
    	}
    	
    	
    	// if n dot step is zero (within tolerence) there's no plane intersection
    	if(Math.abs(nStep)<=1E-14){
    		return false;
    	}

    	// if the distance to the plane is less than a walker radius AFTER THE STEP
    	// we need to do something about it.
    	//
    	// TODO: I have no idea what to do about it!
    	if(Math.abs(nDotPosPlusStep)<walkerRad){
    		throw new TooDamnCloseException("Step takes walker within a radius of the barrier.");
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
    	
    	
    	// if we've got here the step passes through the plane containing
    	// the triangle, now must check if the intersection point is 
    	// contained in the triangle itself.
    	double[] w= new double[D];
    	
    	for(int i=0; i<D; i++){
    		// displacement of interaction point from zeroth vertex
    		w[i]=walkerPos[i]+tInt*step[i]-vertex[0][i];
    	}
    	
    	// evaluate dot products
    	double wv=0.0;
    	double wu=0.0;
    	
    	for(int i=0; i<D; i++){
    		wv+=w[i]*v[i];
    		wu+=w[i]*u[i];
    	}
    	
    	// evaluate triangle coord parameters
    	double denom= (uv*uv) - (uu*vv);
    	
    	double s=(uv*wv - vv*wu)/denom;

    	// if s less than 0 or greater than one, reject
    	if(s<0.0){
    		return false;
    	}
    	if(s>1.0){
    		return false;
    	}
    	
    	double t=(uv*wu - uu*wv)/denom;
    	
    	// check parameter inequalities
		if(t>=0.0){
			if(s+t<=1.0){
				// if we're in here, we're in the triangle
				// so set the geometric quantities...
				
				// d is the distance of the intersection from the
				// origin dotted with the surface normal
				d[0]=0.0;
		    	for(int i=0; i<D; i++){
		    		double Pint_i= walkerPos[i]+tInt*step[i];
		    		d[0]+=Pint_i*this.normal[i];
		    	}
		    	
		    	// int dist is the arclength to the intersection
				intDist[0]=tInt;
				
				// and the the normal is... err... the normal
				for(int i=0; i<D; i++){
					normal[i]=this.normal[i];
				}
				p[0]=this.p;
				
				// ...and leave
				return true;
			}
		}
    	
    	// otherwise we've failed the inequality test
    	return false;
    }
    
    public final BoundingBox getBoundingBox(){
    	logger.warning("attenmpt to retreive bounding box from individual triangle!");
    	logger.warning("returning null.");
    	
    	return null;
    }
    
    public final double getDiffusivityAt(double[] pos){
    	logger.warning("attempt to query individual triangle for diffusivity!");
    	logger.warning("returning extracellular diffusivity");
    	
    	return CL_Initializer.DIFF_CONST;
    }
    
    public final boolean inside(double[] pos){
    	logger.warning("individual triangle queried for inside() calculation!");
    	logger.warning("returning false");
    	
    	return false;
    }
    
    /** 
     * checks if triangle intersects an axis-aligned cubic region
     * 
     * @param bottomLeft lower corner of box
     * @param topRight upper corner of box
     * 
     * @return true or false
     * 
     * TODO: not implemented.
     */
    public final boolean intersectsCubicRegion(double[] bottomLeft, double[] topRight){
    	
    	for(int i=0; i<D; i++){
    		boxcenter[i]=(topRight[i]+bottomLeft[i])/2.0;
    		boxhalfsize[i]=(topRight[i]-bottomLeft[i])/2.0;
    	}
    	
    	return triBoxOverlap();
    }
    
    public final double getPermeability(int i){
    	return p;
    }
    
    public final boolean boundingBoxIntersects(double[] pos, double[] step){
    	return true;
    }
    
    /**
     * returns the details of the triangle as a string
     * @return a string conataining the verts and normal information
     */
    public String toString(){
    	String s = "vertices: (";
    	for(int i =0; i<3;i++){
    		for(int j=0;j<3;j++){
    			s += " " + vertex[i][j];
    		}
    		s += ")\t";
    	}

    	s +="\tnormal: (";
    	s += normal[0] + " " + normal[1] + " " + normal[2] + ")";

    	return s;
    }
    
    
    /**
     * constructs triangles in various planes and prints out the 
     * normals to check the norm cross product routine is working
     * correctly.
     *
     */
    private static void testNormalCalc(){
        	
    	double[] vert1= new double[]{0.0, 0.0, 0.0};
    	double[] vert2= new double[]{0.0, 1.0, 0.0};
    	double[] vert3= new double[]{0.5, 0.5, 0.0};
    	
    	Triangle triangle= new Triangle(vert1, vert2, vert3, 0.0);
    	
    	double[] normal= triangle.getNormal();
    	
    	System.err.println("triangle 1 (xy-plane) normal is: ("+normal[0]+","+normal[1]+","+normal[2]+")");

    	vert3= new double[]{0.0, 0.5, 0.5};
    	
    	triangle= new Triangle(vert1, vert2, vert3, 0.0);
    	
    	normal= triangle.getNormal();

    	System.err.println("triangle 2 (yz-plane) normal is: ("+normal[0]+","+normal[1]+","+normal[2]+")");

    	vert2= new double[]{1.0, 0.0, 0.0};
    	vert3= new double[]{0.5, 0.0, 0.5};
    	
    	triangle= new Triangle(vert1, vert2, vert3, 0.0);
    	
    	normal= triangle.getNormal();

    	System.err.println("triangle 3 (xz-plane) normal is: ("+normal[0]+","+normal[1]+","+normal[2]+")");
	
    }
    
    /**
     *  the following code is a java port of thomas akenine-moller's triangle-box
     *  intersection test. @see http://jgt.akpeters.com/papers/AkenineMoller01/
     *  
     *  port by matt.
     */
    private final int X=0;
    private final int Y=1;
    private final int Z=2;


    private final void CROSS(double[] dest, double[] v1, double[] v2){
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1];
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2];
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
    } 



    private final double DOT(double[] v1, double[] v2){
    	return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
    }



    private final void SUB(double[] dest, double[] v1, double[] v2){
          dest[0]=v1[0]-v2[0];
          dest[1]=v1[1]-v2[1];
          dest[2]=v1[2]-v2[2];
    }


    private final void FINDMINMAX(double x0, double x1, double x2, double min, double max){
    	min = max = x0;
    	if(x1<min) min=x1;
    	if(x1>max) max=x1;
    	if(x2<min) min=x2;
    	if(x2>max) max=x2;
    }


    private final boolean planeBoxOverlap(double[] normal, double d, double[] maxbox){

    	int q;

		double[] vmin= new double[D];
		double[] vmax= new double[D];
	
		for(q=X;q<=Z;q++){
		    if(normal[q]>0.0f)
			{
			    vmin[q]=-maxbox[q];
			    vmax[q]=maxbox[q];
			}
		    else
			{
			    vmin[q]=maxbox[q];
			    vmax[q]=-maxbox[q];
			}
		}
	
		if(DOT(normal,vmin)+d>0.0f) return false;
		if(DOT(normal,vmax)+d>=0.0f) return true;
		
		return false;
    }


    /*======================== X-tests ========================*/
    private final boolean  AXISTEST_X01(double a, double b, double fa, double fb){
	
		double min, max;
														       
		double p0 = a*w[0][Y] - b*w[0][Z];
		double p2 = a*w[2][Y] - b*w[2][Z];
	
	        if(p0<p2) {
		    min=p0; max=p2;
		}
		else {
		    min=p2; max=p0;
		}
	
		double rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];
		
		if(min>rad || max<-rad) 
		    return false;
	
		return true;
    }



    private final boolean AXISTEST_X2(double a, double b, double fa, double fb){
	
    	double min, max;

		double p0 = a*w[0][Y] - b*w[0][Z];
		double p1 = a*w[1][Y] - b*w[1][Z];
	        
	
		if(p0<p1) {
		    min=p0; 
		    max=p1;
		} else {
		    min=p1; 
		    max=p0;
		}
		
		double rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];
		if(min>rad || max<-rad) 
		    return false;
	
		return true;
    }


    /*======================== Y-tests ========================*/
    private final boolean AXISTEST_Y02(double a, double b, double fa, double fb){

    	double min, max;

    	double p0 = -a*w[0][X] + b*w[0][Z];
    	double p2 = -a*w[2][X] + b*w[2][Z];

        if(p0<p2) {
        	min=p0; 
        	max=p2;
        } 
        else {
        	min=p2; 
        	max=p0;
        }

        double rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];

        if(min>rad || max<-rad) 
        	return false;

        return true;
    }


    private final boolean AXISTEST_Y1(double a, double b, double fa, double fb){

		double min, max;
	
		double p0 = -a*w[0][X] + b*w[0][Z];
		double p1 = -a*w[1][X] + b*w[1][Z];
	
        if(p0<p1) {
		    min=p0; 
		    max=p1;
		} 
		else {
		    min=p1; 
		    max=p0;
		}
	
		double rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];
	
		if(min>rad || max<-rad) 
		    return false;
	
		return true;
    }



    /*======================== Z-tests ========================*/

    private final boolean AXISTEST_Z12(double a, double b, double fa, double fb){
	
		double min, max;
	
		double p1 = a*w[1][X] - b*w[1][Y];
		double p2 = a*w[2][X] - b*w[2][Y];
	
	    if(p2<p1) {
		    min=p2; 
		    max=p1;
		} 
		else {
		    min=p1; 
		    max=p2;
		}
	
		double rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];
		
		if(min>rad || max<-rad) 
		    return false;
	
		return true;
    }



    private final boolean AXISTEST_Z0(double a, double b, double fa, double fb){

		double min, max;
	
		double p0 = a*w[0][X] - b*w[0][Y];
		double p1 = a*w[1][X] - b*w[1][Y];
	
	    if(p0<p1) {
		    min=p0; 
		    max=p1;
		} 
		else {
		    min=p1; 
		    max=p0;
		}
	
		double rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];
		
		if(min>rad || max<-rad) 
		    return false;
	
		return true;
    }



    public final boolean triBoxOverlap(){
	
		/*    use separating axis theorem to test overlap between triangle and box */
		/*    need to test for overlap in these directions: */
		/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
		/*       we do not even need to test these) */
		/*    2) normal of the triangle */
		/*    3) crossproduct(edge from tri, {x,y,z}-directin) */
		/*       this gives 3x3=9 more tests */
    	double[] axis= new double[D];
		double min=0;
		double max=0;
		double d,p0,p1,p2,rad,fex,fey,fez;  
		double[] normal= new double[D];
		double[] e0= new double[D];
		double[] e1= new double[D];
		double[] e2= new double[D];
		
		
		
		/* This is the fastest branch on Sun */
		/* move everything so that the boxcenter is in (0,0,0) */
		SUB(w[0], vertex[0],boxcenter);
		SUB(w[1], vertex[1],boxcenter);
		SUB(w[2], vertex[2],boxcenter);
		
		/* compute triangle edges */
		SUB(e0,w[1],w[0]);      /* tri edge 0 */
		SUB(e1,w[2],w[1]);      /* tri edge 1 */
		SUB(e2,w[0],w[2]);      /* tri edge 2 */
		
		/* Bullet 3:  */
		/*  test the 9 tests first (this was faster) */
		fex = Math.abs(e0[X]);
		fey = Math.abs(e0[Y]);
		fez = Math.abs(e0[Z]);
		
		if(!AXISTEST_X01(e0[Z], e0[Y], fez, fey)){
		    return false;
		}
		if(!AXISTEST_Y02(e0[Z], e0[X], fez, fex)){
		    return false;
		}
		if(!AXISTEST_Z12(e0[Y], e0[X], fey, fex)){
		    return false;
		}
		
		fex = Math.abs(e1[X]);
		fey = Math.abs(e1[Y]);
		fez = Math.abs(e1[Z]);
		if(!AXISTEST_X01(e1[Z], e1[Y], fez, fey)){
		    return false;
		}
		if(!AXISTEST_Y02(e1[Z], e1[X], fez, fex)){
		    return false;
		}
		if(!AXISTEST_Z0(e1[Y], e1[X], fey, fex)){
		    return false;
		}
		
		fex = Math.abs(e2[X]);
		fey = Math.abs(e2[Y]);
		fez = Math.abs(e2[Z]);
		if(!AXISTEST_X2(e2[Z], e2[Y], fez, fey)){
		    return false;
		}
		if(!AXISTEST_Y1(e2[Z], e2[X], fez, fex)){
		    return false;
		}
		if(!AXISTEST_Z12(e2[Y], e2[X], fey, fex)){
		    return false;
		}
		
		/* Bullet 1: */
		/*  first test overlap in the {x,y,z}-directions */
		/*  find min, max of the triangle each direction, and test for overlap in */
		/*  that direction -- this is equivalent to testing a minimal AABB around */
		/*  the triangle against the AABB */
		
		/* test in X-direction */
		FINDMINMAX(w[0][X],w[1][X],w[2][X],min,max);
		if(min>boxhalfsize[X] || max< -boxhalfsize[X]) return false;
		
		/* test in Y-direction */
		FINDMINMAX(w[0][Y],w[1][Y],w[2][Y],min,max);
		if(min>boxhalfsize[Y] || max< -boxhalfsize[Y]) return false;
		
		/* test in Z-direction */
		FINDMINMAX(w[0][Z],w[1][Z],w[2][Z],min,max);
		if(min>boxhalfsize[Z] || max< -boxhalfsize[Z]) return false;
		
		/* Bullet 2: */
		/*  test if the box intersects the plane of the triangle */
		/*  compute plane equation of triangle: normal*x+d=0 */
		CROSS(normal,e0,e1);
		d=-DOT(normal,w[0]);  /* plane eq: normal.x+d=0 */
		if(!planeBoxOverlap(normal,d,boxhalfsize)) return false;
		
		return true;   /* box and triangle overlaps */
    }


    public final double getDistanceFrom(double[] r0){
    	
    	// first, get the point projection into the plane
    	final double[] w= new double[D];
    	
    	double qMpDotn= 0.0;
    	
    	for(int i=0; i<D; i++){
    		qMpDotn+=(r0[i]-vertex[0][i])*normal[i];
    	}
    	
    	for(int i=0; i<D; i++){
    		w[i]= r0[i]-qMpDotn*normal[i];
    	}
    	
    	// check if projected point is in the triangle or not
    	// evaluate dot products
    	double wv=0.0;
    	double wu=0.0;
    	
    	for(int i=0; i<D; i++){
    		wv+=w[i]*v[i];
    		wu+=w[i]*u[i];
    	}
    	
    	// evaluate triangle coord parameters
    	double denom= (uv*uv) - (uu*vv);
    	
    	double s=(uv*wv - vv*wu)/denom;    	
    	double t=(uv*wu - uu*wv)/denom;
    	
    	// check parameter inequalities
		if((t>=0.0)&&(s+t<=1.0)){
			// if we're in here, we're in the triangle
			// so the distance is the projected point to the original
			double distSq=0.0;
			for(int i=0; i<D; i++){
				distSq+=(r0[i]-w[i])*(r0[i]-w[i]);
			}

			return Math.sqrt(distSq);
    	}	
		
		// if we're here, the projection isn't in the triangle,
		// so we find the closest point in the triangle by
		// clamping the barycentric coordinates and return
		// the distance from that point to the original
		final double[] rc= new double[D];
		
		// clamp s and t
		if(t<0.0){
			t=0.0;
		}
		if(t>1.0){
			t=1.0;
		}
		if(s<0.0){
			s=0.0;
		}
		if(s>1.0-t){
			s=1.0-t;
		}
		
		// transform barycentric coords of clamped point back to substrate frame
		for(int i=0; i<D; i++){
			rc[i]=vertex[0][i]+s*u[i]+t*v[i];
		}
		
		double sqDist=0.0;
		for(int i=0; i<D; i++){
			sqDist+=(r0[i]-rc[i])*(r0[i]-rc[i]);
		}
    
		return Math.sqrt(sqDist);
    }    
    
    /**
     * returns the triangles associated with this object. In this case that just means
     * a copy of the current triangle in an ArrayList
     * 
     * @return a copy of the current triangle (different instance)
     */
    public ArrayList<Triangle> getTriangles(){
    	
    	ArrayList<Triangle> triangles= new ArrayList<Triangle>(1);
    	
    	Triangle copy= new Triangle(vertex[0], vertex[1], vertex[2], p);
    	
    	triangles.add(copy);
    	
    	return triangles;
    }
    
    /**
     * @return copy of the first vertex (NOT the centroid)
     */
    public final double[] getPosition(){
    	
    	return new double[]{vertex[0][0], vertex[0][1], vertex[0][2]};
    }

    /**
     * checks the triangle intersention code in each of 
     * five cases:
     * 
     * short: ray of step intersects triangle but segment
     *        is too short to make contact
     *        
     * opp: ray intersects triangle but step is in the 
     *      opposite direction (has caused problems in
     *      the past)
     *      
     * para: step is parallel to plane -- no intersection
     * 
     * plane: step intersects the plane of the triangle,
     *        but not the triangle itself
     *        
     * tri: step intersects the triangle.
     *
     * aside from a step entirely in the plane of the 
     * triangle (which is vanishingly rare numerically), 
     * this should cover all the relevent cases.
     *
     */
    private static void testLineTriangleIntersection(){
    	
    	double[] v0= new double[]{1.0, 1.0, 0.0};
    	double[] v1= new double[]{0.0, 1.0, 0.0};
    	double[] v2= new double[]{0.5, 0.5, 0.0};
    	
    	// construct triangle & check normal
    	Triangle triangle = new Triangle(v0, v1, v2, 0.0);
    	
    	double[] normal= triangle.getNormal();
    	
    	System.err.println("triangle has "+triangle);
    	System.err.println("normal is ("+normal[0]+","+normal[1]+","+normal[2]+")");

    	// space for position, step, distance, normal etc
    	double[] walkerpos= new double[D];
    	double[] step= new double[D];
    	double[] d= new double[1];
    	double[] intDist= new double[1];
    	double[] p= new double[1];
    	boolean[] in = new boolean[1];
    	
    	// short case
    	walkerpos[0]=0.5;
    	walkerpos[1]=0.75;
    	walkerpos[2]=0.5;
    	
    	step[0]=0.0;
    	step[1]=0.0;
    	step[2]=-0.4;
    	
    	try{
	    	boolean crosses= triangle.crosses(walkerpos, step, normal, d, false, 0.4, intDist, in, p, 0.004);
	    	
	    	System.err.println("short case: crosses= "+crosses);
	    	
	    	// opposite case
	    	walkerpos[2]=0.1;
	    	
	    	step[2]=0.4;
	    	
	    	crosses= triangle.crosses(walkerpos, step, normal, d, false, 0.4, intDist, in, p, 0.004);
	    	
	    	System.err.println("opposite direction case: crosses= "+crosses);
	    	
	    	// parallel case
	    	walkerpos[2]=1.0;
	    	
	    	step[0]=0.5;
	    	step[2]=0.0;
	
	    	crosses= triangle.crosses(walkerpos, step, normal, d, false, 0.5, intDist, in, p, 0.004);
	
	    	System.err.println("parallel case: crosses= "+crosses);
	    	
	    	// plane case
	    	walkerpos[0]=1.0;
	    	walkerpos[1]=0.5;
	    	walkerpos[2]=0.5;
	    	
	    	step[0]=0.0;
	    	step[1]=0.0;
	    	step[2]=-1.0;
	    	
	    	crosses= triangle.crosses(walkerpos, step, normal, d, false, 1.0, intDist, in, p, 0.004);
	    	System.err.println("plane intersection case: crosses= "+crosses);
	    	
	    	// triangle intersection case
	    	walkerpos[0]=0.5;
	    	walkerpos[1]=0.75;
	    	walkerpos[2]=0.5;
	    	
	    	crosses= triangle.crosses(walkerpos, step, normal, d, false, 1.0, intDist, in, p, 0.004);
	    	System.err.println("triangle intersection case: crosses= "+crosses);
	    	System.err.println("\tnormal= ("+normal[0]+","+normal[1]+","+normal[2]+")");
	    	System.err.println("\td= "+d[0]);
	    	System.err.println("\tintDist= "+intDist[0]);
    	}
    	catch(TooDamnCloseException tdce){
    		throw new LoggedException(tdce);
    	}

    }
    
    public static void testStepIntersection(){

        
        double[] normal = new double[D];
        double[] d= new double[1];
        double[] intDist= new double[1];
        double[] p= new double[1];
        boolean[] in = new boolean[1]; 
        
        double[] v0= new double[]{1.5773502691896257, 0.42264973081037416, 1.5773502691896257};
        double[] v1= new double[]{0.42264973081037416, 1.5773502691896257, 1.5773502691896257};
        double[] v2= new double[]{1.5773502691896257, 1.5773502691896257, 1.5773502691896257};
        
        Triangle triangle= new Triangle(v0, v1, v2, 0.0);
        
        double[] walkerpos= new double[]{1.1, 1.0, 1.0};
        double[] step= new double[]{0.0, 0.0, 2.4};
        
        try{
	        boolean crossing= triangle.crosses(walkerpos, step, normal, d, false, 2.4, intDist, in, p, 0.004);
	        
	        System.err.println("crossing is "+crossing);
	        if(crossing){
	            System.err.println("d= "+d[0]);
	            System.err.println("intDist= "+intDist[0]);
	        }
        }
        catch(TooDamnCloseException tdce){
        	throw new LoggedException(tdce);
        }
        
    }
    
    /**
     * checks intersection between triangle completely contained 
     * in cubic region
     */
    public static void testBoxIntersection1(){
        
        System.err.println("case 1: triangle inside box");
        
        // construct triangle
        double[] v1= new double[]{0.15, 0.15, 0.11};
        double[] v2= new double[]{0.17, 0.17, 0.12};
        double[] v3= new double[]{0.12, 0.12, 0.11};
        Triangle triangle= new Triangle(v1, v2, v3, 0.0);
        
        // construct box
        double[] bottomLeft= new double[]{0.1, 0.1, 0.1};
        double[] topRight= new double[]{0.25, 0.2, 0.2};
        
        
        System.err.println("triangle constructed as:");
        for(int i=0; i<3; i++){
            double[] v=triangle.getVertex(i);
            System.err.print("v"+i+" (");
            for(int j=0; j<D; j++){
                System.err.print(v[j]+" ");
            }
            System.err.println(")");
        }
        double[] n= triangle.getNormal();
        System.err.print("normal = (");
        for(int j=0; j<D; j++){
            System.err.print(n[j]+" ");
        }
        System.err.println(")");
        
        System.err.println();
        
        System.err.println("box constructed as:");
        System.err.print("topRight: (");
        for(int i=0; i<D; i++){
            System.err.print(topRight[i]+" ");
        }
        System.err.println(")");
        System.err.print("bottomLeft: (");
        for(int i=0; i<D; i++){
            System.err.print(bottomLeft[i]+" ");
        }
        System.err.println(")");
        
        
        boolean intersection= triangle.intersectsCubicRegion(bottomLeft, topRight);
        System.err.println("\nintersection is "+intersection+ " (expected true)\n\n\n");
        
    }

    
    
    /**
     * checks intersection between triangle with two
     * vertices inside cubic region
     */
    public static void testBoxIntersection2(){
        
        System.err.println("case 2: triangle has two vertices inside box");
        
        // construct triangle
        double[] v1= new double[]{0.15, 0.15, 0.11};
        double[] v2= new double[]{0.17, 0.17, -0.12};
        double[] v3= new double[]{0.12, 0.12, 0.11};
        Triangle triangle= new Triangle(v1, v2, v3, 0.0);
        
        // construct box
        double[] bottomLeft= new double[]{0.1, 0.1, 0.1};
        double[] topRight= new double[]{0.25, 0.2, 0.2};
        
        
        System.err.println("triangle constructed as:");
        for(int i=0; i<3; i++){
            double[] v=triangle.getVertex(i);
            System.err.print("v"+i+" (");
            for(int j=0; j<D; j++){
                System.err.print(v[j]+" ");
            }
            System.err.println(")");
        }
        double[] n= triangle.getNormal();
        System.err.print("normal = (");
        for(int j=0; j<D; j++){
            System.err.print(n[j]+" ");
        }
        System.err.println(")");
        
        System.err.println();
        
        System.err.println("box constructed as:");
        System.err.print("topRight: (");
        for(int i=0; i<D; i++){
            System.err.print(topRight[i]+" ");
        }
        System.err.println(")");
        System.err.print("bottomLeft: (");
        for(int i=0; i<D; i++){
            System.err.print(bottomLeft[i]+" ");
        }
        System.err.println(")");
        
        
        boolean intersection= triangle.intersectsCubicRegion(bottomLeft, topRight);
        System.err.println("\nintersection is "+intersection+ " (expected true)\n\n\n");
        
    }

    /**
     * checks intersection between triangle with one
     * vertex inside cubic region
     */
    public static void testBoxIntersection3(){
        
        System.err.println("case 3: triangle has one vertex inside box");
        
        // construct triangle
        double[] v1= new double[]{0.15, 0.15, 0.11};
        double[] v2= new double[]{0.17, 0.17, -0.12};
        double[] v3= new double[]{0.12, 0.12, -0.11};
        Triangle triangle= new Triangle(v1, v2, v3, 0.0);
        
        // construct box
        double[] bottomLeft= new double[]{0.1, 0.1, 0.1};
        double[] topRight= new double[]{0.25, 0.2, 0.2};
        
        
        System.err.println("triangle constructed as:");
        for(int i=0; i<3; i++){
            double[] v=triangle.getVertex(i);
            System.err.print("v"+i+" (");
            for(int j=0; j<D; j++){
                System.err.print(v[j]+" ");
            }
            System.err.println(")");
        }
        double[] n= triangle.getNormal();
        System.err.print("normal = (");
        for(int j=0; j<D; j++){
            System.err.print(n[j]+" ");
        }
        System.err.println(")");
        
        System.err.println();
        
        System.err.println("box constructed as:");
        System.err.print("topRight: (");
        for(int i=0; i<D; i++){
            System.err.print(topRight[i]+" ");
        }
        System.err.println(")");
        System.err.print("bottomLeft: (");
        for(int i=0; i<D; i++){
            System.err.print(bottomLeft[i]+" ");
        }
        System.err.println(")");
        
        
        boolean intersection= triangle.intersectsCubicRegion(bottomLeft, topRight);
        System.err.println("\nintersection is "+intersection+ " (expected true)\n\n\n");
        
    }

    
    /**
     * checks intersection between triangle with no vertices
     * inside box that does nonetheless intersect the region
     */
    public static void testBoxIntersection4(){
        
        System.err.println("case 3: intersecting triangle with no vertices inside box");
        
        // construct triangle
        double[] v1= new double[]{0.15, 0.15, 0.25};
        double[] v2= new double[]{0.15, 0.15, 0.0};
        double[] v3= new double[]{0.3, 0.15, 0.15};
        Triangle triangle= new Triangle(v1, v2, v3, 0.0);
        
        // construct box
        double[] bottomLeft= new double[]{0.1, 0.1, 0.1};
        double[] topRight= new double[]{0.25, 0.2, 0.2};
        
        
        System.err.println("triangle constructed as:");
        for(int i=0; i<3; i++){
            double[] v=triangle.getVertex(i);
            System.err.print("v"+i+" (");
            for(int j=0; j<D; j++){
                System.err.print(v[j]+" ");
            }
            System.err.println(")");
        }
        double[] n= triangle.getNormal();
        System.err.print("normal = (");
        for(int j=0; j<D; j++){
            System.err.print(n[j]+" ");
        }
        System.err.println(")");
        
        System.err.println();
        
        System.err.println("box constructed as:");
        System.err.print("topRight: (");
        for(int i=0; i<D; i++){
            System.err.print(topRight[i]+" ");
        }
        System.err.println(")");
        System.err.print("bottomLeft: (");
        for(int i=0; i<D; i++){
            System.err.print(bottomLeft[i]+" ");
        }
        System.err.println(")");
        
        
        boolean intersection= triangle.intersectsCubicRegion(bottomLeft, topRight);
        System.err.println("\nintersection is "+intersection+ " (expected true)\n\n\n");
        
    }

    
    /**
     * checks intersection between triangle with no vertices
     * inside box that does nonetheless intersect the region
     */
    public static void testBoxIntersection5(){
        
        System.err.println("case 3: triangle outside box");
        
        // construct triangle
        double[] v1= new double[]{0.3, 0.25, 0.25};
        double[] v2= new double[]{0.35, 0.25, 0.0};
        double[] v3= new double[]{0.3, 0.25, 0.15};
        Triangle triangle= new Triangle(v1, v2, v3, 0.0);
        
        // construct box
        double[] bottomLeft= new double[]{0.1, 0.1, 0.1};
        double[] topRight= new double[]{0.25, 0.2, 0.2};
        
        
        System.err.println("triangle constructed as:");
        for(int i=0; i<3; i++){
            double[] v=triangle.getVertex(i);
            System.err.print("v"+i+" (");
            for(int j=0; j<D; j++){
                System.err.print(v[j]+" ");
            }
            System.err.println(")");
        }
        double[] n= triangle.getNormal();
        System.err.print("normal = (");
        for(int j=0; j<D; j++){
            System.err.print(n[j]+" ");
        }
        System.err.println(")");
        
        System.err.println();
        
        System.err.println("box constructed as:");
        System.err.print("topRight: (");
        for(int i=0; i<D; i++){
            System.err.print(topRight[i]+" ");
        }
        System.err.println(")");
        System.err.print("bottomLeft: (");
        for(int i=0; i<D; i++){
            System.err.print(bottomLeft[i]+" ");
        }
        System.err.println(")");
        
        
        boolean intersection= triangle.intersectsCubicRegion(bottomLeft, topRight);
        System.err.println("\nintersection is "+intersection+ " (expected false)\n\n\n");
        
    }

    
    
    public static void main(String[] args){
    	
    	//testNormalCalc();
    
    	//testLineTriangleIntersection();
        
        //testStepIntersection();
                
        testBoxIntersection1();
        
        testBoxIntersection2();
        
        testBoxIntersection3();
        
        testBoxIntersection4();
        
        testBoxIntersection5();
    }
    
}

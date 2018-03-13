package numerics;

/**
 * <dl>
 * <dt> Purpose: 
 * <dd> To encapsulate 3D cartesian vector.
 * </dl>
 * @author Philip Cook
 * @version $Id$
 */
public class Vector3D {

    public final double x;
    public final double y;
    public final double z;

    /**
     * Get a Cartesian vector from the spherical polar coordinates. 
     *
     * Coordinate system defined as in Boas, Mathematical Methods in the Physical Sciences.
     *
     * @param r the magnitude of the vector
     * @param theta cos theta is the dot product of this vector.normalized() with the Z axis.
     * @param phi the angle with the positive X axis, between 0 and 2 * PI.
     */
    public static Vector3D vectorFromSPC(double r, double theta, double phi) {

	double rst = r * Math.sin(theta);

	return new Vector3D(rst * Math.cos(phi),
			    rst * Math.sin(phi),
			    r * Math.cos(theta));



    }


    public Vector3D() {
	x = 0.0;
	y = 0.0;
	z = 0.0;
    }

    public Vector3D(double X, double Y, double Z) {
	x = X;
	y = Y;
	z = Z;
    }

    
    public Vector3D(double[] xyz) {
	x = xyz[0];
	y = xyz[1];
	z = xyz[2];
    }

    public Vector3D(Vector3D v) {
	x = v.x;
	y = v.y;
	z = v.z;
    }
    

  
    /** Construct vector (p1 - p0) from two points.
     * @param p1 
     * @param p0 the points to use.
     */
    public Vector3D(Point3D p1, Point3D p0) {
	x = p1.x - p0.x;
	y = p1.y - p0.y;
	z = p1.z - p0.z;
    }
    

    /** Gives scalar (dot) product of this with another vector. 
     * @param v the vector to take the scalar product with.
     * @return the scalar product of this with v.
     */
    public double dot(Vector3D v) {
	return (v.x * x + v.y * y + v.z * z);
    }

    /** Gives vector (cross) product of this with another vector. 
     * @param v the vector to cross with.
     * @return this CROSS v.
     */
    public Vector3D cross(Vector3D v) {

	// v1(x1, y1, z1) CROSS v2(x2, y2, z2) == (y1 * z2 - y2 * z1), (x2 * z1 - x1 * z2), (x1 * y2 - x2 * y1)
	return new Vector3D(y * v.z - v.y * z, v.x * z - x * v.z, x * v.y - v.x * y);

    }

    
    /**
     * Get the negated version of this vector
     * @return the vector (-x, -y, -z).
     */
    public Vector3D negated() {
	return new Vector3D(-x, -y, -z);
    }

    /** Add another vector to this one.
     * @param v the vector to add.
     * @return the vector sum of this and v.
     */
    public Vector3D plus(Vector3D v) {
	return new Vector3D(x + v.x, y + v.y, z + v.z);
    }


    /** Subtract another vector from this one.
     * @param v the vector to subtract from this.
     */
    public Vector3D minus(Vector3D v) {
	return new Vector3D(x - v.x, y - v.y, z - v.z);
    }


    /**
     * @return this vector divided by its modulus.
     */
    public Vector3D normalized() {
	double mod = Math.sqrt(x * x + y * y + z * z);
	return new Vector3D(x / mod, y / mod, z / mod);
    }


    /**
     * @return the modulus of this vector. 
     */
    public double mod() {
	return Math.sqrt(x * x + y * y + z * z);
    }


    /**
     * @return the squared modulus of this vector. 
     */
    public double modSquared() {
	return x * x + y * y + z * z;
    }

    
    /** Scale vector by some constant.
     * @param factor the factor to scale by.
     * @return scaled vector.
     */
    public Vector3D scaled(double factor) {
	return new Vector3D(x * factor, y * factor, z * factor);
    }


    public String toString() {
	return "Vector3D: " + x + "\t" + y + "\t" + z + "\t";
    } 
    
    /**
     * Tests equality to within an epsilon of 1E-10
     */
    public boolean equals(Object o) {

	if (o == this) {
	    return true;
	}
	if (! (o instanceof Vector3D) ) {
	    return false;
	}
	if (o == null) {
	    return false;
	}
	
        Vector3D vec = (Vector3D)o;

        double eps = 1E-10;

	return (Math.abs(x - vec.x) < eps && Math.abs(y - vec.y) < eps && Math.abs(z - vec.z) < eps);
    }


    public int hashCode() {
	return (int)(13.0 * x + 29.0 * y + 31.0 * z);
    }


    /**
     * @return a 3x1 RealMatrix.
     */
    public RealMatrix toRealMatrix() {
	RealMatrix r = new RealMatrix(3,1);
	r.entries[0][0] = x;
	r.entries[1][0] = y;
	r.entries[2][0] = z;
	
	return r;
    }


    /**
     * @return a 3x1 Jama.Matrix.
     */
    public Jama.Matrix toJamaMatrix() {
	
	Jama.Matrix j = new Jama.Matrix(3,1);
	
	j.set(0,0, x);
	j.set(1,0, y);
	j.set(2,0, z);

	return j;

    }


    /**
     * @param v Should be a unit vector but no warning will be given if it is not.
     * @return {theta, phi}. Names of angles are the same as for #VectorFromSPC(double, double, double)
     */
    public static double[] thetaPhi(Vector3D v) {

	double[] tp = new double[2];
	
	tp[0] = Math.acos( v.z );
	
	// phi goes from 0.0 (+x axis) and wraps at 2 * PI
	// theta goes from 0.0 (+z axis) and wraps at PI

	// if x and y are 0.0 or very close, return phi == 0
	if (Math.abs(v.x) + Math.abs(v.y) < 1E-9) {
	    tp[1] = 0.0; 
	}
	else {
	    
	    // ie, if ( x == 0 && y == 0 ) == false

	    if (v.y == 0.0) {

		if (v.x > 0.0) {
		    tp[1] = 0.0;
		}
		else {
		    tp[1] = Math.PI;
		}
		
		
	    }
	    else if (v.x == 0.0) {
		
		// avoid div by zero
		if (v.y > 0) {
		    tp[1] = Math.PI / 2.0;
		}
		else {
		    tp[1] = 1.5 * Math.PI;
		}
	    }
	    else if (v.x > 0.0 && v.y > 0.0) { // first quadrant
		tp[1] = Math.atan(v.y / v.x);
	    }
	    else if (v.x < 0.0 && v.y > 0.0) { // second quadrant
		tp[1] = Math.PI + Math.atan(v.y / v.x);
	    }
	    else if (v.x < 0.0 && v.y < 0.0) { // third quadrant
		tp[1] =  Math.PI + Math.atan(v.y / v.x);
	    }
	    else { // fourth quadrant
		tp[1] = 2.0 * Math.PI + Math.atan(v.y / v.x); 
	    }

	}
	
	return tp;
    }

    /**
     * @param alpha angle of rotation in radians
     * @param axis a unit vector describing an axis of rotation
     * @return rotated vector
     */
    public Vector3D rotate(double alpha, Vector3D axis) {
    	
    	/*
    	 * uses Rodrigues' formula
    	 * v = cos(alpha) * v + sin(alpha) * (axis X v) + (1-cos(alpha)) * (axis.v) * axis
    	 */
    	
    	
    	if (axis.z < 0)
    	{
    		axis=axis.scaled(-1);
    	}
    	else
    	{
    		if ((axis.z == 0) && (axis.x < 0))
    		{
    			axis.scaled(-1);
    		}
    	}
    	
    	Vector3D p1 = this.scaled(Math.cos(alpha));
    	Vector3D p2 = this.cross(axis).scaled(Math.sin(alpha));
    	Vector3D p3 = axis.scaled(axis.dot(this)).scaled(1-Math.cos(alpha));
    	
    	return p1.plus(p2).plus(p3);
    }

}


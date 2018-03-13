package numerics;

/**
 * <dl>
 * <dt> Purpose: To encapsulate a 3D point.
 * 
 * </dl>
 *
 * @author Philip Cook
 * @version $Id$
 */
public final class Point3D {

   
    public final double x;
    public final double y;
    public final double z;
    
    
    /**
     * Creates a point at the origin (0,0,0).
     */
    public Point3D() {
	x = 0.0;
	y = 0.0;
	z = 0.0;
    }

    /**
     * Creates a point at position (xVal, yVal, zVal).
     */
    public Point3D(double xVal, double yVal, double zVal) {
	x = xVal;
	y = yVal;
	z = zVal;
    }

    /**
     * Creates a new <code>Point3D</code> that is a copy of this one.
     */
    public Point3D(Point3D point) {
	x = point.getX();
	y = point.getY();
	z = point.getZ();
    }


    /**
     * Creates a new <code>Point3D</code> at the position {x,y,z}.
     * @param xyz a double array of size 3 containing the x, y and z
     * positions, in that order
     */
    public Point3D(double[] xyz) {
	x = xyz[0];
	y = xyz[1];
	z = xyz[2];
    }


    /**
     * @return the x position of this point.
     */
    public double getX() {
	return x;
    }

    /**
     * @return the y position of this point.
     */
    public double getY() {
	return y;
    }

    /**
     * @return the z position of this point.
     */
    public double getZ() {
	return z;
    }


    /**
     * @return a point displaced from this one by some vector.
     */
    public Point3D displace(Vector3D v) {
	return new Point3D(x + v.x, y + v.y, z + v.z);
    }


    /**
     * @return the distance between p and this point. Assumes the points are specified in the 
     * same space and units. 
     */
    public double distance(Point3D p) {
	return Math.sqrt( (x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z) );
    }


    /**
     * Transform a point by a 4D Affine matrix.
     *
     * @param aff a 4x4 matrix containing an affine transformation. Non-affine elements (the fourth row) are ignored.
     */
    public Point3D transform(RealMatrix aff) {
        double[] abc = new double[3];

        for (int i = 0; i < 3; i++) {
            abc[i] = aff.entries[i][0] * x + aff.entries[i][1] * y + aff.entries[i][2] * z + aff.entries[i][3];
        }

        return new Point3D(abc);
    }


    /**
     * @return a string containing the details of the point.
     */
    public String toString() {
	return "Point3D: " + "(" + x + ", " + y + ", " + z + ")";
    }

    /**
     * @return true if o is a Point3D that represents the same point, to within a delta of 1E-10
     */
    public boolean equals(Object o) {
	
	if (o instanceof Point3D == false) {
	    return false;
	}
	else if (o == this) {
	    return true;
	}
	else if (o == null) {
	    return false;
	}
	
	Point3D p = (Point3D)o;

        double eps = 1E-10;

	return (Math.abs(x - p.x) < eps && Math.abs(y - p.y) < eps && Math.abs(z - p.z) < eps);
	
        
	
    }	


    public int hashCode() {
	return (int)(17.0 * x + 23.0 * y + 31.0 * z);
    }


}

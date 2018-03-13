package numerics;

import java.util.Random;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Various methods for generating and using rotation matrices.
 * 
 * <dt>Description:
 * 
 * <dd>General manipulation of rotations in 3D.
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id: Rotations.java,v 1.3 2005/08/18 11:12:21 ucacmgh
 *         Exp $
 *  
 */
public class Rotations {

    public static final Vector3D X_AXIS = new Vector3D(1.0, 0.0, 0.0);
    public static final Vector3D Y_AXIS = new Vector3D(0.0, 1.0, 0.0);
    public static final Vector3D Z_AXIS = new Vector3D(0.0, 0.0, 1.0);


    /**
     * Computes the rotation matrix that rotates by angle theta about axis
     * (alpha, beta).
     * 
     * @param alpha
     *            Axis angle of colatitude.
     * 
     * @param beta
     *            Axis angle of longitude
     * 
     * @param theta
     *            Angle of rotation.
     * 
     * @return the rotation matrix.
     */
    public static RealMatrix getRotMat(double alpha, double beta, double theta) {

        double[] r = new double[3];
        r[0] = Math.cos(beta) * Math.sin(alpha);
        r[1] = Math.sin(beta) * Math.sin(alpha);
        r[2] = Math.cos(alpha);

        return getRotMat(r, theta);

    }

    /**
     * Computes the rotation matrix that rotates by angle theta about axis r =
     * (rx, ry, rz).
     * 
     * @param r
     *            The axis of rotation.
     * 
     * @param theta
     *            The angle of rotation.
     * 
     * @return The rotation matrix.
     */
    public static RealMatrix getRotMat(double[] r, double theta) {

        RealMatrix rot = new RealMatrix(3, 3);

        double ct = Math.cos(theta);
        double st = Math.sin(theta);

        double modR = Math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        double rx = r[0];
        double ry = r[1];
        double rz = r[2];

        rot.entries[0][0] = rx * rx * (1.0 - ct) + ct;
        rot.entries[1][0] = ry * rx * (1.0 - ct) + rz * st;
        rot.entries[2][0] = rz * rx * (1.0 - ct) - ry * st;
        rot.entries[0][1] = rx * ry * (1.0 - ct) - rz * st;
        rot.entries[1][1] = ry * ry * (1.0 - ct) + ct;
        rot.entries[2][1] = rz * ry * (1.0 - ct) + rx * st;
        rot.entries[0][2] = rx * rz * (1.0 - ct) + ry * st;
        rot.entries[1][2] = ry * rz * (1.0 - ct) - rx * st;
        rot.entries[2][2] = rz * rz * (1.0 - ct) + ct;

        return rot;
    }

    /**
     * Computes the rotation matrix that rotates by angle theta about axis
     * r = (rx, ry, rz).
     *
     * @param r The axis of rotation.
     * 
     * @param theta The angle of rotation.
     *
     * @return The rotation matrix.
     */
    public static RealMatrix getRotMat(Vector3D r, double theta) {
	
	RealMatrix rot = new RealMatrix(3, 3);
	
	double ct = Math.cos(theta);
	double st = Math.sin(theta);
	
	double modR = Math.sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
	double rx = r.x / modR;
	double ry = r.y / modR;
        double rz = r.z / modR;
	
	rot.entries[0][0] = rx*rx*(1.0 - ct) + ct;
	rot.entries[1][0] = ry*rx*(1.0 - ct) + rz*st;
	rot.entries[2][0] = rz*rx*(1.0 - ct) - ry*st;
	rot.entries[0][1] = rx*ry*(1.0 - ct) - rz*st;
	rot.entries[1][1] = ry*ry*(1.0 - ct) + ct;
	rot.entries[2][1] = rz*ry*(1.0 - ct) + rx*st;
	rot.entries[0][2] = rx*rz*(1.0 - ct) + ry*st;
	rot.entries[1][2] = ry*rz*(1.0 - ct) - rx*st;
	rot.entries[2][2] = rz*rz*(1.0 - ct) + ct;
	
	return rot;
    }


    /**
     * Generates a random rotation matrix with random seed taken from the clock
     * time.
     * 
     * @return The rotation matrix.
     */
    public static RealMatrix randomRotMat() {

        Random r = new Random();
        return randomRotMat(r);
    }

    /**
     * Generates a seeded random rotation matrix.
     * 
     * @param seed
     *            The random number generator seed.
     * 
     * @return The rotation matrix.
     */
    public static RealMatrix randomRotMat(int seed) {

        Random r = new Random(seed);
        return randomRotMat(r);
    }

    /**
     * Generates a random rotation matrix drawing using the random number
     * generator provided.
     * 
     * @param r
     *            The random number generator.
     * 
     * @return The rotation matrix.
     */
    public static RealMatrix randomRotMat(Random r) {

        //Select a random rotation axis.
        double[] axis = SphericalPoints.getRandomPoint(r);

        //The rotation angle is uniformly distributed on [0, 2pi).
        double theta = 2.0 * Math.PI * r.nextDouble();

        RealMatrix rot = getRotMat(axis, theta);
        return rot;
    }

    /**
     * Applies a 3D transformation matrix to a point.
     * 
     * @param mat
     *            The rotation matrix.
     * 
     * @param point
     *            The point to be rotated.
     * 
     * @return The rotated point.
     */
    public static double[] transformPoint(RealMatrix mat, double[] point) {
        double[] newPoint = new double[3];

        newPoint[0] = mat.entries[0][0] * point[0] + mat.entries[0][1] * point[1]
                + mat.entries[0][2] * point[2];
        newPoint[1] = mat.entries[1][0] * point[0] + mat.entries[1][1] * point[1]
                + mat.entries[1][2] * point[2];
        newPoint[2] = mat.entries[2][0] * point[0] + mat.entries[2][1] * point[1]
                + mat.entries[2][2] * point[2];

        return newPoint;
    }


    // these methods added by PAC 


    /**
     * Rotate a tensor so that it's e1 is aligned with vector.
     *
     * @param original the tensor to be rotated.
     * @param newE1 a unit vector. The returned tensor's e1 will lie along this vector.
     *
     */
    public static DT rotateTensor(DT original, Vector3D newE1) {
	
	double[][] eSys = original.sortedEigenSystem();

	Vector3D originalE1 = new Vector3D(eSys[1][0], eSys[2][0], eSys[3][0]);

	RealMatrix trans = getRotMat(originalE1, newE1);

	return original.transform(trans);

    }

    /**
     * Rotate original by angle theta about axis.
     *  
     * @param original the vector to be rotated.
     * @param axis the axis of rotation. This should be a unit vector.
     * @param theta the angle of rotation, in radians.
     */
    public static Vector3D rotateVector(Vector3D original, Vector3D axis, double theta) {
	return rotateVector(original, getRotMat(axis, theta));
    }


    /**
     * Apply a rotation matric to a vector. Note: replace with Vector3D.transform?
     *
     */
    public static Vector3D rotateVector(Vector3D original, RealMatrix rot) {
	
	// For speed, do the multiplication manually
	double[] transVector = new double[3];
	
	transVector[0] = 
	    rot.entries[0][0] * original.x + 
	    rot.entries[0][1] * original.y + 
	    rot.entries[0][2] * original.z;
	
	transVector[1] = 
	    rot.entries[1][0] * original.x + 
	    rot.entries[1][1] * original.y + 
	    rot.entries[1][2] * original.z;  

	transVector[2] = 
	    rot.entries[2][0] * original.x + 
	    rot.entries[2][1] * original.y + 
	    rot.entries[2][2] * original.z;
	
	return new Vector3D(transVector[0], transVector[1], transVector[2]);
	


    }



    /**
     * Get rotation matrix of one vector onto another vector.
     *
     * @param vec the vector to be rotated.
     * @param targetAxis the target unit vector. The returned matrix will transform <code>vec</code> 
     * onto <code>targetAxis</code>.
     *
     */
    public static RealMatrix getRotMat(Vector3D vec, Vector3D targetAxis) {
	
	double orPhi, targetPhi;
	double orTheta, targetTheta;

	RealMatrix r1, r2, r3;


	// rotate targetAxis into XZ plane
	if (targetAxis.x == 0.0) {
	    targetPhi = targetAxis.y > 0.0 ? Math.PI / 2.0 : -1.0 * Math.PI / 2.0;
	}
	else {
	    if (targetAxis.x < 0.0) {
		targetPhi = Math.PI + Math.atan(targetAxis.y / targetAxis.x);
	    }
	    else {
		targetPhi = Math.atan(targetAxis.y / targetAxis.x);
	    }
	}

	r1 = Rotations.getRotMat(Z_AXIS, -1.0 * targetPhi);
	
	targetAxis = Rotations.rotateVector(targetAxis, r1);

	// same deal for vec

	if (vec.x == 0.0) {
	    orPhi = vec.y > 0.0 ? Math.PI / 2.0 : -1.0 * Math.PI / 2.0;
	}
	else {
	    if (vec.x < 0.0) {
		orPhi = Math.PI + Math.atan(vec.y / vec.x);
	    }
	    else {
		orPhi = Math.atan(vec.y / vec.x);
	    }
	}
	
	r2 = Rotations.getRotMat(Z_AXIS, -1.0 * orPhi);

	
	// find targetTheta, angle between targetAxis and the z axis
	targetTheta = Math.acos(Z_AXIS.dot(targetAxis));

	// find orTheta, angle between vec and the z axis
	orTheta = Math.acos( Z_AXIS.dot(vec) );
	
	r3 = Rotations.getRotMat(Y_AXIS, targetTheta - orTheta);
	
	RealMatrix invR1 = r1.inverse();

	return invR1.product(r3).product(r2);

    }




}

package misc;

import java.util.Random;
import java.text.MessageFormat;
import tools.*;
import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Library of functions concerning points on the sphere.
 * 
 * <dt>Description:
 * 
 * <dd>Provides methods to generate and manipulate spherical point sets.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class SphericalPoints {

    public static String pointSetDirectory = "/PointSets/";


    /**
     * These points are the vertices of an icosahedron. This point set
     * came from an electrostatic minimization with 6 pairs. The edge
     * length corresponds to that expected for an icosahedron.
     */
    public static final double[][] icosahedron = { { 0.200009, 0.601164, 0.773691 },
            { -0.200009, -0.601164, -0.773691 }, { -0.665299, 0.004440, 0.746564 },
            { 0.665299, -0.004440, -0.746564 }, { -0.567356, 0.818696, 0.088562 },
            { 0.567356, -0.818696, -0.088562 }, { 0.961832, 0.092438, 0.257555 },
            { -0.961832, -0.092438, -0.257555 }, { 0.279792, -0.444412, 0.851008 },
            { -0.279792, 0.444412, -0.851008 }, { 0.438265, 0.873082, -0.213663 },
            { -0.438265, -0.873082, 0.213663 } };


    /**
     * Get a single random point on the unit sphere, drawn from a uniform
     * distribution over its surface.
     * 
     * @param r
     *            A random number generator.
     * 
     * @return The point.
     */
    public static double[] getRandomPoint(Random r) {

        double[] pt = new double[3];

        pt[2] = 2.0 * r.nextDouble() - 1.0;
        double phi = r.nextDouble() * 2.0 * Math.PI;
        double theta = Math.acos(pt[2]);
        pt[0] = Math.sin(theta) * Math.cos(phi);
        pt[1] = Math.sin(theta) * Math.sin(phi);

        return pt;
    }


    /**
     * Reads in a returns the electrostatic pairs point set with the specified
     * number of pairs.
     * 
     * @param n
     *            The number of pairs.
     * 
     * @return The list of points.
     */
    public static double[][] getElecPointSet(int n) {
        Object[] o = new Object[1];
        o[0] = new Integer(n);

        String filename = MessageFormat.format("Elec{0,number,000}.txt", o);

	java.net.URL url = SphericalPoints.class.getResource(pointSetDirectory + filename);
        if(url == null) {
            throw new LoggedException("Could not find file: " + pointSetDirectory + filename);
        }

        return readFromFile(url.getFile());
    }


    /**
     * Returns a list of points on the sphere which is the vertices of
     * a number of icosahedra each having undergone a separate random
     * rotation.
     * 
     * @param density
     *            The number of randomly rotated icosahedra
     * 
     * @param seedOffset
     *            An offset of the seeds used in the random number generator
     *            used to get the random rotations.
     *
     * @return A list of sample points for random sampling in finding
     *         peak directions.
     */
    public static double[][] getIcosahedralPointSet(int density, int seedOffset) {

        double[][] points = new double[SphericalPoints.icosahedron.length * density / 2][SphericalPoints.icosahedron[0].length];

        double[] point = new double[3];

        for (int i = 0; i < density; i++) {
            RealMatrix randomRot = Rotations.randomRotMat(i + seedOffset);

            //We only take one point from every equal and opposite
            //pair in the icosahedron.
            for (int j = 0; j < SphericalPoints.icosahedron.length; j += 2) {
                point[0] = SphericalPoints.icosahedron[j][0];
                point[1] = SphericalPoints.icosahedron[j][1];
                point[2] = SphericalPoints.icosahedron[j][2];

                double[] rotPoint = Rotations.transformPoint(randomRot, point);

                points[i * SphericalPoints.icosahedron.length / 2 + j / 2][0] = rotPoint[0];
                points[i * SphericalPoints.icosahedron.length / 2 + j / 2][1] = rotPoint[1];
                points[i * SphericalPoints.icosahedron.length / 2 + j / 2][2] = rotPoint[2];
            }
        }

        return points;
    }


    /**
     * Rotates a list of points.
     * 
     * @param pointSet This list of points
     *
     * @param rot The rotation matrix
     *
     * @return The list of rotated points.
     */
    public static double[][] rotatePointSet(double[][] pointSet, RealMatrix rot) {

        double[][] rotPoints = new double[pointSet.length][pointSet[0].length];

        for (int i = 0; i < pointSet.length; i++) {
            rotPoints[i] = Rotations.transformPoint(rot, pointSet[i]);
        }

        return rotPoints;
    }


    /**
     * Loads a spherical point set from a file converted by PointSetReformat
     * from the format output by the ElecPairsOnSphere program.
     * 
     * @param filename
     *            The name of the file to read the directions in from.
     * 
     * @return The array of directions.
     */
    public static double[][] readFromFile(String filename) {

        FileInput f = new FileInput(filename);
        int numPairs = f.readInteger();

        // Allocate the array.
        double[][] pointSet = new double[numPairs][3];
        for (int i = 0; i < numPairs; i++) {
            for (int j = 0; j < 3; j++) {
                pointSet[i][j] = f.readDouble();
            }
        }

        return pointSet;
    }

    /**
     * Converts an xyz vector to spherical polar coordinates: r, theta, phi.
     * 
     * @param xyz
     *            Cartesian coordinates {x, y, z}.
     * 
     * @return {r, theta, phi}.
     */
    public static double[] getSphPolars(double[] xyz) {
        double r = Math.sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
        double theta = (r > 0.0) ? Math.acos(xyz[2] / r) : 0.0;
        double phi = 0.0;

        //Need to be careful when theta is very small or close to PI, when
        //phi is not well defined.
        if (r > 0.0 && 1.0 - Math.abs(xyz[2] / r) > 1.0E-10) {
            //Have to be wary of numerical problems here.
            double cphi = xyz[0] / (r * Math.sin(theta));
            cphi = (cphi < -1.0) ? -1.0 : cphi;
            cphi = (cphi > 1.0) ? 1.0 : cphi;
            phi = Math.acos(cphi);

            double sphi = xyz[1] / (r * Math.sin(theta));
            sphi = (sphi < -1.0) ? -1.0 : sphi;
            sphi = (sphi > 1.0) ? 1.0 : sphi;

            if (Math.asin(sphi) < 0.0) {
                phi = 2.0 * Math.PI - phi;
            }
        }

        double[] results = new double[3];
        results[0] = r;
        results[1] = theta;
        results[2] = phi;

        return results;
    }

    /**
     * Returns an array of unit vectors obtained by sampling z and phi
     * evenly at the resolution specified. Archimedian point sets.
     *
     * @param sampRes The sampling resolution
     *
     * @return The pointset.
     */
    public static double[][] getZPhiXs(int sampRes) {
        double[] zs = new double[sampRes];
        double[] phis = new double[sampRes];
        double[][] xvecs = new double[sampRes * sampRes][3];

        //Sample the theta and phi ranges evenly.
        for (int i = 0; i < sampRes; i++) {
            zs[i] = ((double) (2 * i + 1)) / ((double) sampRes) - 1.0;
            phis[i] = ((double) i + 0.5) * 2.0 * Math.PI / ((double) sampRes);
        }

        int x = 0;
        for (int i = 0; i < sampRes; i++) {
            for (int j = 0; j < sampRes; j++) {
                double sinTheta = Math.sqrt(1.0 - zs[i] * zs[i]);
                xvecs[x][0] = sinTheta * Math.cos(phis[j]);
                xvecs[x][1] = sinTheta * Math.sin(phis[j]);
                xvecs[x][2] = zs[i];
                x++;
            }
        }

        return xvecs;
    }

}

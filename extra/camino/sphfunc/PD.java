package sphfunc;

import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Stores a direction vector with proportional magnitude.
 * 
 * <dt>Description:
 * 
 * <dd>This class encapsulates a principal direction of a spherical function.
 * It is used in PDList which contains a list of PDs of a given function.
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id$
 *  
 */
public class PD {

    // Store the PD in spherical polar coordinates.

    /**
     * The angle of colatitude.
     */
    private double pdTheta;

    /**
     * The angle of longitude.
     */
    private double pdPhi;

    /**
     * The value of the function at this PD.
     */
    private double pdProp;

    /**
     * Constructs a PD object from spherical polar coordinates.
     * 
     * @param theta
     *            The angle of colatitude.
     * 
     * @param phi
     *            The angle of longitude.
     * 
     * @param prop
     *            The value of the function at theta, phi.
     */
    public PD(double theta, double phi, double prop) {
        pdTheta = theta;
        pdPhi = phi;
        pdProp = prop;
    }

    /**
     * Constructs a PD object from Cartesian coordinates.
     * 
     * @param xyz
     *            Cartesian coordinates [x, y, z] must have mod 1 (not checked).
     * 
     * @param prop
     *            The value of the function at (x, y, z).
     */
    public PD(double[] xyz, double prop) {
        double[] sps = SphericalPoints.getSphPolars(xyz);

        pdTheta = sps[1];
        pdPhi = sps[2];
        pdProp = prop;
    }

    /**
     * Returns the angle of colatitude.
     * 
     * @return theta
     */
    public double getTheta() {
        return pdTheta;
    }

    /**
     * Returns the angle of longitude.
     * 
     * @return phi
     */
    public double getPhi() {
        return pdPhi;
    }

    /**
     * Returns the value of the function at the PD.
     * 
     * @return f(theta, phi)
     */
    public double getProp() {
        return pdProp;
    }

    /**
     * Returns the x-component of the PD.
     * 
     * @return x
     */
    public double getPDX() {
        return Math.cos(pdPhi) * Math.sin(pdTheta);
    }

    /**
     * Returns the y-component of the PD.
     * 
     * @return y
     */
    public double getPDY() {
        return Math.sin(pdPhi) * Math.sin(pdTheta);
    }

    /**
     * Returns the z-component of the PD.
     * 
     * @return z
     */
    public double getPDZ() {
        return Math.cos(pdTheta);
    }

    public String toString() {
        String s = "Theta: " + pdTheta + ", Phi: " + pdPhi + ", Prop.: " + pdProp;
        return s;
    }

}

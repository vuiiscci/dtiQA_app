package data;

import misc.*;
import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Implements standard test functions.
 * 
 * <dt>Description:
 * 
 * <dd>Provides methods that return the standard test functions as ModelPDF
 * objects, as well as variations on those functions.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: StandardTestFunctions.java,v 1.4 2005/08/18 10:59:37 ucacmgh
 *          Exp $
 *  
 */
public class StandardTestFunctions {

    /**
     * Trace of standard diffusion tensors.
     */
    protected static double traceD = 21.0E-10;

    /**
     * Largest eigenvalue in prolate diffusion tensors.
     */
    protected static double lambda1 = 17.0E-10;

    /**
     * Mixing parameters for p3.
     */
    protected static double[] mix3 = { 0.5, 0.5 };

    /**
     * Mixing parameters for p4.
     */
    protected static double[] mix4 = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

    /**
     * Scaling factor for DTs.
     */
    protected static double scaleFactor = 1.0;

    /**
     * Angle through which to rotate dt2 about the z-axis.
     */
    protected static double dt2rotangle = 0.0;

    /**
     * Linear transformation applied to DTs. Usually a rotation.
     */
    protected static RealMatrix trans = RealMatrix.identity(3);

    /**
     * Isotropic diffusion tensor.
     */
    protected static DT dt0;

    /**
     * Prolate diffusion tensor aligned with x-axis.
     */
    protected static DT dt1;

    /**
     * Prolate diffusion tensor aligned with y-axis.
     */
    protected static DT dt2;

    /**
     * Prolate diffusion tensor aligned with z-axis.
     */
    protected static DT dt3;

    /**
     * Oblate diffusion tensor in the xy-plane.
     */
    protected static DT dt4;

    /**
     * Keeps track of whether anything has changed so that the diffusion tensors
     * need recomputing.
     */
    private static boolean parametersChanged = true;

    /**
     * Returns the indexed test function rotated by the indexed rotation, with
     * specified mixing parameters and scaling factor on the diffusion tensors.
     * 
     * @param index
     *            Test function index.
     * 
     * @return The test function as a ModelPDF object.
     */
    public static ModelPDF getFunction(int index) {

        // Create the diffusion tensors with the current parameter
        // settings.
        if (parametersChanged) {
            constructDTs();
        }

        // Create the test function and return it.
        if (index == 0 || index == 1 || index == 2) {
            DT[] comps = new DT[1];
            if (index == 0) {
                comps[0] = dt0;
            }
            else if (index == 1) {
                comps[0] = dt1;
            }
            else {
                comps[0] = dt4;
            }

            double[] mix = { 1.0 };

            return new GaussianMixture(comps, mix);
        }
        else if (index == 3) {
            DT[] comps = new DT[2];
            comps[0] = dt1;
            comps[1] = dt2;
            return new GaussianMixture(comps, mix3);
        }
        else if (index == 4) {
            DT[] comps = new DT[3];
            comps[0] = dt1;
            comps[1] = dt2;
            comps[2] = dt3;
            return new GaussianMixture(comps, mix4);
        }
        else {
            throw new RuntimeException("Unrecognized test function (index=" + index
                    + ") requested.");
        }

    }

    /**
     * Resets all the functions' components to the defaults.
     */
    public static void reset() {
        traceD = 21.0E-10;
        lambda1 = 17.0E-10;
        double[] newMix3 = { 0.5, 0.5 };
        mix3 = newMix3;
        double[] newMix4 = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
        mix4 = newMix4;
        scaleFactor = 1.0;
        dt2rotangle = 0.0;
        trans = RealMatrix.identity(3);
    }

    /**
     * Sets the scaling factor.
     * 
     * @param scale
     *            The scaling factor.
     */
    public static void setScale(double scale) {
        scaleFactor = scale;
        parametersChanged = true;
    }

    /**
     * Sets the dt2 rotation angle.
     * 
     * @param scale
     *            The new rotation angle.
     */
    public static void setDT2RotationAngle(double angle) {
        dt2rotangle = angle;
        parametersChanged = true;
    }

    /**
     * Sets the trace of the components diffusion tensors.
     * 
     * @param newTrace
     *            The new trace.
     */
    public static void setTraceD(double newTrace) {
        traceD = newTrace;
        parametersChanged = true;
    }

    /**
     * Sets the largest eigenvalue of the prolate diffusion tensors. This can be
     * used to change the anisotropy of the components.
     * 
     * @param newLambda1
     *            The new largest eigenvalue.
     */
    public static void setLambda1(double newLambda1) {
        lambda1 = newLambda1;
        parametersChanged = true;
    }

    /**
     * Sets the mixing parameters for function 3.
     * 
     * @param m1
     *            The mixing parameter for the first component.
     */
    public static void setMix3(double m1) {
        mix3[0] = m1;
        mix3[1] = 1.0 - m1;
    }

    /**
     * Sets the mixing parameters for function 4.
     * 
     * @param m1
     *            The mixing parameter for the first component.
     * 
     * @param m2
     *            The mixing parameter for the second component.
     */
    public static void setMix4(double m1, double m2) {
        mix3[0] = m1;
        mix3[1] = m2;
        mix3[2] = 1.0 - m1 - m2;
    }

    /**
     * Sets the transformation.
     * 
     * @param transformation
     *            The new transformation.
     */
    public static void setTransformation(RealMatrix transformation) {
        trans = transformation;
        parametersChanged = true;
    }

    /**
     * Constructs the DTs the make the test functions using the current
     * transformations, scaling factors, etc.
     */
    protected static void constructDTs() {
        double lambda0 = traceD / 3.0;
        dt0 = new DT(lambda0, 0.0, 0.0, lambda0, 0.0, lambda0);

        double lambda2 = (traceD - lambda1) / 2.0;

        dt1 = new DT(lambda1, 0.0, 0.0, lambda2, 0.0, lambda2);
        dt2 = new DT(lambda2, 0.0, 0.0, lambda1, 0.0, lambda2);
        dt3 = new DT(lambda2, 0.0, 0.0, lambda2, 0.0, lambda1);

        double lambda3 = (lambda1 + lambda2) / 2.0;
        dt4 = new DT(lambda3, 0.0, 0.0, lambda3, 0.0, lambda2);

        // Rotate dt2 about the z-axis if specified.
        if (dt2rotangle != 0.0) {
            double[] zaxis = { 0.0, 0.0, 1.0 };
            RealMatrix dt2trans = Rotations.getRotMat(zaxis, dt2rotangle);
            dt2 = dt2.transform(dt2trans);
        }

        dt0 = dt0.transform(trans).scale(scaleFactor);
        dt1 = dt1.transform(trans).scale(scaleFactor);
        dt2 = dt2.transform(trans).scale(scaleFactor);
        dt3 = dt3.transform(trans).scale(scaleFactor);
        dt4 = dt4.transform(trans).scale(scaleFactor);

    }

    
       
    /** returns the rotation matrix transform
     * 
     * @return transform matrix trans
     */
    public static RealMatrix getRotationMatrix(){
        return trans;
    }
}

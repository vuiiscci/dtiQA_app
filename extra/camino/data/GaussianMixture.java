package data;

import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Simple model for the diffusion displacement density function .
 * 
 * <dt>Description:
 * 
 * <dd>Models the displacement by a mixture of Gaussian densities.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class GaussianMixture implements ModelPDF {

    /**
     * An array of diffusion tensors - one for each mixture model component.
     */
    protected DT[] dt;

    /**
     * Array of inverse diffusion tensors.
     */
    protected DT[] invDT;

    /**
     * Array of square rooted diffusion tensor determinants.
     */
    protected double[] sqrtDeterm;

    /**
     * Array of mixing parameters.
     */
    protected double[] mix;


    /**
     * Constructor initializes the array of diffusion tensors and mixing
     * parameters as well as the diffusion time.
     * 
     * @param components
     *            The diffusion tensors for each component.
     * 
     * @param mixingParameters
     *            The mixing parameters for each component.
     */
    public GaussianMixture(DT[] components, double[] mixingParameters) {

        // Copy the contents of the array parameters.
        dt = new DT[components.length];
        mix = new double[mixingParameters.length];
        for (int i = 0; i < dt.length; i++) {
            dt[i] = (DT) (components[i].clone());
            mix[i] = mixingParameters[i];
        }
        

    }


    public double at(double[] x, double tau) {

        // Only compute the inverses and determinants when required.
        if (invDT == null) {
            computeInverses();
        }

        double val = 0;
        for (int i = 0; i < invDT.length; i++) {
            val += mix[i] * Math.exp(-invDT[i].contractBy(x) / (4 * tau)) / sqrtDeterm[i];
        }

        double z = Math.sqrt(64.0 * Math.PI * Math.PI * Math.PI * tau * tau * tau);
        return val / z;
    }


    public double ftAt(double[] q, double tau) {
        double val = 0;
        for (int i = 0; i < dt.length; i++) {
            val += mix[i] * Math.exp(-tau * dt[i].contractBy(q));
        }

        return val;
    }

    public double ftAtB_Vec(double[] g, double b) {
        double val = 0;
        for (int i = 0; i < dt.length; i++) {
            val += mix[i] * Math.exp(-b * dt[i].contractBy(g));
        }

        return val;
    }


    /**
     * Returns a list of principal directions of the test function.
     * 
     * @return Array of principal directions.
     */
    public double[][] getPDs() {
        double[][] dirs = new double[dt.length][3];
        for (int i = 0; i < dt.length; i++) {
            double[] dir = dt[i].getPD();
            for (int j = 0; j < 3; j++) {
                dirs[i][j] = dir[j];
            }
        }

        return dirs;
    }


    /**
     * 
     * @return the diffusion tensor(s) of this PDF.
     */
    public DT[] getTensors() {

        // since tensors are immutable, we don't have to care about copying them
        // we just need to copy the mutable array
        DT[] copy = new DT[dt.length];

        for (int i = 0; i < dt.length; i++) {
            copy[i] = dt[i];
        }

        return copy;
    }


    /**
     * Computes the inverses and determinants of all the component diffusion
     * tensors.
     */
    private void computeInverses() {
        invDT = new DT[dt.length];
        sqrtDeterm = new double[dt.length];
        for (int i = 0; i < dt.length; i++) {
            invDT[i] = dt[i].inverse();
            sqrtDeterm[i] = Math.sqrt(dt[i].determinant());
        }
    }


    public String toString() {
        String s = "p = \n";
        for (int i = 0; i < dt.length; i++) {
            s += mix[i] + " X " + dt[i].toString();
            if (i < dt.length - 1) {
                s += " +\n";
            }
        }

        return s;
    }
}

package data;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General interface for models for the diffusion displacement density
 * function.
 * 
 * <dt>Description:
 * 
 * <dd>Provides methods to access the value of a test function and its Fourier
 * transform and a specific time and displacement/wavenumber.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public interface ModelPDF {

    /**
     * Returns the value of the function at the specified point and time.
     * 
     * @param x
     *            The point to sample at.
     * 
     * @param tau
     *            The diffusion time.
     * 
     * @return The value of the function at x.
     */
    public double at(double[] x, double tau);


    /**
     * Returns the value of the FT of the function at the specified wavenumber.
     * 
     * @param q
     *            The wavenumber to sample at.
     * 
     * @param tau
     *            The diffusion time.
     * 
     * @return The value of the FT at q.
     */
    public double ftAt(double[] q, double tau);


    /**
     * Returns the value of the FT of the function at the specified gradient direction and b-value.
     * 
     * @param g
     *            The gradient direction, should be a unit vector on the sphere.
     * 
     * @param b
     *            The b-value.
     * 
     * @return The value of the FT at q.
     */
    public double ftAtB_Vec(double[] g, double b);


    /**
     * Returns a list of principal directions of the test function.
     * 
     * @return Array of principal directions. The first index is the principal
     *         direction number. Each principal direction is stored in Cartesian
     *         coordinates: {x, y, z}.
     */
    public double[][] getPDs();

}

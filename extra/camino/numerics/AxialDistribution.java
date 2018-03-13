package numerics;

/**
 * Interface for classes that define an axial PDF on the sphere, giving p(x) == p(-x) for any unit vector
 * x.
 */
public interface AxialDistribution {



    /**
     * Gets the next Vector from the distribution.
     */
    public Vector3D nextVector();


    /**
     * Gets the probability density for a given unit vector.
     *
     * @param v a unit vector.
     *
     */
    public double pdf(Vector3D v);



}

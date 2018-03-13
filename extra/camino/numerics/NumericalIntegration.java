package numerics;

/**
 * Camino fibre reconstruction and tracking toolkit
 * 
 * NumericalIntegration (numerics)
 * 
 * Abstract class providing numerical integration on a sphere via the use of
 * Gaussian Quadrature. Here the integral is approximated by a weighted sum
 * 
 * \int_{S^2}f(\theta, \phi)d\mu(\theta, \phi) = \int_{S^2} f(\cos\theta,
 * \cos\phi) d\mu(\theta, \phi) \simeq
 * \sum_{i=0}^{N}\sum_{j=0}^{N}W_{ij}f(\cos\theta_i, \cos\phi\j)
 * 
 * where W_{ij}=w_{i}v_{j} is the weights matrix
 * 
 * ideal quadrature points correspond to all pairs (\theta_i, \phi_j) where
 * \theta_i are the roots of the Nth chebyshev polynomial, and \phi_j are the
 * roots of the Nth legendre polynomial respectively, but other point sets may
 * be used and weights calculated
 * 
 * currently, spiral points on the sphere are used, as these represent an
 * assymtotically equal spacing. In this case the weights matrix is diagonal
 * 
 * @author ucacmgh
 *  
 */
public abstract class NumericalIntegration {

    /** tolerence on conversion of integral */
    private double TINY = 1E-14;

    /**
     * flag that specifies whether the weights matrix is diagonal or not.
     */
    private final boolean diagonalWeights;

    /** the maximum order of integral sum before termination */
    private int maxN = 10000;

    /** the min number of terms to evaluate in the integration sum */
    private int minN = 2;

    public NumericalIntegration(boolean diagonalWeights) {
        this.diagonalWeights = diagonalWeights;
    }

    public NumericalIntegration(boolean diagonalWeights, int minN, int maxN) {
        this(diagonalWeights);

        setMinN(minN);
        setMaxN(maxN);

    }

    public NumericalIntegration(boolean diagonalWeights, int minN, int maxN,
            double threshold) {
        this(diagonalWeights, minN, maxN);

        setThreshold(threshold);
    }

    /**
     * abstract method to evaluate the function to be integrated
     * 
     * @param u
     *            one coordinate
     * @param v
     *            another
     * 
     * @return some non-pathologial function on the sphere at given coords f:
     *         S^2 -> \mathbb{R}
     */
    public abstract double f(double u, double v);

    public double integrate() {

        double last, current;

        last = integrate(minN);
        current = integrate(minN + 1);

        for (int n = minN + 2; n <= maxN; n++) {
            if (Math.abs(current - last) <= TINY) {
                return current;
            }
            last = current;
            current = integrate(n + 1);
        }

        System.err.println("WARNING: integral conversion not reached in " + maxN
                + " terms with threshold " + TINY);

        return current;

    }

    /**
     * performs an integration with a given index
     * 
     * @param n
     *            the index of the integration
     * 
     * @return the integral of the function f(u,v) over the sphere with points
     *         and weights defined by n
     */
    public double integrate(int n) {

        if (n < 2) {
            throw new RuntimeException(
                    "at least two points on the sphere needed for numerical integration");
        }

        if (diagonalWeights) {
            return diagIntegrate(n);
        }
        else {
            return nondiagIntegrate(n);
        }

    }

    /**
     * performs cubature sum with diagonal weights matrix. i.e. off-diagonal
     * terms in W_{ij} are not included
     * 
     * @param n
     * @return \sum_{i=0}^n W_{ii}f(\cos\theta_i, \cos\phi_i)
     */
    private double diagIntegrate(int n) {

        double[][] point = getSpiralPoints(n);
        double[] W = getDiagWeights(point, n);

        double integral = 0.0;

        for (int i = 0; i < n; i++) {
            integral += W[i] * f(point[i][0], point[i][1]);
        }

        return integral;
    }

    /**
     * returns the array of points on the archimidean spiral around the sphere
     * from pole to pole sampled n times on the way down
     * 
     * @param n
     *            number of points
     * @return array of pairs of points on the spiral
     */
    private double[][] getSpiralPoints(int n) {

        double[][] point = new double[n][2];

        // get the (cos) thetas...
        for (int k = 1; k <= n; k++) {
            point[k - 1][0] = -1.0 + (2.0 * (k - 1) / (double) (n - 1));
        }

        // get the phis (they need the cos thetas)
        point[0][1] = 0.0;
        for (int k = 1; k < n - 1; k++) {
            point[k][1] = point[k - 1][1] + (3.6 / Math.sqrt(n))
                    / (Math.sqrt(1.0 - point[k][0] * point[k][0]));
        }
        point[n - 1][1] = 0.0;

        // take cosines of phis
        for (int k = 0; k < n; k++) {
            point[k][1] = Math.cos(point[k][1]);
        }

        return point;
    }

    /**
     * delivers the array of weights for the set of points given
     * 
     * @param point
     *            array of points on sphere
     * @param n
     *            index of weights (number of points)
     * @return array of weights (diagonal of weights matrix)
     */
    private double[] getDiagWeights(double[][] point, int n) {

        double[] W = new double[n];

        for (int ii = 0; ii < n; ii++) {
            W[ii] = getCubatureWeight(point[ii], n);

        }

        return W;
    }

    /**
     * calculates the weight of the point given in the spherical cubature regime
     * of index n
     * 
     * @param point
     *            the point on the sphere (cos theta, cos phi)
     * @param n
     * @return the weight of the cubtaure sum
     */
    private double getCubatureWeight(double[] point, int n) {

        /*
         * double w_i; double v_j;
         * 
         * 
         * 
         * double P=SphericalHarmonics.plgndr(n+1, 0, point[0]); // nth legendre
         * poly at point[0] double Tn=Math.cos((n)*Math.acos(point[1])); // nth
         * checbyshev poly at point[1] double
         * Tnp1=Math.cos((n+1)*Math.acos(point[1])); // (n+1)th chebyshev poly
         * at point[1] double TnPrimed=
         * n*Math.sqrt((1.0-Tn*Tn)/(1.0-point[1]*point[1])); // deriv of nth
         * Chebyshev poly at point[1]
         * 
         * w_i=-(2*(1.0-point[0]*point[0]))/((n+1)*P*P); // legendre-gauss
         * weight v_j=Math.PI/(TnPrimed*Tnp1); // Chebshev-legendre weight
         *  // the singularities cause odd behaviour in TnPrimed as // x -> +/-
         * 1 so we need to catch the limitting cases
         * if(Math.abs(point[1])>=1.0){ v_j=0.0; }
         * 
         * 
         * //System.out.println("P="+P+" Tn="+Tn+" Tnp1="+Tnp1+"
         * TnPrimed="+TnPrimed);
         * 
         * return w_i*v_j;
         */

        return (4.0 * Math.PI) / n;
    }

    /**
     * performs the cubature summation over the full weights matrix for given
     * index of weights and points
     * 
     * @deprecated 
     * @param n
     *            index of weights and points
     * @return TODO: returns area of sphere, does not integrate.
     */
    private double nondiagIntegrate(int n) {

        return 4.0 * Math.PI;
    }

    /**
     * gets the maximum number of terms in the integration sum (default 10000)
     * 
     * @return Returns the maxN.
     */
    public int getMaxN() {
        return maxN;
    }

    /**
     * sets the max number of terms in the integration sum (default 10000)
     * 
     * @param maxN
     *            The maxN to set.
     */
    public void setMaxN(int maxN) {
        if (maxN < minN + 2) {
            throw new RuntimeException(
                    "integration maximum must be at least 2 greater than minimum. min="
                            + minN + " attempt to set max to " + maxN);
        }
        this.maxN = maxN;
    }

    /**
     * @return Returns the conversion threshold
     */
    public double getThreshold() {
        return TINY;
    }

    /**
     * @param tiny
     *            The new conversion threshold
     */
    public void setThreshold(double tiny) {
        TINY = tiny;
    }

    /**
     * @return Returns the minN.
     */
    public int getMinN() {
        return minN;
    }

    /**
     * @param minN
     *            The minN to set.
     */
    public void setMinN(int minN) {
        if (maxN < minN + 2) {
            throw new RuntimeException(
                    "integration maximum must be at least 2 greater than minimum. max="
                            + maxN + ". attempt to set min to " + minN);
        }
        this.minN = minN;
    }
}
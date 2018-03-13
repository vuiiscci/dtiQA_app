package sphfunc;

import java.util.logging.Logger;

import numerics.*;
import optimizers.*;
import tools.*;
import misc.*;

/**
 * <dl>
 *
 * <dt>Purpose:
 *
 * <dd> General spherical function representation.
 * 
 * <dt>Description:
 *
 * <dd>Abstract class for general spherical functions and defines methods
 * required by all spherical functions.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 *
 * $Id$
 *  
 */
public abstract class SphericalFunction {

    /** logging object */
    Logger logger = Logger.getLogger(this.getClass().getName());
    
    
    
    //These instance variables are parameters of the getPDs method that
    //searches for maxima in the spherical function.

    /**
     * The search radius in the random sampling principal direction
     * extraction algorithm. The default value is empirically good for
     * density 1000.
     */
    private static double searchRadius = 0.4;

    /**
     * Parameter scale used in the Powell routine for finding principal
     * directions. If the optimization proves prone to failure, try reducing
     * this value. It reduces the step sizes in Powell to reduce the chance of
     * missing a sharp peak.
     */
    public static double SFPD_PARAMSCALE = 0.1;

    //Threshold on dot product between PDs to classify them as the same.
    protected double DPTHRESH = 0.001;

    /**
     * Convergence threshold in Powell's method for the direction finding
     * routines.
     */
    protected double MINCONV = 1.0E-12;

    
    /**
     *  Convergence threshold for great circle integrals
     */
    private final double GREAT_CIRCLE_INTEGRAL_THRESHOLD = 1e-5;
    
    /**
     * max number of points to use in great circle integration
     */
    private final int MAX_GREAT_CIRCLE_POINTS = 100;
    
    /**
     * min number of points to use in great circle integral
     */
    private final int MIN_GREAT_CIRCLE_POINTS = 20;
    
    /** 
     * step increment in great circle integral
     */
    private final int GREAT_CIRCLE_POINTS_INC = 5;
    
    /**
     * Return the radius of the function at a specified latitude and longitude.
     * Default is to convert and call unit vector version of the method.
     *
     * @param theta The angle of colatitude.
     *
     * @param phi The angle of longitude.
     *
     * @return The value of the function at (theta, phi).
     */
    public double getRadius(double theta, double phi) {
        double[] uVec = toUnitVec(theta, phi);
        return getRadius(uVec[0], uVec[1], uVec[2]);
    }

    /**
     * Return the radius along a specified vector.
     *
     * @param x The x coordinate of the point on the sphere
     *
     * @param y The y coordinate of the point on the sphere
     *
     * @param z The z coordinate of the point on the sphere
     *
     * @return The value at (x, y, z)
     */
    public abstract double getRadius(double x, double y, double z);

    /**
     * Computes the order-th moment of the function over the
     * sphere. The default method computes the moment numerically.
     *
     * @param order The order of the moment to compute.
     *
     * @return the value of the moment.
     */
    public double moment(int order) {

        try {
            double[][] ksInt = ISCodes.getPointSetForMaxEnt(CL_Initializer.pointSetInd).data;
	    return moment(order, ksInt);
        }
        catch (Exception e) {
	    logger.warning("WARNING: " + e);
        }

	return 0.0;
    }
	

    /**
     * Computes the order-th moment of the function over the      
     * sphere. The default method computes the moment numerically.
     *                                                 
     * @param order The order of the moment to compute.
     *                                 
     * @param ksInt List of points at which to sample the function.
     *                                 
     * @return the value of the moment.
     */
    public double moment(int order, double[][] ksInt) {

	try{
            double sum = 0.0;
            double devSum = 0.0;
            for (int i = 0; i < ksInt.length; i++) {
                double con = getRadius(ksInt[i][0], ksInt[i][1], ksInt[i][2]);
                sum += Math.pow(con, (double) order);
            }

            return sum / (double) ksInt.length;
        }
        catch (Exception e) {
           logger.warning("WARNING: " + e);
        }

        return 0.0;
    }

    /**
     * Computes the normalised order-th moment of the function over
     * the sphere.  The default method computes the moment
     * numerically.  The normalized moment is moment(order) normalized
     * by the order-th moment of p - \bar{p} where \bar{p} is
     * moment(1).
     *
     * @param order The order of the normalized moment.
     *
     * @return the value of the normalized moment.
     */
    public double normMoment(int order) {

        try {
            double[][] ksInt = ISCodes.getPointSetForMaxEnt(CL_Initializer.pointSetInd).data;
	    return normMoment(order, ksInt);
        }
        catch (Exception e) {
	    logger.warning("WARNING: " + e);
        }

	return 0.0;
    }
	

    /**
     * Computes the normalised order-th moment of the function over
     * the sphere.  The default method computes the moment
     * numerically.  The normalized moment is moment(order) normalized
     * by the order-th moment of p - \bar{p} where \bar{p} is
     * moment(1).
     *                                                 
     * @param order The order of the normalized moment to compute.
     *                                 
     * @param ksInt List of points at which to sample the function.
     *                                 
     * @return the value of the normalized moment.
     */
    public double normMoment(int order, double[][] ksInt) {

        // Need the mean of the function.
        double mean = mean();

        //Otherwise compute the moment numerically.
        try {
            double sum = 0.0;
            double devSum = 0.0;
            for (int i = 0; i < ksInt.length; i++) {
                double con = getRadius(ksInt[i][0], ksInt[i][1], ksInt[i][2]);
                sum += Math.pow(con, (double) order);
                devSum += Math.pow(con - mean, (double) order);
            }

            return (sum==0)?0:(devSum / sum);
        }
        catch (Exception e) {
            System.out.println("WARNING: " + e);
        }

        return 0.0;
    }

    /**
     * Computes the order-th central moment of the function over the
     * sphere.  The default method computes the moment numerically.
     * The central moment is the order-th moment of p - \bar{p} where
     * \bar{p} is moment(1).
     *
     * @param order The order of the central moment.
     *
     * @return the value of the central moment.
     */
    public double centralMoment(int order) {

        try {
            double[][] ksInt = ISCodes.getPointSetForMaxEnt(CL_Initializer.pointSetInd).data;
	    return centralMoment(order, ksInt);
        }
        catch (Exception e) {
	    logger.warning("WARNING: " + e);
        }

	return 0.0;
    }
	

    /**
     * Computes the order-th central moment of the function over the
     * sphere.  The default method computes the moment numerically.
     * The central moment is the order-th moment of p - \bar{p} where
     * \bar{p} is moment(1).
     *                                                 
     * @param order The order of the central moment to compute.
     *                                 
     * @param ksInt List of points at which to sample the function.
     *                                 
     * @return the value of the central moment.
     */
    public double centralMoment(int order, double[][] ksInt) {

        // Need the mean of the function.
        double mean = mean(ksInt);

        //Otherwise compute the moment numerically.
        try {
            double devSum = 0.0;
            for (int i = 0; i < ksInt.length; i++) {
                double con = getRadius(ksInt[i][0], ksInt[i][1], ksInt[i][2]);
                devSum += Math.pow(con - mean, (double) order);
            }

            return devSum;
        }
        catch (Exception e) {
            System.out.println("WARNING: " + e);
        }

        return 0.0;
    }

    /**
     * The mean of the function is the order 1 moment.
     *
     * @return The mean of the function.
     */
    public double mean() {
        return moment(1);
    }

    /**
     * The mean of the function is the order 1 moment.
     * Numerical version.
     *
     * @param ksInt List of points at which to sample the function.
     *                                 
     * @return The mean of the function.
     */
    public double mean(double[][] ksInt) {
        return moment(1, ksInt);
    }

    /**
     * The anisotropy is the square root of the order 2 normalised
     * moment.
     *
     * @return The anisotropy
     */
    public double anisotropy() {
        return Math.sqrt(normMoment(2));
    }

    /**
     * The skewness is the cube root of order 3 normalised moment.
     *
     * @return the skewness
     */
    public double skewness() {
        double sk = normMoment(3);
        sk = (sk < 0) ? (-Math.pow(-sk, 1. / 3.)) : (Math.pow(sk, 1. / 3.));

        return sk;
    }

    /**
     * The kurtosis is the 4th root of the order 4 normalised moment
     * normalised by the anisotropy.
     *
     * @return the kurtosis
     */
    public double kurtosis() {
        double fullK = Math.sqrt(Math.sqrt(normMoment(4)));
        double anis = anisotropy();

        //The useful measure normalises the full kurtosis by the
        //anisotropy.
        if (anis != 0.0) {
            return fullK / anis;
        }

        return 1.0;
    }

    /**
     * Computes the mean, variance, max and min of the spherical function and
     * returns them in an array. The default method computes all the statistics
     * numerically.
     * 
     * @return [mean, var, max, min]
     */
    public double[] getStats() {

        double meanme = 0.0;
        double varme = 0.0;
        double maxme = 0.0;
        double minme = 0.0;

        try {

            //Compute integral with higher resolution mid point rule.
            double[][] ksInt = ISCodes.getPointSetForMaxEnt(CL_Initializer.pointSetInd).data;

            double sum = 0.0;
            double sum2 = 0.0;
            for (int i = 0; i < ksInt.length; i++) {
                double con = getRadius(ksInt[i][0], ksInt[i][1], ksInt[i][2]);
                sum += con;
                sum2 += con * con;

                if (i == 0) {
                    maxme = con;
                    minme = con;
                }
                else {
                    maxme = (con > maxme) ? con : maxme;
                    minme = (i == 0 || con < minme) ? con : minme;
                }
            }

            meanme = sum / (double) ksInt.length;
            varme = sum2 / (double) ksInt.length - meanme * meanme;

        }
        catch (Exception e) {
            System.out.println(e);
        }

        double[] stats = new double[4];
        stats[0] = meanme;
        stats[1] = varme;
        stats[2] = maxme;
        stats[3] = minme;

        return stats;

    }

    /**
     * Object for finding the minimum value of the radius of the PAS from some
     * specified starting point on the sphere in spherical polar coordinates.
     */
    protected class SFPD_Minimiser extends PowellMinimiser {

        private double paramScale = 0.1;

        public SFPD_Minimiser() {
            super(2);
        }

        public SFPD_Minimiser(double pScale) {
            super(2);
            paramScale = pScale;
        }

        /**
         * The objective function is minus the value of the function.
         *
         * @param params Values specifying the point on the sphere at
         * which to evaluate the function.
         *
         * @return The value of the function.
         */
        protected double fObj(double[] params) {
            return -getRadius(params[1] * paramScale, params[2] * paramScale);
        }


        /**
         * Scales the parameters before starting the usual
         * minimization routine.
         *
         * @param as Current point to test.
         *
         * @param ct Convergence threshold
         *
         * @return The minimum value of the objective function.
         */
        public double minimise(double[] as, double ct) throws PowellMinimiserException {
            for (int i = 0; i < as.length; i++) {
                as[i] /= paramScale;
            }

            return super.minimise(as, ct);
        }


        /**
         * Returns the parameters after minimization.
         *
         * @return the parameters that minimize the objective
         * function.
         */
        public double[] getMinParams() {
            double[] res = super.getMinParams();
            for (int i = 0; i < res.length; i++) {
                res[i] *= paramScale;
            }

            return res;
        }

    }


    /**
     * Puts a list of peak directions found through sampling
     * followed by local optimization in pdsOpt.
     * 
     * @param samplePoints List of points on the sphere at which to
     * sample the function in the initial search for maxima.
     * 
     * @param pdsOpt
     *            A peak direction list that the method fills with the
     *            directions it finds.
     * 
     * @return An array of the distances moved by each principal direction
     *         during the final optimization phase. The value is -1 if the
     *         optimization failed.
     */
    public double[] getPDs(double[][] samplePoints, PDList pdsOpt) {

        //Get the PDs from the random sampling technique.
        PDList pdsRS = getPDsRS(samplePoints);

        //Now use Powell to refine each one.
        double[] distances = new double[pdsRS.getNoPDs()];
        SFPD_Minimiser p = new SFPD_Minimiser(SFPD_PARAMSCALE);

        boolean allOK = true;

        for (int i = 0; i < pdsRS.getNoPDs(); i++) {
            PD rsPD = pdsRS.getPD(i);
            double[] as = { 0.0, rsPD.getTheta(), rsPD.getPhi() };
            try {
                double minVal = -p.minimise(as, MINCONV);
                double[] minPoint = p.getMinParams();
                PD optPD = new PD(minPoint[1], minPoint[2], minVal);
                pdsOpt.addPD(optPD);

                double xRS = rsPD.getPDX();
                double yRS = rsPD.getPDY();
                double zRS = rsPD.getPDZ();

                double xOpt = optPD.getPDX();
                double yOpt = optPD.getPDY();
                double zOpt = optPD.getPDZ();

                double dotProd = xRS * xOpt + yRS * yOpt + zRS * zOpt;
                double dist = Math.sqrt((xRS - xOpt) * (xRS - xOpt) + (yRS - yOpt)
                        * (yRS - yOpt) + (zRS - zOpt) * (zRS - zOpt));

                distances[i] = dist;

                if (dist > searchRadius) {
                    allOK = false;
                }

            }
            catch (Exception e) {
                pdsOpt.addPD(rsPD);
                distances[i] = -1.0;
                logger.warning(e.toString());
                allOK = false;
            }
        }

        //Print PDs if there were problems.
        if (!allOK) {
            for (int i = 0; i < pdsRS.getNoPDs(); i++) {
                logger.info("Before opt (" + i + "): " + pdsRS.getPD(i));
            }
            for (int i = 0; i < pdsOpt.getNoPDs(); i++) {
                logger.info("After opt  (" + i + "): " + pdsOpt.getPD(i));
            }
        }

        return distances;
    }


    /**
     * Samples the function at each point and returns a list of points
     * that are maximal with a fixed radius.
     * 
     * @param samplePoints List of points on the sphere at which to
     * sample the function
     *
     * @param seedOffset
     *            An offset of the seeds used in the random number generator
     *            used to get the random rotations.
     */
    public PDList getPDsRS(double[][] samplePoints) {

        // Now compute the function value at each sample point.
        double[] rs = new double[samplePoints.length];
        for (int i = 0; i < rs.length; i++) {
            rs[i] = getRadius(samplePoints[i][0], samplePoints[i][1], samplePoints[i][2]);
        }

        // Now thin out the list of points to only those that have
        // maximum value within a fixed radius on the sphere.
        PDList maxima = getMaxima(samplePoints, rs);

        return maxima;
    }


    /**
     * Sets the search radius in the random sampling principal direction
     * extraction algorithm.
     *
     * @param newSearchRadius The new value of the search radius.
     */
    public static void setSearchRadius(double newSearchRadius) {
        searchRadius = newSearchRadius;
    }


    /**
     * Returns the value of the search radius in the random sampling principal
     * direction extraction algorithm.
     *
     * @return The search radius.
     */
    public static double getSearchRadius() {
        return searchRadius;
    }


    /**
     * Finds a list of local maxima from a list of points on the
     * sphere together with a list of the values of the function at
     * each.
     *
     * @param samplePoints The list of points at which the function is
     * sampled.
     *
     * @param rs The values of the function at each sample point.
     *
     * @return The list of local maxima.
     */
    private PDList getMaxima(double[][] samplePoints, double[] rs) {

        boolean[] isMax = new boolean[rs.length];
        for (int j = 0; j < isMax.length; j++) {
            isMax[j] = true;
        }

        // Each point is a local maximum if all points within a fixed
        // radius have smaller PAS value.

        // We ignore the sign here, since the PAS is symmetric on the
        // sphere.
        for (int i = 0; i < rs.length; i++) {
            if (isMax[i]) {
                for (int j = 0; j < rs.length; j++) {

                    double d1 = Math.sqrt((samplePoints[i][0] - samplePoints[j][0])
                            * (samplePoints[i][0] - samplePoints[j][0])
                            + (samplePoints[i][1] - samplePoints[j][1])
                            * (samplePoints[i][1] - samplePoints[j][1])
                            + (samplePoints[i][2] - samplePoints[j][2])
                            * (samplePoints[i][2] - samplePoints[j][2]));

                    double d2 = Math.sqrt((samplePoints[i][0] + samplePoints[j][0])
                            * (samplePoints[i][0] + samplePoints[j][0])
                            + (samplePoints[i][1] + samplePoints[j][1])
                            * (samplePoints[i][1] + samplePoints[j][1])
                            + (samplePoints[i][2] + samplePoints[j][2])
                            * (samplePoints[i][2] + samplePoints[j][2]));

                    double distance = (d1 < d2) ? d1 : d2;

                    if ((distance < searchRadius) && (distance != 0)) {
                        if (rs[j] > rs[i]) {

                            // This point is not a maximum. We can stop
                            // searching.
                            isMax[i] = false;
                            j = rs.length;
                        }
                        else {

                            // We know that the other point cannot
                            // be a maximum.
                            isMax[j] = false;
                        }
                    }
                }
            }

        }

        // Now create the principal direction list containing the
        // maximum points.
        PDList pds = new PDList();
        for (int i = 0; i < isMax.length; i++) {
            if (isMax[i]) {
                // 		System.err.println(i);
                pds.addPD(new PD(samplePoints[i], rs[i]));
            }
        }

        // 	for(int i=0; i<pds.getNoPDs(); i++) {
        // 	    System.err.println("PD " + i + " " + pds.getPD(i));
        // 	}

        return pds;
    }


    /**
     * Sets the convergence threshold for Powell's algorithm used in
     * getPDs.  Default value is 1.0E-12.
     *
     * @param thresh The new convergence threshold.
     */
    public void setConvThresh(double thresh) {
        MINCONV = thresh;
    }


    /**
     * Computes the Hessian of the function at the specified
     * point. Converts the PD object to an array and calls the other
     * version of the method.
     * 
     * @param pd
     *            The principal direction at which to compute the Hessian.
     * 
     * @return The Hessian.
     */
    public RealMatrix getHessian(PD pd) {
        double[] x = new double[3];
        x[0] = pd.getPDX();
        x[1] = pd.getPDY();
        x[2] = pd.getPDZ();

        return getHessian(x);
    }


    /**
     * Computes the Hessian of the function at the specified
     * point. This general method computes the Hessian numerically.
     * 
     * @param x
     *            An array specifying the point on the sphere, in Cartesian
     *            coordinates ([x, y, z]) at which to compute the Hessian
     *            matrix.
     * 
     * @return The Hessian matrix: [[d^2 f/d \mu_1^2, d^2 f/d \mu_1 \mu_2]; [d^2
     *         f/d \mu_2 \mu_1, d^2 f/d \mu_2^2]], where \mu_1 and \mu_2 are
     *         distances along tangent vectors.
     */
    public RealMatrix getHessian(double[] x) {

        // Normalize the input vector.
        double xMag = Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        for (int i = 0; i < 3; i++) {
            x[i] = x[i] / xMag;
        }

        // Compute the value of the function at x.
        double fX = getRadius(x[0], x[1], x[2]);

        // Estimate first and second derivatives numerically.

        // This is the stepsize for numerical derivative estimation.
        double h = 0.00001;

        // These are two orthogonal vectors in the tangent plane at x.
        double[] tangent1 = new double[3];
        double t1Mag = Math.sqrt(x[0] * x[0] + x[1] * x[1]);
        if(t1Mag == 0) {
    	    tangent1[0] = 1;
    	    tangent1[1] = 0;
    	    tangent1[2] = 0;
    	}
    	else {
    	    tangent1[0] = x[1] / t1Mag;
    	    tangent1[1] = -x[0] / t1Mag;
    	    tangent1[2] = 0.0;
    	}

            double[] tangent2 = new double[3];
    	if(t1Mag == 0) {
    	    tangent2[0] = 0;
    	    tangent2[1] = 1;
    	    tangent2[2] = 0;
    	}
    	else {
    	    tangent2[0] = -x[0] * x[2];
    	    tangent2[1] = -x[1] * x[2];
    	    tangent2[2] = t1Mag * t1Mag;
    	}

        double t2Mag = Math.sqrt(tangent2[0] * tangent2[0] + tangent2[1] * tangent2[1]
                + tangent2[2] * tangent2[2]);
        for (int i = 0; i < 3; i++) {
            tangent2[i] = tangent2[i] / t2Mag;
        }

        // These are slightly displaced points at which to compute
        // the function value for estimating the derivatives.
        double[] xpdx = new double[3];
        double[] xpdy = new double[3];
        double[] xmdx = new double[3];
        double[] xmdy = new double[3];
        double[] xpdxpdy = new double[3];
        double[] xmdxpdy = new double[3];
        double[] xmdxmdy = new double[3];
        double[] xpdxmdy = new double[3];
        for (int i = 0; i < 3; i++) {
            xpdx[i] = x[i] + h * tangent1[i];
            xpdy[i] = x[i] + h * tangent2[i];
            xmdx[i] = x[i] - h * tangent1[i];
            xmdy[i] = x[i] - h * tangent2[i];
            xpdxpdy[i] = x[i] + h * tangent1[i] + h * tangent2[i];
            xmdxpdy[i] = x[i] - h * tangent1[i] + h * tangent2[i];
            xmdxmdy[i] = x[i] - h * tangent1[i] - h * tangent2[i];
            xpdxmdy[i] = x[i] + h * tangent1[i] - h * tangent2[i];
        }
        double xpdxMag = Math.sqrt(xpdx[0] * xpdx[0] + xpdx[1] * xpdx[1] + xpdx[2]
                * xpdx[2]);
        double xpdyMag = Math.sqrt(xpdy[0] * xpdy[0] + xpdy[1] * xpdy[1] + xpdy[2]
                * xpdy[2]);
        double xmdxMag = Math.sqrt(xmdx[0] * xmdx[0] + xmdx[1] * xmdx[1] + xmdx[2]
                * xmdx[2]);
        double xmdyMag = Math.sqrt(xmdy[0] * xmdy[0] + xmdy[1] * xmdy[1] + xmdy[2]
                * xmdy[2]);
        double xpdxpdyMag = Math.sqrt(xpdxpdy[0] * xpdxpdy[0] + xpdxpdy[1] * xpdxpdy[1]
                + xpdxpdy[2] * xpdxpdy[2]);
        double xmdxpdyMag = Math.sqrt(xmdxpdy[0] * xmdxpdy[0] + xmdxpdy[1] * xmdxpdy[1]
                + xmdxpdy[2] * xmdxpdy[2]);
        double xmdxmdyMag = Math.sqrt(xmdxmdy[0] * xmdxmdy[0] + xmdxmdy[1] * xmdxmdy[1]
                + xmdxmdy[2] * xmdxmdy[2]);
        double xpdxmdyMag = Math.sqrt(xpdxmdy[0] * xpdxmdy[0] + xpdxmdy[1] * xpdxmdy[1]
                + xpdxmdy[2] * xpdxmdy[2]);
        for (int i = 0; i < 3; i++) {
            xpdx[i] = xpdx[i] / xpdxMag;
            xpdy[i] = xpdy[i] / xpdyMag;
            xmdx[i] = xmdx[i] / xmdxMag;
            xmdy[i] = xmdy[i] / xmdyMag;
            xpdxpdy[i] = xpdxpdy[i] / xpdxpdyMag;
            xmdxpdy[i] = xmdxpdy[i] / xmdxpdyMag;
            xmdxmdy[i] = xmdxmdy[i] / xmdxmdyMag;
            xpdxmdy[i] = xpdxmdy[i] / xpdxmdyMag;
        }

        //Now compute the derivatives
        double fXpdX = getRadius(xpdx[0], xpdx[1], xpdx[2]);
        double fXmdX = getRadius(xmdx[0], xmdx[1], xmdx[2]);
        double fXpdY = getRadius(xpdy[0], xpdy[1], xpdy[2]);
        double fXmdY = getRadius(xmdy[0], xmdy[1], xmdy[2]);
        double fXpdXpdY = getRadius(xpdxpdy[0], xpdxpdy[1], xpdxpdy[2]);
        double fXmdXpdY = getRadius(xmdxpdy[0], xmdxpdy[1], xmdxpdy[2]);
        double fXmdXmdY = getRadius(xmdxmdy[0], xmdxmdy[1], xmdxmdy[2]);
        double fXpdXmdY = getRadius(xpdxmdy[0], xpdxmdy[1], xpdxmdy[2]);

        double dfdxNum = (fXpdX - fXmdX) / (2.0 * h);
        double dfdyNum = (fXpdY - fXmdY) / (2.0 * h);
        double d2fdx2Num = (fXpdX - 2.0 * fX + fXmdX) / (h * h);
        double d2fdy2Num = (fXpdY - 2.0 * fX + fXmdY) / (h * h);
        double d2fdxdyNum = (fXpdXpdY - fXmdXpdY - fXpdXmdY + fXmdXmdY) / (4.0 * h * h);

        // Construct Hessian matrix
        RealMatrix H = new RealMatrix(2, 2);
        H.entries[0][0] = d2fdx2Num;
        H.entries[1][0] = d2fdxdyNum;
        H.entries[0][1] = d2fdxdyNum;
        H.entries[1][1] = d2fdy2Num;

        return H;
    }


    //Some general protected methods used by derived classes.

    /**
     * Convert a pair of angles to a unit vector. The returned array
     * contains the x, y and z components of the unit vector in that
     * order. Theta is the angle with the positive z-axis
     * (latitude). From 0 to pi. Phi is the angle with the positive
     * x-axis (longitude). from 0 to 2pi.
     *
     * @param theta The angle of colatitude.
     *
     * @param phi The angle of longitude
     *
     * @return (x, y, z)
     */
    protected double[] toUnitVec(double theta, double phi) {
        double[] uVec = new double[3];
        uVec[0] = Math.cos(phi) * Math.sin(theta);
        uVec[1] = Math.sin(phi) * Math.sin(theta);
        uVec[2] = Math.cos(theta);

        return uVec;
    }


    /**
     * Convert a unit vector to a pair of angles of colatitude and
     * longitude.
     *
     * @param x The x coordinate of the point on the sphere
     *
     * @param y The y coordinate of the point on the sphere
     *
     * @param z The z coordinate of the point on the sphere
     *
     * @return (theta, phi)
     */
    protected double[] toAngles(double x, double y, double z) {
        double[] angles = new double[2];

        double mag = Math.sqrt(x * x + y * y + z * z);

        //Check against zero vector.
        if (mag == 0.0) {
            angles[0] = 0.0;
            angles[1] = 0.0;

            return angles;
        }

        double ux = x / mag;
        double uy = y / mag;
        double uz = z / mag;

        double theta = Math.acos(uz);
        double phi = 0.0;

        //Careful of tiny values of theta. In such cases phi is ill defined
        // anyway
        //so just set to zero.
        //Also small numerical errors can shift the cos and sin of phi over 1,
        // so
        //we need to explicitly correct.
        double st = Math.sin(theta);
        if (Math.abs(st) >= 1.0E-7) {
            double cosPhi = ux / st;
            cosPhi = (cosPhi > 1.0) ? 1.0 : cosPhi;
            cosPhi = (cosPhi < -1.0) ? -1.0 : cosPhi;
            phi = Math.acos(cosPhi);

            double sinPhi = uy / st;
            sinPhi = (sinPhi > 1.0) ? 1.0 : sinPhi;
            sinPhi = (sinPhi < -1.0) ? -1.0 : sinPhi;
            if (Math.asin(sinPhi) < 0.0) {
                phi = 2.0 * Math.PI - phi;
            }
        }

        angles[0] = theta;
        angles[1] = phi;

        return angles;
    }
    
    
    /**
     * constructs an array of n angles equally spaced
     * from 0 to (n-1)*2.0*pi/n
     * 
     * @param n number of desired angles
     * @return array of n angles
     */
    private final double[] getGreatCircleAngles(int n){
        
        double[] t= new double[n];
        
        for(int i=0; i<n; i++){
            t[i]=i*2.0*Math.PI/n;
        }
        
        return t;
    }
    
    
    /** obtains the points on the surface of a sphere for the great circle intgration
     *  path
     * 
     * @param u axis defining plne of circle
     * @param t angles around circle
     * @param n number of points
     * 
     * @return points on path in spherical (theta, phi) coords
     */
    private final double[][] getGreatCirclePoints(double[] u, double[] t, int n){
        
        
        if(t.length != n){
            String errMess = new String("number of angles for great circle intgration is "+t.length+" but degree is "+n+". no good!");
            logger.severe(errMess);
            throw new RuntimeException(errMess);
        }
        
        double[][] r= new double[n][];
                
        for(int i=1; i<=n; i++){
            double angle = (t[i%n] + t[i-1])/2.0;
            
            r[i%n] = getPointOnCircle(u, angle);
        }
        
        return r;
    }
    
    
    /**
     * returns array of length of contour segments.
     * These are calculated as path lengths along the unit
     * circle for the arc from one point to the next.
     * 
     * sounds grandiose, but these are nothing other than
     * the angle differences.
     * 
     * @param u unit vector defining plane of circle
     * @param t angles used for points
     * 
     * @return array of segment lengths
     */
    private final double[] getGreatCircleSegmentLengths(double[] u, double[] t){
        
        double[] h= new double[t.length];
        int n=t.length;
        
        // first length is difference between 2pi and final point
        h[0]=2.0*Math.PI-t[n-1];
        
        // all the rest are difference between successive points
        for(int i=1; i<n; i++){
            h[i]=t[i]-t[i-1];            
        }
        
        return h;
    }
    
    /**
     * returns a point on the circle defined byt the plane with normal u
     * which is itself defined by the angle t around that circle. The circle
     * is obtained by constructing the spherical unit vectors defined by
     * 
     *          cos th sin ph
     * r= u =   sin th sin ph
     *          cos ph
     * 
     *         -sin th 
     * Th =     cos th
     *             0
     * 
     *          cos th cos ph
     * Ph =     sin th cos ph
     *          -sin th
     * 
     * Th and Ph and constructed from u, and the point on the circle is defined
     * as
     * 
     *  r = Th cos(t) + Ph sin(t)
     * 
     * 
     * @param u unit normal to plane of circle
     * @param t angle around circle
     * 
     * @return cartesian vector defining point on sphere
     */
    private final double[] getPointOnCircle(double[] u, double t){
        
        // get the sine and cosines from the coords of u
        double cosPhi=u[2];
        double sinPhi=Math.sqrt(1.0-cosPhi*cosPhi);
        
        double cosTheta=(sinPhi==0.0)? 0.0 : u[0]/sinPhi;
        double sinTheta=(sinPhi==0.0)? 1.0 : u[1]/sinPhi;
        
        // construct the unit vectors Th and Ph 
        double[] Th= new double[] {-sinTheta, cosTheta, 0.0};
        double[] Ph= new double[] {cosTheta*cosPhi, sinTheta*cosPhi, -sinPhi};
        
        double cosT=Math.cos(t);
        double sinT=Math.sin(t);
        
        // construct point on circle from Th, Ph and t
        double[] pointOnCircle= new double[] {Th[0]*cosT + Ph[0]*sinT,
                							  Th[1]*cosT + Ph[1]*sinT,
                							  Th[2]*cosT + Ph[2]*sinT};
 
        
        return pointOnCircle;
        
    }
    
    
    /**
     * private method which performs a numerical integration
     * around a sphere where the path is divided into n segments
     * 
     * we have n points around the circle (i=0, ... , n-1)
     * 
     * \[
     *   t_i = \frac{2\pi}{n}i
     * \]
     * 
     * and a circle on the sphere defining a set of points
     * \[
     *   \mathbf{r}_i = \mathbf{l}\left(\frac{t_i + u_{t-1}}{2}\right)
     * \]
     * 
     * and a set of segment lengths
     * 
     * \[
     *   h_i = |\mathbf{l}(t_i) - \mathbf{l}(t_{i-1})|
     * \]
     * 
     * then the great circle integral of  integral is approximated as
     * 
     * \[
     *   \int_{\mathbf{l}} A(\mathbf{r})d\mathbf{l} 
     *                  \simeq \sum_{i=0}^{n-1}A(\mathbf{r}_ih_i
     * \]
     * 
     * 
     * @param u unit vector defining plane of great circle
     * @param n number of segments
     * 
     * @return
     */
    private double greatCircleIntegral(double[] u, int n){
        
        /** the angles */
        double[] t= getGreatCircleAngles(n);
        
        /** the points to evaluate the function at */
        double[][] r= getGreatCirclePoints(u, t, n);
        
        /** the line segment lengths */
        double[] h= getGreatCircleSegmentLengths(u, t);
        
        double integral=0.0;
        
        for(int i=0; i<n; i++){
            double radius= getRadius(r[i][0], r[i][1], r[i][2]);
            integral+=radius*h[i];
        }
        
    	return integral;
    }
    
    
    /** 
     * returns the great circle integral given a vector
     *  normal to the plane of the great circle
     * 
     * @param u unit vector defining plane of circle vector
     * 
     * @return result of performing numerical integral 
     */
    public double numGreatCircleIntegral(double[] u){
        
        double result=-1.0;;
        int level=MIN_GREAT_CIRCLE_POINTS;
        
        double last = greatCircleIntegral(u, MIN_GREAT_CIRCLE_POINTS);
        
        
        for(int n=MIN_GREAT_CIRCLE_POINTS+1; ; n+=GREAT_CIRCLE_POINTS_INC){
            double next= greatCircleIntegral(u, n);
            
            result=next;
            level=n; 
            
            if(Math.abs(last-next)<= GREAT_CIRCLE_INTEGRAL_THRESHOLD){
                break;
            }
            
            
            if(n>MAX_GREAT_CIRCLE_POINTS){
                logger.warning("WARNING: great circle integral about u=("+u[0]+","+u[1]+","+u[2]+
                        ") did not converge after max number of steps (max steps = "+MAX_GREAT_CIRCLE_POINTS+
                        ", threshold = "+GREAT_CIRCLE_INTEGRAL_THRESHOLD+", diff = "+Math.abs(last-next));
                logger.warning("returning integral value at n="+n);
                
                break;
            }

            last=next;
            
        }

        
        return result;
        
    }
    
    
    /** 
     * returns the great circle integral given a vector normal to the
     * plane of the great circle.  The default integral calculation is
     * numerical and just calls numGreatCircleIntegral, but we can
     * override this for specific functions where the integral is
     * analytic.
     * 
     * @param u unit vector defining plane of circle vector
     * 
     * @return The integral 
     */
    public double greatCircleIntegral(double[] u){
        return numGreatCircleIntegral(u);
    }


    /**
     * Like greatCircleIntegral, but raises the function to
     * power pow.
     * 
     * @param u unit vector defining plane of great circle
     * @param pow the power to raise the function to
     * @param n number of segments
     * 
     * @return
     */
    private double greatCirclePowerIntegral(double[] u, double pow, int n){
        
        /** the angles */
        double[] t= getGreatCircleAngles(n);
        
        /** the points to evaluate the function at */
        double[][] r= getGreatCirclePoints(u, t, n);
        
        /** the line segment lengths */
        double[] h= getGreatCircleSegmentLengths(u, t);
        
        double integral=0.0;
        
        for(int i=0; i<n; i++){
            double radius= Math.pow(getRadius(r[i][0], r[i][1], r[i][2]), pow);
            integral+=radius*h[i];
        }
        
    	return integral;
    }
    

    /** 
     * returns the great circle integral of the function
     * raised to power pow given a vector
     *  normal to the plane of the great circle
     * 
     * @param u unit vector defining plane of circle vector
     *
     * @param pow the power to raise the function to
     * 
     * @return result of performing numerical integral 
     */
    public double numGreatCirclePowerIntegral(double[] u, double pow){
        
        double result=-1.0;;
        int level=MIN_GREAT_CIRCLE_POINTS;
        
        double last = greatCirclePowerIntegral(u, pow, MIN_GREAT_CIRCLE_POINTS);
        
        
        for(int n=MIN_GREAT_CIRCLE_POINTS+1; ; n+=GREAT_CIRCLE_POINTS_INC){
            double next= greatCirclePowerIntegral(u, pow, n);
            
            result=next;
            level=n; 
            
            if(Math.abs(last-next)<= GREAT_CIRCLE_INTEGRAL_THRESHOLD){
                break;
            }
            
            
            if(n>MAX_GREAT_CIRCLE_POINTS){
                logger.warning("WARNING: great circle integral about u=("+u[0]+","+u[1]+","+u[2]+
                        ") did not converge after max number of steps (max steps = "+MAX_GREAT_CIRCLE_POINTS+
                        ", threshold = "+GREAT_CIRCLE_INTEGRAL_THRESHOLD+", diff = "+Math.abs(last-next));
                logger.warning("returning integral value at n="+n);
                
                break;
            }

            last=next;
            
        }

        
        return result;
        
    }
    
    
    /** 
     * returns the great circle integral of the spherical function 
     * raised to power pow given a vector normal to the
     * plane of the great circle.  The default integral calculation is
     * numerical and just calls numGreatCircleIntegral, but we can
     * override this for specific functions where the integral is
     * analytic.
     * 
     * @param u unit vector defining plane of circle vector
     * 
     * @param pow the power to raise the function to.
     *
     * @return The integral 
     */
    public double greatCirclePowerIntegral(double[] u, double pow){
        return numGreatCirclePowerIntegral(u, pow);
    }


    /**
     * Like greatCirclePowerIntegral, but raises subtracts the mean
     * from the function.
     * 
     * @param u unit vector defining plane of great circle
     * @param pow the power to raise the function to
     * @param n number of segments
     * 
     * @return
     */
    private double greatCircleCentralMoment(double[] u, double pow, int n){
        
        /** the angles */
        double[] t= getGreatCircleAngles(n);
        
        /** the points to evaluate the function at */
        double[][] r= getGreatCirclePoints(u, t, n);
        
        /** the line segment lengths */
        double[] h= getGreatCircleSegmentLengths(u, t);
        
        double integral=0.0;
        
        double mn = greatCircleIntegral(u, n)/(2.0*Math.PI);

        for(int i=0; i<n; i++){
            double radius= Math.pow(getRadius(r[i][0], r[i][1], r[i][2]) - mn, pow);
            integral+=radius*h[i];
        }
        
    	return integral;
    }
    

    /** 
     * returns the great circle central moment of order pow given a
     * vector normal to the plane of the great circle
     * 
     * @param u unit vector defining plane of circle vector
     *
     * @param pow the power to raise the function to
     * 
     * @return result of performing numerical integral 
     */
    public double numGreatCircleCentralMoment(double[] u, double pow){
        
        double result=-1.0;;
        int level=MIN_GREAT_CIRCLE_POINTS;
        
        double last = greatCircleCentralMoment(u, pow, MIN_GREAT_CIRCLE_POINTS);
        
        
        for(int n=MIN_GREAT_CIRCLE_POINTS+1; ; n+=GREAT_CIRCLE_POINTS_INC){
            double next= greatCircleCentralMoment(u, pow, n);
            
            result=next;
            level=n; 
            
            if(Math.abs(last-next)<= GREAT_CIRCLE_INTEGRAL_THRESHOLD){
                break;
            }
            
            
            if(n>MAX_GREAT_CIRCLE_POINTS){
                logger.warning("WARNING: great circle integral about u=("+u[0]+","+u[1]+","+u[2]+
                        ") did not converge after max number of steps (max steps = "+MAX_GREAT_CIRCLE_POINTS+
                        ", threshold = "+GREAT_CIRCLE_INTEGRAL_THRESHOLD+", diff = "+Math.abs(last-next));
                logger.warning("returning integral value at n="+n);
                
                break;
            }

            last=next;
            
        }

        
        return result;
        
    }
    
    
    /** 
     * returns the central moment of order pow of the spherical
     * function given a vector normal to the plane of the great
     * circle.  The default integral calculation is numerical and just
     * calls numGreatCircleIntegral, but we can override this for
     * specific functions where the integral is analytic.
     * 
     * @param u unit vector defining plane of circle vector
     * 
     * @param pow the power to raise the function to.
     *
     * @return The integral 
     */
    public double greatCircleCentralMoment(double[] u, double pow){
        return numGreatCircleCentralMoment(u, pow);
    }


}


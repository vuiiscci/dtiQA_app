package numerics;

import java.util.Random;



/**
 * Bingham Distribution. Some mathematical definitions:
 * <dl>
 * <dd> T - scatter matrix of a series of N samples from the distribution: \sum_{i=1}^N (1 / N) x x^T.
 * <dd> k1 - first concentration parameter.
 * <dd> k2 - second concentration parameter.
 * <dd> x - axis on the sphere, |x| == 1.
 *
 * @see section 9.4.3 of Mardia and Jupp ("Directional Statistics", Wiley, 2000)
 * @version $Id$
 * @author Philip Cook
 */
public final class BinghamDistribution implements AxialDistribution {

    /**
     * The ordering of eigenvectors is reversed...so e1 is associated with the smallest 
     * eigenvector of the scatter matrix you would get from this distribution. 
     * This is the numbering scheme of Watson and Kent, and is preserved for consistency with Kent's code.
     */
    private final Vector3D e1;

    private final Vector3D e2;

    private final Vector3D e3;

    private final double k1;
    private final double k2;

    // k3 == 0

    // normalization constant
    private final double normC;
    
    private final WatsonDistribution wrapper;
    
    private final double wrapperK;

    // multiply wrapper.pdf() by this to get the correct density
    private final double wrapperPDFNormC; 

    private final Random ran;
    private final double logNormC; 

    // don't want anyone to use this
    private BinghamDistribution() {

	e1 = null;
	e2 = null;
	e3 = null;

	k1 = 0.0;
	k2 = 0.0;

	normC = 1.0;
	logNormC = 0.0;
	wrapperK = 0.0;
	wrapperPDFNormC = 0.0;

	ran = new Random();

	wrapper = null;
	
    }

    /**
     * @param axes the axes sorted by their eigenvalues (greatest first). In the limit that k1 == k2, the distribution is Watson with mu == axes[2].
     * @param k1 the first concentration parameter
     * @param k2 the second concentration parameter
     * @param normC the normalization constant
     */
    public BinghamDistribution(Vector3D[] axes, double k1, double k2, double normC, Random random) {
	
	this.k1 = k1;
	this.k2 = k2;

	// sanity check
	if (k2 > 0.0 || k1 > k2) {
	    throw new IllegalArgumentException("Can't have Bingham parameters > 0 or k1 > k2");
	}
        
        for (int i = 0; i < 3; i++) {
            if (Math.abs(1.0 - axes[i].mod()) > 1E-6) {
                throw new java.lang.IllegalArgumentException("Principal axes of Bingham distribution must be a unit vector. Magnitude of axis " + i + " was " + axes[i].mod());
            }
            
        }

        // check orthogonality
        if (Math.abs( axes[0].dot(axes[1]) ) > 1E-6) {
            throw new IllegalArgumentException("Eigenvectors must be orthogonal. Dot product between e3 and e2 was " + axes[0].dot(axes[1]) );
        }
        if (Math.abs( axes[0].dot(axes[2]) ) > 1E-6) {
            throw new IllegalArgumentException("Eigenvectors must be orthogonal. Dot product between e3 and e1 was " + axes[0].dot(axes[2]) );
        }


	e1 = axes[2];
	e2 = axes[1];
	e3 = axes[0];

	this.normC = normC;
	logNormC = Math.log(normC);
	
	ran = random;

	double wrapperNormC;

	if (Math.abs(k2) > 0.02 * Math.abs(k1) + 5) {
	    wrapper = new WatsonDistribution(e3, random);
	    wrapperK = -1.0 * k2;
	    wrapperNormC = normC * Math.exp(-1.0 * k2);
	}
	else {
	    wrapper = new WatsonDistribution(e1, random);
	    wrapperK = k1;
	    wrapperNormC = normC;
	}

	wrapperPDFNormC = Math.exp(WatsonDistribution.logHyper1F1(0.5, 1.5, wrapperK, 1e-9) - Math.log(wrapperNormC));

    }


    /**
     * @param kappas {k1, k2}.
     */
    public static BinghamDistribution getBinghamDistribution(Vector3D[] axes, 
							     double[] kappas, Random random) throws 
								 ConvergenceException {
	
	return getBinghamDistribution(axes, kappas[0], kappas[1], random);
    }


    public static BinghamDistribution getBinghamDistribution(Vector3D[] axes, 
							     double k1, double k2, Random random) throws 
								 ConvergenceException {
	double bingc = 0.0;

	bingc = BinghamFitter.bingc(k1, k2);
	
	return new BinghamDistribution(axes, k1, k2, bingc, random);
    }





    /**
     * Evaluate the pdf of the distribution.
     *
     * @param x a unit vector on the unit sphere.
     */
    public double pdf(Vector3D x) {
	double xDotE1 = x.dot(e1);
	double xDotE2 = x.dot(e2);
	
	return Math.exp( k1 * xDotE1 * xDotE1 + k2 * xDotE2 * xDotE2) / normC;
    }


    /**
     * Evaluate the pdf of the distribution.
     *
     * @param x a unit vector on the unit sphere.
     */
    public double logPDF(Vector3D x) {
	double xDotE1 = x.dot(e1);
	double xDotE2 = x.dot(e2);
	
	return k1 * xDotE1 * xDotE1 + k2 * xDotE2 * xDotE2 - logNormC;
    }

    /**
     * Evaluate the pdf of the distribution.
     *
     * @param x a unit vector on the unit sphere.
     */
    public static double pdf(Vector3D[] axes, double[] bingPars, Vector3D x) {
	double xDotE1 = x.dot(axes[2]);
	double xDotE2 = x.dot(axes[1]);
	
	return Math.exp( bingPars[0] * xDotE1 * xDotE1 + bingPars[1] * xDotE2 * xDotE2) / bingPars[2];
    }


    /**
     * Evaluate the pdf of the distribution.
     *
     * @param x a unit vector on the unit sphere.
     */
    private static double logPDF(Vector3D[] axes, double k1, double k2, double logNormC, Vector3D x) {
	double xDotE1 = x.dot(axes[2]);
	double xDotE2 = x.dot(axes[1]);
	
	return k1 * xDotE1 * xDotE1 + k2 * xDotE2 * xDotE2 - logNormC;
    }



    /**
     *
     * @return the next vector from the distribution.
     * 
     */
    public Vector3D nextVector() {

	int rejectionCounter = 0;

        if (k1 == 0.0) {
            // uniform distribution
	    double phi = Math.PI * 2.0 * ran.nextDouble();
	    double theta = Math.acos(2.0 * ran.nextDouble() - 1.0);
	    return Vector3D.vectorFromSPC(1.0, theta, phi);	    

        } 

	while (true) {
	    
	    Vector3D v = wrapper.nextVector(wrapperK);

	    double x = ran.nextDouble();

	    if (wrapperK < 700.0) {

		double wrapperPDF = 0.0;

		wrapperPDF = wrapper.pdf(v, wrapperK) * wrapperPDFNormC;

		if (x <= pdf(v) / wrapperPDF) {
		    return v;
		}
		else {
		    rejectionCounter++;
		    if (rejectionCounter % 10000 == 0) {
			System.err.println("Still rejecting samples after " + rejectionCounter + " attempts"); 
			throw new IllegalStateException("Couldn't get nextVector.\nConcentration\n\t" + k1 + "\t" + k2 + "\nAxes\n" + e1 + "\n" + e2 + "\n" + e3);
		    }
		}


	    }
	    else {
		double logWrapperPDF = wrapper.logPDF(v, wrapperK) + Math.log(wrapperPDFNormC);

		if (Math.log(x) <= logPDF(v) - logWrapperPDF) {
		    return v;
		}
		else {
		    rejectionCounter++;
		    if (rejectionCounter % 10000 == 0) {
			System.err.println("Still rejecting samples after " + rejectionCounter + " attempts"); 

			throw new IllegalStateException("Couldn't get nextVector.\nConcentration\n\t" + k1 + "\t" + k2 + "\nAxes\n" + e1 + "\n" + e2 + "\n" + e3);
		    }
		}
		    
	    }
	    
	}
    }


    /**
     *
     * @return the next vector from the distribution. 
     * 
     */
    public static Vector3D nextVector(Vector3D[] axes, double k1, double k2, double logNormC, Random ran) {

	int rejectionCounter = 0;

	// sanity check
	if (k1 >= 0.0 || k2 >= 0.0 || k1 > k2) {
	    throw new IllegalArgumentException("Can't have Bingham parameters >= 0 or k1 > k2");
	}
	
	Vector3D e1 = axes[2];
	Vector3D e2 = axes[1];
	Vector3D e3 = axes[0];

	double logWrapperNormC;
	double wrapperK;
	Vector3D wrapperMu;

	if (Math.abs(k2) > 0.02 * Math.abs(k1) + 5) {
	    wrapperMu = e3;
	    wrapperK = -1.0 * k2;
	    logWrapperNormC = logNormC + -1.0 * k2;
	}
	else {
	    wrapperMu = e1;
	    wrapperK = k1;
	    logWrapperNormC = logNormC;
	}

	while (true) {
	    
	    Vector3D v = WatsonDistribution.nextVector(wrapperMu, wrapperK, ran);

	    double x = ran.nextDouble();

	    double logWrapperPDF = logWrapperPDF(v, wrapperMu, wrapperK, logWrapperNormC);

		if (Math.log(x) <= logPDF(axes, k1, k2, logNormC, v) - logWrapperPDF) {
		    return v;
		}
		else {
		    rejectionCounter++;
		    if (rejectionCounter % 10000 == 0) {
			System.out.println("Still rejecting samples after " + rejectionCounter / 1000 + "0000 attempts"); 
		    }
		}
		
	}
	
    }
    
    private static double logWrapperPDF(Vector3D v, Vector3D wrapperMu, double wrapperK, double logWrapperNormC) {

	double vDotMu = v.dot(wrapperMu);
	return vDotMu * vDotMu * wrapperK - logWrapperNormC;
    }

    /**
     * Gets the concentration parameter k1 where k1 <= k2 <= 0. 
     */
    public double k1() {
	return k1;
    }

    /**
     * Gets the concentration parameter k2 where k1 <= k2 <= 0.
     */
    public double k2() {
	return k2;
    }

    /**
     * Eigenvector corresponding to the smallest eigenvalue of T. This is the opposite ordering to the
     * diffusion tensor convention, but we use this ordering to be consistent with Bingham's paper.
     * 
     */
    public Vector3D e1() {
	return e1;
    }

    /**
     *
     * Eigenvector corresponding to the second eigenvalue of T. This is the opposite ordering to the
     * diffusion tensor convention, but we use this ordering to be consistent with Bingham's paper.
     */
    public Vector3D e2() {
	return e2;
    }

    /**
     *
     * Eigenvector corresponding to the largest eigenvalue of T. This is the opposite ordering to the
     * diffusion tensor convention, but we use this ordering to be consistent with Bingham's paper.
     *
     */
    public Vector3D e3() {
	return e3;
    }


    /**
     * @return the normalization constant in the PDF.
     */
    public double normC() {
	return normC;
    }
    
    /**
     * @return the normalization constant in the PDF.
     */
    public double logNormC() {
	return logNormC;
    }

    public String toString() {
	return "Bingham Distribution, k_1 == " + k1 + " k2 == " + k2;
    }


}

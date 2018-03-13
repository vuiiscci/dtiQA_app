package numerics;

import java.util.Random;

/**
 * Watson Distribution, kind of like a Fisher but for axial data, ie p(x) == p(-x). 
 *
 * @see section 9.4.2 of Mardia and Jupp ("Directional Statistics", Wiley, 2000)
 * @version $Id$
 * @author Philip Cook
 */
public final class WatsonDistribution implements AxialDistribution {

    private final Vector3D mu;
    private final double kappa;
    
    private final double muPhi;

    private final double muTheta;

    // for the computation of hyper1F1. 
    private final static double asympThresh = 28.0;

    private final Random ran;


    // don't want anyone to use this
    private WatsonDistribution() {

	mu = null;
	muPhi = 0.0;
	muTheta = 0.0;
	kappa = 0.0;
	ran = null;
	
    }

    /**
     * @param mean the mean direction of a bipolar distribution 
     * (or polar direction of a girdle distribution).
     */
    public WatsonDistribution(Vector3D mean, Random r) {

	mu = mean;
	ran = r;

	if (Math.abs(1.0 - mean.mod()) > 0.000001) {
	    throw new java.lang.IllegalArgumentException("Mean axis of Watson distribution must be a unit vector. Magnitude of mean was " + mean.mod());
	}
	
	double[] thetaPhi = Vector3D.thetaPhi(mu);

	muTheta = thetaPhi[0];
	muPhi = thetaPhi[1];

	kappa = 0.0;

    }


    /**
     * @param mean the mean direction of a bipolar distribution 
     * (or polar direction of a girdle distribution).
     * @param k the concentration that will be used if the no-arg nextVector() is called.
     */
    public WatsonDistribution(Vector3D mean, double k, Random r) {

	mu = mean;
	ran = r;

	if (Math.abs(1.0 - mean.mod()) > 0.000001) {
	    throw new java.lang.IllegalArgumentException("Mean axis of Watson distribution must be a unit vector. Magnitude of mean was " + mean.mod());
	}
	
	double[] thetaPhi = Vector3D.thetaPhi(mu);

	muTheta = thetaPhi[0];
	muPhi = thetaPhi[1];

	kappa = k;

        if (Double.isInfinite(k) || Double.isNaN(k)) {
            throw new IllegalArgumentException("kappa is infinite or NaN");
        }

    }



    /**
     * Evaluate the pdf of the distribution. This integrates to 4.0 * PI over the sphere.
     *
     * @param x a unit vector on the unit sphere.
     */
    public double pdf(Vector3D x) {
	return pdf(mu, x, kappa);
    }

    /**
     * Evaluate the pdf of the distribution.
     *
     * @param x a unit vector on the unit sphere.
     * @return the log of the PDF.
     */
    public double logPDF(Vector3D x) {
	return logPDF(mu, x, kappa);
    }


    /**
     * Evaluate the pdf of the distribution. This integrates to 4.0 * PI over the sphere.
     *
     * @param x a unit vector on the unit sphere.
     * @param k the concentration parameter of the distribution.
     */
    public double pdf(Vector3D x, double k) {
	return pdf(mu, x, k);
    }


    /**
     * Evaluate the pdf of the distribution.
     *
     * @param x a unit vector on the unit sphere.
     * @param k the concentration parameter of the distribution.
     *
     */
    public double logPDF(Vector3D x, double k) {
	return logPDF(mu, x, k);
    }

    /**
     * Evaluate the pdf of any Watson distribution. This integrates to 4.0 * PI over the sphere.
     *
     * @param mean the mean axis of the distribution.
     * @param x a unit vector on the unit sphere.
     * @param k the concentration parameter of the distribution.
     *
     */
    public static double pdf(Vector3D mean, Vector3D x, double k) {

	if (k < 700.0) {
	    double muDotX = mean.dot(x);
	    return Math.exp( k * muDotX * muDotX ) / hyper1F1(0.5, 1.5, k, 1.0e-9);
	}
	else {
	    return Math.exp( logPDF(mean, x, k) ); 
	}

    }


    /**
     * Evaluate the pdf of any Watson distribution. 
     *
     * @param mean the mean axis of the distribution.
     * @param x a unit vector on the unit sphere.
     * @param k the concentration parameter of the distribution.
     *
     */
    public static double logPDF(Vector3D mean, Vector3D x, double k) {

	double muDotX = mean.dot(x);

	return k * muDotX * muDotX - logHyper1F1(0.5, 1.5, k, 1.0e-9);

    }


    /**
     * Calculate the log likelihood of a set of axes.
     *
     * @param mean the mean axis of the distribution.
     * @param x an array of unit vectors.
     * @param k the concentration parameter of the distribution.
     *
     */
    public static double sumLogPDF(Vector3D mean, Vector3D[] x, double k) {
       
	double l = 0.0;
	double logHyper1F1 = logHyper1F1(0.5, 1.5, k, 1.0e-9);
	for (int i = 0; i < x.length; i++) {

	    double muDotX = mean.dot(x[i]);

	    l += k * muDotX * muDotX - logHyper1F1;
	}

	return l;

    }



   /**
     * @return the axis that defines the distribution.
     */
    public Vector3D meanAxis() {
	return mu;
    }


    /**
     * @return the scalar parameter that defines the concentration.
     */
    public double concentration() {
	return kappa;
    }

    
    /**
     * @return the next vector from the distribution.
     */
    public Vector3D nextVector() {
	return nextVector(kappa);
    }


    /**
     * Procedure adapted from Statistical Analysis of Spherical Data, by Fisher, Lewis, Embleton, Cambridge University Press, 1993.
     *
     *
     * @param k the concentration parameter.
     * @return the next vector from the distribution.
     */
    public Vector3D nextVector(double k) {
	//	return nextVector(muTheta, muPhi, k, ran);
	// the extra method call slows the process of getting a vector by about 10%.
	// we can make additional optimization because we know mu

        if (Double.isInfinite(k) || Double.isNaN(k)) {
            throw new IllegalArgumentException("kappa is infinite or NaN");
        }

	double theta = 0.0;
	double phi = 0.0;
	Vector3D sample;
	
	if (k < 0.0) {
	    // girdle distribution

	    // 0

	    // (-1 * k) is +ve
	    double c1 = Math.sqrt(-1.0 * k);
	    double c2 = Math.atan(c1);

	    double s = 0.0;

	    while (true) {
	    
		double u = ran.nextDouble();
		double v = ran.nextDouble();
	    
		s = (1.0 / c1) * Math.tan(c2 * u);

		if ( v <= (1.0 - k * s * s) * Math.exp(k * s * s) ) {
		    break;
		}
                
            }
	
	    theta = muTheta + Math.acos(s);
	    phi = muPhi;

	    // theta wraps at PI
 	    if (theta > Math.PI) {
		
 		theta = theta - Math.PI;
		
 	    }
	
	    sample = Vector3D.vectorFromSPC(1.0, theta, phi);
	
	    if (ran.nextDouble() > 0.5) {
		sample = sample.negated();
	    }

        
	}
	else if (k == 0.0) {
	    // uniform distribution 
	    phi = Math.PI * 2.0 * ran.nextDouble();
	    theta = Math.acos(2.0 * ran.nextDouble() - 1.0);
	    return Vector3D.vectorFromSPC(1.0, theta, phi);	    
	}
	else if (k > 700.0) {
	    
	    double lnc = -1.0 * k;
	    
	    double s = 0.0;
	    
	    while (true) {
		
		// 1
		double u = ran.nextDouble(); 
		double v = ran.nextDouble(); 
		
		// 2
		// was    s = (1.0 / k) * Math.log(1.0 + u / c);
		
		// log (1 + u / c) -----> log(u / c) for (u / c) very large
		s = (1.0 / k) * (Math.log(u) - lnc);
		
		// 3
		if ( v <= Math.exp(k * s * (s - 1.0)) ) {
		    break;
		}
	
	    }
	    
	    theta = muTheta + Math.acos(s);
	    phi = muPhi;

	    // theta wraps at PI
	    if (theta > Math.PI) {
		
 		theta = theta - Math.PI;
		
 	    }

	    sample = Vector3D.vectorFromSPC(1.0, theta, phi);
	    
	    if (ran.nextDouble() > 0.5) {
		sample = sample.negated();
	    }

	    
	}
	else {
	    // bipolar distribution
	    
	    // 0
	    double c = 1.0 / (Math.exp(k) - 1.0); 
	    
	    double s = 0.0;
	    
	    while (true) {
		
		// 1
		double u = ran.nextDouble(); 
		double v = ran.nextDouble(); 
		
		// 2
		s = (1.0 / k) * Math.log(1.0 + u / c);
		
		// 3
		if ( v <= Math.exp(k * s * (s - 1.0)) ) {
		    break;
		}

	    }
	    
	    theta = muTheta + Math.acos(s);
	    phi = muPhi;
	    
	    // theta wraps at PI
	    
 	    if (theta > Math.PI) {
		
 		theta = theta - Math.PI;
		
 	    }

	    sample = Vector3D.vectorFromSPC(1.0, theta, phi);
	    
	    // because p(x) == p(-x)
	    if (ran.nextDouble() > 0.5) {
		sample = sample.negated();
	    }
	    
	}
	 
	// Now rotate sample about the mean axis by random angle gamma, between 0 and 2 * PI
	// this part makes the distribution cylindrically symmetric
	double gamma = Math.PI * 2.0 * ran.nextDouble();
	
	
	// sample
	double cg = Math.cos(gamma);
	double sg = Math.sin(gamma);
	
	double sampX = 
	    (mu.x*mu.x*(1.0 - cg) + cg) * sample.x + 
	    (mu.x*mu.y*(1.0 - cg) - mu.z*sg) * sample.y + 
	    (mu.x*mu.z*(1.0 - cg) + mu.y*sg) * sample.z;
	
	double sampY = 
	    (mu.y*mu.x*(1.0 - cg) + mu.z*sg) * sample.x + 
	    (mu.y*mu.y*(1.0 - cg) + cg) * sample.y + 
	    (mu.y*mu.z*(1.0 - cg) - mu.x*sg) * sample.z;  
	
	double sampZ = 
	    (mu.z*mu.x*(1.0 - cg) - mu.y*sg) * sample.x + 
	    (mu.z*mu.y*(1.0 - cg) + mu.x*sg) * sample.y + 
	    (mu.z*mu.z*(1.0 - cg) + cg) * sample.z;
	
	return new Vector3D(sampX, sampY, sampZ);
    }


    /**
     * @return the next vector from the distribution from arbitrary mean.
     */
    public static Vector3D nextVector(double meanTheta, double meanPhi, double k, Random ran) {

	Vector3D sample;

        if (Double.isInfinite(k) || Double.isNaN(k)) {
            throw new IllegalArgumentException("kappa is infinite or NaN");
        }


	if (k < 0.0) {
	    // girdle distribution

	    // 0
	    // (-1 * k) is +ve
	    double c1 = Math.sqrt(-1.0 * k);
	    double c2 = Math.atan(c1);

	    double s = 0.0;

	    while (true) {
	    
		double u = ran.nextDouble();
		double v = ran.nextDouble();
	    
		s = (1.0 / c1) * Math.tan(c2 * u);

		if ( v <= (1.0 - k * s * s) * Math.exp(k * s * s) ) {
		    break;
		}
	    
	    
	    }
	
	    double theta = meanTheta + Math.acos(s);

	    double phi = meanPhi;

	    // theta wraps at PI
	    if (theta > Math.PI) {
		theta = theta - Math.PI;
	    }
	    
	    sample = Vector3D.vectorFromSPC(1.0, theta, phi);
	
	    if (ran.nextDouble() > 0.5) {
		sample = sample.negated();
	    }


	}
	else if (k == 0.0) {
	    // uniform distribution 
	    double phi = Math.PI * 2.0 * ran.nextDouble();
	    double theta = Math.acos(2.0 * ran.nextDouble() - 1.0);
	    return Vector3D.vectorFromSPC(1.0, theta, phi);	    
	}
	else if (k > 700.0) {
	    return nextVectorLargeK(meanTheta, meanPhi, k, ran);
	}
	else {
	    // bipolar distribution
	    
	    // 0
	    double c = 1.0 / (Math.exp(k) - 1.0); 
	    
	    double s = 0.0;
	    
	    while (true) {
		
		// 1
		double u = ran.nextDouble(); 
		double v = ran.nextDouble(); 
		
		// 2
		s = (1.0 / k) * Math.log(1.0 + u / c);
		
		// 3
	    if ( v <= Math.exp(k * s * (s - 1.0)) ) {
		break;
	    }
	    
	    }
	    
	    double theta = meanTheta + Math.acos(s);
	    
	    double phi = meanPhi;
	    
	    // theta wraps at PI
	    if (theta < 0.0 || theta > Math.PI) {
		phi = Math.PI + phi;
	    }
	    
	    sample = Vector3D.vectorFromSPC(1.0, theta, phi);
	    
	    // because p(x) == p(-x)
	    if (ran.nextDouble() > 0.5) {
		sample = sample.negated();
	    }
	    

	}
	 
	// Now rotate sample about the mean axis by random angle gamma, between 0 and 2 * PI
	// this part makes the distribution cylindrically symmetric

	double gamma = Math.PI * 2.0 * ran.nextDouble();

	double cg = Math.cos(gamma);
	double sg = Math.sin(gamma);

	Vector3D mu = Vector3D.vectorFromSPC(1.0, meanTheta, meanPhi);
	
	double sampX = 
	    (mu.x*mu.x*(1.0 - cg) + cg) * sample.x + 
	    (mu.x*mu.y*(1.0 - cg) - mu.z*sg) * sample.y + 
	    (mu.x*mu.z*(1.0 - cg) + mu.y*sg) * sample.z;
	
	double sampY = 
	    (mu.y*mu.x*(1.0 - cg) + mu.z*sg) * sample.x + 
	    (mu.y*mu.y*(1.0 - cg) + cg) * sample.y + 
	    (mu.y*mu.z*(1.0 - cg) - mu.x*sg) * sample.z;  
	
	double sampZ = 
	    (mu.z*mu.x*(1.0 - cg) - mu.y*sg) * sample.x + 
	    (mu.z*mu.y*(1.0 - cg) + mu.x*sg) * sample.y + 
	    (mu.z*mu.z*(1.0 - cg) + cg) * sample.z;
	
	return new Vector3D(sampX, sampY, sampZ);


    }


    public Vector3D nextVector(Vector3D mean, double k) {
	// interface requires non static method
	double[] thetaPhi = Vector3D.thetaPhi(mean);
	return nextVector(thetaPhi[0], thetaPhi[1], k, ran);

    }


    public static Vector3D nextVector(Vector3D mean, double k, Random r) {
	// interface requires non static method
	double[] thetaPhi = Vector3D.thetaPhi(mean);
	return nextVector(thetaPhi[0], thetaPhi[1], k, r);

    }

    
    /**
     * @return the next vector from the distribution.
     * Approximates a computable solution for k large.
     */
    private static Vector3D nextVectorLargeK(double meanTheta, double meanPhi, double k, Random ran) {
	// was	double c = 1.0 / (Math.exp(k) - 1.0); 
	
	// (\exp{\kappa} - 1.0} ---> \exp{\kappa} for large k

	double lnc = -1.0 * k;

	double s = 0.0;

	while (true) {

	    // 1
	    double u = ran.nextDouble(); 
	    double v = ran.nextDouble(); 
	    
	    // 2
	    // was    s = (1.0 / k) * Math.log(1.0 + u / c);
	    
	    // log (1 + u / c) -----> log(u / c) for (u / c) very large
	    s = (1.0 / k) * (Math.log(u) - lnc);
	    
	    // 3
	    if ( v <= Math.exp(k * s * (s - 1.0)) ) {
		break;
	    }
	    
	}

	double theta = meanTheta + Math.acos(s);

	double phi = meanPhi;

	// theta wraps at PI
	if (theta < 0.0 || theta > Math.PI) {
	    phi = phi + Math.PI;
	}

	Vector3D sample = Vector3D.vectorFromSPC(1.0, theta, phi);

 	if (ran.nextDouble() > 0.5) {
 	    sample = sample.negated();
 	}
    

	double gamma = Math.PI * 2.0 * ran.nextDouble();

	double cg = Math.cos(gamma);
	double sg = Math.sin(gamma);

	Vector3D mu = Vector3D.vectorFromSPC(1.0, meanTheta, meanPhi);
	
	double sampX = 
	    (mu.x*mu.x*(1.0 - cg) + cg) * sample.x + 
	    (mu.x*mu.y*(1.0 - cg) - mu.z*sg) * sample.y + 
	    (mu.x*mu.z*(1.0 - cg) + mu.y*sg) * sample.z;
	
	double sampY = 
	    (mu.y*mu.x*(1.0 - cg) + mu.z*sg) * sample.x + 
	    (mu.y*mu.y*(1.0 - cg) + cg) * sample.y + 
	    (mu.y*mu.z*(1.0 - cg) - mu.x*sg) * sample.z;  
	
	double sampZ = 
	    (mu.z*mu.x*(1.0 - cg) - mu.y*sg) * sample.x + 
	    (mu.z*mu.y*(1.0 - cg) + mu.x*sg) * sample.y + 
	    (mu.z*mu.z*(1.0 - cg) + cg) * sample.z;
	
	return new Vector3D(sampX, sampY, sampZ);
	

    }


    /**
     * Computes the confluent hypergeometric function of the first kind, M(a,b,x). 
     * @param eps the maximum fractional error. The function will return when the series
     * has converged to within this amount.
     *
     * This is taken, with some modification from Thompson, "Atlas for computing mathematical functions:
     * an illustrated guide for practitioners, with programs in C and Mathematica" (Wiley, 1997)
     *
     */
    public static double hyper1F1(double a, double b, double x, double eps) {

	double sum, ratio, aks, bks, kfact, tk, xpow, eps10, term, fk, bMaPkM1, oneMPkM1;
	double termold, termnew, aPkM1, aMbPk, mx;

	int k;

	if (a == b) {
	    return Math.exp(x);
	}

	if (b <= 0.0) {
	    System.err.println("hyper1F1: b <= 0.0, 0 returned");
	    return 0.0;
	}

	eps10 = 0.1 * eps;

	if (Math.abs(x) < asympThresh) {

	    sum = 1.0;
	    ratio = 10.0;
	    aks = a;
	    bks = b;
	    kfact = 1.0;
	    tk = 1.0;
	    xpow = 1.0;

	    while (ratio > eps10) {
		
		tk *= aks / (bks * kfact);
		xpow *= x;
		
		term = tk * xpow;

		sum += term;
		aks += 1.0;
		bks += 1.0;
		kfact += 1.0;

		ratio = Math.abs(term / sum);
		
	    }


	    return sum;
	    
	} // end if |x| < asympThresh
	

	// use asymptotic series
	if (x >= asympThresh) {
	    if (a <= 0.0) {
		System.err.println("a <= 0.0 and x >= " + asympThresh + ". 0 returned");
		return 0.0;
	    }


	    k = 0; 
	    fk = 1.0;
	    bMaPkM1 = b - a;
	    oneMPkM1 = 1.0 - a;
	    term = 1.0;
	    termold = 100.0;
	    sum = 1.0;
	    ratio = 10.0;

	    while (ratio > eps10) {
		k++;
		term *= bMaPkM1 * oneMPkM1 / (fk * x);
		sum += term;
		bMaPkM1 += 1.0;
		oneMPkM1 += 1.0;
		fk += 1.0;
		ratio = Math.abs(term / sum);

		termnew = Math.abs(term);
		
		if (termnew > termold) {
		    System.err.println("Series diverging after " + k + " terms");
		    ratio = 0.0; // returns
		}
		else {
		    termold = termnew;
		}

	    } // end while

	    return Math.exp(GammaFunctions.gammln(b) + x - GammaFunctions.gammln(a)) * Math.pow(x, a - b) * sum;


	} // end if 


	// x <= -asympThresh
	
	if (b <= a) {
	    System.err.println("b <= a, x <= -" + asympThresh + ". 0 returned");
	    return 0.0;
	}

	
	mx = -1.0 * x;
	k = 0;
	fk = 1.0;
	aMbPk = a - b + 1.0;
	aPkM1 = a;
	term = 1.0;
	termold = 100.0;
	sum = 1.0;
	ratio = 10.0;

	while (ratio < eps10) {
	    
	    k++;
	    term *= aPkM1 * aMbPk / (fk * mx);
	    sum += term;
	    aPkM1 += 1.0;
	    aMbPk += 1.0;
	    fk += 1.0;

	    ratio = Math.abs(term / sum);

	    termnew = Math.abs(term);

		
	    if (termnew > termold) {
		System.err.println("Series diverging after " + k + " terms");
		ratio = 0.0; // returns
	    }
	    else {
		termold = termnew;
	    }

	    
	    
	} // end while
	
	return Math.exp( GammaFunctions.gammln(b) - GammaFunctions.gammln(b - a) ) * sum / Math.pow(mx, a);



    }



    /**
     * @return log of the hypergeometric series. Use for large x.
     * 
     *
     */
    public static double logHyper1F1(double a, double b, double x, double eps) {
    	
	double sum, ratio, tk, eps10, term, fk, bMaPkM1, oneMPkM1;
	double termold, termnew, aPkM1, aMbPk, mx;
	
	int k;

	eps10 = 0.1 * eps;
	
	// use asymptotic series
	if (x <= asympThresh) {
	    return Math.log(hyper1F1(a, b, x, eps));
	}
	
	if (a <= 0.0) {
	    System.err.println("a <= 0.0 and x >= " + asympThresh + ". 0 returned");
	    return 0.0;
	}


	k = 0; 
	fk = 1.0;
	bMaPkM1 = b - a;
	oneMPkM1 = 1.0 - a;
	term = 1.0;
	termold = 100.0;
	sum = 1.0;
	ratio = 10.0;

	while (ratio > eps10) {
	    k++;

	    term *= bMaPkM1 * oneMPkM1 / (fk * x);

	    sum += term;
	    bMaPkM1 += 1.0;
	    oneMPkM1 += 1.0;
	    fk += 1.0;

	    ratio = Math.abs(term / sum);
		
	    termnew = Math.abs(term);
		
	    if (termnew > termold) {
		System.err.println("Series diverging after " + k + " terms");
		ratio = 0.0; // returns
	    }
	    else {
		termold = termnew;
	    }

	} // end while
	    
	return (GammaFunctions.gammln(b) + x - GammaFunctions.gammln(a)) + Math.log(Math.pow(x, a - b)) + Math.log(sum);


    }









}

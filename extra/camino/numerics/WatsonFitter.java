package numerics;



import java.util.Arrays;

/**
 * This class implements several routines for fitting and evaluating 
 * goodness-of-fit to a Watson distribution. 
 * <p> The routines come from "Statistical Analysis of spherical data" by Fisher, Lewis, Embleton (1987). 
 *
 * @author Philip Cook
 * @version $Id$
 */
public class WatsonFitter extends SphericalDistributionFitter {


    /**
     * Test of rotational symmetry for a bipolar distribution.
     *
     * @param eig the eigensystem of the normalized sample scatter matrix.
     * @param mu the mean axis
     * @param sampleVecs the sample axes
     * @see Fisher (p. 165)
     */
   public static double testBipolarRotSymm(EigenSystem3D eig, Vector3D mu, Vector3D[] sampleVecs) {
	double gamma = 0.0;
	
	double n = (double)sampleVecs.length;

	for (int i = 0; i < sampleVecs.length; i++) {
	    
	    gamma += (1.0 / n) * Math.pow(mu.dot(sampleVecs[i]), 4.0);
	    
	}
	
	double pn =  2.0 * n * Math.pow(eig.eigenvalues[1] - eig.eigenvalues[2], 2.0) / 
	    ( 1.0 - 2.0 * eig.eigenvalues[0] + gamma);


	return pn;


    }


    /**
     * Get the semi-vertical angle of the 100(1 - \alpha)% confidence cone around \mu. 
     *
     * @param eig the eigensystem of the normalized sample scatter matrix. 
     * The mean axis is presumed to be the first eigenvector, so eig should be sorted.
     * @param sampleVecs the sample axes.
     * @param alpha significance == 100(1 - alpha)%. For a 95% cone, alpha == 0.05.
     * @see Fisher (p. 165)
     */
   public static double getBipolarConfidenceCone(EigenSystem3D eig, Vector3D[] sampleVecs, double alpha) {
	double gamma = 0.0;

	Vector3D mu = eig.eigenvectors[0];

	double n = (double)sampleVecs.length;

	for (int i = 0; i < sampleVecs.length; i++) {
	    
	    gamma += (1.0 / n) * Math.pow(mu.dot(sampleVecs[i]), 4.0);
	    
	}
	
	if ((eig.eigenvalues[0] < gamma)) {
	    gamma = eig.eigenvalues[0];
	}

	double sigma = Math.sqrt( (eig.eigenvalues[0] - gamma) / n ) / 
	    (eig.eigenvalues[0] - eig.eigenvalues[1]);

	double qArg = Math.sqrt(-Math.log(alpha)) * sigma;

	if (Math.abs(qArg) > 1.0) {
	    qArg = 1.0;
	}

	// something wrong here at high concentration
	double q = Math.asin( qArg );

	return q;

    }


    /**
     * Test rotational symmetry for a girdle distribution.
     *
     * @param eig the eigensystem of the normalized sample scatter matrix.
     * @param mu the polar axis.
     * @param sampleVecs the sample axes.
     *
     * @see Fisher (p. 182)
     */
    public static double testGirdleRotSymm(EigenSystem3D eig, Vector3D mu, Vector3D[] sampleVecs) {
	double gamma = 0.0;

	double n = (double)sampleVecs.length;
	
	for (int i = 0; i < sampleVecs.length; i++) {
	    
	    gamma += (1.0 / n) * Math.pow(mu.dot(sampleVecs[i]), 4.0);
	    
	}
	
	return 2.0 * n * Math.pow(eig.eigenvalues[0] - eig.eigenvalues[1], 2.0) / 
	    ( 1.0 - 2.0 * eig.eigenvalues[2] + gamma);
    }


    /**
     * @param mu the polar / mean axis.
     * @param sampleVecs the sample axes.
     * @return modified kuiper test statistic.     
     *
     * @see Fisher (p. 56)
     */
    public static double kuiperTest(Vector3D mu, Vector3D[] sampleVecs) {
	
	double vn = 0.0;

	double n = (double)(sampleVecs.length);
	
	double dnPlus = Double.MIN_VALUE;
	double dnMinus = Double.MIN_VALUE;
	
	double[] xi = new double[sampleVecs.length];
	
	// rotation of mean vector onto z axis
	RealMatrix rot = Rotations.getRotMat(mu, new Vector3D(0.0, 0.0, 1.0));

	for (int i = 0; i < sampleVecs.length; i++) {
	    Vector3D rotated = Rotations.rotateVector(sampleVecs[i], rot);

	    double phi = 0.0; 

	    if (rotated.x == 0.0) {
		phi = rotated.y > 0.0 ? Math.PI / 2.0 : -1.0 * Math.PI / 2.0;
	    }
	    else {
		if (rotated.x < 0.0) {
		    phi = Math.PI + Math.atan(rotated.y / rotated.x);
		}
		else {
		    phi = Math.atan(rotated.y / rotated.x);
		}
	    }

	    xi[i] = (phi / (2.0 * Math.PI));
	    
	    
	}
	
	Arrays.sort(xi);

	double[] fxi = new double[xi.length];

	for (int i = 0; i < sampleVecs.length; i++) {
	    fxi[i] = xi[i];
	}


	for (int i = 0; i < sampleVecs.length; i++) {
	    
	    double dnPn = (i + 1.0) / n - fxi[i];
	    
	    if (dnPn > dnPlus) {
		dnPlus = dnPn; 
	    }

	    double dnMn = fxi[i] - i / n; 
	    
	    if (dnMn > dnMinus) {
		dnMinus = dnMn; 
	    }
	    
	    
	}
	    
	vn = dnPlus + dnMinus;

	return vn * (Math.sqrt(n) + 0.155 + 0.24 / Math.sqrt(n));
    }




    /**
     *
     * Fits kappa for a sample of vectors. The method fits both a girdle and a bipolar 
     * distribution, and then tests which is best.
     *
     * @return the solution that maximises the likelihood of the samples
     *
     * @see Fisher (p.176, 189)
     */
    public static double fitKappa(Vector3D[] sampleVecs) {
	EigenSystem3D eig = tBarEigenSystem(sampleVecs);
	return fitKappa(eig, sampleVecs);
    }



    /**
     *
     * Fits kappa for a sample of vectors. The method fits both a girdle and a bipolar 
     * distribution, and then tests which is best.
     *
     * @return the solution that maximises the likelihood of the samples
     *
     * @see Fisher (p.176, 189)
     */
    public static double fitKappa(EigenSystem3D eig, Vector3D[] sampleVecs) {


	// choose the one that maximises the log likelihood 
	WatsonD d3 = new WatsonD();

	double bipolarKappa = 0.0;
	double girdleKappa = 0.0;
		
	boolean bipolarConverged = false;
	boolean girdleConverged = false;

	
	try {

	    bipolarKappa = NewtonRaphsonSolver.solve(d3, eig.eigenvalues[0], 0.0, 0.00001);
	    bipolarConverged = true;
	}
	catch (ConvergenceException e) {
	    
	}
	try {
	    girdleKappa = NewtonRaphsonSolver.solve(d3, eig.eigenvalues[2], 0.0, 0.00001);
	    girdleConverged = true;
	}
	catch (ConvergenceException e) {
	}

	if (bipolarConverged && girdleConverged) {
	    
	    double logBipolarLikelihood = WatsonDistribution.sumLogPDF(eig.eigenvectors[0], sampleVecs, bipolarKappa);		
	    double logGirdleLikelihood = WatsonDistribution.sumLogPDF(eig.eigenvectors[2], sampleVecs, girdleKappa);		
	    
	    // force bipolar if t_1 > 0.95 (corresponds to k = 20)
	    if (logBipolarLikelihood > logGirdleLikelihood || eig.eigenvalues[0] > 0.95) {
		return bipolarKappa;
	    }
	    else {
		return girdleKappa;
	    }
	    
	}
	else if (bipolarConverged) {
	    return bipolarKappa;
	}
	else if (girdleConverged) {
	    return girdleKappa;
	}
	else {
	    return 0.0;
	}
	
    }


    
    /**
     *
     * Fits kappa to a scatter matrix. The method fits a bipolar 
     * distribution if |t_1 - t_2| > |t_2 - t_3|, and a girdle if
     * |t_1 - t_2| < |t_2 - t_3|
     *
     * @return kappa
     *
     * @see Fisher (p.176, 189)
     */
    public static double fitKappa(EigenSystem3D eig) {


	WatsonD d3 = new WatsonD();

	double kappa = 0.0;
		
	double ev = eig.eigenvalues[0];

	if ( Math.abs(eig.eigenvalues[0] - eig.eigenvalues[1]) < 
	     Math.abs(eig.eigenvalues[1] - eig.eigenvalues[2]) ) {
	    ev = eig.eigenvalues[2];
	}
	
	try {

	    kappa = NewtonRaphsonSolver.solve(d3, ev, 0.0, 0.00001);

	}
	catch (ConvergenceException e) {
	    return 0.0;
	}

	return kappa;
	
    }


    /**
     *
     * Fits kappa to an eigenvalue. 
     *
     * @param t should be t1 for a bipolar distribution about e1, or t3 for a girdle distribution about e3.
     * @see Fisher (p.176, 189)
     */
    public static double fitKappa(double t) throws ConvergenceException {

	WatsonD d3 = new WatsonD();

	double kappa = 0.0;

	kappa = NewtonRaphsonSolver.solve(d3, t, 0.0, 0.00001);
	return kappa;
	
    }

  


    /**
     * Calculate goodness according to the Kolmogorov-Smirnov test. 
     *
     * Doesn't work for kappa < 0. There is a disagreement between the Mardia book
     * and the Fisher book on how to do the test for the girdle distribution. Not sure
     * which is right; neither seems to work.
     *
     * @return the Kolmogorov-Smirnov test statistic. See Fisher for critical values.
     *
     * @see Fisher (p. 169)
     *
     *
     */
    public static double ksTest(Vector3D[] sampleVecs, Vector3D mu, double kappa) {

	// first we need Dn
	double dn = 0.0;
	
	// there are two possible values for dn
	double dnPoss1 = -Double.MAX_VALUE;
	double dnPoss2 = -Double.MAX_VALUE;


	int xiDim = sampleVecs.length;

	double n = (double)xiDim;

	double[] xi = new double[xiDim];

	if (kappa > 0.0) {
	    for (int i = 0; i < xiDim; i++) {
		xi[i] = 1.0 - Math.pow(sampleVecs[i].dot(mu), 2.0);
	    }
	}
	else {
	    for (int i = 0; i < xiDim; i++) {
		xi[i] = Math.pow(sampleVecs[i].dot(mu), 2.0);
	    }
	}
	    
	Arrays.sort(xi);
	

	for (int i = 0; i < xiDim; i++) {
		
	    double fx = 1.0 - Math.exp(-1.0 * kappa * xi[i]);

	    double dn1n = (i + 1.0) / n - fx;

	    if (dn1n > dnPoss1) {
		dnPoss1 = dn1n; 
	    }

	    double dn2n = fx - i / n;  

	    if (dn2n > dnPoss2) {
		dnPoss2 = dn2n; 
	    }
		    

	}
	       
	if (dnPoss1 >= dnPoss2) {
	    dn = dnPoss1;
	}
	else {
	    dn = dnPoss2;	    
	}
	 
	if (kappa > 0.0) {
	    return (dn - 0.2 / n) * (Math.sqrt(n) + 0.26 + 0.5 / Math.sqrt(n));
	}
	else {
	    return dn * (Math.sqrt(n) - 0.04 + 0.7 / Math.sqrt(n));
	}
  
   
    }






}

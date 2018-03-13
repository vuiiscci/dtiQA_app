package numerics;

import optimizers.*;

/**
 * <dl>
 *
 * <dt>Purpose: 
 *
 * <dd> Fits two Bingham distributions to a collection of axes.
 *
 * <dt>Description:
 *
 * <dd> This class takes a collection of axes representing a set of estimates of two fibre directions. 
 * It uses a Conjugate Gradient  optimization routine to find the pair of Bingham distributions that
 *  best fit the axes.
 * </dl>
 *
 * @version $Id$
 * @author Philip Cook
 * 
 */
public class TwoFibreBinghamFitter extends ConjGradMinimizer {


    protected Vector3D[] samples;
    
    protected double[] gi;

    protected double[] mu11DotXSq;
    protected double[] mu12DotXSq;

    protected double[] mu21DotXSq;
    protected double[] mu22DotXSq;
    
    protected double[] normC; // {normC1, d(normC1)/dk11, d(normC1)/dk21, normC2, d(normC2)/dk12, d(normC2)/dk22}
        
    public TwoFibreBinghamFitter(Vector3D[] vecs) {
	samples = vecs;
	gi = new double[samples.length];
	normC = new double[6];
	init(10);

	mu11DotXSq = new double[samples.length];
	mu21DotXSq = new double[samples.length];

	mu12DotXSq = new double[samples.length];
	mu22DotXSq = new double[samples.length];
	

    }

   /**
     * Runs the minimization and returns the answer.
     * @return {0.0, theta11, phi11, psi1, kappa11, kappa21, theta12, phi12, psi2, kappa12, kappa22, normC1, normC2}.
     */
    public double[] fitEstimatedParams(Vector3D[] initAxes, double[] initKs, double ftol) throws ConjGradMinimizerException {

	double[] atry = getParams(initAxes, initKs);

	minimise(atry, ftol);

	atry[4] = -1.0 * atry[4] * atry[4];
	atry[5] = atry[4] * Math.cos(atry[5]) * Math.cos(atry[5]);
	
	atry[9] = -1.0 * atry[9] * atry[9];
	atry[10] = atry[9] * Math.cos(atry[10]) * Math.cos(atry[10]);


	double[] params = new double[13];

	for (int i = 0; i < 11; i++) {
	    params[i] = atry[i];
	}
	
	double normC1 = 0.0;


	double[] der = new double[2];
	double[] hes = new double[4];


	try {
	    normC1 = BinghamFitter.bingc(atry[4], atry[5], der, hes);

	}
	catch (ConvergenceException e) {
	    throw new IllegalStateException("Couldn't calculate normC");
	}
	double normC2 = 0.0;

	try {
	    normC2 = BinghamFitter.bingc(atry[9], atry[10], der, hes);
	}
	catch (ConvergenceException e) {
	    throw new IllegalStateException("Couldn't calculate normC");
	}
	params[11] = normC1;
	params[12] = normC2;
	
	return params;
	
    }

    /**
     * @param params the output of the <code>fitEstimatedParams</code> method.
     *
     * @return {mu31, mu21, mu11, mu32, mu22, mu12}, where mu11 corresponds to the smallest eigenvalue 
     * (Bingham ordering scheme).
     */
    public static Vector3D[] getMus(double[] params) {
	Vector3D fittedE11 = Vector3D.vectorFromSPC(1.0, params[1], params[2]);
	Vector3D fittedE12 = Vector3D.vectorFromSPC(1.0, params[6], params[7]);

	Vector3D fittedE21 = new Vector3D(Math.cos(params[2])*Math.cos(params[3])*Math.cos(params[1]) - Math.sin(params[2])*Math.sin(params[3]), Math.cos(params[3])*Math.cos(params[1])*Math.sin(params[2]) + Math.cos(params[2])*Math.sin(params[3]), -(Math.cos(params[3])*Math.sin(params[1])));
	
	Vector3D fittedE22 = new Vector3D(Math.cos(params[7])*Math.cos(params[8])*Math.cos(params[6]) - Math.sin(params[7])*Math.sin(params[8]), Math.cos(params[8])*Math.cos(params[6])*Math.sin(params[7]) + Math.cos(params[7])*Math.sin(params[8]), -(Math.cos(params[8])*Math.sin(params[6])));
        
	Vector3D fittedE31 = fittedE11.cross(fittedE21);
	Vector3D fittedE32 = fittedE12.cross(fittedE22);

        return new Vector3D[] {fittedE31, fittedE21, fittedE11, fittedE32, fittedE22, fittedE12};

    }


    /**
     * Gets the parameter array for fObj.
     */
    public double[] getParams(Vector3D[] initAxes, double[] initKs) {

	double[] atry = new double[11];

	double[] tp = Vector3D.thetaPhi(initAxes[0]);
	
	atry[1] = tp[0];
	atry[2] = tp[1];

	Vector3D psiZero = Rotations.rotateVector(Rotations.X_AXIS, Rotations.Y_AXIS, atry[1]);
	psiZero = Rotations.rotateVector(psiZero, Rotations.Z_AXIS, atry[2]);

	atry[3] = Math.acos(-psiZero.dot(initAxes[1]));
	
	if (Math.abs( initAxes[1].dot(Rotations.rotateVector(psiZero, initAxes[0], atry[3])) ) < 0.99) {
	    atry[3] = -atry[3];
	}


	atry[4] = Math.sqrt(Math.abs(initKs[0]));
	atry[5] = Math.acos( Math.sqrt(Math.abs(initKs[1])) / atry[4] );

	tp = Vector3D.thetaPhi(initAxes[2]);
	
	atry[6] = tp[0];
	atry[7] = tp[1];

	psiZero = Rotations.rotateVector(Rotations.X_AXIS, Rotations.Y_AXIS, atry[6]);
	psiZero = Rotations.rotateVector(psiZero, Rotations.Z_AXIS, atry[7]);

	atry[8] = Math.acos(-psiZero.dot(initAxes[3]));

	if (Math.abs( initAxes[3].dot(Rotations.rotateVector(psiZero, initAxes[2], atry[8])) ) < 0.99) {
	    atry[8] = -atry[8];
	}

	atry[9] = Math.sqrt(Math.abs(initKs[2]));

	atry[10] = Math.acos( Math.sqrt(Math.abs(initKs[3])) / atry[9] );

	return atry;
	
    }
    


    protected void init(int noParams) {		
	super.init(noParams);
    }
    

    /**
     * @param atry should be {0.0, theta1_1, phi1_1, psi1, kappa1_1, alpha1, theta1_2, phi1_2, psi2, kappa1_2, alpha2}.
     */
    protected double fObj(double[] atry) {
	// set gi

	double sinT11 = Math.sin(atry[1]);
	double cosT11 = Math.cos(atry[1]);
	double sinP11 = Math.sin(atry[2]);
	double cosP11 = Math.cos(atry[2]);

	double sinPSI1 = Math.sin(atry[3]);
	double cosPSI1 = Math.cos(atry[3]);

	double sinT12 = Math.sin(atry[6]);
	double cosT12 = Math.cos(atry[6]);
	double sinP12 = Math.sin(atry[7]);
	double cosP12 = Math.cos(atry[7]);

	double sinPSI2 = Math.sin(atry[8]);
	double cosPSI2 = Math.cos(atry[8]);

	double cosAlpha1Sq = Math.cos(atry[5]) * Math.cos(atry[5]);
	double cosAlpha2Sq = Math.cos(atry[10]) * Math.cos(atry[10]);

	double minusK11Sq = -1.0 * atry[4] * atry[4];
	double minusK21Sq = -1.0 * atry[4] * atry[4] * cosAlpha1Sq;
	
	double normC1 = 0.0;

	try {

	    double[] der = new double[2];
	    double[] hes = new double[4];

	    normC1 = 4.0 * Math.PI * BinghamFitter.bingc(minusK11Sq, minusK21Sq, der, hes);
	    normC[0] = normC1;

	    normC[1] = 4.0 * Math.PI * der[0] * normC1;
	    normC[2] = 4.0 * Math.PI * der[1] * normC1;

	}
	catch (ConvergenceException e) {
	    throw new IllegalStateException("Couldn't calculate normC");
	}

	double minusK12Sq = -1.0 * atry[9] * atry[9];
	double minusK22Sq = -1.0 * atry[9] * atry[9] * cosAlpha2Sq;

	double normC2 = 0.0;

	try {

	    double[] der = new double[2];
	    double[] hes = new double[4];

		normC2 = 4.0 * Math.PI * BinghamFitter.bingc(minusK12Sq, minusK22Sq, der, hes);
		normC[3] = normC2 ;
		normC[4] = 4.0 * Math.PI * der[0] * normC2;
		normC[5] = 4.0 * Math.PI * der[1] * normC2;
	}
	catch (ConvergenceException e) {
	    throw new IllegalStateException("Couldn't calculate normC");
	}


	double f = 0.0;

	for (int i = 0; i < samples.length; i++) {

	    double mu11DotX = samples[i].z*cosT11 + samples[i].x*cosP11*sinT11 + samples[i].y*sinP11*sinT11;
 
	    double mu21DotX = samples[i].y*(cosPSI1*cosT11*sinP11 + cosP11*sinPSI1) + samples[i].x*(cosP11*cosPSI1*cosT11 - sinP11*sinPSI1) - samples[i].z*cosPSI1*sinT11;

	    double mu12DotX = samples[i].z*cosT12 + samples[i].x*cosP12*sinT12 + samples[i].y*sinP12*sinT12;

	    double mu22DotX = samples[i].y*(cosPSI2*cosT12*sinP12 + cosP12*sinPSI2) + samples[i].x*(cosP12*cosPSI2*cosT12 - sinP12*sinPSI2) - samples[i].z*cosPSI2*sinT12;

	    mu11DotXSq[i] = mu11DotX*mu11DotX;
	    mu21DotXSq[i] = mu21DotX*mu21DotX;

	    mu12DotXSq[i] = mu12DotX*mu12DotX;
	    mu22DotXSq[i] = mu22DotX*mu22DotX;
	    
	    gi[i] = Math.exp(minusK21Sq*mu21DotXSq[i] + minusK11Sq*mu11DotXSq[i])/ (2.*normC[0]) + Math.exp(minusK22Sq*mu22DotXSq[i] + minusK12Sq*mu12DotXSq[i]) / (2.*normC[3]);

	    f -= Math.log(gi[i]);
	    
	}
	
	return f;

    }


    /**
     * @param atry should be {0.0, theta1, phi1, psi1, kappa1_1, alpha1, 
     * theta2, phi2, psi2, kappa1_2, alpha2}.
     */
    protected double[] dfObj(double[] atry) {
	
	double sinT11 = Math.sin(atry[1]);
	double cosT11 = Math.cos(atry[1]);
	double sinP11 = Math.sin(atry[2]);
	double cosP11 = Math.cos(atry[2]);

	double sinPSI1 = Math.sin(atry[3]);
	double cosPSI1 = Math.cos(atry[3]);

	double sinT12 = Math.sin(atry[6]);
	double cosT12 = Math.cos(atry[6]);
	double sinP12 = Math.sin(atry[7]);
	double cosP12 = Math.cos(atry[7]);

	double sinPSI2 = Math.sin(atry[8]);
	double cosPSI2 = Math.cos(atry[8]);

	double k11 = atry[4];
	double k12 = atry[9];

	double alpha1 = atry[5];
	double alpha2 = atry[10];

	double cosAlpha1 = Math.cos(alpha1);
	double sinAlpha1 = Math.sin(alpha1);

	double cosAlpha2 = Math.cos(alpha2);
	double sinAlpha2 = Math.sin(alpha2);

	double cosAlpha1Sq = cosAlpha1 * cosAlpha1;
	double cosAlpha2Sq = cosAlpha2 * cosAlpha2;

	double[] dfda = new double[11];

	for (int i = 0; i < samples.length; i++) {

	    double[] dgda = new double[11];
 
	    double xx = samples[i].x;
	    double yy = samples[i].y;
	    double zz = samples[i].z;
	    
	    double aaa = Math.exp(-((k11*k11)*cosAlpha1Sq* mu21DotXSq[i]) - (k11*k11)*mu11DotXSq[i]);
	    double bbb = Math.exp(-((k12*k12)*cosAlpha2Sq* mu22DotXSq[i]) - (k12*k12)*mu12DotXSq[i]);

	    dgda[1] = (aaa* (-2*(k11*k11)*(xx*cosP11*cosT11 + yy*cosT11*sinP11 - zz*sinT11)* (zz*cosT11 + xx*cosP11*sinT11 + yy*sinP11*sinT11) - 2*(k11*k11)*cosAlpha1Sq* (yy*(cosPSI1*cosT11*sinP11 + cosP11*sinPSI1) + xx*(cosP11*cosPSI1*cosT11 - sinP11*sinPSI1) - zz*cosPSI1*sinT11)*(-(zz*cosPSI1*cosT11) - xx*cosP11*cosPSI1*sinT11 - yy*cosPSI1*sinP11*sinT11)))/ (2.*normC[0]);
	  
	  dfda[1] -= dgda[1] / gi[i];

	  dgda[2] = (aaa* (-2*(k11*k11)*cosAlpha1Sq* (xx*(-(cosPSI1*cosT11*sinP11) - cosP11*sinPSI1) + yy*(cosP11*cosPSI1*cosT11 - sinP11*sinPSI1))* (yy*(cosPSI1*cosT11*sinP11 + cosP11*sinPSI1) + xx*(cosP11*cosPSI1*cosT11 - sinP11*sinPSI1) - zz*cosPSI1*sinT11) - 2*(k11*k11)*(yy*cosP11*sinT11 - xx*sinP11*sinT11)* (zz*cosT11 + xx*cosP11*sinT11 + yy*sinP11*sinT11)))/(2.*normC[0]);

	 dfda[2] -= dgda[2] / gi[i];

	 dgda[3] = -((aaa*(k11*k11)*cosAlpha1Sq* (yy*(cosPSI1*cosT11*sinP11 + cosP11*sinPSI1) + xx*(cosP11*cosPSI1*cosT11 - sinP11*sinPSI1) - zz*cosPSI1*sinT11)*(xx* (-(cosPSI1*sinP11) - cosP11*cosT11*sinPSI1) + yy*(cosP11*cosPSI1 - cosT11*sinP11*sinPSI1) + zz*sinPSI1*sinT11))/ normC[0]);


	 dfda[3] -= dgda[3] / gi[i];


	 dgda[4] = (aaa* (-2*k11*cosAlpha1Sq* mu21DotXSq[i] - 2*k11*mu11DotXSq[i]))/ (2.*normC[0]) - (aaa* (-2*k11*cosAlpha1Sq* normC[2] - 2*k11*normC[1]))/(2.*(normC[0]*normC[0]));

	 dfda[4] -= dgda[4] / gi[i];


	 dgda[5] = (aaa*(k11*k11)*cosAlpha1*sinAlpha1* mu21DotXSq[i])/ normC[0] - (aaa*(k11*k11)*cosAlpha1*sinAlpha1* normC[2])/ (normC[0]*normC[0]);

	 dfda[5] -= dgda[5] / gi[i];
	 
	 
	 dgda[6] = (bbb* (-2*(k12*k12)*(xx*cosP12*cosT12 + yy*cosT12*sinP12 - zz*sinT12)* (zz*cosT12 + xx*cosP12*sinT12 + yy*sinP12*sinT12) - 2*(k12*k12)*cosAlpha2Sq* (yy*(cosPSI2*cosT12*sinP12 + cosP12*sinPSI2) + xx*(cosP12*cosPSI2*cosT12 - sinP12*sinPSI2) - zz*cosPSI2*sinT12)*(-(zz*cosPSI2*cosT12) - xx*cosP12*cosPSI2*sinT12 - yy*cosPSI2*sinP12*sinT12)))/(2.*normC[3]);
		
	 dfda[6] -= dgda[6] / gi[i];


	 dgda[7] = (bbb* (-2*(k12*k12)*cosAlpha2Sq* (xx*(-(cosPSI2*cosT12*sinP12) - cosP12*sinPSI2) + yy*(cosP12*cosPSI2*cosT12 - sinP12*sinPSI2))* (yy*(cosPSI2*cosT12*sinP12 + cosP12*sinPSI2) + xx*(cosP12*cosPSI2*cosT12 - sinP12*sinPSI2) - zz*cosPSI2*sinT12) - 2*(k12*k12)*(yy*cosP12*sinT12 - xx*sinP12*sinT12)* (zz*cosT12 + xx*cosP12*sinT12 + yy*sinP12*sinT12)))/(2.*normC[3]);


	 dfda[7] -= dgda[7] / gi[i];


	 dgda[8] = -((bbb*(k12*k12)*cosAlpha2Sq* (yy*(cosPSI2*cosT12*sinP12 + cosP12*sinPSI2) + xx*(cosP12*cosPSI2*cosT12 - sinP12*sinPSI2) - zz*cosPSI2*sinT12)*(xx* (-(cosPSI2*sinP12) - cosP12*cosT12*sinPSI2) + yy*(cosP12*cosPSI2 - cosT12*sinP12*sinPSI2) + zz*sinPSI2*sinT12))/normC[3]);

	 dfda[8] -= dgda[8] / gi[i];

	 dgda[9] = (bbb* (-2*k12*cosAlpha2Sq* mu22DotXSq[i] - 2*k12*mu12DotXSq[i]))/ (2.*normC[3]) - (bbb* (-2*k12*cosAlpha2Sq* normC[5] - 2*k12*normC[4]))/(2.*(normC[3]*normC[3]));

	 dfda[9] -= dgda[9] / gi[i];


	 dgda[10] = (bbb*(k12*k12)*cosAlpha2*sinAlpha2* mu22DotXSq[i])/ normC[3] - (bbb*(k12*k12)*cosAlpha2*sinAlpha2* normC[5])/(normC[3]*normC[3]);


	 dfda[10] -= dgda[10] / gi[i]; 

	}
	
	return dfda;

	} 
 

}

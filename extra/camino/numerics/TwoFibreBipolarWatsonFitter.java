package numerics;

import tools.*;

import optimizers.*;
import java.util.ArrayList;
import Jama.Matrix;

/**
 * <dl>
 *
 * <dt>Purpose: 
 *
 * <dd> Fits two Watson distributions to a collection of axes.
 *
 * <dt>Description:
 *
 * <dd> This class takes a collection of axes representing a set of estimates of two fibre directions. It uses the Marquardt optimization routine to find the pair of Watson distributions that best fit the axes, subject to the constraint that both distributions must be bipolar.
 * <dd>
 * </dl>
 *
 * @version $Id$
 * @author Philip Cook
 * 
 */
public class TwoFibreBipolarWatsonFitter extends TwoFibreFixedPropWatsonFitter {
	
  	
    // a[1]     a[2]   a[3]    a[4]    a[5]   a[6]   
    // theta1   phi1   kappa1  theta2  phi2   kappa2 
  

    double[] k1Mu1DotXSq;
    double[] k2Mu2DotXSq;

    /**
     * Returns the value of the objective function with parameters
     * atry.  Array dfda is filled with values of the first derivative
     * of fObj at atry wrt each of the parameters and d2fda2 is filled
     * with the second derivatives.
     *
     * @param atry The point at which to evaluate the objective
     * function. The parameters are as follows: {theta_1, phi_1, kappa_1, theta_2, phi_2, kappa_2, alpha}, which are the coordinates of the mean axis (1 and 2), the concentration parameters, and the mixing parameter. The mixing parameter is \sin^2 \alpha.
     *
     * @param dfda The array to be filled with the first
     * derivatives.
     *
     * @param d2fda2 The array to be filled with the second
     * derivatives.
     *
     * @return The value of the objective function.
     */
    protected double fObj(double[] atry, double[] dfda, double[][] d2fda2) {
		
	// fObj = -\sum_i \log(g) 
	// where g = \left( sin^2(\alpha) W(\mu_1, \kappa_1, x_i) + (1 - sin^2(\alpha))W(\mu_2, \kappa_2, x_i) \right)
	
	// want to MAXimise sumLogPDF, ie minimise this
	double minusSumLogPDF = 0.0;
		
	Vector3D mu1 = Vector3D.vectorFromSPC(1.0, atry[1], atry[2]);
	Vector3D mu2 = Vector3D.vectorFromSPC(1.0, atry[4], atry[5]);
		
	double[] g = new double[samples.length];
	double[][] dgda = new double[ma+1][samples.length];

	if (atry[3] > 15.0 || atry[6] > 15.0) {

	    	double logHyper1F1k1 = WatsonDistribution.logHyper1F1(0.5, 1.5, atry[3]*atry[3], 1e-9);
		double logHyper1F1k2 = WatsonDistribution.logHyper1F1(0.5, 1.5, atry[6]*atry[6], 1e-9);
		
		for (int i = 0; i < samples.length; i++) {
		    
		    mu1DotX[i] = mu1.x * samples[i].x + mu1.y * samples[i].y + mu1.z * samples[i].z;
		    mu2DotX[i] = mu2.x * samples[i].x + mu2.y * samples[i].y + mu2.z * samples[i].z;
		    
		    k1Mu1DotXSq[i] = atry[3]*atry[3] * mu1DotX[i] * mu1DotX[i];
		    k2Mu2DotXSq[i] = atry[6]*atry[6] * mu2DotX[i] * mu2DotX[i];
		    
		    //		    eK1Mu1DotXSq[i] = Math.exp(k1Mu1DotXSq[i]);
		    //		    eK2Mu2DotXSq[i] = Math.exp(k2Mu2DotXSq[i]);
		    
		    g[i] = Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1) / 2.0 + Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2) / 2.0; 
		    
		    minusSumLogPDF -= Math.log(g[i]);
		}
		
		double logDHyperDk1 = WatsonDistribution.logHyper1F1(1.5, 2.5, atry[3]*atry[3], 1e-9);
		double logDHyperDk2 = WatsonDistribution.logHyper1F1(1.5, 2.5, atry[6]*atry[6], 1e-9);

		firstDerivLog(atry,g, dgda, dfda, mu1DotX, mu2DotX, k1Mu1DotXSq, k2Mu2DotXSq, logHyper1F1k1, logHyper1F1k2, logDHyperDk1, logDHyperDk2);
		
		secondDerivLog(atry,g, dgda, d2fda2, mu1DotX, mu2DotX, k1Mu1DotXSq, k2Mu2DotXSq, logHyper1F1k1, logHyper1F1k2, logDHyperDk1, logDHyperDk2);


	}
	else {
	    double hyper1F1k1 = WatsonDistribution.hyper1F1(0.5, 1.5, atry[3]*atry[3], 1e-9);
	    double hyper1F1k2 = WatsonDistribution.hyper1F1(0.5, 1.5, atry[6]*atry[6], 1e-9);
	    
	    
	    for (int i = 0; i < samples.length; i++) {
		
		mu1DotX[i] = mu1.x * samples[i].x + mu1.y * samples[i].y + mu1.z * samples[i].z;
		mu2DotX[i] = mu2.x * samples[i].x + mu2.y * samples[i].y + mu2.z * samples[i].z;
		
		eK1Mu1DotXSq[i] = Math.exp(atry[3]*atry[3] * mu1DotX[i] * mu1DotX[i]);
		eK2Mu2DotXSq[i] = Math.exp(atry[6]*atry[6] * mu2DotX[i] * mu2DotX[i]);
		
		g[i] = ( (1.0 / hyper1F1k1) * eK1Mu1DotXSq[i] + (1.0 / hyper1F1k2) * eK2Mu2DotXSq[i] ) / 2.0;
		
		minusSumLogPDF -= Math.log(g[i]);
	    }
	    
	    double dHyperDk1 = WatsonDistribution.hyper1F1(1.5, 2.5, atry[3]*atry[3], 1e-9);
	    double dHyperDk2 = WatsonDistribution.hyper1F1(1.5, 2.5, atry[6]*atry[6], 1e-9);
	    
	    firstDeriv(atry,g, dgda, dfda, mu1DotX, mu2DotX, eK1Mu1DotXSq, eK2Mu2DotXSq, hyper1F1k1, hyper1F1k2, dHyperDk1, dHyperDk2);
	    
	    secondDeriv(atry,g, dgda, d2fda2, mu1DotX, mu2DotX, eK1Mu1DotXSq, eK2Mu2DotXSq, hyper1F1k1, hyper1F1k2, dHyperDk1, dHyperDk2);
	    
	}
	    
	return minusSumLogPDF;
    }
	
	
	
    /**
     * @param atry the parameter values
     * @param g the values of the function g = \left( \alpha W(\mu_1, \kappa_1, x_i) + (1 - \alpha) W(\mu_2, \kappa_2, x_i) \right) for all i
     * @param dgda Will be filled with the derivatives of g_i wrt a
     * @param dfda Will be filled with the derivatives of f wrt a
     */
    protected void firstDeriv(double[] atry, double[] g, double[][] dgda, double[] dfda, double[] mu1DotX, double[] mu2DotX, double[] eK1Mu1DotXSq, double[] eK2Mu2DotXSq, double hyper1F1k1, double hyper1F1k2, double dHyperDk1, double dHyperDk2) {
		
	// a[1]    a[2]   a[3]    a[4]    a[5]   a[6]   
	// theta1  phi1   kappa1  theta2  phi2   kappa2 
		
	// g = \alpha W_1(\theta_1, \phi_1, \kappa_1) + \alpha W_2(\theta_2, \phi_2, \kappa_2)
		
	// (d/dx)(log g) = (1/g)(dg/dx)
	// (d/dx)(-\sum_i log g) = -\sum_i (1/g_i)(dg_i/dx)

	double cosAtry1 = Math.cos(atry[1]);
	double sinAtry1 = Math.sin(atry[1]);
			
	double cosAtry2 = Math.cos(atry[2]);
	double sinAtry2 = Math.sin(atry[2]);
		
	double cosAtry4 = Math.cos(atry[4]);
	double sinAtry4 = Math.sin(atry[4]);
		
	double cosAtry5 = Math.cos(atry[5]);
	double sinAtry5 = Math.sin(atry[5]);
		
	// code pasted from Mathematica
	for (int i = 0; i < samples.length; i++) {
			
	    double xx = samples[i].x;
	    double yy = samples[i].y;
	    double zz = samples[i].z;

	    double mu1DotXSq = mu1DotX[i]*mu1DotX[i];
	    double mu2DotXSq = mu2DotX[i]*mu2DotX[i];

	    dgda[1][i] = (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1;
			
	    dfda[1] -= (1.0 / g[i]) * dgda[1][i];
			
	    dgda[2][i] = (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(yy*cosAtry2*sinAtry1 - xx*sinAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1;

			
	    dfda[2] -= (1.0 / g[i]) * dgda[2][i];
			
	    dgda[3][i] = (-0.3333333333333333*eK1Mu1DotXSq[i]*atry[3]*dHyperDk1)/(hyper1F1k1*hyper1F1k1) + (eK1Mu1DotXSq[i]*atry[3]* mu1DotXSq)/hyper1F1k1;
			 
	    dfda[3] -= (1.0 / g[i]) * dgda[3][i];
			
	    dgda[4][i] = (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(xx*cosAtry4*cosAtry5 - zz*sinAtry4 + yy*cosAtry4*sinAtry5)* (mu2DotX[i]))/hyper1F1k2;
			
	    dfda[4] -= (1.0 / g[i]) * dgda[4][i];
			
			
	    dgda[5][i] = (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(yy*cosAtry5*sinAtry4 - xx*sinAtry4*sinAtry5)*(mu2DotX[i]))/hyper1F1k2;
			
	    dfda[5] -= (1.0 / g[i]) * dgda[5][i];
			
	    dgda[6][i] = (-0.3333333333333333*eK2Mu2DotXSq[i]*atry[6]*dHyperDk2)/(hyper1F1k2*hyper1F1k2) + (eK2Mu2DotXSq[i]*atry[6]* mu2DotXSq)/hyper1F1k2;

	    dfda[6] -= (1.0 / g[i]) * dgda[6][i];
			
		
	}
		
    }
	



    /**
     *
     * Computes the first derivative. Takes a logarithmic form of some parameters, which is useful for high concentration.
     * @param atry the parameter values
     * @param g the values of the function g = \left( \alpha W(\mu_1, \kappa_1, x_i) + (1 - \alpha) W(\mu_2, \kappa_2, x_i) \right) for all i
     * @param dgda Will be filled with the derivatives of g_i wrt a
     * @param dfda Will be filled with the derivatives of f wrt a
     */
    protected void firstDerivLog(double[] atry, double[] g, double[][] dgda, double[] dfda, double[] mu1DotX, double[] mu2DotX, double[] k1Mu1DotXSq, double[] k2Mu2DotXSq, double logHyper1F1k1, double logHyper1F1k2, double logDHyperDk1, double logDHyperDk2) {

	// a[1]    a[2]   a[3]    a[4]    a[5]   a[6]   
	// theta1  phi1   kappa1  theta2  phi2   kappa2 
		
	// g = \alpha W_1(\theta_1, \phi_1, \kappa_1) + \alpha W_2(\theta_2, \phi_2, \kappa_2)
		
	// (d/dx)(log g) = (1/g)(dg/dx)
	// (d/dx)(-\sum_i log g) = -\sum_i (1/g_i)(dg_i/dx)

	double cosAtry1 = Math.cos(atry[1]);
	double sinAtry1 = Math.sin(atry[1]);
			
	double cosAtry2 = Math.cos(atry[2]);
	double sinAtry2 = Math.sin(atry[2]);
		
	double cosAtry4 = Math.cos(atry[4]);
	double sinAtry4 = Math.sin(atry[4]);
		
	double cosAtry5 = Math.cos(atry[5]);
	double sinAtry5 = Math.sin(atry[5]);
		
	// code pasted from Mathematica
	for (int i = 0; i < samples.length; i++) {
			
	    double xx = samples[i].x;
	    double yy = samples[i].y;
	    double zz = samples[i].z;

	    double mu1DotXSq = mu1DotX[i]*mu1DotX[i];
	    double mu2DotXSq = mu2DotX[i]*mu2DotX[i];

	    //	    dgda[1][i] = (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1;
		
	    dgda[1][i] = (atry[3]*atry[3])*(xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2)* (mu1DotX[i]) * Math.exp( k1Mu1DotXSq[i] - logHyper1F1k1 );
		
	
	    dfda[1] -= (1.0 / g[i]) * dgda[1][i];
			
	    //	    dgda[2][i] = (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(yy*cosAtry2*sinAtry1 - xx*sinAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1;

	    dgda[2][i] = 
		(atry[3]*atry[3])*(yy*cosAtry2*sinAtry1 - xx*sinAtry1*sinAtry2)*(mu1DotX[i])*
		Math.exp( k1Mu1DotXSq[i] - logHyper1F1k1 );
			
	    dfda[2] -= (1.0 / g[i]) * dgda[2][i];
			
	    //	    dgda[3][i] = (-0.3333333333333333*eK1Mu1DotXSq[i]*atry[3]*dHyperDk1)/(hyper1F1k1*hyper1F1k1) + (eK1Mu1DotXSq[i]*atry[3]* mu1DotXSq)/hyper1F1k1;
	    
	    dgda[3][i] = -0.3333333333333333*atry[3] * Math.exp(
				  k1Mu1DotXSq[i] + logDHyperDk1 - (logHyper1F1k1 + logHyper1F1k1)
				  ) 
		+ atry[3]* mu1DotXSq * Math.exp(
						k1Mu1DotXSq[i] - logHyper1F1k1
						);
			 
	    dfda[3] -= (1.0 / g[i]) * dgda[3][i];
			
	    //	    dgda[4][i] = (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(xx*cosAtry4*cosAtry5 - zz*sinAtry4 + yy*cosAtry4*sinAtry5)* (mu2DotX[i]))/hyper1F1k2;

	    dgda[4][i] = (atry[6]*atry[6])*(xx*cosAtry4*cosAtry5 - zz*sinAtry4 + yy*cosAtry4*sinAtry5)* (mu2DotX[i])*Math.exp(
															      k2Mu2DotXSq[i] - logHyper1F1k2
															      );
				
	
	    dfda[4] -= (1.0 / g[i]) * dgda[4][i];


	    //	dgda[5][i] = (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(yy*cosAtry5*sinAtry4 - xx*sinAtry4*sinAtry5)*(mu2DotX[i]))/hyper1F1k2;		
			
	    dgda[5][i] = (atry[6]*atry[6])*(yy*cosAtry5*sinAtry4 - xx*sinAtry4*sinAtry5)*(mu2DotX[i]) * 
		Math.exp( 
			 k2Mu2DotXSq[i] +  - logHyper1F1k2
			 );
	    
	    dfda[5] -= (1.0 / g[i]) * dgda[5][i];
			

	    // dgda[6][i] = (-0.3333333333333333*eK2Mu2DotXSq[i]*atry[6]*dHyperDk2)/(hyper1F1k2*hyper1F1k2) + (eK2Mu2DotXSq[i]*atry[6]* mu2DotXSq)/hyper1F1k2;
  
	    dgda[6][i] = -0.3333333333333333*atry[6] * Math.exp(
				  k2Mu2DotXSq[i] + logDHyperDk2 - (logHyper1F1k2 + logHyper1F1k2)
				  ) + atry[6]*mu2DotXSq*Math.exp(
					       k2Mu2DotXSq[i] - logHyper1F1k2
					       );

	    dfda[6] -= (1.0 / g[i]) * dgda[6][i];
			
		
	}
		
    }
	
	
    /** 
     * @param g the value of g(x) for each sample x. 
     * @param d2fda2 the second derivative matrix of fObj wrt each parameter. 
     * @param dgda the derivatives of g(x) wrt each parameter. 
     */ 
    protected void secondDeriv(double[] atry, double[] g, double[][] dgda, double[][] d2fda2, double[] mu1DotX, double[] mu2DotX, double[] eK1Mu1DotXSq, double[] eK2Mu2DotXSq, double hyper1F1k1, double hyper1F1k2, double dHyperDk1, double dHyperDk2) {
	
	
	// (d^2/dxdz)(Log g) = - (1 / g^2)(dg/dx)(dg/dz) + (d^2g/dzdx)(1/g)
	// -\sum_i (d^2/dxdz)(Log g) = -\sum_i -(1 / g_i^2)(dg/dx)(dg/dz) + (d^2g/dzdx)(1/g_i)

	double d2Hyper1F1k1 = WatsonDistribution.hyper1F1(2.5,3.5,atry[3]*atry[3], 1e-9);
	double d2Hyper1F1k2 = WatsonDistribution.hyper1F1(2.5,3.5,atry[6]*atry[6], 1e-9);
		
	double cosAtry1 = Math.cos(atry[1]);
	double sinAtry1 = Math.sin(atry[1]);
		
	double cosAtry2 = Math.cos(atry[2]);
	double sinAtry2 = Math.sin(atry[2]);
		
	double cosAtry4 = Math.cos(atry[4]);
	double sinAtry4 = Math.sin(atry[4]);
		
	double cosAtry5 = Math.cos(atry[5]);
	double sinAtry5 = Math.sin(atry[5]);

	for (int i = 0; i < samples.length; i++) {
			
	    double xx = samples[i].x;
	    double yy = samples[i].y;
	    double zz = samples[i].z;
			
	    double invgi = 1.0 / g[i];
	    double invgiSq = invgi * invgi;
	 
	    double mu1DotXSq = mu1DotX[i]*mu1DotX[i];
	    double mu2DotXSq = mu2DotX[i]*mu2DotX[i];

	    double aa = xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2;
	    double bb = yy*cosAtry2*sinAtry1 - xx*sinAtry1*sinAtry2;
	    double cc = xx*cosAtry4*cosAtry5 - zz*sinAtry4 + yy*cosAtry4*sinAtry5;
	    double dd = yy*cosAtry5*sinAtry4 - xx*sinAtry4*sinAtry5;

	    // In[48]:= CForm[D[dfda1, atry[1]]] 
	    // Out[48]//CForm=
	    d2fda2[1][1] -= invgi * (
	
				     (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(aa*aa))/ hyper1F1k1 + (eK1Mu1DotXSq[i]*(atry[3]*atry[3])* (-(zz*cosAtry1) - xx*cosAtry2*sinAtry1 - yy*sinAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1 + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3]*atry[3])*(aa*aa)* mu1DotXSq)/ hyper1F1k1

			   
				     ) - invgiSq * dgda[1][i] * dgda[1][i];
			
	    // In[49]:= CForm[D[dfda1, atry[2]]]
	    // Out[49]//CForm=
	    d2fda2[1][2] -= invgi * (
				 
				     (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(aa)* (bb))/hyper1F1k1 + (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(yy*cosAtry1*cosAtry2 - xx*cosAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1 + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3]*atry[3])*(aa)* (bb)* mu1DotXSq)/ hyper1F1k1
				     ) - invgiSq * dgda[1][i] * dgda[2][i];
			
			
	    // In[50]:=
	    // CForm[D[dyda1, atry[3]]]
	    // Out[50]//CForm=
	    d2fda2[1][3] -= invgi * ( 
				 
				     (2*eK1Mu1DotXSq[i]*atry[3]* (aa)* (mu1DotX[i]))/hyper1F1k1 - (0.6666666666666666*eK1Mu1DotXSq[i]*(atry[3]*atry[3]*atry[3])*dHyperDk1* (aa)* (mu1DotX[i]))/ (hyper1F1k1*hyper1F1k1) + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3])*(aa)* (mu1DotXSq*mu1DotX[i]))/ hyper1F1k1
 
				     ) - invgiSq * dgda[1][i] * dgda[3][i];
			
			
	    // 	  // In[52]:=
	    // 	  // CForm[D[dfda1, atry[4]]]
	    // 	  // Out[52]//CForm=
	    // 	  d2fda2[1][4] -= 0.0;
			
	    // 	  // In[53]:=
	    // 	  // CForm[D[dfda1, atry[5]]]
	    // 	  // Out[53]//CForm=
	    // 	  d2fda2[1][5] -= 0.0;
			
	    // 	  // In[54]:=
	    // 	  // CForm[D[dfda1, atry[6]]]
	    // 	  // Out[54]//CForm=
	    // 	  d2fda2[1][6] -= 0.0;
			
			
	    //In[56]:=
	    // CForm[D[dfda2, atry[2]]]
	    // Out[56]//CForm=
	    d2fda2[2][2] -= invgi * (
			
				     (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(bb*bb))/ hyper1F1k1 + (eK1Mu1DotXSq[i]*(atry[3]*atry[3])* (-(xx*cosAtry2*sinAtry1) - yy*sinAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1 + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3]*atry[3])*(bb*bb)* mu1DotXSq)/ hyper1F1k1


				     ) - invgiSq * dgda[2][i] * dgda[2][i]; 
			
	    // In[57]:=
	    // CForm[D[dfda2, atry[3]]]
	    // Out[57]//CForm=
	    d2fda2[2][3] -= invgi * (

				     (2*eK1Mu1DotXSq[i]*atry[3]* (bb)* (mu1DotX[i]))/hyper1F1k1 - (0.6666666666666666*eK1Mu1DotXSq[i]*(atry[3]*atry[3]*atry[3])*dHyperDk1* (bb)* (mu1DotX[i]))/ (hyper1F1k1*hyper1F1k1) + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3])*(bb)* (mu1DotXSq*mu1DotX[i]))/ hyper1F1k1


				     ) - invgiSq * dgda[2][i] * dgda[3][i]; 
			
	    // 	  // In[58]:=
	    // 	  // CForm[D[dfda2, atry[4]]]
	    // 	  // Out[58]//CForm=
	    // 	  d2fda2[2][4] -= 0.0;
			
	    // 	  // In[59]:=
	    // 	  // CForm[D[dfda2, atry[5]]]
	    // 	  // Out[59]//CForm=
	    // 	  d2fda2[2][5] -= 0.0;
			
	    // 	  // In[60]:=
	    // 	  // CForm[D[dfda2, atry[6]]]
	    // 	  // Out[60]//CForm=
	    // 	  d2fda2[2][6] -= 0.0;
			
			
	    // In[62]:=
	    // CForm[D[dfda3, atry[3]]]
	    // Out[62]//CForm=
	    d2fda2[3][3] -= invgi * (

				     (-0.3333333333333333*eK1Mu1DotXSq[i]*dHyperDk1)/(hyper1F1k1*hyper1F1k1) + (0.4444444444444444*eK1Mu1DotXSq[i]*(atry[3]*atry[3])*(dHyperDk1*dHyperDk1))/ (hyper1F1k1*hyper1F1k1*hyper1F1k1) - (0.4*eK1Mu1DotXSq[i]* (atry[3]*atry[3])*d2Hyper1F1k1)/(hyper1F1k1*hyper1F1k1) + (eK1Mu1DotXSq[i]* mu1DotXSq)/ hyper1F1k1 - (1.3333333333333333* eK1Mu1DotXSq[i]* (atry[3]*atry[3])*dHyperDk1* mu1DotXSq)/ (hyper1F1k1*hyper1F1k1) + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(mu1DotXSq*mu1DotXSq))/ hyper1F1k1


				     ) - invgiSq * dgda[3][i] * dgda[3][i]; 
			
	    // 	  // In[63]:=
	    // 	  // CForm[D[dfda3, atry[4]]]
	    // 	  // Out[63]//CForm=
	    // 	  d2fda2[3][4] -= 0.0;
			
	    // 	  // In[64]:=
	    // 	  // CForm[D[dfda3, atry[5]]]
	    // 	  // Out[64]//CForm=
	    // 	  d2fda2[3][5] -= 0.0;
			
	    // 	  // In[65]:=
	    // 	  // CForm[D[dfda3, atry[6]]]
	    // 	  // Out[65]//CForm=
	    // 	  d2fda2[3][6] -= 0.0;
			
			
	    // In[67]:=
	    // CForm[D[dfda4, atry[4]]]
	    // Out[67]//CForm=
	    d2fda2[4][4] -= invgi * (

				     (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(cc*cc))/ hyper1F1k2 + (eK2Mu2DotXSq[i]*(atry[6]*atry[6])* (-(zz*cosAtry4) - xx*cosAtry5*sinAtry4 - yy*sinAtry4*sinAtry5)* (mu2DotX[i]))/hyper1F1k2 + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6]*atry[6])*(cc*cc)* mu2DotXSq)/ hyper1F1k2


				     ) - invgiSq * dgda[4][i] * dgda[4][i]; 
			
			
	    // In[68]:=
	    // CForm[D[dfda4, atry[5]]]
	    // Out[68]//CForm=
	    d2fda2[4][5] -= invgi * (


				     (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(cc)* (dd))/hyper1F1k2 + (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(yy*cosAtry4*cosAtry5 - xx*cosAtry4*sinAtry5)* (mu2DotX[i]))/hyper1F1k2 + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6]*atry[6])*(cc)* (dd)* mu2DotXSq)/ hyper1F1k2


			
				     ) - invgiSq * dgda[4][i] * dgda[5][i]; 
			
	    // In[69]:=
	    // CForm[D[dfda4, atry[6]]]
	    // Out[69]//CForm=
	    d2fda2[4][6] -= invgi * (

				     (2*eK2Mu2DotXSq[i]*atry[6]* (cc)* (mu2DotX[i]))/hyper1F1k2 - (0.6666666666666666*eK2Mu2DotXSq[i]*(atry[6]*atry[6]*atry[6])*dHyperDk2* (cc)* (mu2DotX[i]))/ (hyper1F1k2*hyper1F1k2) + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6])*(cc)* (mu2DotXSq*mu2DotX[i]))/ hyper1F1k2



				     ) - invgiSq * dgda[4][i] * dgda[6][i]; 
			
	    // In[70]:=
	    // CForm[D[dfda4, atry[7]]]
			
	    // In[71]:=
	    // CForm[D[dfda5, atry[5]]]
	    // Out[71]//CForm=
	    d2fda2[5][5] -= invgi * (

				     (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(dd*dd))/ hyper1F1k2 + (eK2Mu2DotXSq[i]*(atry[6]*atry[6])* (-(xx*cosAtry5*sinAtry4) - yy*sinAtry4*sinAtry5)* (mu2DotX[i]))/hyper1F1k2 + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6]*atry[6])*(dd*dd)* mu2DotXSq)/ hyper1F1k2

				     ) - invgiSq * dgda[5][i] * dgda[5][i]; 
			
	    // In[72]:=
	    // CForm[D[dfda5, atry[6]]]
	    // Out[72]//CForm=
	    d2fda2[5][6] -= invgi * (

				     (2*eK2Mu2DotXSq[i]*atry[6]* (dd)* (mu2DotX[i]))/hyper1F1k2 - (0.6666666666666666*eK2Mu2DotXSq[i]*(atry[6]*atry[6]*atry[6])*dHyperDk2* (dd)* (mu2DotX[i]))/ (hyper1F1k2*hyper1F1k2) + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6])*(dd)* (mu2DotXSq*mu2DotX[i]))/ hyper1F1k2
				     ) - invgiSq * dgda[5][i] * dgda[6][i]; 
			
	    // In[73]:=
	    // CForm[D[dfda5, atry[7]]]
	    // Out[73]//CForm=
			
	    // In[74]:=
	    // CForm[D[dfda6, atry[6]]]
	    // Out[74]//CForm=
	    d2fda2[6][6] -= invgi * ( 


				     (-0.3333333333333333*eK2Mu2DotXSq[i]*dHyperDk2)/(hyper1F1k2*hyper1F1k2) + (0.4444444444444444*eK2Mu2DotXSq[i]*(atry[6]*atry[6])*(dHyperDk2*dHyperDk2))/ (hyper1F1k2*hyper1F1k2*hyper1F1k2) - (0.4*eK2Mu2DotXSq[i]* (atry[6]*atry[6])*d2Hyper1F1k2)/(hyper1F1k2*hyper1F1k2) + (eK2Mu2DotXSq[i]* mu2DotXSq)/ hyper1F1k2 - (1.3333333333333333* eK2Mu2DotXSq[i]* (atry[6]*atry[6])*dHyperDk2* mu2DotXSq)/ (hyper1F1k2*hyper1F1k2) + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(mu2DotXSq*mu2DotXSq))/ hyper1F1k2

			
				     ) - invgiSq * dgda[6][i] * dgda[6][i]; 
			
			
			
		
	}
		
	for (int i = 1; i <= ma; i++) {
	    for (int j = 1; j <= i; j++) {
		d2fda2[i][j] = d2fda2[j][i]; 
	    }
	}
		
    }


	
	
    /** 
     * Computes the second derivative. Takes a logarithmic form of some parameters, which is useful for high concentration.
     *
     * @param g the value of g(x) for each sample x. 
     * @param d2fda2 the second derivative matrix of fObj wrt each parameter. 
     * @param dgda the derivatives of g(x) wrt each parameter. 
     */ 
    protected void secondDerivLog(double[] atry, double[] g, double[][] dgda, double[][] d2fda2, double[] mu1DotX, double[] mu2DotX, double[] k1Mu1DotXSq, double[] k2Mu2DotXSq, double logHyper1F1k1, double logHyper1F1k2, double logDHyperDk1, double logDHyperDk2) {
	
	// (d^2/dxdz)(Log g) = - (1 / g^2)(dg/dx)(dg/dz) + (d^2g/dzdx)(1/g)
	// -\sum_i (d^2/dxdz)(Log g) = -\sum_i -(1 / g_i^2)(dg/dx)(dg/dz) + (d^2g/dzdx)(1/g_i)

	double logD2Hyper1F1k1 = WatsonDistribution.logHyper1F1(2.5,3.5,atry[3]*atry[3], 1e-9);
	double logD2Hyper1F1k2 = WatsonDistribution.logHyper1F1(2.5,3.5,atry[6]*atry[6], 1e-9);
		
	double cosAtry1 = Math.cos(atry[1]);
	double sinAtry1 = Math.sin(atry[1]);
		
	double cosAtry2 = Math.cos(atry[2]);
	double sinAtry2 = Math.sin(atry[2]);
		
	double cosAtry4 = Math.cos(atry[4]);
	double sinAtry4 = Math.sin(atry[4]);
		
	double cosAtry5 = Math.cos(atry[5]);
	double sinAtry5 = Math.sin(atry[5]);

	for (int i = 0; i < samples.length; i++) {
			
	    double xx = samples[i].x;
	    double yy = samples[i].y;
	    double zz = samples[i].z;
			
	    double invgi = 1.0 / g[i];
	    double invgiSq = invgi * invgi;
	 
	    double mu1DotXSq = mu1DotX[i]*mu1DotX[i];
	    double mu2DotXSq = mu2DotX[i]*mu2DotX[i];

	    double aa = xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2;
	    double bb = yy*cosAtry2*sinAtry1 - xx*sinAtry1*sinAtry2;
	    double cc = xx*cosAtry4*cosAtry5 - zz*sinAtry4 + yy*cosAtry4*sinAtry5;
	    double dd = yy*cosAtry5*sinAtry4 - xx*sinAtry4*sinAtry5;


	    // In[48]:= CForm[D[dfda1, atry[1]]] 
	    // Out[48]//CForm=
	    d2fda2[1][1] -= invgi * (
	
				     // (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*Math.pow(xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2,2))/ hyper1F1k1 + (eK1Mu1DotXSq[i]*(atry[3]*atry[3])* (-(zz*cosAtry1) - xx*cosAtry2*sinAtry1 - yy*sinAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1 + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3]*atry[3])*Math.pow(xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2,2)* mu1DotXSq)/ hyper1F1k1
				     
				     
				     (atry[3]*atry[3])*(aa*aa)*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1) + (atry[3]*atry[3])* (-(zz*cosAtry1) - xx*cosAtry2*sinAtry1 - yy*sinAtry1*sinAtry2)* (mu1DotX[i])*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1) + 2.0 * (atry[3]*atry[3]*atry[3]*atry[3])*(aa*aa)* mu1DotXSq * Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1)

			   
				     ) - invgiSq * dgda[1][i] * dgda[1][i];
			
	    // In[49]:= CForm[D[dfda1, atry[2]]]
	    // Out[49]//CForm=
	    d2fda2[1][2] -= invgi * (

// 				     (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2)* (bb))/hyper1F1k1 + (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(yy*cosAtry1*cosAtry2 - xx*cosAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1 + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3]*atry[3])*(xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2)* (bb)* mu1DotXSq)/ hyper1F1k1
				 
				     (atry[3]*atry[3])*(aa)*(bb)*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1) + (atry[3]*atry[3])*(yy*cosAtry1*cosAtry2 - xx*cosAtry1*sinAtry2)* (mu1DotX[i])*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1) + (2*(atry[3]*atry[3]*atry[3]*atry[3])*(aa)* (bb)* mu1DotXSq)*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1)
				     ) - invgiSq * dgda[1][i] * dgda[2][i];
			
			
	    // In[50]:=
	    // CForm[D[dyda1, atry[3]]]
	    // Out[50]//CForm=
	    d2fda2[1][3] -= invgi * ( 
				 
// 				     (2*eK1Mu1DotXSq[i]*atry[3]* (xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1 - (0.6666666666666666*eK1Mu1DotXSq[i]*(atry[3]*atry[3]*atry[3])*dHyperDk1* (xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2)* (mu1DotX[i]))/ (hyper1F1k1*hyper1F1k1) + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3])*(xx*cosAtry1*cosAtry2 - zz*sinAtry1 + yy*cosAtry1*sinAtry2)* (mu1DotXSq*mu1DotX[i]))/ hyper1F1k1

				     (2.0*atry[3]* (aa)* (mu1DotX[i]))*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1) - 0.6666666666666666*(atry[3]*atry[3]*atry[3])* (aa)* (mu1DotX[i]) * Math.exp(k1Mu1DotXSq[i] + logDHyperDk1 - (logHyper1F1k1+logHyper1F1k1)) + 2.0*(atry[3]*atry[3]*atry[3])*(aa)* (mu1DotXSq*mu1DotX[i])*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1)
 
				     ) - invgiSq * dgda[1][i] * dgda[3][i];
			
			
	    // 	  // In[52]:=
	    // 	  // CForm[D[dfda1, atry[4]]]
	    // 	  // Out[52]//CForm=
	    // 	  d2fda2[1][4] -= 0.0;
			
	    // 	  // In[53]:=
	    // 	  // CForm[D[dfda1, atry[5]]]
	    // 	  // Out[53]//CForm=
	    // 	  d2fda2[1][5] -= 0.0;
			
	    // 	  // In[54]:=
	    // 	  // CForm[D[dfda1, atry[6]]]
	    // 	  // Out[54]//CForm=
	    // 	  d2fda2[1][6] -= 0.0;
			
			
	    //In[56]:=
	    // CForm[D[dfda2, atry[2]]]
	    // Out[56]//CForm=
	    d2fda2[2][2] -= invgi * (
			
// 				     (eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(bb*bb))/ hyper1F1k1 + (eK1Mu1DotXSq[i]*(atry[3]*atry[3])* (-(xx*cosAtry2*sinAtry1) - yy*sinAtry1*sinAtry2)* (mu1DotX[i]))/hyper1F1k1 + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3]*atry[3])*(bb*bb)* mu1DotXSq)/ hyper1F1k1

				     (atry[3]*atry[3])*(bb*bb)*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1) + (atry[3]*atry[3])* (-(xx*cosAtry2*sinAtry1) - yy*sinAtry1*sinAtry2)* (mu1DotX[i])*Math.exp(k1Mu1DotXSq[i]- logHyper1F1k1) + 2.0*(atry[3]*atry[3]*atry[3]*atry[3])*(bb*bb)*mu1DotXSq *Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1)

				     ) - invgiSq * dgda[2][i] * dgda[2][i]; 
			
	    // In[57]:=
	    // CForm[D[dfda2, atry[3]]]
	    // Out[57]//CForm=
	    d2fda2[2][3] -= invgi * (

// 				     (2*eK1Mu1DotXSq[i]*atry[3]* (bb)* (mu1DotX[i]))/hyper1F1k1 - (0.6666666666666666*eK1Mu1DotXSq[i]*(atry[3]*atry[3]*atry[3])*dHyperDk1* (bb)* (mu1DotX[i]))/ (hyper1F1k1*hyper1F1k1) + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3]*atry[3])*(bb)* (mu1DotXSq*mu1DotX[i]))/ hyper1F1k1

				     2.0*atry[3]* (bb)* (mu1DotX[i])*Math.exp(k1Mu1DotXSq[i]-logHyper1F1k1) - 0.6666666666666666*(atry[3]*atry[3]*atry[3])*(bb)*(mu1DotX[i])*Math.exp(k1Mu1DotXSq[i] + logDHyperDk1 - (logHyper1F1k1+logHyper1F1k1)) + 2.0*(atry[3]*atry[3]*atry[3])*(bb)* (mu1DotXSq*mu1DotX[i])*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1)


				     ) - invgiSq * dgda[2][i] * dgda[3][i]; 
			
	    // 	  // In[58]:=
	    // 	  // CForm[D[dfda2, atry[4]]]
	    // 	  // Out[58]//CForm=
	    // 	  d2fda2[2][4] -= 0.0;
			
	    // 	  // In[59]:=
	    // 	  // CForm[D[dfda2, atry[5]]]
	    // 	  // Out[59]//CForm=
	    // 	  d2fda2[2][5] -= 0.0;
			
	    // 	  // In[60]:=
	    // 	  // CForm[D[dfda2, atry[6]]]
	    // 	  // Out[60]//CForm=
	    // 	  d2fda2[2][6] -= 0.0;
			
			
	    // In[62]:=
	    // CForm[D[dfda3, atry[3]]]
	    // Out[62]//CForm=
	    d2fda2[3][3] -= invgi * (

// 				     (-0.3333333333333333*eK1Mu1DotXSq[i]*dHyperDk1)/(hyper1F1k1*hyper1F1k1) + (0.4444444444444444*eK1Mu1DotXSq[i]*(atry[3]*atry[3])*(dHyperDk1*dHyperDk1))/ (hyper1F1k1*hyper1F1k1*hyper1F1k1) - (0.4*eK1Mu1DotXSq[i]* (atry[3]*atry[3])*d2Hyper1F1k1)/(hyper1F1k1*hyper1F1k1) + (eK1Mu1DotXSq[i]* mu1DotXSq)/ hyper1F1k1 - (1.3333333333333333* eK1Mu1DotXSq[i]* (atry[3]*atry[3])*dHyperDk1* mu1DotXSq)/ (hyper1F1k1*hyper1F1k1) + (2*eK1Mu1DotXSq[i]* (atry[3]*atry[3])*(mu1DotXSq*mu1DotXSq))/ hyper1F1k1

				     (-0.3333333333333333)*Math.exp(k1Mu1DotXSq[i]+logDHyperDk1-(logHyper1F1k1+logHyper1F1k1)) + 0.4444444444444444*(atry[3]*atry[3])*Math.exp(k1Mu1DotXSq[i] + logDHyperDk1 + logDHyperDk1 - (logHyper1F1k1+logHyper1F1k1+logHyper1F1k1)) - 0.4*(atry[3]*atry[3])*Math.exp(k1Mu1DotXSq[i]+logD2Hyper1F1k1 -(logHyper1F1k1+logHyper1F1k1)) + mu1DotXSq*Math.exp(k1Mu1DotXSq[i] - logHyper1F1k1) - 1.3333333333333333*(atry[3]*atry[3])*mu1DotXSq*Math.exp(k1Mu1DotXSq[i] + logDHyperDk1 - (logHyper1F1k1+logHyper1F1k1)) + 2.0*(atry[3]*atry[3])*(mu1DotXSq*mu1DotXSq)*Math.exp(k1Mu1DotXSq[i]- logHyper1F1k1)


				     ) - invgiSq * dgda[3][i] * dgda[3][i]; 
			
	    // 	  // In[63]:=
	    // 	  // CForm[D[dfda3, atry[4]]]
	    // 	  // Out[63]//CForm=
	    // 	  d2fda2[3][4] -= 0.0;
			
	    // 	  // In[64]:=
	    // 	  // CForm[D[dfda3, atry[5]]]
	    // 	  // Out[64]//CForm=
	    // 	  d2fda2[3][5] -= 0.0;
			
	    // 	  // In[65]:=
	    // 	  // CForm[D[dfda3, atry[6]]]
	    // 	  // Out[65]//CForm=
	    // 	  d2fda2[3][6] -= 0.0;
			
			
	    // In[67]:=
	    // CForm[D[dfda4, atry[4]]]
	    // Out[67]//CForm=
	    d2fda2[4][4] -= invgi * (

// 				     (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(cc*cc))/ hyper1F1k2 + (eK2Mu2DotXSq[i]*(atry[6]*atry[6])* (-(zz*cosAtry4) - xx*cosAtry5*sinAtry4 - yy*sinAtry4*sinAtry5)* (mu2DotX[i]))/hyper1F1k2 + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6]*atry[6])*(cc*cc)* mu2DotXSq)/ hyper1F1k2

				     (atry[6]*atry[6])*(cc*cc)*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2) + (atry[6]*atry[6])* (-(zz*cosAtry4) - xx*cosAtry5*sinAtry4 - yy*sinAtry4*sinAtry5)* (mu2DotX[i]) * Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2) + (2.0* (atry[6]*atry[6]*atry[6]*atry[6])*(cc*cc)* mu2DotXSq)*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2)


				     ) - invgiSq * dgda[4][i] * dgda[4][i]; 
			
			
	    // In[68]:=
	    // CForm[D[dfda4, atry[5]]]
	    // Out[68]//CForm=
	    d2fda2[4][5] -= invgi * (


// 				     (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(cc)* (dd))/hyper1F1k2 + (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(yy*cosAtry4*cosAtry5 - xx*cosAtry4*sinAtry5)* (mu2DotX[i]))/hyper1F1k2 + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6]*atry[6])*(cc)* (dd)* mu2DotXSq)/ hyper1F1k2

				     (atry[6]*atry[6])*(cc)* (dd)*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2) + (atry[6]*atry[6])*(yy*cosAtry4*cosAtry5 - xx*cosAtry4*sinAtry5)*(mu2DotX[i])*Math.exp(k2Mu2DotXSq[i]-logHyper1F1k2) + 2.0*(atry[6]*atry[6]*atry[6]*atry[6])*(cc)* (dd)*mu2DotXSq*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2)
			
				     ) - invgiSq * dgda[4][i] * dgda[5][i]; 
			
	    // In[69]:=
	    // CForm[D[dfda4, atry[6]]]
	    // Out[69]//CForm=
	    d2fda2[4][6] -= invgi * (


// 				     (2*eK2Mu2DotXSq[i]*atry[6]* (cc)* (mu2DotX[i]))/hyper1F1k2 - (0.6666666666666666*eK2Mu2DotXSq[i]*(atry[6]*atry[6]*atry[6])*dHyperDk2* (cc)* (mu2DotX[i]))/ (hyper1F1k2*hyper1F1k2) + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6])*(cc)* (mu2DotXSq*mu2DotX[i]))/ hyper1F1k2


				     2.0*atry[6]* (cc)* (mu2DotX[i])*Math.exp(k2Mu2DotXSq[i]-logHyper1F1k2) - 0.6666666666666666*(atry[6]*atry[6]*atry[6])*(cc)* (mu2DotX[i])* Math.exp(k2Mu2DotXSq[i] + logDHyperDk2 - (logHyper1F1k2+logHyper1F1k2)) + 2.0*(atry[6]*atry[6]*atry[6])*(cc)* (mu2DotXSq*mu2DotX[i])*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2)

				     ) - invgiSq * dgda[4][i] * dgda[6][i]; 
			
	    // In[71]:=
	    // CForm[D[dfda5, atry[5]]]
	    // Out[71]//CForm=
	    d2fda2[5][5] -= invgi * (

// 				     (eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(dd*dd))/ hyper1F1k2 + (eK2Mu2DotXSq[i]*(atry[6]*atry[6])* (-(xx*cosAtry5*sinAtry4) - yy*sinAtry4*sinAtry5)* (mu2DotX[i]))/hyper1F1k2 + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6]*atry[6])*(dd*dd)* mu2DotXSq)/ hyper1F1k2


				     (atry[6]*atry[6])*(dd*dd)*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2) + (atry[6]*atry[6])* (-(xx*cosAtry5*sinAtry4) - yy*sinAtry4*sinAtry5)* (mu2DotX[i])*Math.exp(k2Mu2DotXSq[i]-logHyper1F1k2) + (2.0*(atry[6]*atry[6]*atry[6]*atry[6])*(dd*dd)* mu2DotXSq)*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2)

				     ) - invgiSq * dgda[5][i] * dgda[5][i]; 
			
	    // In[72]:=
	    // CForm[D[dfda5, atry[6]]]
	    // Out[72]//CForm=
	    d2fda2[5][6] -= invgi * (

// 				     (2*eK2Mu2DotXSq[i]*atry[6]* (dd)* (mu2DotX[i]))/hyper1F1k2 - (0.6666666666666666*eK2Mu2DotXSq[i]*(atry[6]*atry[6]*atry[6])*dHyperDk2* (dd)* (mu2DotX[i]))/ (hyper1F1k2*hyper1F1k2) + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6]*atry[6])*(dd)* (mu2DotXSq*mu2DotX[i]))/ hyper1F1k2

				     2.0*atry[6]* (dd)* (mu2DotX[i])*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2) - 0.6666666666666666*(atry[6]*atry[6]*atry[6])* (dd)* (mu2DotX[i])*Math.exp(k2Mu2DotXSq[i] + logDHyperDk2 - (logHyper1F1k2+logHyper1F1k2)) + 2.0* (atry[6]*atry[6]*atry[6])*(dd)* (mu2DotXSq*mu2DotX[i])*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2)


				     ) - invgiSq * dgda[5][i] * dgda[6][i]; 
			
	    // In[73]:=
	    // CForm[D[dfda5, atry[7]]]
	    // Out[73]//CForm=
			
	    // In[74]:=
	    // CForm[D[dfda6, atry[6]]]
	    // Out[74]//CForm=
	    d2fda2[6][6] -= invgi * ( 


// 				     (-0.3333333333333333*eK2Mu2DotXSq[i]*dHyperDk2)/(hyper1F1k2*hyper1F1k2) + (0.4444444444444444*eK2Mu2DotXSq[i]*(atry[6]*atry[6])*(dHyperDk2*dHyperDk2))/ (hyper1F1k2*hyper1F1k2*hyper1F1k2) - (0.4*eK2Mu2DotXSq[i]* (atry[6]*atry[6])*d2Hyper1F1k2)/(hyper1F1k2*hyper1F1k2) + (eK2Mu2DotXSq[i]* mu2DotXSq)/ hyper1F1k2 - (1.3333333333333333* eK2Mu2DotXSq[i]* (atry[6]*atry[6])*dHyperDk2* mu2DotXSq)/ (hyper1F1k2*hyper1F1k2) + (2*eK2Mu2DotXSq[i]* (atry[6]*atry[6])*(mu2DotXSq*mu2DotXSq))/ hyper1F1k2


				     -0.3333333333333333*Math.exp(k2Mu2DotXSq[i]+logDHyperDk2 - (logHyper1F1k2+logHyper1F1k2)) + 0.4444444444444444*(atry[6]*atry[6])*Math.exp(k2Mu2DotXSq[i]+(logDHyperDk2+logDHyperDk2) - (logHyper1F1k2+logHyper1F1k2+logHyper1F1k2)) - 0.4* (atry[6]*atry[6])*Math.exp(k2Mu2DotXSq[i]+logD2Hyper1F1k2 - (logHyper1F1k2+logHyper1F1k2)) + mu2DotXSq*Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2) - 1.3333333333333333*mu2DotXSq * (atry[6]*atry[6]) * Math.exp(k2Mu2DotXSq[i] + logDHyperDk2 - (logHyper1F1k2+logHyper1F1k2)) + 2.0* (atry[6]*atry[6])*(mu2DotXSq*mu2DotXSq) * Math.exp(k2Mu2DotXSq[i] - logHyper1F1k2)

			
				     ) - invgiSq * dgda[6][i] * dgda[6][i]; 
			
			
			
		
	}
		
	for (int i = 1; i <= ma; i++) {
	    for (int j = 1; j <= i; j++) {
		d2fda2[i][j] = d2fda2[j][i]; 
	    }
	}
		
    }

	
	 
    public TwoFibreBipolarWatsonFitter(Vector3D[] axes) {
	super(axes, 0.5);
	k1Mu1DotXSq = new double[axes.length];
	k2Mu2DotXSq = new double[axes.length];
    }


    /**
     * Sets the initial values of the parameters. The values in aInit
     * start counting from zero.       
     * 
     * @param aInit Array containing the new parameter values starting
     * from index 0: {theta1, phi1, kappa1, theta2, phi2, kappa2}. (theta, phi) are the spherical
     * polar coordinates of the mean axes of the distributions. 
     */
    public void setInitParams(double[] aInit) throws MarquardtMinimiserException {
	double[] realA = new double[aInit.length];
	for (int i = 0; i < realA.length; i++) {
	    realA[i] = aInit[i];
	}

	realA[2] = Math.sqrt(realA[2]);
	realA[5] = Math.sqrt(realA[5]);

	super.setInitParams(realA);
    }

    /**
     * Set an initial estimate of the parameters. 
     *
     * @param k1 >= 0.0
     * @param k2 >= 0.0
     */
    public void setInitParams(Vector3D mu1, Vector3D mu2, double k1, double k2) throws MarquardtMinimiserException {
	double[] tp = Vector3D.thetaPhi(mu1);

	double[] aInit = new double[6];

	aInit[0] = tp[0];
	aInit[1] = tp[1];
	aInit[2] = Math.sqrt(k1);

	tp = Vector3D.thetaPhi(mu2);
	aInit[3] = tp[0];
	aInit[4] = tp[1];
	aInit[5] = Math.sqrt(k2);

	setInitParams(aInit);
	
    }

    public double[] getParameters() {
		
	double[] param = new double[ma];
		
	for (int i = 0;  i < param.length; i++) {
	    param[i] = a[i+1];
	}
	
	param[2] = param[2] * param[2];
	param[5] = param[5] * param[5];
	
	return param;
		
    }



    /**
     * @return the fitted Watson Distributions, with a specified random number source.
     *
     */
    public WatsonDistribution[] getDistributions(java.util.Random r) {
		
	Vector3D mu1 = Vector3D.vectorFromSPC(1.0, a[1], a[2]);
	Vector3D mu2 = Vector3D.vectorFromSPC(1.0, a[4], a[5]);
		
	return new WatsonDistribution[] {new WatsonDistribution(mu1, a[3]*a[3], r), 
					 new WatsonDistribution(mu2, a[6]*a[6], r)};
		
    }


    /** 
     * @return {k1, k2}, which correspond to W1 and W2, such that the minimized fObj is alpha * W1 + (1-alpha) * W2.
     */
    public double[] getKappas() {
	return new double[] {a[3] * a[3], a[6] * a[6]};
    }	

    public static void main(String[] args) {
		
	if (args.length == 0) {
	    System.out.println("Usage: TwoFibreBipolarWatsonFitter -sampleAxes [filename] -initMu1 [theta phi] -initMu2 [theta phi] -initK1 [value] -initK2 [value]\n\nThe axes file should contain a list of unit vectors with x y z delimited by spaces and each vector on a separate line. Parameters default to zero if arguments are missing.");
	    System.exit(0);
	}
		
	Vector3D[] axes = null;
		
	double[] initParams = new double[6];
		
	double mixingPar = 0.5; 
	for (int i = 0; i < args.length; i++) {
			
	    if (args[i].toLowerCase().equals("-sampleaxes")) {
				
		FileInput in = new FileInput(args[i+1]);
				
		ArrayList<Vector3D> a = new ArrayList<Vector3D>();
						
		while (!in.eof()) {
					
		    String[] axisString = in.readString().split(" ");
		    double x = Double.parseDouble(axisString[0]);
		    double y = Double.parseDouble(axisString[1]);
		    double z = Double.parseDouble(axisString[2]);
					
		    a.add(new Vector3D(x,y,z));
					
		}
				
	    }
	    else if (args[i].toLowerCase().equals("-initmu1")) {
		initParams[0] = Double.parseDouble(args[i+1]);
		initParams[1] = Double.parseDouble(args[i+2]);
	    }
	    else if (args[i].toLowerCase().equals("-initk1")) {
		initParams[2] = Double.parseDouble(args[i+1]);
	    }
	    else if (args[i].toLowerCase().equals("-initmu2")) {
		initParams[3] = Double.parseDouble(args[i+1]);
		initParams[4] = Double.parseDouble(args[i+2]);
	    }
	    else if (args[i].toLowerCase().equals("-initk2")) {
		initParams[5] = Double.parseDouble(args[i+1]);
	    }
	    else if (args[i].toLowerCase().equals("-alpha")) {
		mixingPar = Double.parseDouble(args[i+1]);
	    }
			
			
			
	}
		
		
	TwoFibreBipolarWatsonFitter fitter = new TwoFibreBipolarWatsonFitter(axes);
		
	try {
	    fitter.setInitParams(initParams);
	    fitter.minimise();
	}
	catch (MarquardtMinimiserException e) {
	    System.out.println("Minimization failed with exception: " + e);
	    System.exit(1);
	}
		
	System.out.println("Fitted parameters");
		
	Vector3D[] mus = fitter.getMus();
	double[] kappas = fitter.getKappas();
	double alpha = fitter.getMixingParameter();
		
	System.out.println("mu1 " + mus[0]);
	System.out.println("kappa1 " + kappas[0]);
		
		
	System.out.println("\nmu2 " + mus[1]);
	System.out.println("kappa2 " + kappas[1]);
		
	System.out.println("alpha " + alpha);




	
    }



}

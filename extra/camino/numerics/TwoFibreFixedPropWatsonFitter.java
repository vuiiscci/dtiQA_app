package numerics;

import tools.*;

import misc.LoggedException;
import optimizers.*;
import java.util.ArrayList;
import java.util.logging.Logger;

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
* <dd> This class takes a collection of axes representing a set of estimates of two fibre directions. It uses the Marquardt optimization routine to find the pair of Watson distributions that best fit the axes.
* </dl>
* <p>
* <dl>
* <dd> Some mathematical definitions used in the code:
* <dd> g(x) = \alpha W_1(\theta_1, \phi_1, \kappa_1, x) + (1-\alpha)W_2(\theta_2, \phi_2, \kappa_2, x)
* <dd> x = a vector on the unit sphere
* <dd> f = the objective function -\sum_i \log(g(x_i)) 
* <dd> W = Watson distribution
* <dd> alpha = a mixing parameter. 0 <= alpha <= 1
* </dl>
*
* @version $Id$
* @author Philip Cook
* 
*/
public class TwoFibreFixedPropWatsonFitter extends MarquardtMinimiser {
	
   	/** logging object */
    private static final Logger logger = Logger.getLogger("numerics.TwoFibreFixedPropWatsonFitter");
    
    
    // sample vectors
    Vector3D[] samples;


    final double alpha;
	
    // a[1]     a[2]   a[3]    a[4]    a[5]   a[6]   
    // theta1   phi1   kappa1  theta2  phi2   kappa2 
	

    protected final double[] mu1DotX;
    protected final double[] mu2DotX;

    protected final double[] eK1Mu1DotXSq;
    protected final double[] eK2Mu2DotXSq;


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

	double hyper1F1k1 = WatsonDistribution.hyper1F1(0.5, 1.5, atry[3], 1e-9);
	double hyper1F1k2 = WatsonDistribution.hyper1F1(0.5, 1.5, atry[6], 1e-9);
	

	for (int i = 0; i < samples.length; i++) {
			
	    mu1DotX[i] = mu1.x * samples[i].x + mu1.y * samples[i].y + mu1.z * samples[i].z;
	    mu2DotX[i] = mu2.x * samples[i].x + mu2.y * samples[i].y + mu2.z * samples[i].z;

	    eK1Mu1DotXSq[i] = Math.exp(atry[3] * mu1DotX[i] * mu1DotX[i]);
	    eK2Mu2DotXSq[i] = Math.exp(atry[6] * mu2DotX[i] * mu2DotX[i]);
	    
	    g[i] = alpha * (1.0 / hyper1F1k1) * eK1Mu1DotXSq[i] + (1.0 - alpha) * (1.0 / hyper1F1k2) * eK2Mu2DotXSq[i];
	    
	    minusSumLogPDF -= Math.log(g[i]);
	}

	double dHyperDk1 = WatsonDistribution.hyper1F1(1.5, 2.5, atry[3], 1e-9);
	double dHyperDk2 = WatsonDistribution.hyper1F1(1.5, 2.5, atry[6], 1e-9);

	firstDeriv(atry,g, dgda, dfda, mu1DotX, mu2DotX, eK1Mu1DotXSq, eK2Mu2DotXSq, hyper1F1k1, hyper1F1k2, dHyperDk1, dHyperDk2);
		
	secondDeriv(atry,g, dgda, d2fda2, mu1DotX, mu2DotX, eK1Mu1DotXSq, eK2Mu2DotXSq, hyper1F1k1, hyper1F1k2, dHyperDk1, dHyperDk2);
		
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
			
	    double x = samples[i].x;
	    double y = samples[i].y;
	    double z = samples[i].z;


	    //	    double eK1Mu1DotXSq = Math.exp(atry[3]*mu1DotX[i]*mu1DotX[i]);
	    //	    double eK2Mu2DotXSq = Math.exp(atry[6]*mu2DotX[i]*mu2DotX[i]);

	    dgda[1][i] = (2.0*eK1Mu1DotXSq[i]*atry[3]*(x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2)*mu1DotX[i]*alpha)/hyper1F1k1;
			
			
	    dfda[1] -= (1.0 / g[i]) * dgda[1][i];
			
	    dgda[2][i] = (2.0*eK1Mu1DotXSq[i]*atry[3]*(y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2)*mu1DotX[i]*alpha)/hyper1F1k1;
			
	    dfda[2] -= (1.0 / g[i]) * dgda[2][i];
			
	    dgda[3][i] = (-(1.0 / 3.0) *eK1Mu1DotXSq[i]*dHyperDk1*alpha)/(hyper1F1k1 * hyper1F1k1) + (eK1Mu1DotXSq[i]*mu1DotX[i]*mu1DotX[i]*alpha)/hyper1F1k1;
			 
	    dfda[3] -= (1.0 / g[i]) * dgda[3][i];
			
	    dgda[4][i] = (2.0*eK2Mu2DotXSq[i]*atry[6]*(x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5)*mu2DotX[i]*(1.0 - alpha))/hyper1F1k2;
			
	    dfda[4] -= (1.0 / g[i]) * dgda[4][i];
			
			
	    dgda[5][i] = (2.0*eK2Mu2DotXSq[i]*atry[6]*(y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*(1.0 - alpha))/hyper1F1k2;
			
			
	    dfda[5] -= (1.0 / g[i]) * dgda[5][i];
			
	    dgda[6][i] = (-(1.0 / 3.0)*eK2Mu2DotXSq[i]*dHyperDk2*(1.0 - alpha))/(hyper1F1k2 * hyper1F1k2) + (eK2Mu2DotXSq[i]*mu2DotX[i]*mu2DotX[i]*(1.0 - alpha))/hyper1F1k2;
			
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

	double d2Hyper1F1k1 = WatsonDistribution.hyper1F1(2.5,3.5,atry[3], 1e-9);
	double d2Hyper1F1k2 = WatsonDistribution.hyper1F1(2.5,3.5,atry[6], 1e-9);
		
	double cosAtry1 = Math.cos(atry[1]);
	double sinAtry1 = Math.sin(atry[1]);
		
	double cosAtry2 = Math.cos(atry[2]);
	double sinAtry2 = Math.sin(atry[2]);
		
	double cosAtry4 = Math.cos(atry[4]);
	double sinAtry4 = Math.sin(atry[4]);
		
	double cosAtry5 = Math.cos(atry[5]);
	double sinAtry5 = Math.sin(atry[5]);

	double k1Sq = atry[3] * atry[3];
	double k2Sq = atry[6] * atry[6];
		
	for (int i = 0; i < samples.length; i++) {
			
	    double x = samples[i].x;
	    double y = samples[i].y;
	    double z = samples[i].z;
			
	    double invgi = 1.0 / g[i];
	    double invgiSq = invgi * invgi;
	
	    double a = x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2;
	    double b = mu1DotX[i];
	    double c = x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5;
	    double d = mu2DotX[i];
	    double e = y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2;
	    double f = y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5;

	    double mu1DotXiSq = b * b;
	    double mu2DotXiSq = d * d;

	    //	    double eK1Mu1DotXSq = Math.exp(atry[3]*mu1DotXiSq);
	    //	    double eK2Mu2DotXSq = Math.exp(atry[6]*mu2DotXiSq);

	    double eK1Mu1DotXSqAtry3x2 = 2.0*eK1Mu1DotXSq[i]*atry[3];
	   
	    // In[48]:= CForm[D[dfda1, atry[1]]] 
	    // Out[48]//CForm=
	    d2fda2[1][1] -= invgi * (
				     (eK1Mu1DotXSqAtry3x2 * a * a * alpha)/hyper1F1k1 + (eK1Mu1DotXSqAtry3x2 *
				      (-1.0 * b)*b* alpha)/hyper1F1k1 + (4.0 * eK1Mu1DotXSq[i]*k1Sq* a * a * mu1DotXiSq*alpha)/hyper1F1k1
				     ) - invgiSq * dgda[1][i] * dgda[1][i];
			
	    // In[49]:= CForm[D[dfda1, atry[2]]]
	    // Out[49]//CForm=
	    d2fda2[1][2] -= invgi * (
				     (eK1Mu1DotXSqAtry3x2 *
				      a*e*alpha)/
				     hyper1F1k1 + (eK1Mu1DotXSqAtry3x2*
						   (y*cosAtry1*cosAtry2 - x*cosAtry1*sinAtry2)*b*alpha)/
				     hyper1F1k1 + (4.0*eK1Mu1DotXSq[i]*k1Sq*
						   a*e*
						   mu1DotXiSq*alpha)/hyper1F1k1
				     ) - invgiSq * dgda[1][i] * dgda[2][i];
			
			
	    // In[50]:=
	    // CForm[D[dyda1, atry[3]]]
	    // Out[50]//CForm=
	    d2fda2[1][3] -= invgi * ( 
				     (2.0*eK1Mu1DotXSq[i] *
				      a*b*
				      alpha)/hyper1F1k1 - (0.6666666666666666*
							   eK1Mu1DotXSq[i]*atry[3]*dHyperDk1*
							   a*b* alpha)/(hyper1F1k1 * hyper1F1k1) + 
				     (eK1Mu1DotXSqAtry3x2 *a*mu1DotXiSq*b*alpha)/hyper1F1k1
				     
				     ) - invgiSq * dgda[1][i] * dgda[3][i];
			
			
// 	    // In[52]:=
// 	    // CForm[D[dfda1, atry[4]]]
// 	    // Out[52]//CForm=
// 	    d2fda2[1][4] -= 0.0;
			
// 	    // In[53]:=
// 	    // CForm[D[dfda1, atry[5]]]
// 	    // Out[53]//CForm=
// 	    d2fda2[1][5] -= 0.0;
			
// 	    // In[54]:=
// 	    // CForm[D[dfda1, atry[6]]]
// 	    // Out[54]//CForm=
// 	    d2fda2[1][6] -= 0.0;
			
			
	    //In[56]:=
	    // CForm[D[dfda2, atry[2]]]
	    // Out[56]//CForm=
	    d2fda2[2][2] -= invgi * (
				     (eK1Mu1DotXSqAtry3x2*
				      e*e*alpha)/hyper1F1k1 + 
				     (eK1Mu1DotXSqAtry3x2*
				      (-(x*cosAtry2*sinAtry1) - y*sinAtry1*sinAtry2)*b*alpha)/
				     hyper1F1k1 + (4*eK1Mu1DotXSq[i]*k1Sq*
						   e*e*mu1DotXiSq*
						   alpha)/hyper1F1k1
				     ) - invgiSq * dgda[2][i] * dgda[2][i]; 
			
	    // In[57]:=
	    // CForm[D[dfda2, atry[3]]]
	    // Out[57]//CForm=
	    d2fda2[2][3] -= invgi * (
				     (2.0*eK1Mu1DotXSq[i]*e*
				      b*alpha)/hyper1F1k1 - 
				     ((2.0 / 3.0)*eK1Mu1DotXSq[i]*atry[3]*dHyperDk1*
				      e*b*alpha)/
				     (hyper1F1k1 * hyper1F1k1) + (eK1Mu1DotXSqAtry3x2*
							       e*mu1DotXiSq*b*
							       alpha)/hyper1F1k1
				     ) - invgiSq * dgda[2][i] * dgda[3][i]; 
			
// 	    // In[58]:=
// 	    // CForm[D[dfda2, atry[4]]]
// 	    // Out[58]//CForm=
// 	    d2fda2[2][4] -= 0.0;
			
// 	    // In[59]:=
// 	    // CForm[D[dfda2, atry[5]]]
// 	    // Out[59]//CForm=
// 	    d2fda2[2][5] -= 0.0;
			
// 	    // In[60]:=
// 	    // CForm[D[dfda2, atry[6]]]
// 	    // Out[60]//CForm=
// 	    d2fda2[2][6] -= 0.0;
			
			
	    // In[62]:=
	    // CForm[D[dfda3, atry[3]]]
	    // Out[62]//CForm=
	    d2fda2[3][3] -= invgi * (
				     ( (2.0 / 9.0)*eK1Mu1DotXSq[i]*(dHyperDk1*dHyperDk1)*
				       alpha)/(hyper1F1k1*hyper1F1k1*hyper1F1k1) - 
				     (0.2*eK1Mu1DotXSq[i]*d2Hyper1F1k1*alpha)/
				     (hyper1F1k1 * hyper1F1k1) - ((2.0 / 3.0)*
							       eK1Mu1DotXSq[i]*dHyperDk1*
							       mu1DotXiSq*alpha)/(hyper1F1k1 * hyper1F1k1) + 
				     (eK1Mu1DotXSq[i]*
				      mu1DotXiSq*mu1DotXiSq*alpha)/hyper1F1k1
				     ) - invgiSq * dgda[3][i] * dgda[3][i]; 
			
// 	    // In[63]:=
// 	    // CForm[D[dfda3, atry[4]]]
// 	    // Out[63]//CForm=
// 	    d2fda2[3][4] -= 0.0;
			
// 	    // In[64]:=
// 	    // CForm[D[dfda3, atry[5]]]
// 	    // Out[64]//CForm=
// 	    d2fda2[3][5] -= 0.0;
			
// 	    // In[65]:=
// 	    // CForm[D[dfda3, atry[6]]]
// 	    // Out[65]//CForm=
// 	    d2fda2[3][6] -= 0.0;
			
			
	    // In[67]:=
	    // CForm[D[dfda4, atry[4]]]
	    // Out[67]//CForm=
	    d2fda2[4][4] -= invgi * (
				     (2.0*eK2Mu2DotXSq[i]*atry[6]*
				      (c*c)*(1.0 - alpha))/hyper1F1k2 + 
				     (2*eK2Mu2DotXSq[i]*atry[6]*
				      (-(z*cosAtry4) - x*cosAtry5*sinAtry4 - y*sinAtry4*sinAtry5)*d*
				      (1.0 - alpha))/hyper1F1k2 + 
				     (4.0*eK2Mu2DotXSq[i]*k2Sq*
				      (c*c)*
				      mu2DotXiSq*(1.0 - alpha))/hyper1F1k2
				     ) - invgiSq * dgda[4][i] * dgda[4][i]; 
			
			
	    // In[68]:=
	    // CForm[D[dfda4, atry[5]]]
	    // Out[68]//CForm=
	    d2fda2[4][5] -= invgi * (
				     (2.0*eK2Mu2DotXSq[i]*atry[6]*
				      (c)*f*
				      (1.0 - alpha))/hyper1F1k2 + 
				     (2.0*eK2Mu2DotXSq[i]*atry[6]*
				      (y*cosAtry4*cosAtry5 - x*cosAtry4*sinAtry5)*d*
				      (1.0 - alpha))/hyper1F1k2 + 
				     (4.0*eK2Mu2DotXSq[i]*k2Sq*
				      c*f*
				      mu2DotXiSq*(1.0 - alpha))/hyper1F1k2
			
				     ) - invgiSq * dgda[4][i] * dgda[5][i]; 
			
	    // In[69]:=
	    // CForm[D[dfda4, atry[6]]]
	    // Out[69]//CForm=
	    d2fda2[4][6] -= invgi * (
				     (2*eK2Mu2DotXSq[i]*
				      (c)*d*
				      (1.0 - alpha))/hyper1F1k2 - 
				     ((2.0 / 3.0)*eK2Mu2DotXSq[i]*atry[6]*dHyperDk2*
				      (c)*d*
				      (1.0 - alpha))/(hyper1F1k2 * hyper1F1k2) + 
				     (2*eK2Mu2DotXSq[i]*atry[6]*
				      c*mu2DotXiSq*d*
				      (1.0 - alpha))/hyper1F1k2
				     ) - invgiSq * dgda[4][i] * dgda[6][i]; 
			
	    // In[70]:=
	    // CForm[D[dfda4, atry[7]]]
			
	    // In[71]:=
	    // CForm[D[dfda5, atry[5]]]
	    // Out[71]//CForm=
	    d2fda2[5][5] -= invgi * (
				     (2.0*eK2Mu2DotXSq[i]*atry[6]*
				      f*f*(1.0 - alpha))/hyper1F1k2 + 
				     (2.0*eK2Mu2DotXSq[i]*atry[6]*
				      (-(x*cosAtry5*sinAtry4) - y*sinAtry4*sinAtry5)*d*
				      (1.0 - alpha))/hyper1F1k2 + 
				     (4.0*eK2Mu2DotXSq[i]*k2Sq*
				      f*f*mu2DotXiSq*
				      (1.0 - alpha))/hyper1F1k2
				     ) - invgiSq * dgda[5][i] * dgda[5][i]; 
			
	    // In[72]:=
	    // CForm[D[dfda5, atry[6]]]
	    // Out[72]//CForm=
	    d2fda2[5][6] -= invgi * (
				     (2.0*eK2Mu2DotXSq[i]*f*
				      d*(1.0 - alpha))/hyper1F1k2 - 
				     ((2.0 / 3.0)*eK2Mu2DotXSq[i]*atry[6]*dHyperDk2*
				      f*d*
				      (1.0 - alpha))/(hyper1F1k2*hyper1F1k2) + 
				     (2*eK2Mu2DotXSq[i]*atry[6]*
				      f*mu2DotXiSq*d*
				      (1.0 - alpha))/hyper1F1k2
				     ) - invgiSq * dgda[5][i] * dgda[6][i]; 
			
	    // In[73]:=
	    // CForm[D[dfda5, atry[7]]]
	    // Out[73]//CForm=
			
	    // In[74]:=
	    // CForm[D[dfda6, atry[6]]]
	    // Out[74]//CForm=
	    d2fda2[6][6] -= invgi * ( 
				     ((2.0 / 9.0)*eK2Mu2DotXSq[i]*(dHyperDk2*dHyperDk2)*
				      (1.0 - alpha))/(hyper1F1k2*hyper1F1k2*hyper1F1k2) - 
				     (0.2*eK2Mu2DotXSq[i]*d2Hyper1F1k2*
				      (1.0 - alpha))/(hyper1F1k2*hyper1F1k2) - 
				     ((2.0 / 3.0)*eK2Mu2DotXSq[i]*dHyperDk2*
				      mu2DotXiSq*(1.0 - alpha))/(hyper1F1k2*hyper1F1k2) + 
				     (eK2Mu2DotXSq[i]*
				      mu2DotXiSq*mu2DotXiSq*(1.0 - alpha))/hyper1F1k2
			
				     ) - invgiSq * dgda[6][i] * dgda[6][i]; 
			
			
			
		
	}
		
	for (int i = 1; i <= ma; i++) {
	    for (int j = 1; j <= i; j++) {
		d2fda2[i][j] = d2fda2[j][i]; 
	    }
	}
		
		
    }
	
	
    /**
     * @return the fitted Watson Distributions
     *
     */
    public WatsonDistribution[] getDistributions() {
		
	java.util.Random r = new java.util.Random();

	return getDistributions(r);
    }


    /**
     * @return the fitted Watson Distributions, with a specified random number source.
     *
     */
    public WatsonDistribution[] getDistributions(java.util.Random r) {
		
	Vector3D mu1 = Vector3D.vectorFromSPC(1.0, a[1], a[2]);
	Vector3D mu2 = Vector3D.vectorFromSPC(1.0, a[4], a[5]);
		
	return new WatsonDistribution[] {new WatsonDistribution(mu1, a[3], r), 
					 new WatsonDistribution(mu2, a[6], r)};
		
    }
	
    /**
     * @return {mu1, mu2}, which correspond to W1 and W2, such that the minimized fObj is alpha * W1 + (1-alpha) * W2.
     *
     */
    public Vector3D[] getMus() {
		
	Vector3D mu1 = Vector3D.vectorFromSPC(1.0, a[1], a[2]);
	Vector3D mu2 = Vector3D.vectorFromSPC(1.0, a[4], a[5]);
		
	return new Vector3D[] {mu1, mu2};
		
    }
	
	
    /** 
     * @return alpha, such that the minimized fObj is alpha * W1 + (1-alpha) * W2
     */
    public double getMixingParameter() {
	return alpha;
    }
	
	
    /** 
     * @return {k1, k2}, which correspond to W1 and W2, such that the minimized fObj is alpha * W1 + (1-alpha) * W2.
     */
    public double[] getKappas() {
	return new double[] {a[3], a[6]};
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
	// interception
	double[] realA = new double[aInit.length];
	for (int i = 0; i < realA.length; i++) {
	    realA[i] = aInit[i];
	}
	super.setInitParams(realA);
    }
	
    
    /**
     * Set an initial estimate of the parameters.
     *
     */
    public void setInitParams(Vector3D mu1, Vector3D mu2, double k1, double k2) throws MarquardtMinimiserException {
	double[] tp = Vector3D.thetaPhi(mu1);

	double[] aInit = new double[6];

	aInit[0] = tp[0];
	aInit[1] = tp[1];
	aInit[2] = k1;

	tp = Vector3D.thetaPhi(mu2);
	aInit[3] = tp[0];
	aInit[4] = tp[1];
	aInit[5] = k2;

	setInitParams(aInit);
	
    }



    /**
     * Make estimates of the initial parameters and attempt a minimisation with each. The fitter will retain the best fObj obtained with the estimated initial parameters.
     * <p> 
     * Two cones are defined for each mean axis, centered on the mean axis. The fitter will attempt a minimization with the mean axes oriented on the radius of each of the cones. 
     *
     * @param attemptsOnCone the number of sets of parameters to try on the inner and outer cone.
 
     */
    public void fitEstimatedParams(Vector3D mu1, Vector3D mu2, double innerCone1, double outerCone1, double innerCone2, double outerCone2, double initK1, double initK2, int attemptsOnCone) throws MarquardtMinimiserException {


 	double[] initParams = new double[6];
	
 	initParams[2] = initK1;
 	initParams[5] = initK2;
	
	double bestFObjVal = Double.MAX_VALUE;
	double[] bestFittedParams = new double[6];
	boolean convergedAtLeastOnce = false;


	double[] tpMu1 = Vector3D.thetaPhi(mu1);
	double[] tpMu2 = Vector3D.thetaPhi(mu2);

	RealMatrix rot1 = Rotations.getRotMat(mu1, 2.0 * Math.PI / attemptsOnCone);
	RealMatrix rot2 = Rotations.getRotMat(mu2, 2.0 * Math.PI / attemptsOnCone);

	// try initial guess
	initParams[0] = tpMu1[0];
	initParams[1] = tpMu1[1];
	
	initParams[3] = tpMu2[0];
	initParams[4] = tpMu2[1];

	for(int i = 0; i < 6; i++) {
	    bestFittedParams[i] = initParams[i];
	}

	try {
	    setInitParams(initParams); // important to do this for subclass compatibility
	    minimise();
	    convergedAtLeastOnce = true;
	    
	}
	catch (MarquardtMinimiserException e) {
	    //  System.out.println(e);
	}
	
	if (fObjVal < bestFObjVal) {
	    bestFObjVal = fObjVal;
	    bestFittedParams = getParameters();
	}


	// try inner cone
	double[] tp1 = new double[2];
	double[] tp2 = new double[2];

	tp1[0] = tpMu1[0] + innerCone1;
	tp1[1] = tpMu1[1];

	if (tp1[0] > Math.PI) {
	    tp1[0] -= Math.PI; 
	}

	Vector3D initMu1 = Vector3D.vectorFromSPC(1.0,tp1[0], tp1[1]); 	


	tp2[0] = tpMu2[0] + innerCone2;
	tp2[1] = tpMu2[1];

	if (tp2[0] > Math.PI) {
	    tp2[0] -= Math.PI; 
	}
	
	Vector3D initMu2 = Vector3D.vectorFromSPC(1.0,tp2[0], tp2[1]); 	
	
		
	for (int i = 0; i < attemptsOnCone; i++) {
	    
	    if (i > 0) {

		initMu1 = Vector3D.vectorFromSPC(1.0,tp1[0], tp1[1]); 	
		initMu1 = Rotations.rotateVector(initMu1, rot1);
		tp1 = Vector3D.thetaPhi(initMu1);


		initMu2 = Vector3D.vectorFromSPC(1.0,tp2[0], tp2[1]); 	
		initMu2 = Rotations.rotateVector(initMu2, rot2);
		tp2 = Vector3D.thetaPhi(initMu2);
	    }
	    
	    initParams[0] = tp1[0];
	    initParams[1] = tp1[1];

	    initParams[3] = tp2[0];
	    initParams[4] = tp2[1];
	    
	    
	    try {
		setInitParams(initParams);
		minimise();
		convergedAtLeastOnce = true;
		
	    }
	    catch (MarquardtMinimiserException e) {
		//		System.out.println(e);
	    }
	    
	    if (fObjVal < bestFObjVal) {
		bestFObjVal = fObjVal;
		bestFittedParams = getParameters();
	    }
	    
	}


	// try outer cone
	tp1[0] = tpMu1[0] + outerCone1;
	tp1[1] = tpMu1[1];

	if (tp1[0] > Math.PI) {
	    tp1[0] -= Math.PI; 
	}

	initMu1 = Vector3D.vectorFromSPC(1.0, tp1[0], tp1[1]); 	

	// offset by half a rotation
	initMu1 = Rotations.rotateVector(initMu1, mu1, Math.PI / attemptsOnCone);


	tp2[0] = tpMu2[0] + outerCone2;
	tp2[1] = tpMu2[1];

	if (tp2[0] > Math.PI) {
	    tp2[0] -= Math.PI; 
	}

	initMu2 = Vector3D.vectorFromSPC(1.0, tp2[0], tp2[1]); 	
	
	initMu2 = Rotations.rotateVector(initMu2, mu2, Math.PI / attemptsOnCone);
		
	for (int i = 0; i < attemptsOnCone; i++) {
	    
	    if (i > 0) {

		initMu1 = Vector3D.vectorFromSPC(1.0,tp1[0], tp1[1]); 	
		initMu1 = Rotations.rotateVector(initMu1, rot1);
		tp1 = Vector3D.thetaPhi(initMu1);


		initMu2 = Vector3D.vectorFromSPC(1.0,tp2[0], tp2[1]); 	
		initMu2 = Rotations.rotateVector(initMu2, rot2);
		tp2 = Vector3D.thetaPhi(initMu2);
	    }
	    
	    initParams[0] = tp1[0];
	    initParams[1] = tp1[1];

	    initParams[3] = tp2[0];
	    initParams[4] = tp2[1];
	    
	    
	    try {
		setInitParams(initParams);
		minimise();
		convergedAtLeastOnce = true;
		
	    }
	    catch (MarquardtMinimiserException e) {
		//		System.out.println(e);
	    }
	    
	    if (fObjVal < bestFObjVal) {
		bestFObjVal = fObjVal;
		bestFittedParams = getParameters();
	    }
	    
	}


	setInitParams(bestFittedParams);	

	
	if (convergedAtLeastOnce == false) {
	    throw new MarquardtMinimiserException("None of the parameters could converge");
	}
	
	
    }

    /**
     * Two cones are defined for each mean axis, centered on the mean axis. The fitter will attempt a minimization with the mean axes oriented on the radius of each of the cones. 
     * <p>
     * This method attmempts to guess the initial values of kappa, and the semi-angle of the cones. The cone radii come from the spherical standard error of the mean axes when the vectors are sorted into two separate distributions (the vectors are sorted to maximise concentration about mu1 and mu2).
     * 
     * @param mu1 an estimate of the first principal axis.
     * @param mu2 an estimate of the second principal axis.
     *
     * @param attemptsOnCone the number of mean axes to try one each cone. Starting points will be distributed equally around the radius of the cone.
 
     */
    public void fitEstimatedParams(Vector3D mu1, Vector3D mu2, int attemptsOnCone) throws MarquardtMinimiserException {
	
	
	Vector3D[] dt1Samples = new Vector3D[samples.length];
	Vector3D[] dt2Samples = new Vector3D[samples.length];

	int dt1Counter = 0;
	int dt2Counter = 0;
	
	for (int i = 0; i < samples.length; i++) {
	    Vector3D sample = samples[i];
	    
	    // Given mu1 and mu2, we can just assign to one or the other, 
	    // no need to keep an even division. This probably leads to bias, 
	    // but the optimization should sort that out.

	    // If we assign the samples in pairs, one to dt1 and one to dt2, 
	    // we need to know the way the samples are ordered. 
	    // This can easily lead to problems if someone passes axes ordered in the wrong way.
	
	    if (Math.abs(mu1.dot(sample)) > Math.abs(mu2.dot(sample))) {
		dt1Samples[dt1Counter++] = sample;
	    }
	    else {
		dt2Samples[dt2Counter++] = sample;
	    }

	}

	Vector3D[] tmpDT1 = new Vector3D[dt1Counter];
	Vector3D[] tmpDT2 = new Vector3D[dt2Counter];


	for (int i = 0; i < dt1Counter; i++) {
	    tmpDT1[i] = dt1Samples[i];
	}
	for (int i = 0; i < dt2Counter; i++) {
	    tmpDT2[i] = dt2Samples[i];
	}
	
	dt1Samples = tmpDT1;
	dt2Samples = tmpDT2;
	
	EigenSystem3D dt1Eig = WatsonFitter.tBarEigenSystem(dt1Samples);
	EigenSystem3D dt2Eig = WatsonFitter.tBarEigenSystem(dt2Samples);
	
	double innerConeRadius1 = WatsonFitter.getBipolarConfidenceCone(dt1Eig, dt1Samples, 0.25) * 180.0 / Math.PI;
	
	double outerConeRadius1 = WatsonFitter.getBipolarConfidenceCone(dt1Eig, dt1Samples, 0.05) * 180.0 / Math.PI;
	
	double innerConeRadius2 = WatsonFitter.getBipolarConfidenceCone(dt2Eig, dt2Samples, 0.25) * 180.0 / Math.PI;

	double outerConeRadius2 = WatsonFitter.getBipolarConfidenceCone(dt2Eig, dt2Samples, 0.05) * 180.0 / Math.PI;

	// guess concentration
	double initK1 = 20.0;

	double initK2 = initK1;


        try {
            fitEstimatedParams(dt1Eig.eigenvectors[0], dt2Eig.eigenvectors[0], innerConeRadius1, 
                               outerConeRadius1, innerConeRadius2, outerConeRadius2, 
                               initK1, initK2, attemptsOnCone);
        }
        catch (MarquardtMinimiserException e) {
            double k1 = WatsonFitter.fitKappa(dt1Eig, dt1Samples); 
            double k2 = WatsonFitter.fitKappa(dt2Eig, dt2Samples); 
            
            if (k1 < 0.0) {
                k1 = 0.0;
            }
            if (k2 < 0.0) {
                k2 = 0.0;
            }
            
            setInitParams(dt1Eig.eigenvectors[0], dt2Eig.eigenvectors[0], k1, k2);
            throw e;
        }

        double[] kappas = getKappas();
     
        // check for infinity
        if (Double.isInfinite(kappas[0]) || Double.isInfinite(kappas[1])) {

            double k1 = WatsonFitter.fitKappa(dt1Eig, dt1Samples); 
            double k2 = WatsonFitter.fitKappa(dt2Eig, dt2Samples); 
            
            if (k1 < 0.0) {
                k1 = 0.0;
            }
            if (k2 < 0.0) {
                k2 = 0.0;
            }
            
            setInitParams(dt1Eig.eigenvectors[0], dt2Eig.eigenvectors[0], k1, k2);

            throw new MarquardtMinimiserException("Converged but concentration was infinite");
        }
    }

    

    public double[] getParameters() {
		
	double[] param = new double[ma];
		
	for (int i = 0;  i < param.length; i++) {
	    param[i] = a[i+1];
	}
		
	return param;
		
    }


    /**
     * @param ai an array of boolean. Parameter i is fixed if ai[i] == false.
     * 
     */
    public void setFixedParams(boolean[] ai) {
	ia = new boolean[ai.length + 1];
	ia[0] = true;
	for (int i = 0; i < ai.length; i++) {
	    ia[i+1] = ai[i];
	}
    }

  
    /**
     * @param variable true if the principal axis mu1 is to be variable (default is variable)
     *
     */
    public void setMu1Variable(boolean variable) {
	ia[1] = variable;
	ia[2] = variable;
    }

    /**
     * @param variable true if the principal axis mu2 is to be variable (default is variable)
     *
     */
    public void setMu2Variable(boolean variable) {
	ia[4] = variable;
	ia[5] = variable;
    }

    /**
     * @param variable true if kappa1 is to be variable (default is variable)
     *
     */
    public void setKappa1Variable(boolean variable) {
	ia[3] = variable;
    }

    /**
     * @param variable true if kappa2 is to be variable (default is variable)
     *
     */
    public void setKappa2Variable(boolean variable) {
	ia[6] = variable;
    }


    public TwoFibreFixedPropWatsonFitter(Vector3D[] axes, double prop) {
	super();
	alpha = prop;

	init(6);
	samples = axes;
	//	CONVERGETHRESH = 1.0E-12;
	mu1DotX = new double[samples.length];
	mu2DotX = new double[samples.length];

	eK1Mu1DotXSq = new double[samples.length];
	eK2Mu2DotXSq = new double[samples.length];
	
    }
	

	
    public static void main(String[] args) {
		
	if (args.length == 0) {
	    System.out.println("Usage: TwoFibreFixedPropWatsonFitter -sampleAxes [filename] -initMu1 [theta phi] -initMu2 [theta phi] -initK1 [value] -initK2 [value] -initAlpha [value]\n\nThe axes file should contain a list of unit vectors with x y z delimited by spaces and each vector on a separate line. Parameters default to zero if arguments are missing.");
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
		
		
	TwoFibreFixedPropWatsonFitter fitter = new TwoFibreFixedPropWatsonFitter(axes, mixingPar);
		
	try {
	    fitter.setInitParams(initParams);
	    fitter.minimise();
	}
	catch (MarquardtMinimiserException e) {
	    logger.severe("Minimization failed with exception: ");
	    LoggedException.logExceptionSevere(e, Thread.currentThread().getName());
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

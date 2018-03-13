package numerics;

import tools.FileInput;

import misc.LoggedException;
import optimizers.*;

import java.util.ArrayList;
import java.util.logging.Logger;

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
* @version $Id $
* @author Philip Cook
* 
*/
public class TwoFibreWatsonFitter extends MarquardtMinimiser {
	
    /** logging object */
    private static final Logger logger= Logger.getLogger("numerics.TwoFibreWtsonFitter");
   	   
    // sample vectors
    Vector3D[] samples;
	
    // a[1]     a[2]   a[3]    a[4]    a[5]   a[6]    a[7]
    // theta1   phi1   kappa1  theta2  phi2   kappa2  alpha
	
	
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
		
	// want 0 <= alpha <= 1.0 in final solution
	// So use sin^2(alpha) in here
	double sinSqAlpha = Math.sin(atry[7]);
	sinSqAlpha = sinSqAlpha * sinSqAlpha;
		
	double[] g = new double[samples.length];
	double[][] dgda = new double[ma+1][samples.length];

	double hyper1F1k1 = WatsonDistribution.hyper1F1(0.5, 1.5, atry[3], 1e-9);
	double hyper1F1k2 = WatsonDistribution.hyper1F1(0.5, 1.5, atry[6], 1e-9);
	
	for (int i = 0; i < samples.length; i++) {
			
	    g[i] = sinSqAlpha * (1.0 / hyper1F1k1) * Math.exp(atry[3] * Math.pow(Math.sin(atry[1])*Math.cos(atry[2])*samples[i].x + Math.sin(atry[1]) * Math.sin(atry[2])*samples[i].y + Math.cos(atry[1]) * samples[i].z, 2.0)) + (1 - sinSqAlpha) * (1.0 / hyper1F1k2) * Math.exp(atry[6] * Math.pow( Math.sin(atry[4])*Math.cos(atry[5])*samples[i].x + Math.sin(atry[4])*Math.sin(atry[5])*samples[i].y + Math.cos(atry[4])*samples[i].z, 2.0) );
	    
	    //	    g[i] = sinSqAlpha * WatsonDistribution.pdf(mu1, samples[i], atry[3]) + (1.0 - sinSqAlpha) * WatsonDistribution.pdf(mu2, samples[i], atry[6]);
	    minusSumLogPDF -= Math.log(g[i]);
	}

	double dHyperDk1 = WatsonDistribution.hyper1F1(1.5, 2.5, atry[3], 1e-9);
	double dHyperDk2 = WatsonDistribution.hyper1F1(1.5, 2.5, atry[6], 1e-9);

	//	System.out.println("kappa1 == " + atry[3]);
	//	System.out.println("kappa2 == " + atry[6]);
		
	firstDeriv(atry,g, dgda, dfda, hyper1F1k1, hyper1F1k2, dHyperDk1, dHyperDk2);
		
	secondDeriv(atry,g, dgda, d2fda2, hyper1F1k1, hyper1F1k2, dHyperDk1, dHyperDk2);
		
	//	System.out.println("fObj == " + minusSumLogPDF);
		
	return minusSumLogPDF;
    }
	
	
	
    /**
     * @param atry the parameter values
     * @param g the values of the function g = \left( \alpha W(\mu_1, \kappa_1, x_i) + (1 - \alpha)W(\mu_2, \kappa_2, x_i) \right) for all i
     * @param dgda Will be filled with the derivatives of g_i wrt a
     * @param dfda Will be filled with the derivatives of f wrt a
     */
    public void firstDeriv(double[] atry, double[] g, double[][] dgda, double[] dfda, double hyper1F1k1, double hyper1F1k2, double dHyperDk1, double dHyperDk2) {
		
	// a[1]    a[2]   a[3]    a[4]    a[5]   a[6]    a[7]
	// theta1  phi1   kappa1  theta2  phi2   kappa2  alpha
		
	// g = sin^2(\alpha) W_1(\theta_1, \phi_1, \kappa_1) + (1-sin^2(\alpha))W_2(\theta_2, \phi_2, \kappa_2)
		
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
		
	double sinSqAlpha = Math.sin(atry[7]);
	sinSqAlpha = sinSqAlpha * sinSqAlpha;


		
	// code pasted from Mathematica
	for (int i = 0; i < samples.length; i++) {
			
	    double x = samples[i].x;
	    double y = samples[i].y;
	    double z = samples[i].z;

	    double mu1DotXSq = Math.pow(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2,2);
	    double mu2DotXSq = Math.pow(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5,2);

	    double eK1Mu1DotXSq = Math.exp(atry[3]*mu1DotXSq);
	    double eK2Mu2DotXSq = Math.exp(atry[6]*mu2DotXSq);

			
	    dgda[1][i] = (2*eK1Mu1DotXSq*atry[3]*(x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*sinSqAlpha)/hyper1F1k1;
			
			
	    dfda[1] -= (1.0 / g[i]) * dgda[1][i];
			
	    dgda[2][i] = (2*eK1Mu1DotXSq*atry[3]*(y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*sinSqAlpha)/hyper1F1k1;
			
	    dfda[2] -= (1.0 / g[i]) * dgda[2][i];
			
	    dgda[3][i] = (-(1.0 / 3.0) *eK1Mu1DotXSq*dHyperDk1*sinSqAlpha)/Math.pow(hyper1F1k1,2) + (eK1Mu1DotXSq*mu1DotXSq*sinSqAlpha)/hyper1F1k1;
			
	    dfda[3] -= (1.0 / g[i]) * dgda[3][i];
			
	    dgda[4][i] = (2*eK2Mu2DotXSq*atry[6]*(x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*(1 - sinSqAlpha))/hyper1F1k2;
			
	    dfda[4] -= (1.0 / g[i]) * dgda[4][i];
			
			
	    dgda[5][i] = (2*eK2Mu2DotXSq*atry[6]*(y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*(1 - sinSqAlpha))/hyper1F1k2;
			
			
	    dfda[5] -= (1.0 / g[i]) * dgda[5][i];
			
	    dgda[6][i] = (-(1.0 / 3.0)*eK2Mu2DotXSq*dHyperDk2*(1 - sinSqAlpha))/Math.pow(hyper1F1k2,2) + (eK2Mu2DotXSq*mu2DotXSq*(1 - sinSqAlpha))/hyper1F1k2;
			
	    dfda[6] -= (1.0 / g[i]) * dgda[6][i];
			
	    dgda[7][i] = (2*eK1Mu1DotXSq*Math.cos(atry[7])*Math.sin(atry[7]))/hyper1F1k1 - (2*eK2Mu2DotXSq*Math.cos(atry[7])*Math.sin(atry[7]))/hyper1F1k2;
			
	    dfda[7] -= (1.0 / g[i]) * dgda[7][i];
			
	}
		
		
    }
	
	
    /**
     * Compute the first derivative numerically.
     */
    protected double[] firstDerivNumerical(final double[] atry) {
	double[] derivs = new double[ma + 1];
		
	final double[] a = new double[atry.length];
	for(int k=1; k<atry.length; k++) {
	    a[k] = atry[k];
	}
		
	for(int j=1; j<=ma; j++) {
			
	    final int index = j;
			
	    NumDeriv df = new NumDeriv() {
		    private double[] dfda = new double[ma+1];
		    private double[][] d2fda2 = new double[ma+1][ma+1];
		    protected float func(float arg) {
					
			a[index] = arg;
			float t = (float)fObj(a, dfda, d2fda2);
			return t; 
		    }
		};
			
	    float[] err = new float[1];
	    derivs[j] = df.dfridr((float)atry[j], 0.1f*(float)atry[j], err);
			
	    //Replace the original value in a[index] ready to compute
	    //the next derivative.
	    a[index] = atry[index];
	}
		
	return derivs;
    }
	
	
    /**
     * @param g the value of g(x) for each sample x.
     * @param d2fda2 the second derivative matrix of fObj wrt each parameter.
     * @param dgda the derivatives of g(x) wrt each parameter.
     */
    protected void secondDeriv(double[] atry, double[] g, double[][] dgda, double[][] d2fda2, double hyper1F1k1, double hyper1F1k2, double dHyperDk1, double dHyperDk2) {
	

	double sinSqAlpha = Math.sin(atry[7]);
	sinSqAlpha = sinSqAlpha * sinSqAlpha;
		
	// (d^2/dxdz)(Log g) = - (1 / g^2)(dg/dx)(dg/dz) + (d^2g/dzdx)(1/g)
	// -\sum_i (d^2/dxdz)(Log g) = -\sum_i -(1 / g_i^2)(dg/dx)(dg/dz) + (d^2g/dzdx)(1/g_i)
		
	double cosAtry1 = Math.cos(atry[1]);
	double sinAtry1 = Math.sin(atry[1]);
		
	double cosAtry2 = Math.cos(atry[2]);
	double sinAtry2 = Math.sin(atry[2]);
		
	double cosAtry4 = Math.cos(atry[4]);
	double sinAtry4 = Math.sin(atry[4]);
		
	double cosAtry5 = Math.cos(atry[5]);
	double sinAtry5 = Math.sin(atry[5]);

		
	for (int i = 0; i < samples.length; i++) {
			
	    double x = samples[i].x;
	    double y = samples[i].y;
	    double z = samples[i].z;
			
	    double invgi = 1.0 / g[i];
	    double invgiSq = invgi * invgi;
	
	    double mu1DotXSq = Math.pow(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2,2.0);
	    double mu2DotXSq = Math.pow(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5,2.0);

	    double eK1Mu1DotXSq = Math.exp(atry[3]*mu1DotXSq);
	    double eK2Mu2DotXSq = Math.exp(atry[6]*mu2DotXSq);

	    // In[48]:= CForm[D[dfda1, atry[1]]] 
	    // Out[48]//CForm=
	    d2fda2[1][1] -= invgi * (
				     (2*eK1Mu1DotXSq*atry[3]*
				      Math.pow(x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2,2)*sinSqAlpha)/hyper1F1k1 + 
				     (2*eK1Mu1DotXSq*atry[3]*
				      (-(z*cosAtry1) - x*cosAtry2*sinAtry1 - y*sinAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*
				      sinSqAlpha)/hyper1F1k1 + (4*
										   eK1Mu1DotXSq*Math.pow(atry[3],2)*
										   Math.pow(x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2,2)*
										   mu1DotXSq*sinSqAlpha)/hyper1F1k1
				     ) - invgiSq * dgda[1][i] * dgda[1][i];
			
	    // In[49]:= CForm[D[dfda1, atry[2]]]
	    // Out[49]//CForm=
	    d2fda2[1][2] -= invgi * (
				     (2*eK1Mu1DotXSq*atry[3]*
				      (x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2)*(y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2)*sinSqAlpha)/
				     hyper1F1k1 + (2*eK1Mu1DotXSq*atry[3]*
						   (y*cosAtry1*cosAtry2 - x*cosAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*sinSqAlpha)/
				     hyper1F1k1 + (4*eK1Mu1DotXSq*Math.pow(atry[3],2)*
						   (x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2)*(y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2)*
						   mu1DotXSq*sinSqAlpha)/hyper1F1k1
				     ) - invgiSq * dgda[1][i] * dgda[2][i];
			
			
	    // In[50]:=
	    // CForm[D[dyda1, atry[3]]]
	    // Out[50]//CForm=
	    d2fda2[1][3] -= invgi * ( 
				     (2*eK1Mu1DotXSq*
				      (x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*
				      sinSqAlpha)/hyper1F1k1 - (0.6666666666666666*
										   eK1Mu1DotXSq*atry[3]*dHyperDk1*
										   (x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*
										   sinSqAlpha)/Math.pow(hyper1F1k1,2) + 
				     (2*eK1Mu1DotXSq*atry[3]*
				      (x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2)*Math.pow(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2,3)*
				      sinSqAlpha)/hyper1F1k1
			
				     ) - invgiSq * dgda[1][i] * dgda[3][i];
			
			
	    // In[52]:=
	    // CForm[D[dfda1, atry[4]]]
	    // Out[52]//CForm=
	    d2fda2[1][4] -= 0.0;
			
	    // In[53]:=
	    // CForm[D[dfda1, atry[5]]]
	    // Out[53]//CForm=
	    d2fda2[1][5] -= 0.0;
			
	    // In[54]:=
	    // CForm[D[dfda1, atry[6]]]
	    // Out[54]//CForm=
	    d2fda2[1][6] -= 0.0;
			
	    // In[55]:=
	    // CForm[D[dfda1, atry[7]]]
	    // Out[55]//CForm=
	    d2fda2[1][7] -= invgi * (
				     (4*eK1Mu1DotXSq*atry[3]*Math.cos(atry[7])*
				      (x*cosAtry1*cosAtry2 - z*sinAtry1 + y*cosAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*
				      Math.sin(atry[7]))/hyper1F1k1
				     ) - invgiSq * dgda[1][i] * dgda[7][i];
			
			
	    //In[56]:=
	    // CForm[D[dfda2, atry[2]]]
	    // Out[56]//CForm=
	    d2fda2[2][2] -= invgi * (
				     (2*eK1Mu1DotXSq*atry[3]*
				      Math.pow(y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2,2)*sinSqAlpha)/hyper1F1k1 + 
				     (2*eK1Mu1DotXSq*atry[3]*
				      (-(x*cosAtry2*sinAtry1) - y*sinAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*sinSqAlpha)/
				     hyper1F1k1 + (4*eK1Mu1DotXSq*Math.pow(atry[3],2)*
						   Math.pow(y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2,2)*mu1DotXSq*
						   sinSqAlpha)/hyper1F1k1
				     ) - invgiSq * dgda[2][i] * dgda[2][i]; 
			
	    // In[57]:=
	    // CForm[D[dfda2, atry[3]]]
	    // Out[57]//CForm=
	    d2fda2[2][3] -= invgi * (
				     (2*eK1Mu1DotXSq*(y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2)*
				      (z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*sinSqAlpha)/hyper1F1k1 - 
				     ((2.0 / 3.0)*eK1Mu1DotXSq*atry[3]*dHyperDk1*
				      (y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*sinSqAlpha)/
				     Math.pow(hyper1F1k1,2) + (2*eK1Mu1DotXSq*atry[3]*
							       (y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2)*Math.pow(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2,3)*
							       sinSqAlpha)/hyper1F1k1
				     ) - invgiSq * dgda[2][i] * dgda[3][i]; 
			
	    // In[58]:=
	    // CForm[D[dfda2, atry[4]]]
	    // Out[58]//CForm=
	    d2fda2[2][4] -= 0.0;
			
	    // In[59]:=
	    // CForm[D[dfda2, atry[5]]]
	    // Out[59]//CForm=
	    d2fda2[2][5] -= 0.0;
			
	    // In[60]:=
	    // CForm[D[dfda2, atry[6]]]
	    // Out[60]//CForm=
	    d2fda2[2][6] -= 0.0;
			
	    // In[61]:=
	    // CForm[D[dfda2, atry[7]]]
	    // Out[61]//CForm=
	    d2fda2[2][7] -= invgi * ((4*eK1Mu1DotXSq*atry[3]*Math.cos(atry[7])*
				      (y*cosAtry2*sinAtry1 - x*sinAtry1*sinAtry2)*(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2)*Math.sin(atry[7]))/
				     hyper1F1k1
				     ) - invgiSq * dgda[2][i] * dgda[7][i]; 
			
	    // In[62]:=
	    // CForm[D[dfda3, atry[3]]]
	    // Out[62]//CForm=
	    d2fda2[3][3] -= invgi * (
				     ( (2.0 / 9.0)*eK1Mu1DotXSq*Math.pow(dHyperDk1,2)*
				       sinSqAlpha)/Math.pow(hyper1F1k1,3) - 
				     (0.2*eK1Mu1DotXSq*WatsonDistribution.hyper1F1(2.5,3.5,atry[3], 1e-7)*sinSqAlpha)/
				     Math.pow(hyper1F1k1,2) - ((2.0 / 3.0)*
							       eK1Mu1DotXSq*dHyperDk1*
							       mu1DotXSq*sinSqAlpha)/Math.pow(hyper1F1k1,2) + 
				     (eK1Mu1DotXSq*
				      Math.pow(z*cosAtry1 + x*cosAtry2*sinAtry1 + y*sinAtry1*sinAtry2,4)*sinSqAlpha)/hyper1F1k1
				     ) - invgiSq * dgda[3][i] * dgda[3][i]; 
			
	    // In[63]:=
	    // CForm[D[dfda3, atry[4]]]
	    // Out[63]//CForm=
	    d2fda2[3][4] -= 0.0;
			
	    // In[64]:=
	    // CForm[D[dfda3, atry[5]]]
	    // Out[64]//CForm=
	    d2fda2[3][5] -= 0.0;
			
	    // In[65]:=
	    // CForm[D[dfda3, atry[6]]]
	    // Out[65]//CForm=
	    d2fda2[3][6] -= 0.0;
			
	    // In[66]:=
	    // CForm[D[dfda3, atry[7]]]
	    // Out[66]//CForm=
	    d2fda2[3][7] -= invgi * (
				     (-(2.0 / 3.0)*eK1Mu1DotXSq*Math.cos(atry[7])*
				      dHyperDk1*Math.sin(atry[7]))/Math.pow(hyper1F1k1,2) + 
				     (2*eK1Mu1DotXSq*Math.cos(atry[7])*
				      mu1DotXSq*Math.sin(atry[7]))/hyper1F1k1
				     ) - invgiSq * dgda[3][i] * dgda[7][i]; 
			
	    // In[67]:=
	    // CForm[D[dfda4, atry[4]]]
	    // Out[67]//CForm=
	    d2fda2[4][4] -= invgi * (
				     (2*eK2Mu2DotXSq*atry[6]*
				      Math.pow(x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5,2)*(1 - sinSqAlpha))/hyper1F1k2 + 
				     (2*eK2Mu2DotXSq*atry[6]*
				      (-(z*cosAtry4) - x*cosAtry5*sinAtry4 - y*sinAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*
				      (1 - sinSqAlpha))/hyper1F1k2 + 
				     (4*eK2Mu2DotXSq*Math.pow(atry[6],2)*
				      Math.pow(x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5,2)*
				      mu2DotXSq*(1 - sinSqAlpha))/hyper1F1k2
				     ) - invgiSq * dgda[4][i] * dgda[4][i]; 
			
			
	    // In[68]:=
	    // CForm[D[dfda4, atry[5]]]
	    // Out[68]//CForm=
	    d2fda2[4][5] -= invgi * (
				     (2*eK2Mu2DotXSq*atry[6]*
				      (x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5)*(y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5)*
				      (1 - sinSqAlpha))/hyper1F1k2 + 
				     (2*eK2Mu2DotXSq*atry[6]*
				      (y*cosAtry4*cosAtry5 - x*cosAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*
				      (1 - sinSqAlpha))/hyper1F1k2 + 
				     (4*eK2Mu2DotXSq*Math.pow(atry[6],2)*
				      (x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5)*(y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5)*
				      mu2DotXSq*(1 - sinSqAlpha))/hyper1F1k2
			
				     ) - invgiSq * dgda[4][i] * dgda[5][i]; 
			
	    // In[69]:=
	    // CForm[D[dfda4, atry[6]]]
	    // Out[69]//CForm=
	    d2fda2[4][6] -= invgi * (
				     (2*eK2Mu2DotXSq*
				      (x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*
				      (1 - sinSqAlpha))/hyper1F1k2 - 
				     ((2.0 / 3.0)*eK2Mu2DotXSq*atry[6]*dHyperDk2*
				      (x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*
				      (1 - sinSqAlpha))/Math.pow(hyper1F1k2,2) + 
				     (2*eK2Mu2DotXSq*atry[6]*
				      (x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5)*Math.pow(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5,3)*
				      (1 - sinSqAlpha))/hyper1F1k2
				     ) - invgiSq * dgda[4][i] * dgda[6][i]; 
			
	    // In[70]:=
	    // CForm[D[dfda4, atry[7]]]
	    // Out[70]//CForm=
	    d2fda2[4][7] -= invgi * (
				     (-4*eK2Mu2DotXSq*atry[6]*Math.cos(atry[7])*
				      (x*cosAtry4*cosAtry5 - z*sinAtry4 + y*cosAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*
				      Math.sin(atry[7]))/hyper1F1k2
				     ) - invgiSq * dgda[4][i] * dgda[7][i]; 
			
	    // In[71]:=
	    // CForm[D[dfda5, atry[5]]]
	    // Out[71]//CForm=
	    d2fda2[5][5] -= invgi * (
				     (2*eK2Mu2DotXSq*atry[6]*
				      Math.pow(y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5,2)*(1 - sinSqAlpha))/hyper1F1k2 + 
				     (2*eK2Mu2DotXSq*atry[6]*
				      (-(x*cosAtry5*sinAtry4) - y*sinAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*
				      (1 - sinSqAlpha))/hyper1F1k2 + 
				     (4*eK2Mu2DotXSq*Math.pow(atry[6],2)*
				      Math.pow(y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5,2)*mu2DotXSq*
				      (1 - sinSqAlpha))/hyper1F1k2
				     ) - invgiSq * dgda[5][i] * dgda[5][i]; 
			
	    // In[72]:=
	    // CForm[D[dfda5, atry[6]]]
	    // Out[72]//CForm=
	    d2fda2[5][6] -= invgi * (
				     (2*eK2Mu2DotXSq*(y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5)*
				      (z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*(1 - sinSqAlpha))/hyper1F1k2 - 
				     ((2.0 / 3.0)*eK2Mu2DotXSq*atry[6]*dHyperDk2*
				      (y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*
				      (1 - sinSqAlpha))/Math.pow(hyper1F1k2,2) + 
				     (2*eK2Mu2DotXSq*atry[6]*
				      (y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5)*Math.pow(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5,3)*
				      (1 - sinSqAlpha))/hyper1F1k2
				     ) - invgiSq * dgda[5][i] * dgda[6][i]; 
			
	    // In[73]:=
	    // CForm[D[dfda5, atry[7]]]
	    // Out[73]//CForm=
			
	    d2fda2[5][7] -= invgi * (
				     (-4*eK2Mu2DotXSq*atry[6]*Math.cos(atry[7])*
				      (y*cosAtry5*sinAtry4 - x*sinAtry4*sinAtry5)*(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5)*Math.sin(atry[7]))/
				     hyper1F1k2
				     ) - invgiSq * dgda[5][i] * dgda[7][i]; 
			
	    // In[74]:=
	    // CForm[D[dfda6, atry[6]]]
	    // Out[74]//CForm=
	    d2fda2[6][6] -= invgi * ( 
				     ((2.0 / 9.0)*eK2Mu2DotXSq*Math.pow(dHyperDk2,2)*
				      (1 - sinSqAlpha))/Math.pow(hyper1F1k2,3) - 
				     (0.2*eK2Mu2DotXSq*WatsonDistribution.hyper1F1(2.5,3.5,atry[6], 1e-7)*
				      (1 - sinSqAlpha))/Math.pow(hyper1F1k2,2) - 
				     ((2.0 / 3.0)*eK2Mu2DotXSq*dHyperDk2*
				      mu2DotXSq*(1 - sinSqAlpha))/Math.pow(hyper1F1k2,2) + 
				     (eK2Mu2DotXSq*
				      Math.pow(z*cosAtry4 + x*cosAtry5*sinAtry4 + y*sinAtry4*sinAtry5,4)*(1 - sinSqAlpha))/hyper1F1k2
			
				     ) - invgiSq * dgda[6][i] * dgda[6][i]; 
			
			
	    // In[75]:=
	    // CForm[D[dfda6, atry[7]]]
	    // Out[75]//CForm=
	    d2fda2[6][7] -= invgi * (
				     ((2.0 / 3.0)*eK2Mu2DotXSq*Math.cos(atry[7])*dHyperDk2*
				      Math.sin(atry[7]))/Math.pow(hyper1F1k2,2) - (2*
										   eK2Mu2DotXSq*Math.cos(atry[7])*
										   mu2DotXSq*Math.sin(atry[7]))/hyper1F1k2
				     ) - invgiSq * dgda[6][i] * dgda[7][i]; 
			
	    // In[76]:=
	    // CForm[D[dfda7, atry[7]]]
	    // Out[76]//CForm=
	    d2fda2[7][7] -= invgi * (
				     (2*eK1Mu1DotXSq*Math.pow(Math.cos(atry[7]),2))/hyper1F1k1 - 
				     (2*eK2Mu2DotXSq*Math.pow(Math.cos(atry[7]),2))/hyper1F1k2 - 
				     (2*eK1Mu1DotXSq*sinSqAlpha)/hyper1F1k1 + 
				     (2*eK2Mu2DotXSq*sinSqAlpha)/hyper1F1k2
			
				     ) - invgiSq * dgda[7][i] * dgda[7][i]; 
			
			
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
	return Math.sin(a[7]) * Math.sin(a[7]);
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
     * from index 0: {theta1, phi1, kappa1, theta2, phi2, kappa2, alpha}. (theta, phi) are the spherical
     * polar coordinates of the mean axes of the distributions. Alpha is a mixing parameter, and should be
     * between 0.0 and 1.0.
     */
    public void setInitParams(double[] aInit) throws MarquardtMinimiserException {
	// interception
	double[] realA = new double[aInit.length];
	for (int i = 0; i < realA.length; i++) {
	    realA[i] = aInit[i];
	}
	realA[6] = Math.asin(Math.sqrt(realA[6]));
	super.setInitParams(realA);
    }
	
	

    /**
     * Make estimates of the initial parameters and attempt a minimisation with each. The fitter will retain the best fObj obtained with the estimated initial parameters.
     * <p>This method works well when the data comes from two orthogonal bipolar distributions of moderate concentration (kappa approx. 30-250). It fits a girdle distribution to the combined set of samples. The polar axis of the girdle distribution is taken as a normal to mu1 and mu2. The fitter will try sets of axes rotated at evenly spaced angles about the polar axis. The values of kappa1 and kappa2 are initialised to the negated kappa of the girdle distribution.
     *
     *
     * @param initAlpha the mixing parameter. This cannot be estimated automatically from the sample. Set to 0.5 or give a better guess if you have one.
     * @param attempts the number of sets of parameters to try.
     */
    public void fitEstimatedParams(double initAlpha, int attempts) throws MarquardtMinimiserException {

	EigenSystem3D wEig = WatsonFitter.tBarEigenSystem(samples);
	
	double[] initParams = new double[7];
	
	initParams[2] = -1.0 * WatsonFitter.fitKappa(wEig, samples);
	initParams[5] = initParams[2];
	initParams[6] = initAlpha;
	
	double bestFObjVal = Double.MAX_VALUE;
	double[] bestFittedParams = new double[7];
	boolean convergedAtLeastOnce = false;

	double[] tp = new double[2];

	Vector3D e1 = wEig.eigenvectors[0];
	Vector3D e2 = wEig.eigenvectors[1];
	Vector3D e3 = wEig.eigenvectors[2];

	
	for (int i = 0; i < attempts; i++) {
	    	    		    
	    // try to find best solution

	    e1 = Rotations.rotateVector(e1, e3, (90.0 / attempts) * i * Math.PI / 180.0);
	    e2 = Rotations.rotateVector(e2, e3, (90.0 / attempts) * i * Math.PI / 180.0);
	    
	    tp = getThetaPhi(e1);
	    
	    initParams[0] = tp[0];
	    initParams[1] = tp[1];
	    
	    tp = getThetaPhi(e2);
	    
	    initParams[3] = tp[0];
	    initParams[4] = tp[1];
	    
	    try {
		setInitParams(initParams);
		minimise();
		convergedAtLeastOnce = true;
		
	    }
	    catch (MarquardtMinimiserException e) {
		//				System.out.println(e);
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
    

    public double[] getParameters() {
		
	double[] param = new double[ma];
		
	for (int i = 0;  i < param.length; i++) {
	    param[i] = a[i+1];
	}
		
	param[6] = Math.pow(Math.sin(param[6]),2.0);
		
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
     * @param variable true if alpha is to be variable (default is variable)
     *
     */
    public void setMixingParameterVariable(boolean variable) {
	ia[7] = variable;
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


    public TwoFibreWatsonFitter(Vector3D[] axes) {
	super();
	//	MarquardtMinimiser.CONVERGETHRESH = 1.0E-12;
	init(7);
	samples = axes;
    }
	

    private static double[] getThetaPhi(Vector3D v) {

	double[] tp = new double[2];

	tp[0] = Math.acos( v.z );

	// phi goes from 0.0 (+x axis) and wraps at 2 * PI
	// theta goes from 0.0 (+z axis) and wraps at PI

	// if x and y are 0.0
	if (v.x == 0.0 && v.y == 0.0) {
	    tp[1] = 0.0; 
	}
	else {
	    
	    // ie, if ( x == 0 && y == 0 ) == false

	    if (v.y == 0.0) {

		if (v.x > 0.0) {
		    tp[1] = 0.0;
		}
		else {
		    tp[1] = Math.PI;
		}
		
		
	    }
	    else if (v.x == 0.0) {
		
		// avoid div by zero
		if (v.y > 0) {
		    tp[1] = Math.PI / 2.0;
		}
		else {
		    tp[1] = 1.5 * Math.PI;
		}
	    }
	    else if (v.x > 0.0 && v.y > 0.0) { // first quadrant
		tp[1] = Math.atan(v.y / v.x);
	    }
	    else if (v.x < 0.0 && v.y > 0.0) { // second quadrant
		tp[1] = Math.PI + Math.atan(v.y / v.x);
	    }
	    else if (v.x < 0.0 && v.y < 0.0) { // third quadrant
		tp[1] =  Math.PI + Math.atan(v.y / v.x);
	    }
	    else { // fourth quadrant
		tp[1] = 2.0 * Math.PI + Math.atan(v.y / v.x); 
	    }

	}
	
	return tp;
    }

	
    public static void main(String[] args) {
		
	if (args.length == 0) {
	    System.out.println("TwoFibreWatsonFitter -sampleAxes [filename] -initMu1 [theta phi] -initMu2 [theta phi] -initK1 [value] -initK2 [value] -initAlpha [value]\n\nThe axes file should contain a list of unit vectors with x y z delimited by spaces and each vector on a separate line. Parameters default to zero if arguments are missing.");
	    System.exit(0);
	}
		
	Vector3D[] axes = null;
		
	double[] initParams = new double[7];
		
	for (int i = 0; i < args.length; i++) {
			
	    if (args[i].toLowerCase().equals("-sampleaxes")) {
				
		FileInput in = new FileInput(args[i+1]);
				
		ArrayList<Vector3D> a = new ArrayList<Vector3D>();
		//J1.4 ArrayList a = new ArrayList();
				
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
	    else if (args[i].toLowerCase().equals("-initalpha")) {
		initParams[6] = Double.parseDouble(args[i+1]);
	    }
			
			
			
	}
		
		
	TwoFibreWatsonFitter fitter = new TwoFibreWatsonFitter(axes);
		
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

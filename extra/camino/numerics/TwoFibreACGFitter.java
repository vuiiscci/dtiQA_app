package numerics;

import optimizers.*;

import java.util.Random;

/**
 * <dl>
 *
 * <dt>Purpose: 
 *
 * <dd> Fits two ACG distributions to a collection of axes.
 *
 * <dt>Description:
 *
 * <dd> This class takes a collection of axes representing a set of estimates of two fibre directions. 
 * It uses a Conjugate Gradient optimization routine to find the pair of ACG distributions that
 * best fit the axes.
 * </dl>
 *
 * @version $Id$
 * @author Philip Cook
 * 
 */
public class TwoFibreACGFitter extends ConjGradMinimizer {

    protected Vector3D[] samples;
    
    protected double[] gi;
    
    public TwoFibreACGFitter(Vector3D[] vecs) {
	samples = vecs;
	init(12);
    }


    /**
     *
     * Gets the parameters of the matrices A and B given the Cholesky parameters. 
     * Note that this returns array with data starting from index 0, while <code>getCholParams</code>
     * expects and returns an array with the data stating at index 1.
     *
     */
   protected static double[] getParams_Chol(double[] a) {
       
       // {a11, a12, a13, a22, a23, a33, b11, b12, b13, b22, b23, b33}
       double[] params = new double[12];
       
       params[0] = a[1]*a[1];
       params[1] = a[1]*a[2]; 
       params[2] = a[1]*a[3];
       params[3] = a[2]*a[2] + a[4]*a[4];
       params[4] = a[2]*a[3] + a[4]*a[5];
       params[5] = a[3]*a[3] + a[5]*a[5] + a[6]*a[6];

       params[6] = a[7]*a[7];
       params[7] = a[7]*a[8]; 
       params[8] = a[7]*a[9];
       params[9] = a[8]*a[8] + a[10]*a[10];
       params[10] = a[8]*a[9] + a[10]*a[11];
       params[11] = a[9]*a[9] + a[11]*a[11] + a[12]*a[12];

    
       return params;
    }


    /**
     * Returns the elements of the Cholesky decomposition U of the
     * matrices A and B. Feed these parameters to the #getParams_Chol(double[]) to recover A and B.
     *
     */
    public static double[] getCholParams(double[] a) {
	
	double a11 = a[1];
	double a12 = a[2];
	double a13 = a[3];
	double a22 = a[4];
	double a23 = a[5];
	double a33 = a[6];

	double b11 = a[7];
	double b12 = a[8];
	double b13 = a[9];
	double b22 = a[10];
	double b23 = a[11];
	double b33 = a[12];

	double[] u = new double[13];

	u[1] = (a11>0.0)?Math.sqrt(a11):0.0;
	u[2] = (u[1]>0)?a12/u[1]:0.0;
	u[3] = (u[1]>0)?a13/u[1]:0.0;
	double eq2 = a22 - u[2]*u[2];
	u[4] = (eq2>0.0)?Math.sqrt(eq2):0.0;
	u[5] = (u[4]>0)?(a23 - u[2]*u[3])/u[4]:0.0;
	double eq3 = a33 - u[3]*u[3] - u[5]*u[5];
	u[6] = (eq3>0.0)?Math.sqrt(eq3):0.0;

	u[7] = (b11>0.0)?Math.sqrt(b11):0.0;
	u[8] = (u[7]>0)?b12/u[7]:0.0;
	u[9] = (u[7]>0)?b13/u[7]:0.0;
	eq2 = b22 - u[8]*u[8];
	u[10] = (eq2>0.0)?Math.sqrt(eq2):0.0;
	u[11] = (u[10]>0)?(b23 - u[8]*u[9])/u[10]:0.0;
	eq3 = b33 - u[9]*u[9] - u[11]*u[11];
	u[12] = (eq3>0.0)?Math.sqrt(eq3):0.0;

	return u;
    }

    /**
     * Runs the minimization and returns the answer.
     *
     * @param dt1E1Vec an estimate of the first principal direction
     * @param dt2E1Vec an estimate of the second principal direction
     * 
     */
    public double[] fitEstimatedParams(Vector3D dt1E1Vec, Vector3D dt2E1Vec, double ftol) {


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
	
	    if (Math.abs(dt1E1Vec.dot(sample)) > Math.abs(dt2E1Vec.dot(sample))) {
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
	
	RealMatrix a = ACG_Fitter.findA(dt1Samples);
	RealMatrix b = ACG_Fitter.findA(dt2Samples);

	double[] aParams = new double[6];
	
	aParams[0] = a.entries[0][0];
	aParams[1] = a.entries[0][1];
	aParams[2] = a.entries[0][2];
	aParams[3] = a.entries[1][1];
	aParams[4] = a.entries[1][2];		
	aParams[5] = a.entries[2][2];
	
	double[] bParams = new double[6];
	
	bParams[0] = b.entries[0][0];
	bParams[1] = b.entries[0][1];
	bParams[2] = b.entries[0][2];
	bParams[3] = b.entries[1][1];
	bParams[4] = b.entries[1][2];		
	bParams[5] = b.entries[2][2];
        
        try {
            return fitEstimatedParams(aParams, bParams, ftol);
        }
        catch (ConjGradMinimizerException e) {

            double[] ab = new double[12];
            
            for (int i = 0; i < 6; i++) {
                ab[i] = aParams[i];
                ab[i+6] = bParams[i];
            }

            return ab;
        }
    }


    /**
     * Runs the minimization and returns the answer.
     *
     * @param a the elements of matrix a: {axx, axy, axz, ayy, ayz, axx}
     * @param b the elements of matrix b: {bxx, bxy, bxz, byy, byz, bxx}
     * 
     */
    public double[] fitEstimatedParams(double[] a, double[] b, double ftol) throws ConjGradMinimizerException {

	// write method to: take initial guess, do minimization, get output and return it
	double[] ab = new double[13];

	for (int i = 0; i < 6; i++) {
	    ab[i+1] = a[i];
	    ab[i+7] = b[i];
	}

	double[] atry = getCholParams(ab);

	minimise(atry, ftol);

	double[] params = getParams_Chol(atry);

	return params;
	
    }


    /**
     * Converts output of <code>fitEstimatedParams</code>.
     */
    public static ACG_Distribution[] getDistributions(double[] params, Random r) {

	RealMatrix A = new RealMatrix(3,3);
	RealMatrix B = new RealMatrix(3,3);

	A.entries[0][0] = params[0];
	A.entries[0][1] = params[1];
	A.entries[0][2] = params[2];

	A.entries[1][0] = params[1];
	A.entries[1][1] = params[3];
	A.entries[1][2] = params[4];

	A.entries[2][0] = params[2];
	A.entries[2][1] = params[4];
	A.entries[2][2] = params[5];
	
	
	B.entries[0][0] = params[6];
	B.entries[0][1] = params[7];
	B.entries[0][2] = params[8];

	B.entries[1][0] = params[7];
	B.entries[1][1] = params[9];
	B.entries[1][2] = params[10];

	B.entries[2][0] = params[8];
	B.entries[2][1] = params[10];
	B.entries[2][2] = params[11];

	
	ACG_Distribution[] d = new ACG_Distribution[2];

	d[0] = new ACG_Distribution(A, r);	
	d[1] = new ACG_Distribution(B, r);

	return d;

    }

    


    protected void init(int noParams) {		
	super.init(noParams);
	gi = new double[samples.length];
    }
    

    protected double fObj(double[] atry) {
	// set gi

	double f = 0.0;

	final double ua11 = atry[1];
	final double ua12 = atry[2];
	final double ua13 = atry[3];
	final double ua22 = atry[4];
	final double ua23 = atry[5];
	final double ua33 = atry[6];

	final double ub11 = atry[7];
	final double ub12 = atry[8];
	final double ub13 = atry[9];
	final double ub22 = atry[10];
	final double ub23 = atry[11];
	final double ub33 = atry[12];
	

	for (int i = 0; i < samples.length; i++) {

	    gi[i] = 0.5/(Math.sqrt(ua11*ua11*ua22*ua22*ua33*ua33)*
        Math.pow((ua13*ua13*ua22*ua22*samples[i].x * samples[i].x + ua22*ua22*ua33*ua33*samples[i].x * samples[i].x +
            ua12*ua12*(ua23*ua23 + ua33*ua33)*samples[i].x * samples[i].x +
            ua11*ua11*ua23*ua23*samples[i].y * samples[i].y + ua11*ua11*ua33*ua33*samples[i].y * samples[i].y -
            2*ua11*ua11*ua22*ua23*samples[i].y*samples[i].z + ua11*ua11*ua22*ua22*samples[i].z * samples[i].z -
            2*ua13*ua22*samples[i].x*(ua12*ua23*samples[i].x - ua11*ua23*samples[i].y + ua11*ua22*samples[i].z) -
            2*ua11*ua12*samples[i].x*(ua23*ua23*samples[i].y + ua33*ua33*samples[i].y - ua22*ua23*samples[i].z))/
          (ua11*ua11*ua22*ua22*ua33*ua33),1.5)) +
     0.5/(Math.sqrt(ub11*ub11*ub22*ub22*ub33*ub33)*
        Math.pow((ub13*ub13*ub22*ub22*samples[i].x * samples[i].x + ub22*ub22*ub33*ub33*samples[i].x * samples[i].x +
            ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x +
            ub11*ub11*ub23*ub23*samples[i].y * samples[i].y + ub11*ub11*ub33*ub33*samples[i].y * samples[i].y -
            2*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z + ub11*ub11*ub22*ub22*samples[i].z * samples[i].z -
            2*ub13*ub22*samples[i].x*(ub12*ub23*samples[i].x - ub11*ub23*samples[i].y + ub11*ub22*samples[i].z) -
            2*ub11*ub12*samples[i].x*(ub23*ub23*samples[i].y + ub33*ub33*samples[i].y - ub22*ub23*samples[i].z))/
	      (ub11*ub11*ub22*ub22*ub33*ub33),1.5));


// 	    if (gi[i] < 0.0 || Double.isNaN(gi[i]) || Double.isInfinite(gi[i])) {
// 		System.out.println("(after ub) gi[" + i + "] = " + gi[i]);
// 		System.out.println(samples[i]);
// 		System.exit(0);
// 	    }
	    

	    f -= Math.log(gi[i]);

	}
    
	return f;

    }


    /**
     * @param atry should be {0.0, ua11, ua12, ua13, ua22, ua23, ua33, ub11, ub12, ub13, ub22, ub23, ub33}.
     */
    protected double[] dfObj(double[] atry) {
	
	double[] dfda = new double[13];

	double ua11 = atry[1];
	double ua12 = atry[2];
	double ua13 = atry[3];
	double ua22 = atry[4];
	double ua23 = atry[5];
	double ua33 = atry[6];

	double ub11 = atry[7];
	double ub12 = atry[8];
	double ub13 = atry[9];
	double ub22 = atry[10];
	double ub23 = atry[11];
	double ub33 = atry[12];


	// atry are Cholesky params of matrices A and B


	for (int i = 0; i < samples.length; i++) {

	    double alpha = ua13*ua13*ua22*ua22*samples[i].x * samples[i].x + ua22*ua22*ua33*ua33*samples[i].x * samples[i].x;
	    
	    double beta = ua12*ua12*(ua23*ua23 + ua33*ua33)*samples[i].x * samples[i].x;
	    
	    
	    double gamma = Math.pow((alpha + beta + ua11*ua11*ua23*ua23*samples[i].y * samples[i].y + ua11*ua11*ua33*ua33*samples[i].y * samples[i].y - 2*ua11*ua11*ua22*ua23*samples[i].y*samples[i].z + ua11*ua11*ua22*ua22*samples[i].z * samples[i].z - 2*ua13*ua22*samples[i].x*(ua12*ua23*samples[i].x - ua11*ua23*samples[i].y + ua11*ua22*samples[i].z) - 2*ua11*ua12*samples[i].x*(ua23*ua23*samples[i].y + ua33*ua33*samples[i].y - ua22*ua23*samples[i].z))/ (ua11*ua11*ua22*ua22*ua33*ua33),2.5);
	    
	    double delta = Math.pow(alpha + beta + ua11*ua11*ua23*ua23*samples[i].y * samples[i].y + ua11*ua11*ua33*ua33*samples[i].y * samples[i].y - 2.*ua11*ua11*ua22*ua23*samples[i].y*samples[i].z + ua11*ua11*ua22*ua22*samples[i].z * samples[i].z + ua13*ua22*samples[i].x*(-2.*ua12*ua23*samples[i].x + 2.*ua11*ua23*samples[i].y - 2.*ua11*ua22*samples[i].z) + ua11*ua12*samples[i].x*(-2.*ua23*ua23*samples[i].y - 2.*ua33*ua33*samples[i].y + 2.*ua22*ua23*samples[i].z),2);
    
	    double epsilon = Math.pow(ub13*ub13*ub22*ub22*samples[i].x * samples[i].x + ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x + ub11*ub11*ub23*ub23*samples[i].y * samples[i].y + ub11*ub11*ub33*ub33*samples[i].y * samples[i].y - 2.*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z + ub11*ub11*ub22*ub22*samples[i].z * samples[i].z + ub13*ub22*samples[i].x*(-2.*ub12*ub23*samples[i].x + 2.*ub11*ub23*samples[i].y - 2.*ub11*ub22*samples[i].z) + ub11*ub12*samples[i].x*(-2.*ub23*ub23*samples[i].y - 2.*ub33*ub33*samples[i].y + 2.*ub22*ub23*samples[i].z),2);

	    double[] dgda = new double[13];

	    //	    In[22]:=
	    //dfdaxx=CForm[D[f,ua11]]
	    // Out[22]//CForm=
	    dgda[1] = (Math.sqrt(ua11*ua11*ua22*ua22*ua33*ua33)* (alpha + beta - 0.5*ua11*ua11*ua23*ua23*samples[i].y * samples[i].y - 0.5*ua11*ua11*ua33*ua33*samples[i].y * samples[i].y + ua11*ua11*ua22*ua23*samples[i].y*samples[i].z - 0.5*ua11*ua11*ua22*ua22*samples[i].z * samples[i].z + ua13*ua22*samples[i].x*(-2.*ua12*ua23*samples[i].x + 0.5*ua11*ua23*samples[i].y - 0.5*ua11*ua22*samples[i].z) + ua11*ua12*samples[i].x*(-0.5*ua23*ua23*samples[i].y - 0.5*ua33*ua33*samples[i].y + 0.5*ua22*ua23*samples[i].z)))/ (ua11*Math.sqrt((alpha + beta + ua11*ua11*ua23*ua23*samples[i].y * samples[i].y + ua11*ua11*ua33*ua33*samples[i].y * samples[i].y - 2*ua11*ua11*ua22*ua23*samples[i].y*samples[i].z + ua11*ua11*ua22*ua22*samples[i].z * samples[i].z - 2*ua13*ua22*samples[i].x*(ua12*ua23*samples[i].x - ua11*ua23*samples[i].y + ua11*ua22*samples[i].z) - 2*ua11*ua12*samples[i].x*(ua23*ua23*samples[i].y + ua33*ua33*samples[i].y - ua22*ua23*samples[i].z))/ (ua11*ua11*ua22*ua22*ua33*ua33))* delta);
				
	    dfda[1] -= dgda[1] / gi[i];
				
	    // In[23]:=
	    
	    // dfdaxy = CForm[D[f,ua12]]
	    // Out[23]//CForm=
	    dgda[2] = (-1.5*samples[i].x*(-(ua13*ua22*ua23*samples[i].x) + ua12*(ua23*ua23 + ua33*ua33)*samples[i].x - ua11*(ua23*ua23*samples[i].y + ua33*ua33*samples[i].y - ua22*ua23*samples[i].z)))/ (Math.pow(ua11*ua11*ua22*ua22*ua33*ua33,1.5)* gamma);
	    
	    dfda[2] -= dgda[2] / gi[i];
	    
	    
	    // In[24]:=

	    // dfdaxz = CForm[D[f,ua13]]
	    // Out[24]//CForm=
	    
	    dgda[3] = (-1.5*ua22*samples[i].x*(ua13*ua22*samples[i].x - ua12*ua23*samples[i].x + ua11*ua23*samples[i].y - ua11*ua22*samples[i].z))/ (Math.pow(ua11*ua11*ua22*ua22*ua33*ua33,1.5)* gamma);

		dfda[3] -= dgda[3] / gi[i];

	    // In[25]:=
	    
	    // dfdayy = CForm[D[f,ua22]]
	    // Out[25]//CForm=
	    
	    dgda[4] = (Math.sqrt(ua11*ua11*ua22*ua22*ua33*ua33)* (-0.5*ua13*ua13*ua22*ua22*samples[i].x * samples[i].x - 0.5*ua22*ua22*ua33*ua33*samples[i].x * samples[i].x + beta + ua11*ua11*ua23*ua23*samples[i].y * samples[i].y + ua11*ua11*ua33*ua33*samples[i].y * samples[i].y - 0.5*ua11*ua11*ua22*ua23*samples[i].y*samples[i].z - 0.5*ua11*ua11*ua22*ua22*samples[i].z * samples[i].z + ua13*ua22*samples[i].x*(-0.5*ua12*ua23*samples[i].x + 0.5*ua11*ua23*samples[i].y + ua11*ua22*samples[i].z) + ua11*ua12*samples[i].x*(-2.*ua23*ua23*samples[i].y - 2.*ua33*ua33*samples[i].y + 0.5*ua22*ua23*samples[i].z)))/ (ua22*Math.sqrt((alpha + beta + ua11*ua11*ua23*ua23*samples[i].y * samples[i].y + ua11*ua11*ua33*ua33*samples[i].y * samples[i].y - 2*ua11*ua11*ua22*ua23*samples[i].y*samples[i].z + ua11*ua11*ua22*ua22*samples[i].z * samples[i].z - 2*ua13*ua22*samples[i].x*(ua12*ua23*samples[i].x - ua11*ua23*samples[i].y + ua11*ua22*samples[i].z) - 2*ua11*ua12*samples[i].x*(ua23*ua23*samples[i].y + ua33*ua33*samples[i].y - ua22*ua23*samples[i].z))/ (ua11*ua11*ua22*ua22*ua33*ua33))* delta);


	    dfda[4] -= dgda[4] / gi[i];

	    // In[26]:=
	    
	    // dfdayz = CForm[D[f, ua23]]
	    // Out[26]//CForm=
	    dgda[5] = (-1.5*(ua12*samples[i].x - ua11*samples[i].y)*(-(ua13*ua22*samples[i].x) + ua12*ua23*samples[i].x - ua11*ua23*samples[i].y + ua11*ua22*samples[i].z))/ (Math.pow(ua11*ua11*ua22*ua22*ua33*ua33,1.5)* gamma);


	    dfda[5] -= dgda[5] / gi[i];

	    
// 	    In[27]:=
		
// 		dfdazz = CForm[D[f,ua33]]
// 		Out[27]//CForm=

	    dgda[6] = (Math.sqrt(ua11*ua11*ua22*ua22*ua33*ua33)* (ua13*ua13*ua22*ua22*samples[i].x * samples[i].x - 0.5*ua22*ua22*ua33*ua33*samples[i].x * samples[i].x + ua12*ua12*(ua23*ua23 - 0.5*ua33*ua33)*samples[i].x * samples[i].x + ua11*ua11*ua23*ua23*samples[i].y * samples[i].y - 0.5*ua11*ua11*ua33*ua33*samples[i].y * samples[i].y - 2.*ua11*ua11*ua22*ua23*samples[i].y*samples[i].z + ua11*ua11*ua22*ua22*samples[i].z * samples[i].z + ua13*ua22*samples[i].x*(-2.*ua12*ua23*samples[i].x + 2.*ua11*ua23*samples[i].y - 2.*ua11*ua22*samples[i].z) + ua11*ua12*samples[i].x*(-2.*ua23*ua23*samples[i].y + ua33*ua33*samples[i].y + 2.*ua22*ua23*samples[i].z)))/ (ua33*Math.sqrt((alpha + beta + ua11*ua11*ua23*ua23*samples[i].y * samples[i].y + ua11*ua11*ua33*ua33*samples[i].y * samples[i].y - 2*ua11*ua11*ua22*ua23*samples[i].y*samples[i].z + ua11*ua11*ua22*ua22*samples[i].z * samples[i].z - 2*ua13*ua22*samples[i].x*(ua12*ua23*samples[i].x - ua11*ua23*samples[i].y + ua11*ua22*samples[i].z) - 2*ua11*ua12*samples[i].x*(ua23*ua23*samples[i].y + ua33*ua33*samples[i].y - ua22*ua23*samples[i].z))/ (ua11*ua11*ua22*ua22*ua33*ua33))* delta);

	    dfda[6] -= dgda[6] / gi[i];

	    // In[28]:=

	    // dfdaxx = CForm[D[f,ub11]]
	    // Out[28]//CForm=
	    dgda[7] = (Math.sqrt(ub11*ub11*ub22*ub22*ub33*ub33)* (ub13*ub13*ub22*ub22*samples[i].x * samples[i].x + ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x - 0.5*ub11*ub11*ub23*ub23*samples[i].y * samples[i].y - 0.5*ub11*ub11*ub33*ub33*samples[i].y * samples[i].y + ub11*ub11*ub22*ub23*samples[i].y*samples[i].z - 0.5*ub11*ub11*ub22*ub22*samples[i].z * samples[i].z + ub13*ub22*samples[i].x*(-2.*ub12*ub23*samples[i].x + 0.5*ub11*ub23*samples[i].y - 0.5*ub11*ub22*samples[i].z) + ub11*ub12*samples[i].x*(-0.5*ub23*ub23*samples[i].y - 0.5*ub33*ub33*samples[i].y + 0.5*ub22*ub23*samples[i].z)))/ (ub11*Math.sqrt((ub13*ub13*ub22*ub22*samples[i].x * samples[i].x + ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x + ub11*ub11*ub23*ub23*samples[i].y * samples[i].y + ub11*ub11*ub33*ub33*samples[i].y * samples[i].y - 2*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z + ub11*ub11*ub22*ub22*samples[i].z * samples[i].z - 2*ub13*ub22*samples[i].x*(ub12*ub23*samples[i].x - ub11*ub23*samples[i].y + ub11*ub22*samples[i].z) - 2*ub11*ub12*samples[i].x*(ub23*ub23*samples[i].y + ub33*ub33*samples[i].y - ub22*ub23*samples[i].z))/ (ub11*ub11*ub22*ub22*ub33*ub33))* epsilon);

	    dfda[7] -= dgda[7] / gi[i];

	    // In[29]:=
	    
	    // dfdaxy = CForm[D[f,ub12]]
	    // Out[29]//CForm=

	    dgda[8] = (-1.5*samples[i].x*(-(ub13*ub22*ub23*samples[i].x) + ub12*(ub23*ub23 + ub33*ub33)*samples[i].x - ub11*(ub23*ub23*samples[i].y + ub33*ub33*samples[i].y - ub22*ub23*samples[i].z)))/ (Math.pow(ub11*ub11*ub22*ub22*ub33*ub33,1.5)* Math.pow((ub13*ub13*ub22*ub22*samples[i].x * samples[i].x + ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x + ub11*ub11*ub23*ub23*samples[i].y * samples[i].y + ub11*ub11*ub33*ub33*samples[i].y * samples[i].y - 2*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z + ub11*ub11*ub22*ub22*samples[i].z * samples[i].z - 2*ub13*ub22*samples[i].x*(ub12*ub23*samples[i].x - ub11*ub23*samples[i].y + ub11*ub22*samples[i].z) - 2*ub11*ub12*samples[i].x*(ub23*ub23*samples[i].y + ub33*ub33*samples[i].y - ub22*ub23*samples[i].z))/ (ub11*ub11*ub22*ub22*ub33*ub33),2.5));

	    dfda[8] -= dgda[8] / gi[i];


// In[30]:=

// dfdaxz = CForm[D[f,ub13]]
// Out[30]//CForm=
	    dgda[9] = (-1.5*ub22*samples[i].x*(ub13*ub22*samples[i].x - ub12*ub23*samples[i].x + ub11*ub23*samples[i].y - ub11*ub22*samples[i].z))/ (Math.pow(ub11*ub11*ub22*ub22*ub33*ub33,1.5)* Math.pow((ub13*ub13*ub22*ub22*samples[i].x * samples[i].x + ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x + ub11*ub11*ub23*ub23*samples[i].y * samples[i].y + ub11*ub11*ub33*ub33*samples[i].y * samples[i].y - 2*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z + ub11*ub11*ub22*ub22*samples[i].z * samples[i].z - 2*ub13*ub22*samples[i].x*(ub12*ub23*samples[i].x - ub11*ub23*samples[i].y + ub11*ub22*samples[i].z) - 2*ub11*ub12*samples[i].x*(ub23*ub23*samples[i].y + ub33*ub33*samples[i].y - ub22*ub23*samples[i].z))/ (ub11*ub11*ub22*ub22*ub33*ub33),2.5));

	    dfda[9] -= dgda[9] / gi[i];


	    // In[31]:=
	    
	    // dfdayy = CForm[D[f,ub22]]
	    // Out[31]//CForm=
	    dgda[10] = (Math.sqrt(ub11*ub11*ub22*ub22*ub33*ub33)* (-0.5*ub13*ub13*ub22*ub22*samples[i].x * samples[i].x - 0.5*ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x + ub11*ub11*ub23*ub23*samples[i].y * samples[i].y + ub11*ub11*ub33*ub33*samples[i].y * samples[i].y - 0.5*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z - 0.5*ub11*ub11*ub22*ub22*samples[i].z * samples[i].z + ub13*ub22*samples[i].x*(-0.5*ub12*ub23*samples[i].x + 0.5*ub11*ub23*samples[i].y + ub11*ub22*samples[i].z) + ub11*ub12*samples[i].x*(-2.*ub23*ub23*samples[i].y - 2.*ub33*ub33*samples[i].y + 0.5*ub22*ub23*samples[i].z)))/ (ub22*Math.sqrt((ub13*ub13*ub22*ub22*samples[i].x * samples[i].x + ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x + ub11*ub11*ub23*ub23*samples[i].y * samples[i].y + ub11*ub11*ub33*ub33*samples[i].y * samples[i].y - 2*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z + ub11*ub11*ub22*ub22*samples[i].z * samples[i].z - 2*ub13*ub22*samples[i].x*(ub12*ub23*samples[i].x - ub11*ub23*samples[i].y + ub11*ub22*samples[i].z) - 2*ub11*ub12*samples[i].x*(ub23*ub23*samples[i].y + ub33*ub33*samples[i].y - ub22*ub23*samples[i].z))/ (ub11*ub11*ub22*ub22*ub33*ub33))* epsilon);

	    // In[32]:=
		dfda[10] -= dgda[10] / gi[i];

	    // dfdayz = CForm[D[f, ub23]]
	    // Out[32]//CForm=
	    dgda[11] = (-1.5*(ub12*samples[i].x - ub11*samples[i].y)*(-(ub13*ub22*samples[i].x) + ub12*ub23*samples[i].x - ub11*ub23*samples[i].y + ub11*ub22*samples[i].z))/ (Math.pow(ub11*ub11*ub22*ub22*ub33*ub33,1.5)* Math.pow((ub13*ub13*ub22*ub22*samples[i].x * samples[i].x + ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x + ub11*ub11*ub23*ub23*samples[i].y * samples[i].y + ub11*ub11*ub33*ub33*samples[i].y * samples[i].y - 2*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z + ub11*ub11*ub22*ub22*samples[i].z * samples[i].z - 2*ub13*ub22*samples[i].x*(ub12*ub23*samples[i].x - ub11*ub23*samples[i].y + ub11*ub22*samples[i].z) - 2*ub11*ub12*samples[i].x*(ub23*ub23*samples[i].y + ub33*ub33*samples[i].y - ub22*ub23*samples[i].z))/ (ub11*ub11*ub22*ub22*ub33*ub33),2.5));


	    // In[33]:=
	      dfda[11] -= dgda[11] / gi[i];

	    // dfdazz = CForm[D[f,ub33]]
	    // Out[33]//CForm=

	    dgda[12] = (Math.sqrt(ub11*ub11*ub22*ub22*ub33*ub33)* (ub13*ub13*ub22*ub22*samples[i].x * samples[i].x - 0.5*ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 - 0.5*ub33*ub33)*samples[i].x * samples[i].x + ub11*ub11*ub23*ub23*samples[i].y * samples[i].y - 0.5*ub11*ub11*ub33*ub33*samples[i].y * samples[i].y - 2.*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z + ub11*ub11*ub22*ub22*samples[i].z * samples[i].z + ub13*ub22*samples[i].x*(-2.*ub12*ub23*samples[i].x + 2.*ub11*ub23*samples[i].y - 2.*ub11*ub22*samples[i].z) + ub11*ub12*samples[i].x*(-2.*ub23*ub23*samples[i].y + ub33*ub33*samples[i].y + 2.*ub22*ub23*samples[i].z)))/ (ub33*Math.sqrt((ub13*ub13*ub22*ub22*samples[i].x * samples[i].x + ub22*ub22*ub33*ub33*samples[i].x * samples[i].x + ub12*ub12*(ub23*ub23 + ub33*ub33)*samples[i].x * samples[i].x + ub11*ub11*ub23*ub23*samples[i].y * samples[i].y + ub11*ub11*ub33*ub33*samples[i].y * samples[i].y - 2*ub11*ub11*ub22*ub23*samples[i].y*samples[i].z + ub11*ub11*ub22*ub22*samples[i].z * samples[i].z - 2*ub13*ub22*samples[i].x*(ub12*ub23*samples[i].x - ub11*ub23*samples[i].y + ub11*ub22*samples[i].z) - 2*ub11*ub12*samples[i].x*(ub23*ub23*samples[i].y + ub33*ub33*samples[i].y - ub22*ub23*samples[i].z))/ (ub11*ub11*ub22*ub22*ub33*ub33))* epsilon);

		dfda[12] -= dgda[12] / gi[i];


	}

	return dfda;
	
    }


}

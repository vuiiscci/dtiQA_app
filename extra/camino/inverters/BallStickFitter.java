package inverters;

import optimizers.*;
import numerics.*;
import misc.*;
import imaging.*;

/**
 * Fits the Behrens ball and stick model to diffusion-weighted data.
 * <BR>
 * The model is S(g, b) = S_0 (f \exp[-b d (-g^T v)^2] + [1-f] \exp[-b d]), where S(g,b) is the
 * DW signal along gradient direction g with b-value b, d is a diffusion coefficient, v is the
 * orientation of anisotropic diffusion, and f is a mixing parameter (0 <= f <= 1).
 *
 * @see data.BallStick
 *
 * @author  Philip Cook
 * @version $Id
 * 
 */
public class BallStickFitter extends MarquardtChiSqFitter {

    // param  name
    // a[1]   sqrt(S_0) [unweighted signal]
    // a[2]   sqrt(d) [diffusion coefficient d]
    // a[3]   f [mixing fraction is cos^2(f)]
    // a[4]   theta
    // a[5]   phi

    // following Behrens, we start with a tensor D = R^T A R
    // where A = diag[d,0,0], and R is a composite rotation of theta radians about the Y-axis
    // followed by phi radians about the Z-axis. We find the theta and phi that rotate the
    // X-axis to the orientation of the stick.

    private DW_Scheme ip; 

    private LinearDT_Inversion inv;

  
    /**
     * The b-values, indexed from 1.
     */
    protected double[] bValues;
    

    /**
     * The number of parameters in the model
     */
    protected final int noParams = 5; 

    
    


    /**
     *
     */
    public BallStickFitter(DW_Scheme ip) throws MarquardtMinimiserException {
	
	this.ip = ip;
	
	// this value in ndata after construction
	int numMeas = ip.numMeasurements();

	double[][] g = new double[numMeas][];

	double[] b = new double[numMeas];
	
	for (int i = 0; i < numMeas; i++) {
	    g[i] = ip.getG_Dir(i);
	    b[i] = ip.getB_Value(i);
	}
	
	inv = new LinearDT_Inversion(ip);

        initialize(g, b);

    }



    /**
     * Initializes the instance variables and data structures. Call newDepVals to set
     * the data.
     * 
     * @param indepVals
     *            The matrix of wavenumbers q.
     * 
     * @param bVals
     *            The array of b-values.
     * 
     *
     */
    protected void initialize(double[][] indepVals, double[] bVals) 
	throws MarquardtMinimiserException {

	// index diffusion times from 1 like y and q
        bValues = new double[bVals.length + 1];

	for (int i = 0; i < bVals.length; i++) {
	    bValues[i+1] = bVals[i];
	}

        // dummy dependent values (ie data)
	double[] depVals = new double[bVals.length];

        // Pass sampled points to initialiser, as well as the
        // number of parameters of the model.
        initData(indepVals, depVals, noParams);

    }


    /**
     * Initializes the fitting procedure with a new set of measurements (dependent variables).
     * 
     * @param depVals
     *            The new set of measurements.
     */
    public void newDepVals(double[] depVals) throws MarquardtMinimiserException {
        if (depVals.length != y.length - 1) {
            System.err.println(depVals.length + " " + (y.length - 1));
            throw new MarquardtMinimiserException(
                    "New data contains the wrong number of values.");
        }
        for (int i = 0; i < ndata; i++) {
            y[i + 1] = depVals[i];
        }

        // Reinitialize the starting parameters and expected variances
        // of the measurements.
        initASigs();

    }

    /**
     * Initialises the parameter values and their standard deviations.
     */
    protected void initASigs() {

 	// Do a linear DT fit

	double[] data = new double[ndata];

	for (int i = 0; i < ndata; i++) {
	    data[i] = y[i+1];
	}

	double[] params = inv.invert(data);

	DT dt = new DT(params[2], params[3], params[4], params[5], params[6], params[7]);

	double s0 = Math.exp(params[1]);

	double[][] seig = dt.sortedEigenSystem();

	Vector3D v = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);
	
	double[] tp = Vector3D.thetaPhi(v);
	
	// find angle that rotates this to the XZ plane (phi), then onto x axis (theta)

	double theta = Math.PI / 2.0 - tp[0];
	double phi = -tp[1];


	// check FA between 0 and 1

	// set mixing parameter in range 0.1 - 0.9
	
	double fa = dt.fa();

 	a[1] = Math.sqrt(s0);
 	a[2] = Math.sqrt(dt.trace() / 3.0);

	if (fa < 0.1) {
	    a[3] = Math.acos(Math.sqrt(0.1));
	}
	else if (fa > 0.9) {
	    a[3] = Math.acos(Math.sqrt(0.9));
	}
	else {
	    a[3] = Math.acos(Math.sqrt(fa));
	}

 	a[4] = theta;
 	a[5] = phi;

//         // Initialize standard deviations from residuals of DT fit

// 	double var = 0.0;

// 	for (int i = 1; i <= ndata; i++) {
// 	    double sDT = s0 * Math.exp(-bValues[i] * dt.contractBy(x[i]));

// 	    var += (y[i] - sDT) * (y[i] - sDT) / (ndata - 5.0);
// 	}

// 	double std = Math.sqrt(var);
	
//         for(int i = 1; i <= ndata; i++) {
//             sig[i] = std;
//         }
    }


    /** 
     * Initializes the standard deviations of the samples. The SD of each 
     * measurement y is sigma / y, where sigma is the noise SD.
     * We ignore the constant scaling factor of the (unknown) sigma.
     *
     */
    protected void initSigs() {
	
        for (int i = 1; i <= ndata; i++) {
	    if (y[i] > 0.0) { 
		sig[i] = 1.0 / (y[i]);
	    }
	    else {
		// lower the weight of bad data
		sig[i] = 1.0;
	    }

        }
    }




    /**
     * @return [0.0, ln(S_0), d, f, theta, phi].
     *
     */
    public double[] getParameters() {
	return new double[] {0.0, Math.log(a[1] * a[1]), a[2] * a[2], 
			     Math.cos(a[3]) * Math.cos(a[3]), a[4], a[5]};
    }


    /**
     * Sets the initial values of the parameters. The values in aInit start
     * counting from zero.
     * 
     * @param aInit
     *            Array containing the new parameter values starting from index
     *            0: [ln(S_0), d, f, theta, phi]  
     */
    public void setInitParams(double[] aInit) throws MarquardtMinimiserException {
	a[1] = Math.sqrt(Math.exp(aInit[0]));
	a[2] = Math.sqrt(aInit[1]);
	a[3] = Math.acos(Math.sqrt(aInit[2]));
	a[4] = aInit[3];
        a[5] = aInit[4];
    }




    /**
     * Returns the value of the model function at x[i] with parameters
     * atry.
     *
     * @param atry The parameters of the model to use.
     *
     * @param i The index of the data point to use.
     *
     * @return The error between the fit and the data point.
     */
    protected double yfit(double[] atry, int i) {


	double S = atry[1];
	double d = atry[2];
	double f = atry[3];
	double theta = atry[4];
	double phi = atry[5];

	double sinT = Math.sin(theta);
	double cosT = Math.cos(theta);

	double sinP = Math.sin(phi);
	double cosP = Math.cos(phi);

	double qx = x[i][0];
	double qy = x[i][1];
	double qz = x[i][2];

	double cosF = Math.cos(f);
	
	double g = (S*S)*((cosF*cosF)/ Math.exp((d*d)*bValues[i]*(qx*(qx*(cosP*cosP)*(cosT*cosT) + qy*cosP*(cosT*cosT)*sinP - qz*cosP*cosT*sinT) + qy* (qx*cosP*(cosT*cosT)*sinP + qy*(cosT*cosT)*(sinP*sinP) - qz*cosT*sinP*sinT) + qz* (-(qx*cosP*cosT*sinT) - qy*cosT*sinP*sinT + qz*(sinT*sinT)))) + (1 - (cosF*cosF))/Math.exp((d*d)*((qx*qx) + (qy*qy) + (qz*qz))*bValues[i]));

	return g;
 
    }
 

    /**
     * Returns an array of values of the derivatives of the model at
     * x[i] with respect to each of the parameters in a, using values
     * in atry. 
     *
     * @param atry The parameters of the model to use.
     *
     * @param i The index of the data point to use.
     *
     * @return The derivates of the model at atry for data point i. 
     */
    protected double[] dydas(double[] atry, int i) {

	double S = atry[1];
	double d = atry[2];
	double f = atry[3];
	double theta = atry[4];
	double phi = atry[5];


	double[] dyda = new double[noParams+1];

	double sinT = Math.sin(theta);
	double cosT = Math.cos(theta);

	double sinP = Math.sin(phi);
	double cosP = Math.cos(phi);


	double sinF = Math.sin(f);
	double cosF = Math.cos(f);

 	double qx = x[i][0];
	double qy = x[i][1];
	double qz = x[i][2];
	
	double alpha = (qx*cosP*cosT + qy*cosT*sinP - qz*sinT)*(qx*cosP*cosT + qy*cosT*sinP - qz*sinT);
	double beta = Math.exp((d*d)*bValues[i]*alpha);
	
	// In[195]:=
	// CForm[dg3]
	// Out[195]CForm=

	dyda[1] = 2*S*((cosF*cosF)/beta + (1.0 - cosF*cosF)/Math.exp((d*d)*((qx*qx) + (qy*qy) + (qz*qz))*bValues[i]));

	 
	// In[196]:=
	// CForm[dg4]
	// Out[196]//CForm=


	dyda[2] = (S*S)*((-2*d*((qx*qx) + (qy*qy) + (qz*qz))*bValues[i]*(1.0 - cosF*cosF))/ Math.exp((d*d)*((qx*qx) + (qy*qy) + (qz*qz))*bValues[i]) - (2*d*bValues[i]*(cosF*cosF)*alpha)/ beta);


	// In[197]:=
	// CForm[dg5]
	// Out[197]//CForm=

	dyda[3] = (-2*(Math.exp((d*d)*((qx*qx) + (qy*qy) + (qz*qz))*bValues[i]) - beta)*(S*S)*cosF*sinF )/Math.exp((d*d)*bValues[i]*((qx*qx) + (qy*qy) + (qz*qz) + (qx*qx)*(cosP*cosP)*(cosT*cosT) + qy*(cosT*cosT)*sinP*(2*qx*cosP + qy*sinP) - 2*qx*qz*cosP*cosT*sinT - 2*qy*qz*cosT*sinP*sinT + (qz*qz)*(sinT*sinT)));

	// 	CForm[dg1]
	// 	Out[193]//CForm=


	dyda[4] = (2*(d*d)*(S*S)*bValues[i]*(cosF*cosF)*(qy*qz*(cosT*cosT)*sinP + (qx*qx)*(cosP*cosP)*cosT*sinT + cosT*(-(qz*qz) + (qy*qy)*(sinP*sinP))*sinT - qy*qz*sinP*(sinT*sinT) + qx*cosP*(qz*(cosT*cosT) - qz*(sinT*sinT) + qy*sinP*Math.sin(2.0*theta))))/ beta;


	// 	In[194]:=
	// 	 CForm[dg2]
	// 	 Out[194]//CForm=

	dyda[5] = (-2*(d*d)*(S*S)*bValues[i]*(cosF*cosF)*cosT*(qy*cosP - qx*sinP)* (qx*cosP*cosT + qy*cosT*sinP - qz*sinT))/ beta;


	return dyda;

 }


 
}

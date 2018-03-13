package inverters;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.CL_Initializer;

import java.io.*;
import Jama.Matrix;

/**
 *
 * Uses a weighted-linear least-squares algorithm to fit the diffusion tensor. Initializes the 
 * solution to the ordinary least squares estimate, then iteratively solves for the weighted solution.
 * <p> 
 * We treat the problem as ordinary least squares solving Y' = X'B, where:
 * <BR> Y' = WY 
 * <BR> X' = WX
 * <BR> Y = log(A) ie the log diffusion-weighted data
 * <BR> X is the design matrix
 * <BR> B are the seven parameters (the unweighted log signal + 6 DT components) we solve for
 * <BR> W is the weight matrix containing weights W_{ii} = A(i)
 * <BR> A(i) is the predicted diffusion-weighted signal for measurement i. 
 * The predicted signal is initialized from an OLS tensor fit, then updated iteratively.
 * <p>
 * Bad data (A(i) <= 0) is given a zero weight.
 *
 * @author Philip Cook
 * @version $Id$
 * @see inverters.LinearDT_Inversion
 */
public class WeightedLinearDT_Inversion extends DT_Inversion {

    /**
     * This is the threshold on the singular values. If a singular value has
     * value less that the largest singular value divided by this value, it is
     * inverted to zero.
     */
    protected static final double SVTHRESH = 1.0E12;

    private final LinearDT_Inversion olsInv;

    // number of measurements per voxel
    private final int numMeas;

    // b-matrix
    private RealMatrix X;
    private RealMatrix XT;

    // optionally write noise variance to this output stream
    private DataOutputStream noiseMap = null;
    private DataOutputStream residualVectorMap = null;

    private static final int MAX_ITERATIONS = 10;


    /**
     * The constructor requires the details of the sequence used to generate the
     * data that will be processed.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     */
    public WeightedLinearDT_Inversion(DW_Scheme imParams) {
    	ip = imParams;

    	numMeas = ip.numMeasurements();

        // Set up the matrix X, which is the b-matrix

    	X = ip.getB_Matrix();

    	XT = X.transpose();

    	olsInv = new LinearDT_Inversion(imParams);

        // Initialize the outlier map output stream.
    	if (CL_Initializer.noiseVarianceMapFile != null) {
    		try {
    			FileOutputStream fout = new FileOutputStream(CL_Initializer.noiseVarianceMapFile);

    			noiseMap = 
    			new DataOutputStream(new BufferedOutputStream(fout, OutputManager.FILEBUFFERSIZE));
    		}
    		catch (java.io.IOException e) {
    			throw new LoggedException(e);
    		}
    	}

    	if (CL_Initializer.residualVectorMapFile != null) {
    		try {
    			FileOutputStream fout = new FileOutputStream(CL_Initializer.residualVectorMapFile);

    			residualVectorMap = 
    			new DataOutputStream(new BufferedOutputStream(fout, OutputManager.FILEBUFFERSIZE));
    		}
    		catch (java.io.IOException e) {
    			throw new LoggedException(e);
    		}
    	}

    }


    /**
     *
     * Fits the diffusion tensor.
     * <p> 
     * If any of the data is bad (log of the data cannot be computed), the inversion attempts to compensate by 
     * giving bad measurements zero weight. If this fails, the unweighted linear fit is returned.  
     *
     * @param data
     *            The MRI data.
     * 
     * @return {exitcode, ln A^\star(0), Dxx, Dxy, Dxz, Dyy, Dyz, Dzz}
     *
     * 
     */
    public double[] invert(double[] data) {

    	double[] olsSolution = olsInv.invert(data);

    	if ( !(olsSolution[0] >= 0) ) {
            // could not get linear fit
            return badData();
    	}

        return computeInversion(data, olsSolution);

    }


    /**
     *
     * Fits the diffusion tensor.
     * <p> 
     * Same as invert(double[] data) method except that gradient corrections are applied 
     * to the inverse matrix before it is used.  
     *
     * @param data
     *            The MRI data.
     * @param  gradAdj specifies the gradient adjustment per voxel.  
     * 
     * @return {exitcode, ln A^\star(0), Dxx, Dxy, Dxz, Dyy, Dyz, Dzz}
     *
     * 
     */
    public double[] invert(double[] data, double[] gradAdj) {

        // Applying the gradient adjustment to the matrix X.
        Matrix g = new Matrix(gradAdj,3);
        Matrix eye = Matrix.identity(3,3);
        Matrix iPlusG = eye.plus(g);
        // Loop through number of Measurements to correct each row in the B matrix.     
        for ( int i = 0; i < numMeas; i++ ) {
            double[] bvec = ip.getG_Dir(i);
            double bval = ip.getB_Value(i);
            Matrix bmat = new Matrix(bvec,bvec.length);
            // compute the new bvec
            Matrix v = iPlusG.times(bmat);
            double[][] a = v.getArray();
            // magnitude of new bvec
            double n = Math.sqrt((a[0][0]*a[0][0])+(a[1][0]*a[1][0])+(a[2][0]*a[2][0]));
            // normalise and set new bvecs and bvals
            double[] new_gDir = new double[]{a[0][0]/n, a[1][0]/n, a[2][0]/n};
            double new_bValue = (n*n)*bval;
            // reset the matrix entries. 
            X.setEntry(i, 0, 1.0);
            X.setEntry(i, 1, -new_bValue * new_gDir[0] * new_gDir[0]);
            X.setEntry(i, 2, -2.0 * new_bValue * new_gDir[0] * new_gDir[1]);
            X.setEntry(i, 3, -2.0 * new_bValue * new_gDir[0] * new_gDir[2]);
            X.setEntry(i, 4, -new_bValue * new_gDir[1] * new_gDir[1]);
            X.setEntry(i, 5, -2.0 * new_bValue * new_gDir[1] * new_gDir[2]);
            X.setEntry(i, 6, -new_bValue * new_gDir[2] * new_gDir[2]);
        }
        // recalculate the transpose of X to override the one from the constructor.
        XT = X.transpose();

        // initialize OLS solution with the gradient adjusment applied. 
        double[] olsSolution = olsInv.invert(data, gradAdj);

        if ( !(olsSolution[0] >= 0) ) {
            // could not get linear fit
            return badData();
        }

        return computeInversion(data, olsSolution);
    } 

        
    // Does the inversion and formats output
    private double[] computeInversion(double[] data, double[] olsSolution) {
        
    	RealMatrix B = null;
        
    	RealMatrix[] weightedFit = null;
        
        // Always have to output these, even if weighted fit fails
        // If fitting fails, we output zeros for these values
    	double sigmaSq = 0.0;
        
    	double[] residuals = new double[numMeas];
        
    	double[] result = new double[8];

    	result[0] = olsSolution[0];
        
        
    	try {
            weightedFit = computeWeightedFit(data, olsSolution);
            B = weightedFit[3];
            
            if (noiseMap != null) {
                sigmaSq = computeNoiseVariance(weightedFit);
            }
            if (residualVectorMap != null) {
                residuals = computeResidualVector(weightedFit, true, Math.sqrt(sigmaSq));
            }
            
    	}
    	catch(ConvergenceException e) {
            B = new RealMatrix(7,1);
            
            for (int i = 0; i < 7; i++) {
                B.entries[i][0] = olsSolution[i+1];
            }
            
            result[0] = 2.0;
    	}
    	catch(SVD_Exception e) {
            B = new RealMatrix(7,1);
            
            for (int i = 0; i < 7; i++) {
                B.entries[i][0] = olsSolution[i+1];
            }
            
            result[0] = 2.0;
    	}
        
    	if (noiseMap != null) {
            try { 
                noiseMap.writeDouble(sigmaSq);
            }
            catch (IOException e) {
                throw new LoggedException(e);
            }
    	}
        
    	if (residualVectorMap != null) {
            try { 
                for (int i = 0; i < numMeas; i++) {
                    residualVectorMap.writeDouble(residuals[i]);
                }
            }
            catch (IOException e) {
                throw new LoggedException(e);
            }
    	}
        
        
    	for (int i = 0; i < 7; i++) {
            result[i+1] = B.entries[i][0];
    	}
        
    	return result;
        
    }

    
    private double[] badData() {
        double[] badData = new double[8];
        badData[0] = -100.0;

        // writes zeros for additional output if needed
        background();

        return badData;
        
    }


    /**
     * Closes the noise map output stream, if open.
     */
    public void close() {

    	if (noiseMap != null) {
    		try { 
    			noiseMap.close();
    		}
    		catch (IOException e) {
    			throw new LoggedException(e);
    		}
    	}
    	if (residualVectorMap != null) {
    		try { 
    			residualVectorMap.close();
    		}
    		catch (IOException e) {
    			throw new LoggedException(e);
    		}
    	}

    }


    /**
     * In background voxels this inverter needs to output zero if there is a noise map.
     */
    public void background() {
    	if (noiseMap != null) {

    		try { 
    			noiseMap.writeDouble(0.0);  
    		}
    		catch (IOException e) {
    			throw new LoggedException(e);
    		}

    	}
    	if (residualVectorMap != null) {

    		try { 
    			for (int i = 0; i < numMeas; i++) {
    				residualVectorMap.writeDouble(0.0);
    			}
    		}
    		catch (IOException e) {
    			throw new LoggedException(e);
    		}

    	}
    }


    /**
     * 
     * @param weightedSolution array array [Y, W, XB, B, H], returned from the
     * computeWeightedFit method.
     *
     * The first three elements are used to compute the weighted residuals, where
     * Y = log of the DWI measurements, W = diag[exp(XB(i))], ie diagonal matrix containing 
     * predicted DWI measurements, and XB = product of design matrix with parameter matrix.
     *
     * @return WY - WXB
     *
     */
     public double[] computeResidualVector(RealMatrix[] weightedSolution) {
         return computeResidualVector(weightedSolution, false, 0.0);
     }


    /**
     *
     * @param weightedSolution array [Y, W, XB, B, H].
     * @param studentize if true, studentize residuals
     * @param sigma estimated noise SD, used only for studentized residuals
     *
     * @return Residuals for each measurement.
     */
    public double[] computeResidualVector(RealMatrix[] weightedSolution, boolean studentize, double sigma) {
    	RealMatrix Y = weightedSolution[0];
    	RealMatrix W = weightedSolution[1];
    	RealMatrix XB = weightedSolution[2];
	RealMatrix H = weightedSolution[4];

    	RealMatrix E = W.product(Y).sub(W.product(XB));

    	double[] result = new double[numMeas];

        if (studentize) {
    	    for (int i = 0; i < numMeas; i++) {
    		result[i] = E.entries[i][0] / (sigma * Math.sqrt(1.0 - H.entries[i][i]));
            }
	}
        else {
            for (int i = 0; i < numMeas; i++) {
                result[i] = E.entries[i][0];
            }
        }

    	return result;
    }
    

    /**
     * Computes the noise variance \sigma^2 = (Y - XB)^T W^2 (Y - XB) / (DOF), where 
     * DOF is the number of degrees of freedom in the model. This is  
     * (number of measurements used to fit model) - 7. Bad data is excluded during the fitting 
     * process, so the DOF is reduced accordingly.
     *
     * @param weightedSolution array [Y, W, XB]. The array [Y, W, XB, B], returned from the
     * computeWeightedFit method, may also be used as the parameter (but B is not required). 
     * The first three elements are used to compute the noise variance, where
     * Y = log of the DWI measurements, W = diag[exp(XB(i))], ie diagonal matrix containing 
     * predicted DWI measurements, and XB = product of design matrix with parameter matrix.
     * 
     *
     * @return \sigma^2. Assumes that the residuals have zero mean and equal variance. 
     *
     */
    public double computeNoiseVariance(RealMatrix[] weightedSolution) {

    	RealMatrix Y = weightedSolution[0];
    	RealMatrix W = weightedSolution[1];
    	RealMatrix XB = weightedSolution[2];

    	RealMatrix WSQ = new RealMatrix(numMeas, numMeas);

    	int numExcludedMeas = 0;

    	for (int i = 0; i < numMeas; i++) {

    		WSQ.entries[i][i] = W.entries[i][i] * W.entries[i][i];

    		if (WSQ.entries[i][i] == 0.0) {
    			numExcludedMeas++;
    		}

    	}

    	RealMatrix sumResiduals = (Y.sub(XB).transpose()).product(WSQ).product(Y.sub(XB));


	// degrees of freedom in data is numMeas - 7 - numExcludedMeas
	// ie number of measurements actually used in the calculation of B, minus the number of
	// parameters in B

    	double sigmaSq = sumResiduals.entries[0][0] / (numMeas - 7 - numExcludedMeas);

    	return sigmaSq;

    }


    /**
     * Solves WY = WXB for B. 
     * Y = log(data)
     * W = weights, 
     * B = [ln A^\star(0), Dxx, Dxy, Dxz, Dyy, Dyz, Dzz]^T.
     * H = projection matrix = (WX)(WX)^*, where (WX)^* is the pseudo inverse used to compute B = (WX)^*(WY)
     * The fit is initialized with the results of a linear inversion.
     * 
     *
     * @return {Y, W, XB, B, H}.
     * 
     */
    public RealMatrix[] computeWeightedFit(double[] data) throws SVD_Exception, ConvergenceException {
    	return computeWeightedFit(data, olsInv.invert(data));
    }


    /**
     * @param olsSolution array in the format returned by LinearDT_Inversion.invert(data).  This solution is used to 
     * initialize B, and hence W, where W(i,i) = exp(XB(i)) for each measurement i.
     *
     * @return {Y, W, XB, B, H}.
     * 
     */
    public RealMatrix[] computeWeightedFit(double[] data, double[] olsSolution) throws SVD_Exception, ConvergenceException {

    	double[] logData = new double[numMeas];

    	RealMatrix Y = new RealMatrix(numMeas, 1);

	// B = [ln[S_0], dxx, dxy, dxz, dyy, dyz, dzz]^T
	// initialize from result of OLS solution
    	RealMatrix B = new RealMatrix(7,1);

	// weights
    	RealMatrix W = new RealMatrix(numMeas, numMeas);


    	for (int i = 0; i < 7; i++) {
    		B.entries[i][0] = olsSolution[i+1];
    	}	

    	RealMatrix XB = X.product(B);

    	double sumResidualSq = 0.0;

    	boolean[] badData = new boolean[numMeas];

    	int numExcludedMeas = 0;

	// if exit code from linear fit is zero, there are no problems with bad data
    	if (olsSolution[0] == 0.0) {

    		for (int i = 0; i < numMeas; i++) {

    			logData[i] = Math.log(data[i]);

    			Y.entries[i][0] = logData[i];

    			W.entries[i][i] = Math.exp(XB.entries[i][0]);

    			sumResidualSq += (logData[i] - XB.entries[i][0]) * (logData[i] - XB.entries[i][0]); 
    		}

    	}
    	else {

	    // initialize weights with 0.0 for bad data

    		for (int i = 0; i < numMeas; i++) {

    			badData[i] = !(data[i] > 0.0);

    			if (badData[i]) {
    				numExcludedMeas++;
    				logData[i] = 0.0;
    			}
    			else {
    				logData[i] = Math.log(data[i]);
    				sumResidualSq += (logData[i] - XB.entries[i][0]) * (logData[i] - XB.entries[i][0]); 
    			}

    			Y.entries[i][0] = logData[i];
    			W.entries[i][i] = badData[i] ? 0.0 : Math.exp(XB.entries[i][0]);

    		}

    	}

    	if (numMeas - numExcludedMeas < 8) {

	    // not enough data to estimate DT + S_0 + \sigma
    		throw new SVD_Exception("Insufficient data to invert WX");
    	}


    	int iter = 0;

	// difference in residuals from last iteration, | \epsilon_{i+1}^2 - \epsilon_i^2 | / \epsilon_{i}^2
    	double relDiffResidual = Double.MAX_VALUE;

	// quit iterating if change in residuals is less than this
    	double convergeThresh = 1E-2;

    	RealMatrix inv = null;

    	while (relDiffResidual > convergeThresh && iter < MAX_ITERATIONS) {

    		RealMatrix[] svd = null;

    		svd = W.product(X).svd();

	    // Find the maximum singular value.
    		double maxSV = svd[1].entry(0, 0);
    		for (int i = 0; (i < svd[1].rows() && i < svd[1].columns()); i++) {
    			if (svd[1].entry(i, i) > maxSV) {
    				maxSV = svd[1].entry(i, i);
    			}
    		}

	    // Invert the singular values.
    		for (int i = 0; (i < svd[1].rows() && i < svd[1].columns()); i++) {
    			double curSV = svd[1].entry(i, i);
    			svd[1].setEntry(i, i, (curSV > maxSV / SVTHRESH) ? (1.0 / curSV) : 0.0);
    		}

	    // Get the pseudo inverse.
    		inv = svd[2].product(svd[1].transpose()).product(svd[0].transpose());

    		B = inv.product(W.product(Y));

    		double nextSumResidualSq = 0.0;	    

    		XB = X.product(B);

    		for (int i = 0; i < numMeas; i++) {
		// continue to set bad data to zero weight, but weight good measurements
		// to the model parameters
    			if (!badData[i]) {
    				nextSumResidualSq += (logData[i] - XB.entries[i][0]) * (logData[i] - XB.entries[i][0]); 
    			}

    			W.entries[i][i] = badData[i] ? 0.0 : Math.exp(XB.entries[i][0]);

    		}

    		if (sumResidualSq > 0.0) {
    			relDiffResidual = Math.abs(nextSumResidualSq - sumResidualSq) / sumResidualSq;
    		}
    		else {
    			relDiffResidual = 0.0;
    		}

    		sumResidualSq = nextSumResidualSq;

    		iter++;
    	}	


    	if (relDiffResidual > convergeThresh * 10) {
    		throw new ConvergenceException("relDiffResidual == " + relDiffResidual);
    	}

    	return new RealMatrix[] {Y, W, XB, B, W.product(X).product(inv)};

    }

}

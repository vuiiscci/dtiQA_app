package inverters;

import misc.*;
import numerics.*;
import imaging.*;
import Jama.Matrix;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fits the diffusion tensor.
 * 
 * <dt>Description:
 * 
 * <dd>Uses the standard least squares linear approach using singular value
 * decomposition.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class LinearDT_Inversion extends DT_Inversion {

    /**
     * This matrix is the linear inversion matrix. It is B^-1, where the rows of
     * B are {1, -b g_1^2, -2 b g_1 g_2, -2 b g_1 g_3, -b g_2^2, -2 b g_2 g_3,
     * -b g_3^2}.
     */
    protected RealMatrix linearInv;

    /**
     * This is the threshold on the singular values. If a singular value has
     * value less that the largest singular value divided by this value, it is
     * inverted to zero.
     */
    protected static final double SVTHRESH = 1.0E12;


    /** 
     * B, the B-matrix, the inverse of linearInv
     */
    protected RealMatrix B;

    /**
     * Number of measurements we expect
     */
    private int numMeas;

    /**
     * The constructor requires the details of the sequence used to generate the
     * data that will be processed.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     */
    public LinearDT_Inversion(DW_Scheme imParams) {
        ip = imParams;

        numMeas = ip.numMeasurements();

        // Set up the matrix X.
        B = ip.getB_Matrix();

	// Now get the singular value decomposition.
        RealMatrix[] svd = null;
        try {
            svd = B.svd();
        }
        catch (Exception e) {
            throw new LoggedException(e);
        }

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
        linearInv = svd[2].product(svd[1].transpose()).product(svd[0].transpose());
    }



    /**
     * Fits the diffusion tensor using the inverse matrix.
     * 
     * @param data
     *            The MRI data.
     * 
     * @return {exitcode, ln A^\star(0), Dxx, Dxy, Dxz, Dyy, Dyz, Dzz}
     */
    public double[] invert(double[] data) {

        // The failure free exit code is 0.
        double exitCode = 0.0;

        if (data.length != numMeas) {
            throw new LoggedException("Wrong number of measurements, expected " + numMeas);
        }

        // Make the data matrix.
        RealMatrix A = new RealMatrix(numMeas, 1);
        for (int i = 0; i < numMeas; i++) {
            if (data[i] > 0.0) {
                A.setEntry(i, 0, Math.log(data[i]));
            }
            else {
                A.setEntry(i, 0, 0.0);

                // The exit code is 6 if the data is bad and has to
                // be changed to perform the inversion.
                exitCode = 6.0;
            }
        }

        RealMatrix D;

        // Do the inversion.
        if (exitCode == 0.0) {
            D = linearInv.product(A);
        }
        else {
            // some bad data, remove it

            // basically do a weighted fit with one iterations, and binary weights
            // 0 for bad data

            // Data == A, already defined
                    
            // Weights
            RealMatrix W = new RealMatrix(numMeas, numMeas);
            
            boolean[] badData = new boolean[numMeas];

            int numExcludedMeas = 0;
            
	    for (int i = 0; i < numMeas; i++) {
		
		badData[i] = !(data[i] > 0.0);

		if (badData[i]) {
		    numExcludedMeas++;
		}

		W.entries[i][i] = badData[i] ? 0.0 : 1.0;
		
	    }

	
            if (numMeas - numExcludedMeas < 7) {
                
                // Solution might not exist even if there are 7 measurements, because they might not have the necessary
                // 6 non-collinear gradients
                
                // but definitely nothing to do if there are < 7 measurements
                return fitFailure(data);
            }

            // The linearInv for this voxel, with bad measurements zeroed out
	    RealMatrix inv = null;
	    
	    RealMatrix[] svd = null;

            try {
                svd = W.product(B).svd();
            }
            catch (Exception e) {
                return fitFailure(data);
            }

	    
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
	    
	    D = inv.product(W.product(A));
            
        }

        // Create the return array.
        double[] res = new double[8];
        res[0] = exitCode;
        for (int i = 0; i < 7; i++) {
            res[i + 1] = D.entry(i, 0);
        }

        return res;
    }

     /**
     * Fits the diffusion tensor after applying gradient corrections to the inverse matrix.
     * 
     * @param data
     *            The MRI data.
     * @param gradAdj specifies the gradient adjustment per voxel. 
     * 
     * @return {exitcode, ln A^\star(0), Dxx, Dxy, Dxz, Dyy, Dyz, Dzz}
     */
    public double[] invert(double[] data, double[] gradAdj){
        //RealMatrix B_Adj = ip.getB_Matrix. not sure if need to get new copy of B 
        //matrix or is ok to use current one as function only called with 2 args. 
         
        // Applying the grad adj to B matrix.
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
            // reset the B Matrix entries. 
            B.setEntry(i, 0, 1.0);
            B.setEntry(i, 1, -new_bValue * new_gDir[0] * new_gDir[0]);
            B.setEntry(i, 2, -2.0 * new_bValue * new_gDir[0] * new_gDir[1]);
            B.setEntry(i, 3, -2.0 * new_bValue * new_gDir[0] * new_gDir[2]);
            B.setEntry(i, 4, -new_bValue * new_gDir[1] * new_gDir[1]);
            B.setEntry(i, 5, -2.0 * new_bValue * new_gDir[1] * new_gDir[2]);
            B.setEntry(i, 6, -new_bValue * new_gDir[2] * new_gDir[2]);
               
        }
        // Copy stuff from constructor to repeat singular value decomp to get a new linearInv
        RealMatrix[] svd = null;
        try {
            svd = B.svd();
        }
        catch (Exception e) {
            throw new LoggedException(e);
        }

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
        linearInv = svd[2].product(svd[1].transpose()).product(svd[0].transpose());

        // Now do original invert method.
        return invert(data);
       
    }


    /**
     * Returns the linear inversion matrix.
     * 
     * @return the linear inversion matrix for fitting the DT.
     */
    public RealMatrix getMatrix() {
        return (RealMatrix)linearInv.clone();
    }


    // could do some debugging here, but for now just return -100
    private double[] fitFailure(double[] data) {
        double exitCode = -100.0;
        
        double[] res = new double[8];
        
        res[0] = exitCode;
        
        return res;
    }

}

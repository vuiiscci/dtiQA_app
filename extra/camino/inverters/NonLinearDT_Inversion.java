package inverters;

import imaging.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fit the diffusion tensor.
 * 
 * <dt>Description:
 * 
 * <dd>Fits the diffusion tensor using Levenberg-Marquardt and constrains the
 * diffusion tensor to be positive definite by using the Cholesky decomposition.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: NonLinearDT_Inversion.java,v 1.4 2005/08/18 11:07:44 ucacmgh
 *          Exp $
 *  
 */
public class NonLinearDT_Inversion extends DT_Inversion {

    /**
     * The diffusion tensor fitter.
     */
    protected DiffTensorFitter fitter;


    /**
     * The inversion uses a linear inversion to provide a starting
     * point for the non-linear fitting procedure.
     */
    protected LinearDT_Inversion ldti;


    /**
     * This constructor takes only the imaging sequence parameters and creates a
     * fitter that does not use the Hessian of the objective function.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     */
    public NonLinearDT_Inversion(DW_Scheme imParams) {
        init(imParams, ModelIndex.NLDT_POS);
    }

    /**
     * This constructor takes the imaging sequence parameters and an index
     * specifying the type of fitting to use.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param index
     *            The fitter index.
     */
    public NonLinearDT_Inversion(DW_Scheme imParams, ModelIndex index) {
        init(imParams, index);
    }

    /**
     * Sets up the fitter to use for inversion.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param index
     *            The fitter index.
     */
    protected void init(DW_Scheme imParams, ModelIndex index) {

        ip = imParams;

        // Set up the linear fitter
        ldti = new LinearDT_Inversion(ip);

        // Set up the non-linear fitter.

	double[][] g = ip.getNonZeroG_Dirs();
	
	double[] b = ip.getNonZeroB_Values();

	int M = ip.numZeroMeasurements();

	// dummy data
	double[] dummy = new double[ip.numMeasurements() - M];

        try {
	    if (index == ModelIndex.NLDT_POS) {
                fitter = new DiffTensorFitter(g, dummy, b, M);
            }
	    else if (index == ModelIndex.NLDT) {        
                fitter = new DiffTensorUnConFitter(g, dummy, b, M);
            }
	    else if (index == ModelIndex.NLDT_MINPACK) {        
                fitter = new DiffTensorMinPackFitter(g, dummy, b, M);
            }
	    else {
		throw new LoggedException("Unknown nonlinear DT inversion " + index);
	    }
        }
        catch (Exception e) {
            throw new LoggedException(e);
        }
    }

    int testindex = 0;

    /**
     * Fits the diffusion tensor.
     * 
     * @param data
     *            The MRI data.
     * 
     * @return {exitcode, log A*(0), Dxx, Dxy, Dxz, Dyy, Dyz, Dzz}
     */
    public double[] invert(double[] data) {

        // Do the linear inversion to get a starting point
        double[] linearInv = ldti.invert(data);

        // Copy the exit code from the linear inversion as a default.
        double exitCode = linearInv[0];

        // Get the normalized data.
        double[] normData = ip.normalizeData(data);

        try {
            fitter.newDepVals(normData);
            fitter.setStartPoint(linearInv);
            fitter.minimise();
        }
        catch (Exception e) {
            LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
	    
            // If the fitter fails, the exit code is 2.
            exitCode = 2.0;
        }

        DT d = fitter.getDiffTensor();
        //System.err.println("Test index: " + testindex);
        testindex++;
        //System.err.println(d);
        double[] dArr = d.getComponents();

        // Construct the array to return. The first element is
        // the exit code. The second element is ln A^\star(0).
        double[] res = new double[8];
        res[0] = exitCode;
        res[1] = linearInv[1];
        for (int i = 0; i < 6; i++) {
            res[i + 2] = dArr[i];
        }

        return res;
    }

}

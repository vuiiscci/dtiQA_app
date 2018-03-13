package inverters;

import imaging.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fit the three-tensor model.
 * 
 * <dt>Description:
 * 
 * <dd>Fits the three-tensor using Levenberg-Marquardt. The starting point
 * comes from a single DT fit.
 * 
 * The exact fitter, ie constraints and parametrization, can be specified using
 * the three tensor index and one tensor index in the constructor.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class ThreeTensorInversion extends DiffusionInversion {

    /**
     * The diffusion tensor inversion.
     */
    protected DT_Inversion dtInverter;

    /**
     * The three-tensor fitter.
     */
    protected ThreeTensorFitter fitter;

    /**
     * The number of elements of the output array.
     */
    public static final int ITEMSPERVOX = 24;

    /**
     * This constructor requires the details of the sequence used to generate
     * the data that will be processed. It creates an object with default
     * choices for the ThreeTensorFitter (unconstrained) and tensor fitter (linear).
     * 
     * @param imParams
     *            The imaging sequence parameters.
     */
    public ThreeTensorInversion(DW_Scheme imParams) {
        init(imParams, ModelIndex.DTDTDT, ModelIndex.LDT);
    }

    /**
     * This constructor requires the details of the sequence and the indexes of
     * the choices of ThreeTensorFitter and DT_Inverter.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param twoTensorIndex
     *            Index of the ThreeTensorFitter.
     * 
     * @param dtInverterIndex
     *            Index of the DT_Inverter.
     */
    public ThreeTensorInversion(DW_Scheme imParams, ModelIndex threeTensorIndex,
            ModelIndex dtInverterIndex) {
        init(imParams, threeTensorIndex, dtInverterIndex);
    }

    /**
     * Called by the constructors to create a ThreeTensorInverter with specified
     * parameters.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param threeTensorIndex
     *            Index of the ThreeTensorFitter.
     * 
     * @param dtInverterIndex
     *            Index of the DT_Inverter. 
     */
    protected void init(DW_Scheme imParams, ModelIndex threeTensorIndex,
			ModelIndex dtInverterIndex) {
        ip = imParams;

        // Set up the single DT inversion.
        dtInverter = DT_Inversion.getIndexedDT_Inversion(dtInverterIndex, imParams);

        // Set up the three-tensor fitter.
        int M = ip.numZeroMeasurements();
        double[] bValues = ip.getNonZeroB_Values();

        double[][] gNoZeros = ip.getNonZeroG_Dirs();

        // Initialize it with dummy data.
        double[] dummyData = new double[ip.numMeasurements() - M];

        try {
            fitter = ThreeTensorFitter.getIndexedThreeTensorFitter(gNoZeros, dummyData,
                    bValues, M, threeTensorIndex);
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
    }


    /**
     * Fits the model.
     * 
     * @param data
     *            The MRI data.
     * 
     * @return {exitcode, ln A^\star(0), 3, mix1, D1xx, D1xy, D1xz, D1yy,
     *         D1yz, D1zz, D2xx, mix2, D2xy, D2xz, D2yy, D2yz, D2zz, mix3, 
     *         D3xx, D3xy, D3xz, D3yy, D3yz, D3zz}
     */
    public double[] invert(double[] data) {

        // Get the normalized data.
        double[] normData = ip.normalizeData(data);

        // Get the single DT.
        double[] dtInversion = dtInverter.invert(data);

        // Need the array without the exit code and ln zero measurement.
        double errorCode = dtInversion[0];
        double logQ0 = dtInversion[1];
        double[] dt = new double[6];
        for (int i = 0; i < 6; i++) {
            dt[i] = dtInversion[i + 2];
        }
        DT singleDT = new DT(dt);

        // Initialise the two-tensor fitter.
        try {
            fitter.newDepVals(normData);
            fitter.setStartFromSingleDT(singleDT);
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }

        // Initialize the fitted parameters from the single DT
        // in case the fitting fails.
        double[] d1Comps = singleDT.getComponents();
        double[] d2Comps = singleDT.getComponents();
        double[] d3Comps = singleDT.getComponents();
        double mix1 = 1.0 / 3.0;
        double mix2 = 1.0 / 3.0;

        // Now do the fitting.
        try {
            fitter.minimise();
            d1Comps = fitter.getDT1().getComponents();
            d2Comps = fitter.getDT2().getComponents();
            d3Comps = fitter.getDT3().getComponents();
            mix1 = fitter.getMix1();
            mix2 = fitter.getMix2();
        }
        catch (Exception e) {
	    logger.info(e.toString() + "Fitting failed.  Outputting single DT.");
            errorCode = 2.0;
        }

        double[] params = new double[ITEMSPERVOX];
        params[0] = errorCode;
        params[1] = logQ0;
        params[2] = 3.0;
        params[3] = mix1;
        params[10] = mix2;
        params[17] = 1.0-mix1-mix2;
        for (int i = 0; i < 6; i++) {
            params[i + 4] = d1Comps[i];
            params[i + 11] = d2Comps[i];
            params[i + 18] = d3Comps[i];
        }

        return params;
    }

    
    public int itemsPerVoxel() {
        return ITEMSPERVOX;
    }

}

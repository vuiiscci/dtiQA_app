package inverters;

import imaging.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fit the two-tensor model.
 * 
 * <dt>Description:
 * 
 * <dd>Fits the two-tensor using Levenberg-Marquardt. The starting point comes
 * from a single DT fit.
 * 
 * The exact fitter, ie constraints and parametrization, can be specified using
 * the two tensor index and one tensor index in the constructor.
 * 
 * </dl>
 * 
 * @see inverters.ModelIndex 
 *
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class TwoTensorInversion extends DiffusionInversion {

    /**
     * The diffusion tensor inversion.
     */
    protected DT_Inversion dtInverter;

    /**
     * The two-tensor fitter.
     */
    protected TwoTensorFitter fitter;

    /**
     * The number of elements of the output array.
     */
    public static final int ITEMSPERVOX = 17;

    /**
     * This constructor requires the details of the sequence used to generate
     * the data that will be processed. It creates an object with default
     * choices for the TwoTensorFitter (unconstrained) and DT inverter (linear).
     * 
     * @param imParams
     *            The imaging sequence parameters.
     */
    public TwoTensorInversion(DW_Scheme imParams) {
        init(imParams, ModelIndex.DTDT, ModelIndex.LDT);
    }

    /**
     * This constructor requires the details of the sequence and the indexes of
     * the choices of TwoTensorFitter and DT_Inverter.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param twoTensorIndex
     *            Index of the TwoTensorFitter.
     * 
     * @param dtInverterIndex
     *            Index of the DT_Inverter.
     */
    public TwoTensorInversion(DW_Scheme imParams, ModelIndex twoTensorIndex,
            ModelIndex dtInverterIndex) {
        init(imParams, twoTensorIndex, dtInverterIndex);
    }

    /**
     * Called by the constructors to create a TwoTensorInverter with specified
     * parameters.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     * 
     * @param twoTensorIndex
     *            Index of the TwoTensorFitter.
     * 
     * @param dtInverterIndex
     *            Index of the DT_Inverter. 
     */
    protected void init(DW_Scheme imParams, ModelIndex twoTensorIndex,
			ModelIndex dtInverterIndex) {
        ip = imParams;

        // Set up the single DT inversion.
        dtInverter = DT_Inversion.getIndexedDT_Inversion(dtInverterIndex, imParams);

        // Set up the two-tensor fitter.
        int M = ip.numZeroMeasurements();

        double[] bValues = ip.getNonZeroB_Values();

        double[][] gNoZeros = ip.getNonZeroG_Dirs();

        // Initialize it with dummy data.
        double[] dummyData = new double[ip.numMeasurements() - M];

        try {
            fitter = TwoTensorFitter.getIndexedTwoTensorFitter(gNoZeros, dummyData, bValues,
                    M, twoTensorIndex);
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
     * @return {exitcode, ln A^\star(0), 2, mix1, D1xx, D1xy, D1xz, D1yy, D1yz,
     *         D1zz, mix2, D2xx, D2xy, D2xz, D2yy, D2yz, D2zz}
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
        double mix = 0.5;

        // Now do the fitting.
        try {
            fitter.minimise();
            d1Comps = fitter.getDT1().getComponents();
            d2Comps = fitter.getDT2().getComponents();
            mix = fitter.getMix();
        }
        catch (Exception e) {
	    logger.info(e.toString() + "Fitting failed.  Outputting single DT.");
            errorCode = 2.0;
        }

        double[] params = new double[ITEMSPERVOX];
        params[0] = errorCode;
        params[1] = logQ0;
        params[2] = 2.0;
        params[3] = mix;
        params[10] = 1.0-mix;
        for (int i = 0; i < 6; i++) {
            params[i + 4] = d1Comps[i];
            params[i + 11] = d2Comps[i];
        }

        return params;
    }


    public int itemsPerVoxel() {
        return ITEMSPERVOX;
    }

}

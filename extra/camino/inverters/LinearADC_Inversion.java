package inverters;

import numerics.*;
import imaging.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fits the apparent diffusion coefficient on the assumption of isotropic
 * diffusion.
 * 
 * <dt>Description:
 * 
 * <dd>Uses the standard least squares linear approach using singular value
 * decomposition. The output is {exit code, ln A^\star(0), ADC}.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 */
public class LinearADC_Inversion extends DiffusionInversion {

    /**
     * This matrix is the linear inversion matrix. It is X^-1, where the rows of
     * X are {1, -b}.
     */
    protected RealMatrix linearInv;

    /**
     * This is the threshold on the singular values. If a singular value has
     * value less that the largest singular value divided by this value, it is
     * inverted to zero.
     */
    protected static final double SVTHRESH = 1.0E12;

    /**
     * The number of data items returned by the inversion.
     */
    public static final int ITEMSPERVOX = 3;


    /**
     * The constructor requires the details of the sequence used to generate the
     * data that will be processed.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     */
    public LinearADC_Inversion(DW_Scheme imParams) {
        ip = imParams;

        // Set up the matrix X.
        RealMatrix X = new RealMatrix(ip.numMeasurements(), 2);

        for (int i = 0; i < X.rows(); i++) {
            double[] g = ip.getG_Dir(i);

            X.setEntry(i, 0, 1.0);
            X.setEntry(i, 1, -ip.getB_Value(i));

        }

        // Now get the singular value decomposition.
        RealMatrix[] svd = null;
        try {
            svd = X.svd();
        }
        catch (Exception e) {
            throw new RuntimeException(e);
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
     * Fits the apparent diffusion coefficient using the inverse matrix.
     * 
     * @param data
     *            The MRI data.
     * 
     * @return {exitcode, ln A^\star(0), ADC}
     */
    public double[] invert(double[] data) {

        // The failure free exit code is 0.
        double exitCode = 0.0;

        // Make the data matrix.
        RealMatrix A = new RealMatrix(data.length, 1);
        for (int i = 0; i < data.length; i++) {
            if (data[i] > 0) {
                A.setEntry(i, 0, Math.log(data[i]));
            }
            else {
                A.setEntry(i, 0, 0.0);

                // The exit code is 6 if the data is bad and has to
                // be changed to perform the inversion.
                exitCode = 6.0;
            }
        }

        // Do the inversion.
        RealMatrix D = linearInv.product(A);

        // Create the return array.
        double[] res = new double[3];
        res[0] = exitCode;
        for (int i = 0; i < 2; i++) {
            res[i + 1] = D.entry(i, 0);
        }

        return res;
    }

    public int itemsPerVoxel() {
        return ITEMSPERVOX;
    }

}

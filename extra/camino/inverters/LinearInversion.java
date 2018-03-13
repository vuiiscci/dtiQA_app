package inverters;

import numerics.*;
import imaging.*;
import data.*;
import java.util.Vector;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Performs a linear inversion
 * 
 * <dt>Description:
 * 
 * <dd>Reads any linear transformation matrix in from a file and uses it to
 * transform input data.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class LinearInversion extends DiffusionInversion {

    /**
     * This matrix is the linear inversion matrix.
     */
    protected RealMatrix linearInv;

    /**
     * Specifies whether or not to normalize the data before reconstruction.
     */
    protected boolean normalize;

    /**
     * Specifies whether or not to use the log measurements.
     */
    protected boolean useLog;

    /**
     * The constructor requires the details of the sequence used to generate the
     * data that will be processed and the name of the file containing the
     * linear transformation matrix.
     * 
     * @param imParams
     *            The imaging sequence parameters.
     */
    public LinearInversion(DW_Scheme imParams, String matrixFile, boolean norm,
            boolean log) {

        ip = imParams;
        normalize = norm;
        useLog = log;

        // The number of columns in the matrix must equal the number of
        // measurements, which we know from the schemefile.
        int mxCols = ip.numMeasurements();
        if (normalize) {
            mxCols = mxCols - ip.numZeroMeasurements();
        }

        // Read in the matrix chunk by chunk
        VoxelOrderDataSource vods = new VoxelOrderDataSource(matrixFile, mxCols, "double");
        Vector<Object> v = new Vector<Object>();
        while (vods.more())
            try {
                v.addElement(vods.nextVoxel());
            }
            catch (Exception e) {
                System.err.println("Could not read the matrix.  Expecting " + mxCols
                        + " elements per row.");
                throw new RuntimeException(e);
            }

        // Construct the matrix
        int mxRows = v.size();
        linearInv = new RealMatrix(mxRows, mxCols);
        for (int i = 0; i < mxRows; i++) {
            double[] row = (double[]) v.elementAt(i);
            for (int j = 0; j < mxCols; j++) {
                linearInv.entries[i][j] = row[j];
            }
        }
    }

    /**
     * Fits the diffusion tensor using the inverse matrix.
     * 
     * @param data
     *            The MRI data.
     * 
     * @return {exitcode, ln A^\star(0), p1, ..., pN}
     */
    public double[] invert(double[] data) {

        // The failure free exit code is 0.
        double exitCode = 0.0;

        // Do normalization if required.
        double[] pData = data;
        if (normalize) {
            pData = ip.normalizeData(data);
        }

        // Make the data matrix.
        RealMatrix A = new RealMatrix(pData.length, 1);
        for (int i = 0; i < pData.length; i++) {
            if (useLog) {
                if (data[i] > 0) {
                    A.setEntry(i, 0, Math.log(pData[i]));
                }
                else {
                    A.setEntry(i, 0, 0.0);

                    // The exit code is 6 if the data is bad and has to
                    // be changed to perform the inversion.
                    exitCode = 6.0;
                }
            }
            else {
                A.setEntry(i, 0, pData[i]);
            }
        }

        // Do the inversion.
        RealMatrix recon = linearInv.product(A);

        // Create the return array.
        int outputExtraSlots = normalize ? 2 : 1;
        double[] res = new double[linearInv.rows() + outputExtraSlots];
        res[0] = exitCode;
        if (normalize) {
            double meanZ = ip.geoMeanZeroMeas(data);
            res[1] = (meanZ > 0.0) ? Math.log(meanZ) : 0.0;
        }
        for (int i = 0; i < linearInv.rows(); i++) {
            res[i + outputExtraSlots] = recon.entry(i, 0);
        }

        return res;
    }

    /**
     * Returns the length of the arrays returned by invert.
     */
    public int itemsPerVoxel() {
        int items = linearInv.rows() + (normalize ? 2 : 1);
        return items;
    }

    /**
     * Returns the linear inversion matrix.
     * 
     * @return the linear inversion matrix for fitting the DT.
     */
    public RealMatrix getMatrix() {
        return linearInv;
    }

}

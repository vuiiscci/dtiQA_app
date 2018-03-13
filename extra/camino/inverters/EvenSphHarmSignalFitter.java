package inverters;

import numerics.*;
import imaging.*;
import misc.*;
import tools.ArrayOps;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Fits an even spherical harmonic series to diffusion MRI signal as
 * a function of direction.  
 * 
 * <dt>Description:
 * 
 * <dd> Differs from EvenSphHarmFitter, which fits an ADC profile.  This
 * class simply fits the spherical harmonic series to the raw diffusion
 * weighted signal.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: EvenSphHarmFitter.java 728 2009-07-15 21:16:01Z ucacpco $  
 */
public class EvenSphHarmSignalFitter extends EvenSphHarmFitter {


    /**
     * The constructor computes the inversion matrix from the imaging parameters
     * for the sequence used in the data acquisition.
     * 
     * @param imParams
     *            The imaging parameters of the acquisition sequence.
     * 
     * @param order
     *            The maximum order of the spherical harmonic series to fit.
     */
    public EvenSphHarmSignalFitter(DW_Scheme imParams, int order) {

        maxOrder = order;
        ip = imParams;

        meanNonZeroB = ArrayOps.mean(ip.getNonZeroB_Values());

        // Compute the number of free parameters in the spherical
        // harmonic coefficients for a real even series up to the
        // specified order (see MRM02).
        numParams = (maxOrder + 1) * (maxOrder / 2 + 1);

        // Set up the matrix X.
        X = new RealMatrix(ip.numMeasurements(), numParams + 1);

        for (int i = 0; i < X.rows(); i++) {
            Vector3D g = new Vector3D(ip.getG_Dir(i));
            
 	    double[] tp = Vector3D.thetaPhi(g);

 	    double theta = tp[0];
 	    double phi = tp[1];

            X.setEntry(i, 0, 1.0);

            int nextEntry = 1;

            for (int l = 0; l <= maxOrder; l += 2)
                try {
                    Complex c = SphericalHarmonics.Y(l, 0, theta, phi);
                    X.setEntry(i, nextEntry, c.real());
                    nextEntry += 1;

                    for (int m = 1; m <= l; m++) {
                        c = SphericalHarmonics.Y(l, m, theta, phi);
                        X.setEntry(i, nextEntry, -2.0 * c.real());
                        X.setEntry(i, nextEntry + 1, 2.0 * c.imag());
                        nextEntry += 2;
                    }
                }
                catch (Exception e) {
                    // This should never happen.
                    throw new RuntimeException(e);
                }

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
        invX = svd[2].product(svd[1].transpose()).product(svd[0].transpose());

    }

    /**
     * Fits the spherical harmonic series using the inverse matrix.
     * 
     * @param data
     *            The MRI data.
     * 
     * @return The fitted series. The array contains [exitCode, log A^\star(0),
     *         c00, c20, Re(c21), Im(c21), Re(c22), Im(c22), c40, Re(c41),
     *         Im(c41), ...]
     */
    public double[] fit(double[] data) {

        // The failure free exit code is 0.
        double exitCode = 0.0;

        // Make the data matrix.
        RealMatrix A = new RealMatrix(data.length, 1);
        for (int i = 0; i < data.length; i++) {
            if (data[i] > 0) {
                A.setEntry(i, 0, data[i]);
            }
            else {
                A.setEntry(i, 0, 0.0);

                // The exit code is 6 if the data is bad and has to
                // be changed to perform the inversion.
                exitCode = 6.0;
            }
        }

        // Do the inversion.
        RealMatrix S = invX.product(A);

        // Create the return array.
        double[] res = new double[numParams + 2];
        res[0] = exitCode;
        for (int i = 0; i < numParams + 1; i++) {
            res[i + 1] = S.entry(i, 0);
        }

        return res;
    }


}

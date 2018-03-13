package apps;

import java.io.*;
import java.text.*;
import java.util.logging.Logger;

import tools.*;
import data.*;
import misc.LoggedException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Computes simple statistics of data files.
 * 
 * <dt>Description:
 * 
 * <dd>Computes the mean, mean squared and variance of the signal in
 * each component of a diffusion-weighted MRI data set within the
 * foreground region specified by a mask or background threshold.
 * Useful for estimating the noise level from background regions.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: ModelFit.java 58 2006-10-26 11:38:18Z ucacmgh $
 *  
 */
public class DataStats {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.ModelFit");

    public static void main(String[] args) {

        // Parse the command line arguments
        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
        CL_Initializer.initImagingScheme();
        CL_Initializer.initDataSynthesizer();


        // Loop over the data
        int fgVox = 0;
        double[] means = new double[CL_Initializer.imPars.numMeasurements()];
        double[] meanSqs = new double[CL_Initializer.imPars.numMeasurements()];

        while (CL_Initializer.data.more())
            try {

                double[] nextVoxel = CL_Initializer.data.nextVoxel();

                // Fit or output background default.
                double backgroundB0 = CL_Initializer.imPars.geoMeanZeroMeas(nextVoxel);
                if (!ModelFit.isBG(backgroundB0)) {
                    fgVox++;
                    for(int i=0; i<nextVoxel.length; i++) {
                        means[i] += nextVoxel[i];
                        meanSqs[i] += nextVoxel[i]*nextVoxel[i];
                    }
                }

            }
            catch (DataSourceException e) {
                throw new LoggedException("The data file does not contain a whole number of voxels." +
                                          "Check the scheme file. Got Exception " + e);
            }


        // Output the results
        DecimalFormat sciForm = new DecimalFormat("0.000000E00;-0.000000E00");

        System.out.println("Foreground voxel count: " + fgVox);
        System.out.println("Component   E(S)           E(S^2)         Var(S)         Std(S)");
        Object[] o = new Object[1];
        for(int i=0; i<means.length; i++) {

            means[i] = means[i]/(double)fgVox;
            meanSqs[i] = meanSqs[i]/(double)fgVox;

            o[0] = new Integer(i+1);
            System.out.print(MessageFormat.format("{0,number, 0;-0}", o));

            o[0] = new Double(means[i]);
            System.out.print("          " + sciForm.format(((Number) o[0]).doubleValue()));
            
            o[0] = new Double(meanSqs[i]);
            System.out.print("    " + sciForm.format(((Number) o[0]).doubleValue()));
            
            o[0] = new Double(meanSqs[i] - means[i]*means[i]);
            System.out.print("    " + sciForm.format(((Number) o[0]).doubleValue()));
            
            o[0] = new Double(Math.sqrt(meanSqs[i] - means[i]*means[i]));
            System.out.println("    " + sciForm.format(((Number) o[0]).doubleValue()));
        }


    }



}

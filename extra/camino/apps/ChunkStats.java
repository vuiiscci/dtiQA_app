package apps;

import java.util.logging.Logger;

import data.*;
import misc.*;
import tools.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Computes statistics of chunks of data.
 * 
 * <dt>Description:
 * 
 * <dd>Reads in a data stream consisting of repeated sets of samples,
 * eg the output of multidimensional MCMC sampling, and computes
 * statistics, eg mean, std, max, min, for each element of the set
 * over a chunk of repeats.
 *
 * For example for a series of input values:
 *
 * S11
 * S12
 * S13
 * S21
 * S22
 * S23
 * :
 * :
 * SN1
 * SN2
 * SN3
 *
 * ChunkStats -chunksize 3 -samples 2 -mean
 *
 * returns
 * 
 * (S11+S21)/2
 * (S12+S22)/2
 * (S13+S23)/2
 * :
 * :
 * (S(N-1)1+SN1)/2
 * (S(N-1)2+SN2)/2
 * (S(N-1)3+SN3)/2
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: SphFuncPD_Stats.java 889 2010-11-21 01:04:29Z ucacpco $
 *  
 */
public class ChunkStats extends Executable {

    private static Logger logger = Logger.getLogger("camino.apps.ChunkStats");



    /**
     * Size of each chunk.
     */
    private int chunkSize;

    /**
     * Switches for different statistics.
     */
    private boolean mean;
    private boolean std;
    private int max;
    private int min;


    public ChunkStats(String[] args) {
        super(args);
    }


    public void initDefaultVals() {
        
        chunkSize = 1;

        mean = false;
        std = false;
        max = -1;
        min = -1;

    }


    public void initOptions(String[] args) {

        // Parse the command line arguments
        CL_Initializer.inputDataType = "double";
        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-chunksize")) {
                chunkSize = (Integer.parseInt(args[i + 1]));
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-mean")) {
                mean = true;
                CL_Initializer.markAsParsed(i,1);
            }
            if (args[i].equals("-std")) {
                std = true;
                CL_Initializer.markAsParsed(i,1);
            }
            if (args[i].equals("-max")) {
                max = (Integer.parseInt(args[i + 1]));
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-min")) {
                min = (Integer.parseInt(args[i + 1]));
                CL_Initializer.markAsParsed(i,2);
            }
            
        }

        CL_Initializer.checkParsing(args);

        if(max >= chunkSize || min >= chunkSize || max < -1 || min < -1)
            throw new LoggedException("Index to max or min out of range.");

    }


    public void initVariables() {
        CL_Initializer.data = new VoxelOrderDataSource(CL_Initializer.inputFile, chunkSize, CL_Initializer.inputDataType);
    }

    
    public void execute(OutputManager om) {
            
        // Loop over the data
        while (CL_Initializer.data.more()) try {
                
            double[] setMean = null;
            double[] setStd = null;
            double[] setMax = null;
            double[] setMin = null;
            double maxSoFar = 0;
            double minSoFar = 0;

            for(int i=0; i<CL_Initializer.samples; i++) {

                // Read in the coefficients.
                double[] coeffs = CL_Initializer.data.nextVoxel();

                // If first chunk, initialize.
                if(i==0) {
                    setMean = new double[coeffs.length];
                    setStd = new double[coeffs.length];
                    setMax = new double[coeffs.length];
                    setMin = new double[coeffs.length];
                    for(int j=0; j<coeffs.length; j++) {
                        setMean[j] = coeffs[j];
                        setStd[j] = coeffs[j]*coeffs[j];
                        setMax[j] = coeffs[j];
                        setMin[j] = coeffs[j];
                    }
                    maxSoFar = (max>=0)?coeffs[max]:0;
                    minSoFar = (min>=0)?coeffs[min]:0;
                }

                else {
                    boolean bigger = max>=0 && coeffs[max]>maxSoFar;
                    maxSoFar = bigger?coeffs[max]:maxSoFar;
                    boolean smaller = min>=0 && coeffs[min]<minSoFar;
                    minSoFar = smaller?coeffs[min]:minSoFar;

                    for(int j=0; j<coeffs.length; j++) {
                        setMean[j] += coeffs[j];
                        setStd[j] += coeffs[j]*coeffs[j];
                        setMax[j] = bigger?coeffs[j]:setMax[j];
                        setMin[j] = smaller?coeffs[j]:setMin[j];
                    }

                }

            }

            // Normalize etc.
            for(int j=0; j<setMean.length; j++) {
                setMean[j] = setMean[j]/CL_Initializer.samples;
                setStd[j] = Math.sqrt(setStd[j]/(CL_Initializer.samples-1) - setMean[j]*setMean[j]);
            }

            // Output it.
            if(mean)
                om.output(setMean);
            if(std)
                om.output(setStd);
            if(max>=0)
                om.output(setMax);
            if(min>=0)
                om.output(setMin);

        }
        catch (Exception e) {
            throw new LoggedException(e);
        }

        // Tidy up.
        om.close();
    }


}

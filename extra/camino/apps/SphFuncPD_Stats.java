package apps;

import java.util.logging.Logger;

import data.*;
import misc.*;
import tools.*;
import numerics.*;
import sphfunc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Computes the peak directions of spherical functions.
 * 
 * <dt>Description:
 * 
 * <dd>Reads coefficients of the spherical functions and uses random sampling
 * and subsequent optimization to find the list of peak directions of the
 * function.
 * 
 * Uses standard input and output streams for input and output data.
 * 
 * The output for each voxel is:
 *
 *  - exit code inherited from reconstruction.
 *
 *  - ln(A(0))
 *
 *  - number of directions after pruning. 
 *
 *  - flag for consistency with repeated run (number of directions is
 * the same and the directions are the same to within a threshold.)
 *
 *  - mean(f). 
 *
 *  - std(f). 
 *
 *  - direction 1 (x, y, z, f, H00, H01, H10, H11). 
 *
 *  - direction 2 (x, y, z, f, H00, H01, H10, H11).
 *
 *  - up to max directions
 * 
 * H is the Hessian of f. It is the matrix:
 * 
 * [d^2f/ds^2 d^2f/dsdt] [d^2f/dtds d^2f/dt^2]
 * 
 * where s and t are orthogonal coordinates local to the peak.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class SphFuncPD_Stats extends Executable {

    private static Logger logger = Logger.getLogger("camino.apps.SphFuncPD_Stats");



    /**
     * Base threshold on the actual pd strength divided by the mean of the
     * function.
     */
    private double pdThresh;

    /**
     * This is the number of standard deviations of the function to be added to
     * the pdThresh in the pd pruning.
     */
    private double stdsFromMeanF;

    /**
     * The sampling density in the principal direction extraction algorithm.
     */
    private int density;

    /**
     * The total number global statistics of the input function output
     * before the list of PDs
     */
    public static final int GLOBALSTATS = 6;

    /**
     * The total number of values output for each direction.
     */
    public static final int STATSPERPD = 8;

    /**
     * Can be used to turn off the binary output for debugging.
     */
    private boolean noOutput;

    /**
     * Turn on to skip the consistency check and speed things up a bit.
     */
    private boolean doConsistencyCheck;


    public SphFuncPD_Stats(String[] args) {
        super(args);
    }


    public void initDefaultVals() {
        
        pdThresh = 1.0;
        stdsFromMeanF = 0.0;
        density = 1000;
        noOutput = false;
        doConsistencyCheck = true;
    }


    public void initOptions(String[] args) {

        // Parse the command line arguments
        CL_Initializer.inputDataType = "double";
        CL_Initializer.CL_init(args);
        
	if(CL_Initializer.inputModel.equals("maxent")) {
	    CL_Initializer.initImagingScheme();
	    CL_Initializer.initMaxEnt();
	}
        
        CL_Initializer.initSphFuncDataSource();

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-noconsistencycheck")) {
                doConsistencyCheck = false;
                CL_Initializer.markAsParsed(i);
            }
	    else if (args[i].equals("-pdthresh") || args[i].equals("-meanthreshscale")) {
                pdThresh = Double.parseDouble(args[i + 1]);
                CL_Initializer.markAsParsed(i, 2);
            }
	    else if (args[i].equals("-stdsfrommean") || args[i].equals("-stdthresh")) {
                stdsFromMeanF = Double.parseDouble(args[i + 1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-density")) {
                density = (Integer.parseInt(args[i + 1]));
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-nooutput")) {
                noOutput = true;
                CL_Initializer.markAsParsed(i);
            }
            
        }

        CL_Initializer.checkParsing(args);

    }


    public void initVariables() {

    }

    
    public void execute(OutputManager om) {
            
        // Initialize the point set for sampling in the peak finder.
        double[][] samplePoints = null;
        double[][] samplePointsRot = null;
        if(CL_Initializer.pointSetIndSet) try {
                samplePoints = ISCodes.getPointSetForMaxEnt(CL_Initializer.pointSetInd).data;
                RealMatrix randomRot = Rotations.randomRotMat(0);
                samplePointsRot = SphericalPoints.rotatePointSet(samplePoints, randomRot);
            } catch(Exception e) {
                throw new LoggedException(e);
            }
        else {
            samplePoints = SphericalPoints.getIcosahedralPointSet(density, 0);
            samplePointsRot = SphericalPoints.getIcosahedralPointSet(density, density);
        }
        
        // Loop over the data
        while (CL_Initializer.data.more())
            
            try {
                
                // Read in the coefficients.
                double[] coeffs = CL_Initializer.data.nextVoxel();

                // Test for background.
                if(coeffs[0]>=0 && Math.exp(coeffs[1]) > CL_Initializer.BACKGROUNDTHRESHOLD) {

                    // Construct a spherical function using the
                    // coefficients.
                    SphericalFunction sf = null;
                    if (CL_Initializer.inputModel.equals("sh")) {
                        sf = new EvenSHS(coeffs, CL_Initializer.maxOrder);
                    }
                    else if (CL_Initializer.inputModel.equals("maxent")) {
                        sf = new MaxEntProfile(coeffs, CL_Initializer.kernelParams);
                    }
                    else {
                        sf = new TuchRBF_Sum(coeffs);
                    }
                    
                    
                    // Get the list of peak directions.
                    PDList pds = new PDList(SphericalFunction.getSearchRadius());
                    double[] distances = sf.getPDs(samplePoints, pds);

                    //Check that the algorithm did not do anything
                    //strange.
                    boolean allOK = true;
                    for (int i = 0; i < distances.length; i++) {
                        
                        if (distances[i] > SphericalFunction.getSearchRadius()) {
                            allOK = false;
                        }
                        else if (distances[i] < 0.0) {
                            allOK = false;
                        }
                    }
                    
                    int noPDsPrePrune = pds.getNoPDs();
                    
                    //Prune any principal direction whose value is less
                    //than the a threshold, which by default is the mean
                    //of the function.
                    double[] fStats = sf.getStats();
                    double meanF = fStats[0];
                    double stdF = fStats[1];
                    pds.pruneByValue(pdThresh * meanF + stdsFromMeanF * stdF);
                    
                    
                    //Get PDs again using a separate set of starting
                    //location.  This time, we don't worry about the fine
                    //optimization.  We don't need to find the maxima that
                    //accurately to check consistency.
                    boolean pdsOK = true;
                    if(doConsistencyCheck) {
                        PDList pdtest = sf.getPDsRS(samplePointsRot);
                        int noPDsPrePruneTest = pds.getNoPDs();
                        pdtest.pruneByValue(pdThresh * meanF + stdsFromMeanF * stdF);
                                        
                        //Compare the sets of pds.
                        pdsOK = pds.equivalent(pdtest, pdThresh);
                    }
                    
                    // Collate all the information to output.
                    double[] results = new double[GLOBALSTATS + CL_Initializer.numPDsIO * STATSPERPD];

                    // Inherit the exitcode and ln(A(0)) from the input data.
                    results[0] = coeffs[0];
                    results[1] = coeffs[1];
                    
                    results[2] = (double) pds.getNoPDs();
                    results[3] = allOK ? 1.0 : 0.0;
                    
                    results[4] = meanF;
                    results[5] = stdF;
                    
                    for (int i = 0; i < CL_Initializer.numPDsIO; i++) {
                        if (i < pds.getNoPDs()) {
                            results[GLOBALSTATS + i*STATSPERPD] = pds.getPD(i).getPDX();
                            results[GLOBALSTATS + i*STATSPERPD + 1] = pds.getPD(i).getPDY();
                            results[GLOBALSTATS + i*STATSPERPD + 2] = pds.getPD(i).getPDZ();
                            results[GLOBALSTATS + i*STATSPERPD + 3] = pds.getPD(i).getProp();
                            
                            RealMatrix hess = sf.getHessian(pds.getPD(i));
                            //RealMatrix[] eig = hess.jacobi();
                            
                            results[GLOBALSTATS + i*STATSPERPD + 4] = hess.entries[0][0];
                            results[GLOBALSTATS + i*STATSPERPD + 5] = hess.entries[0][1];
                            results[GLOBALSTATS + i*STATSPERPD + 6] = hess.entries[1][0];
                            results[GLOBALSTATS + i*STATSPERPD + 7] = hess.entries[1][1];
                        }
                        else {
                            for (int j = 0; j < STATSPERPD; j++) {
                                results[GLOBALSTATS + i*STATSPERPD + j] = 0.0;
                            }
                        }
                    } //for i

                    // Output it.
                    om.output(results);
                }

                // Otherwise output some default values to indicate
                // background.
                else {
                    double[] results = new double[GLOBALSTATS + CL_Initializer.numPDsIO * STATSPERPD];
                    results[0] = coeffs[0];
                    results[1] = coeffs[1];

                     om.output(results);
                }

            }
            catch (Exception e) {
                throw new LoggedException(e);
            }

        // Tidy up.
        om.close();
    }


    /**
     * Set the sampling density in the peak search.
     *
     * @param density
     *            The density to set.
     */
    public void setDensity(int density) {
        this.density = density;
    }


    /**
     * Set the maximum number of peaks to output per voxel.
     *
     * @param noOutput
     *            The noOutput to set.
     */
    public void setNoOutput(boolean noOutput) {
        this.noOutput = noOutput;
    }


    /**
     * Set the threshold on difference between peak directions in the
     * consistency test.
     *
     * @param pdThresh
     *            The pdThresh to set.
     */
    public void setPdThresh(double pdThresh) {
        this.pdThresh = pdThresh;
    }


    /**
     * Adjust the threshold on small peaks to ignore.
     *
     * @param stdsFromMeanF
     *            The stdsFromMeanF to set.
     */
    public void setStdsFromMeanF(double stdsFromMeanF) {
       this.stdsFromMeanF = stdsFromMeanF;
    }
}

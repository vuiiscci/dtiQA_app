package apps;

import imaging.*;
import misc.*;
import numerics.*;
import tools.*;
import data.*; //added for OutputManager

import java.io.*;
import java.util.logging.*;
import java.util.Random;
import java.util.Scanner;
import java.text.*;


/**
 * Produces information about the quality of electrostatic subsets, and produces scheme
 * files for those subsets given the scheme file for the full set.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class SubsetScheme extends Executable{
    
    /**
     * logging object
     */
    private static final Logger logger = Logger.getLogger("apps/SubsetScheme");
    
    
    public SubsetScheme(String[] args){
        super(args);
    }
    
    // contains non-zero normalized qs
    int fullSetSize;
    double[][] fullSetPoints;
    int numSubsets;
    // non zero qs in subsets. Subset schemes will contain these plus
    // all q=0 measurements
    int[] subsetSizes;
    // subsetPoints[s] is the set of points for set s
    double[][][] subsetPoints;
    // total measurements in full scheme
    int fullMeas;
    String imageListFile;
    String[] imageList;
    String outputRoot;
    String subsetFile;
    
    public void initDefaultVals(){
        fullSetSize = 0;
        fullSetPoints = null;
        numSubsets = 0;
        subsetSizes = null;
        subsetPoints = null;
        fullMeas = 0;
        imageListFile = null;
        imageList = null;
        outputRoot = "subset";
        subsetFile = null;
    }
    
    public void initOptions(String[] args){
        //    public static void main(String[] args) {
        
        if (args.length == 0) {
            System.err.println("SubsetScheme -schemefile <schemefile> -subsetpoints <points> [-imagelist <list> -outputroot <root>]");
            System.exit(0);
        }
        
        CL_Initializer.CL_init(args);
        
        CL_Initializer.initImagingScheme();
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-subsetpoints")) {
                subsetFile = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-imagelist")) {
                imageListFile = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-outputroot")) {
                outputRoot = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
        }
        
        CL_Initializer.checkParsing(args);
    }
    
    public void initVariables(){
        fullMeas = CL_Initializer.imPars.numMeasurements();
        
        fullSetPoints = CL_Initializer.imPars.getNonZeroG_Dirs();
        
        fullSetSize = fullSetPoints.length;
    }
    
    
    public void execute(OutputManager om){
        // get subset points
        try {
            Scanner subsetHeader = new Scanner(new File(subsetFile)).useDelimiter("\\s+");
            
            // total number of points in subset scheme
            int totalSubsetPoints = subsetHeader.nextInt();
            
            if (fullSetSize != totalSubsetPoints) {
                throw new LoggedException("Full scheme has " + fullSetSize + " non-zero g, but " +
                "subset is from a set of " + totalSubsetPoints + " points");
            }
            
            numSubsets = subsetHeader.nextInt();
            
            subsetSizes = new int[numSubsets];
            
            for (int s = 0; s < numSubsets; s++) {
                subsetSizes[s] = subsetHeader.nextInt();
            }
            
            // file contains all points in scheme
            double[][] tmp = PointSetToScheme.readPoints(subsetFile, false, false, false);
            
            int pointCounter = 0;
            
            subsetPoints = new double[numSubsets][][];
            
            for (int s = 0; s < numSubsets; s++) {
                
                subsetPoints[s] = new double[subsetSizes[s]][];
                
                for (int p = 0; p < subsetSizes[s]; p++) {
                    subsetPoints[s][p] = tmp[pointCounter++];
                }
                
            }
            
            
            
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }
        
        
        
        
        if (imageListFile != null) {
            try {
                Scanner images = new Scanner(new File(imageListFile));
                
                images.useDelimiter("\\n");
                
                imageList = new String[fullMeas];
                
                for (int i = 0; i < fullMeas; i++) {
                    imageList[i] = images.next();
                }
                
            }
            catch (IOException e) {
                logger.warning("Can't process image list " + imageListFile + ", image list will " +
                "not be written");
                
                imageList = null;
            }
        }
        
        for (int s = 0; s < numSubsets; s++) {
            
            int numB0 = CL_Initializer.imPars.numZeroMeasurements();
            
            int subsetMeas = subsetSizes[s] + numB0;
            
            // which measurements are in the subset
            int[] subsetIndices = new int[subsetMeas];
            
            int indexCounter = 0;
            
            // iterate over measurements in scheme, not points in fullSet (want to get
            // the indices right including the b=0).
            
            // first get all of the b=0 measurements
            for (int m = 0; m < fullMeas; m++) {
                if (CL_Initializer.imPars.zero(m)) {
                    subsetIndices[indexCounter++] = m;
                }
                
            }
            
            for (int i = 0; i < subsetSizes[s]; i++) {
                
                // subset measurement is the closest measurement to the subset direction
                
                int maxInd = 0;
                
                double maxDot = 0.0;
                
                for (int m = 0; m < fullMeas; m++) {
                    
                    double[] g = CL_Initializer.imPars.getG_Dir(m);
                    
                    // no Math.abs here because we are looking for a subset of the directions,
                    // not a common axis
                    // precision errors prevent an exact match
                    double dotProd =
                    g[0] * subsetPoints[s][i][0] +
                    g[1] * subsetPoints[s][i][1] +
                    g[2] * subsetPoints[s][i][2];
                    
                    // fails if we can't resolve all points, but this should not happen with
                    // a decent angular separation (> 1 degree should be fine)
                    if (dotProd > maxDot) {
                        maxInd = m;
                        maxDot = dotProd;
                    }
                }
                
                
                // sanity check
                if (maxDot < 0.99984769516) {
                    logger.warning("Couldn't match point " + i + " in subset " + s +
                    " to a scheme direction within 1 degree. Closest angle was " +
                    Math.acos(maxDot) * 180.0 / Math.PI);
                }
                
                subsetIndices[indexCounter++] = maxInd;
                
            }
            
            if (indexCounter != subsetIndices.length) {
                throw new
                LoggedException("Could not match all points in subset to points " +
                "in the full set");
            }
            
            
            DW_Scheme subsetScheme = CL_Initializer.imPars.getSubsetScheme(subsetIndices);
            FileOutput out = new FileOutput(outputRoot +  "_" + s + ".scheme");
            out.writeString(subsetScheme.toString());
            out.close();
            
            out = new FileOutput(outputRoot +  "_" + s + ".indices");
            
            for (int i = 0; i < subsetMeas; i++) {
                out.writeString(subsetIndices[i] + "\n");
            }
            
            out.close();
            
            if (imageList != null) {
                
                out = new FileOutput(outputRoot + "_" + s + ".imagelist");
                
                for (int i = 0; i < subsetMeas; i++) {
                    out.writeString(imageList[subsetIndices[i]] + "\n");
                }
                
                out.close();
            }
            
            // now test how this subset ranks against a randomly chosen sample
            
            DecimalFormat df = new DecimalFormat("0.000");
            
            double[][] optimalSet = SphericalPoints.getElecPointSet(subsetSizes[s]);
            
            double[] randomEnergy = getRandomSubsetEnergy(fullSetPoints, subsetSizes[s], 100,
            new MTRandom(CL_Initializer.seed));
            
            double meanRandomEnergy = ArrayOps.mean(randomEnergy);
            double medianRandomEnergy = ArrayOps.median(randomEnergy);
            
            System.out.println("Energy of subset " + s + "\t\t" + df.format(pointsetEnergy(subsetPoints[s])));
            System.out.println("Energy of optimal " + subsetSizes[s] + "-point set\t" +
            df.format(pointsetEnergy(optimalSet)));
            System.out.println("Mean energy of random subset\t" +
            df.format(meanRandomEnergy));
            System.out.println("Median energy of random subset\t" +
            df.format(medianRandomEnergy));
            System.out.println();
        }
        
    }
    
    
    
    /**
     * Gets a randomly chosen subset of points.
     *
     *
     */
    public static double[][] getRandomSubset(double[][] points, int size, Random ran) {
        
        int fullSize = points.length;
        
        boolean[] chosen = new boolean[fullSize];
        
        int numChosen = 0;
        
        while (numChosen < size) {
            
            int next = ran.nextInt(fullSize);
            
            if (!chosen[next]) {
                chosen[next] = true;
                numChosen++;
            }
        }
        
        int indexCounter = 0;
        
        double[][] subsetPoints = new double[size][3];
        
        for (int i = 0; i < fullSize; i++) {
            if (chosen[i]) {
                subsetPoints[indexCounter][0] = points[i][0];
                subsetPoints[indexCounter][1] = points[i][1];
                subsetPoints[indexCounter][2] = points[i][2];
                indexCounter++;
            }
        }
        
        return subsetPoints;
        
    }
    
    
    /**
     * Gets a sample of random subsets, and returns their energy.
     *
     *
     */
    public static double[] getRandomSubsetEnergy(double[][] points, int size, int samples,
    Random ran) {
        
        double[] subsetEnergy = new double[samples];
        
        for (int i = 0; i < samples; i++) {
            double[][] subset = getRandomSubset(points, size, ran);
            subsetEnergy[i] = pointsetEnergy(subset);
        }
        
        return subsetEnergy;
        
    }
    
    
    /**
     * Calculates the total energy of a point set.
     *
     */
    public static double pointsetEnergy(double[][] points) {
        double[][] energy = OrderedAcqMinimizer.pairEnergyMatrix(points);
        
        double subsetEnergy = 0.0;
        
        for (int n = 0; n < points.length; n++) {
            for (int m = 0; m < n; m++) {
                subsetEnergy += energy[n][m];
            }
        }
        
        return subsetEnergy;
    }
    
}

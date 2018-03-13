package apps;

import data.*;
import inverters.*;
import misc.*;
import numerics.*;
import tools.*;

import java.text.*;
import java.util.*;
import java.util.logging.Logger;

/** 
 * Iterates over the hemisphere, fitting synthetic data at each point and evaluating its error statistics.
 * Useful for evaluating scheme files / subset schemes. For the sake of simplicity, we sample theta and 
 * phi evenly, which oversamples at the pole relative to the equator. An alternative would be to use
 * a dense point set and then interpolate the results onto a grid somehow. However, this is difficult to 
 * visualize because you have to project the results onto a 2D grid in a principled way. 
 * <p>
 * This app could be improved and extended by making use of InversionStats methods.
 * 
 * @author Philip Cook
 *
 */
 public class OrientationBiasMap extends Executable {
     
     private static Logger logger = Logger.getLogger("camino.apps.OrientationBiasMap");
     
     private int numSamples;

     private double testFuncFA;

     private double testFuncTr;

     // degrees
     private int thetaStep;
     private int phiStep;

     private double alpha;

     private String outputRoot;

 
     public OrientationBiasMap(String[] args){
         super(args);
     }
     
     
     public void initDefaultVals() {

         testFuncFA = 0.8;
         testFuncTr = 2.1E-9;
         
         numSamples = 500;

         thetaStep = 2;
         phiStep = 4;

         alpha = 0.05;

         outputRoot = "orientationMap_";
         
     }
     
     /**
      * default commandline parsing and initialisation
      */
     public void initOptions(String[] args){
         
         // Parse the command line arguments
         CL_Initializer.CL_init(args);

         for (int i = 0; i < args.length; i++) {
             
             if (args[i].equals("-samples")) {
                 numSamples = Integer.parseInt(args[i+1]);
                 CL_Initializer.markAsParsed(i, 2);
             }
             if (args[i].equals("-step")) {
                 thetaStep =  Integer.parseInt(args[i+1]);
                 phiStep =  Integer.parseInt(args[i+2]);
                 CL_Initializer.markAsParsed(i, 3);
             }
             if (args[i].equals("-fa")) {
                 testFuncFA = Double.parseDouble(args[i+1]);
                 CL_Initializer.markAsParsed(i, 2);
             }
             if (args[i].equals("-trace")) {
                 testFuncTr = Double.parseDouble(args[i+1]);
                 CL_Initializer.markAsParsed(i, 2);
             }
             if (args[i].equals("-outputroot")) {
                 outputRoot = args[i+1];
                 CL_Initializer.markAsParsed(i, 2);
             }
             if (args[i].equals("-alpha")) {
                 alpha = Double.parseDouble(args[i+1]);
                 CL_Initializer.markAsParsed(i, 2);
             }
     
         }

        
         CL_Initializer.checkParsing(args);

         CL_Initializer.initImagingScheme();

         if (CL_Initializer.inversionIndices.length > 1) {
             throw new LoggedException("This program fits single DT models only");
         }
         if (CL_Initializer.inversionIndices[0].numDTs > 1) {
             throw new LoggedException("This program fits single DT models only");
         }

     }
     

     public void initVariables() {


     }
	 
     
     public void execute(OutputManager om) {

         int thetaGridSize = 1 + (int)Math.floor((90.0 - thetaStep / 2.0) / thetaStep);
         int phiGridSize = 1 + (int)Math.floor((360.0 - phiStep / 2.0) / phiStep);


         double[][] faMean = new double[thetaGridSize][phiGridSize];

         double[][] faStd = new double[thetaGridSize][phiGridSize];

         double[][] trStd = new double[thetaGridSize][phiGridSize];

         double[][] trMean = new double[thetaGridSize][phiGridSize];

         double[][] rdMean = new double[thetaGridSize][phiGridSize];

         double[][] rdStd = new double[thetaGridSize][phiGridSize];

         double[][] adMean = new double[thetaGridSize][phiGridSize];
                  
         double[][] adStd = new double[thetaGridSize][phiGridSize];

         
         double[][] cu = new double[thetaGridSize][phiGridSize];

         // angle (degrees) between mean dyad E1 and the true test function orientation
         double[][] e1BiasAngle = new double[thetaGridSize][phiGridSize];

         // L1 of the mean dyad, 1 if the samples are completely aligned
         double[][] dyadL1 = new double[thetaGridSize][phiGridSize];

         DT tensor = getTensor(testFuncFA, testFuncTr);

         double[][] testFuncEig = tensor.sortedEigenSystem();

         double[] testFuncEigenvalues = new double[] {testFuncEig[0][0], testFuncEig[0][1], testFuncEig[0][2]};

         DT_Inversion inv = DT_Inversion.getIndexedDT_Inversion(CL_Initializer.inversionIndices[0], CL_Initializer.imPars);

         DecimalFormat df = new DecimalFormat("0.000E00");

         MTRandom ran = new MTRandom(CL_Initializer.seed);


         for (int i = 0; i < thetaGridSize; i++) {

             // degrees
             double theta = (i + 0.5) * thetaStep;

             //    logger.info("theta == " + df.format(theta));

             for (int j = 0; j < phiGridSize; j++) {

                 double phi = (j + 0.5) * phiStep;

                 Vector3D orientation = Vector3D.vectorFromSPC(1.0, theta * Math.PI / 180.0, phi * Math.PI / 180.0);

                 DT testFuncD = Rotations.rotateTensor(tensor, orientation);

                 GaussianMixture pdf = new GaussianMixture(new DT[] {testFuncD}, new double [] {1.0});

                 DataSynthesizer ds = new DataSynthesizer(pdf, CL_Initializer.imPars, CL_Initializer.SNR, numSamples, ran);

                 double[] fa = new double[numSamples];
                 double[] tr = new double[numSamples];
                 double[] rd = new double[numSamples];
                 double[] ad = new double[numSamples];

                 Vector3D[] e1 = new Vector3D[numSamples];

                 for (int v = 0; v < numSamples; v++) {

                     double[] fit = inv.invert(ds.nextVoxel());
                     
                     DT noisyTensor = new DT(fit[2], fit[3], fit[4], fit[5], fit[6], fit[7]);

                     fa[v] = noisyTensor.fa();
                     
                     double[][] seig = noisyTensor.sortedEigenSystem();

                     e1[v] = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);

                     tr[v] = seig[0][0] + seig[0][1] + seig[0][2];

                     ad[v] = seig[0][0];

                     rd[v] = (seig[0][1] + seig[0][2]) / 2.0;

                     
                 }
                 
                 faMean[i][j] = ArrayOps.mean(fa);
                 faStd[i][j] = Math.sqrt(ArrayOps.var(fa, faMean[i][j]));

                 trMean[i][j] = ArrayOps.mean(tr);
                 trStd[i][j] = Math.sqrt(ArrayOps.var(tr, trMean[i][j]));

                 rdMean[i][j] = ArrayOps.mean(rd);
                 rdStd[i][j] = Math.sqrt(ArrayOps.var(rd, rdMean[i][j]));

                 adMean[i][j] = ArrayOps.mean(ad);
                 adStd[i][j] = Math.sqrt(ArrayOps.var(ad, adMean[i][j]));

                 EigenSystem3D dyad = SphericalDistributionFitter.tBarEigenSystem(e1);
                 
                 dyadL1[i][j] = dyad.eigenvalues[0];

                 e1BiasAngle[i][j] = Math.acos(Math.abs(dyad.eigenvectors[0].dot(orientation))) * 180.0 / Math.PI;

                 cu[i][j] = coneOfUncertainty(e1, dyad.eigenvectors[0], alpha);
             }
         }

         logger.info("Writing results");
         
         logger.info( "\n    Theta range: (min step max):  " + df.format(thetaStep / 2.0) + " " + df.format(thetaStep) + 
                     " " + df.format( (thetaGridSize - 0.5) * thetaStep ) + "\n");

         logger.info( "\n    Phi range: (min step max):  " + df.format(phiStep / 2.0) + " " + df.format(phiStep) + 
                     " " + df.format( (phiGridSize - 0.5) * phiStep ) + "\n");

         logger.info( "\n    Test function FA:  " + df.format(testFuncFA) + "\n");
         logger.info( "\n    Test function Tr:  " + df.format(testFuncTr) + "\n");
         logger.info( "\n    Test function eigenvalues:  " + df.format(testFuncEigenvalues[0]) + " " + df.format(testFuncEigenvalues[1]) + " " + df.format(testFuncEigenvalues[2])); 

         DecimalFormat dataDF =  new DecimalFormat("0.00000E00");
         
         // Write results to CSV files
         writeMatrixCSV(faMean, outputRoot + "meanfa.csv", dataDF);
         writeMatrixCSV(faStd, outputRoot + "stdfa.csv", dataDF);

         writeMatrixCSV(trMean, outputRoot + "meantr.csv", dataDF);
         writeMatrixCSV(trStd, outputRoot + "stdtr.csv", dataDF);

         writeMatrixCSV(rdMean, outputRoot + "meanrd.csv", dataDF);
         writeMatrixCSV(rdStd, outputRoot + "stdrd.csv", dataDF);

         writeMatrixCSV(adMean, outputRoot + "meanad.csv", dataDF);
         writeMatrixCSV(adStd, outputRoot + "stdad.csv", dataDF);

         writeMatrixCSV(dyadL1, outputRoot + "dyadl1.csv", dataDF);
         writeMatrixCSV(e1BiasAngle, outputRoot + "dyadbias.csv", dataDF);

         writeMatrixCSV(cu, outputRoot + "coneUncertainty.csv", dataDF); 
     }


     protected void writeMatrixCSV(double[][] matrix, String fileName, DecimalFormat df) {

         int cols = matrix[0].length;

         StringBuffer buffer = new StringBuffer();

          for (int i = 0; i < matrix.length; i++) {
              for (int j = 0; j < cols - 1; j++) {
                  buffer.append(df.format(matrix[i][j]));
                  buffer.append(",");
              }
              buffer.append(df.format(matrix[i][cols - 1]));
              buffer.append("\n");
          }

          
          FileOutput out = new FileOutput(fileName); 
          
          out.writeString(buffer.toString());
          
          out.close();
          
     }

     
     /**
      * Get a prolate axi-symmetric tensor with e1 == x, e2 == y, e3 == z.
      *
      * @param fa the fractional anisotropy.
      * @param trace the trace of the tensor.
      */
     protected DT getTensor(double fa, double trace) {
         
        if (fa < 0.0 || fa > 1.0) {
            throw new IllegalArgumentException("Require 0.0 <= FA <= 1.0");
        }
        
        
        // from Mathematica, solving the FA equation for L1, using L2 = L3 = (Trace - L1) / 2
        
        double t = trace;
        
        double l1 = (-3.0 * t + 2 * fa * fa * t - 2.0 *
                     Math.sqrt(3.0 * fa * fa * t * t - 2 * fa * fa * fa * fa * t * t))
            /(3.0 * (-3.0 + 2.0 * fa * fa));

        double l2 = (t - l1) / 2.0;
        double l3 = l2;

        return new DT(l1, 0.0, 0.0, l2, 0.0, l3);
     }

    /**
     * Empirically determine the aperture (or opening angle) of the cone of uncertainty for some sample vectors, which are treated as axes.
     *
     * @param vecs the sample vectors, assumed to be unit vectors 
     * @param meanAxis the sample mean axis, should also be a unit vector.
     *
     * @return theta, in degrees. 100 * (1-alpha) % of the samples are contained within a cone centered on the mean, with aperture theta.
     * 
     */	
    protected double coneOfUncertainty(Vector3D[] vecs, Vector3D meanAxis, double alpha) {

        int numVecs = vecs.length;

        double[] dots = new double[numVecs];

        for (int i = 0; i < numVecs; i++) {
            dots[i] = Math.abs(vecs[i].dot(meanAxis));
        }

        // Ascending order
        Arrays.sort(dots);

        int alphaInd = (int)(Math.round(numVecs * alpha));

        return 2.0 * 180.0 * Math.acos(dots[alphaInd]) / Math.PI;
    } 
}

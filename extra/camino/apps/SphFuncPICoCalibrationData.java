package apps;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;

import java.util.Random;

/**
 * Generates data for Hessian-based PICo calibration.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class SphFuncPICoCalibrationData {

    // variables here refer to DT1 and DT2 eigenvectors BEFORE any rotations

    // radians 
    public static final double DT1_E1_THETA = 0.0;
    public static final double DT1_E1_PHI = 0.0;

    public static final double DT1_E2_THETA = Math.PI / 2.0;
    public static final double DT1_E2_PHI = 0.0;

    public static final double DT1_E3_THETA = Math.PI / 2.0;
    public static final double DT1_E3_PHI = Math.PI / 2.0;


    // DT1_E1 must be orthogonal to DT1_E1 and DT1_E2
    public static final double DT2_E1_THETA = DT1_E3_THETA;
    public static final double DT2_E1_PHI = DT1_E3_PHI;

    public static final double DT2_E2_THETA = DT1_E2_THETA;
    public static final double DT2_E2_PHI = DT1_E2_PHI;

    public static final double DT2_E3_THETA = DT1_E1_THETA;
    public static final double DT2_E3_PHI = DT1_E1_PHI;


    public static final Vector3D DT1_E1 = Vector3D.vectorFromSPC(1.0,  DT1_E1_THETA, DT1_E1_PHI);
    public static final Vector3D DT1_E2 = Vector3D.vectorFromSPC(1.0,  DT1_E2_THETA, DT1_E2_PHI);
    public static final Vector3D DT1_E3 = Vector3D.vectorFromSPC(1.0,  DT1_E3_THETA, DT1_E3_PHI);
    
    public static final Vector3D DT2_E1 = Vector3D.vectorFromSPC(1.0,  DT2_E1_THETA, DT2_E1_PHI);
    public static final Vector3D DT2_E2 = Vector3D.vectorFromSPC(1.0,  DT2_E2_THETA, DT2_E2_PHI);
    public static final Vector3D DT2_E3 = Vector3D.vectorFromSPC(1.0,  DT2_E3_THETA, DT2_E3_PHI);
    

    private static RealMatrix dt1_e1_e1T = null;
    private static RealMatrix dt1_e2_e2T = null;
    private static RealMatrix dt1_e3_e3T = null;

    private static RealMatrix dt2_e1_e1T = null;
    private static RealMatrix dt2_e2_e2T = null;
    private static RealMatrix dt2_e3_e3T = null;
	

    static {

        RealMatrix dt1_e1m = DT1_E1.toRealMatrix();
        RealMatrix dt1_e2m = DT1_E2.toRealMatrix();
        RealMatrix dt1_e3m = DT1_E3.toRealMatrix();

        RealMatrix dt2_e1m = DT2_E1.toRealMatrix();
        RealMatrix dt2_e2m = DT2_E2.toRealMatrix();
        RealMatrix dt2_e3m = DT2_E3.toRealMatrix();
        
        dt1_e1_e1T = dt1_e1m.product(dt1_e1m.transpose());
        dt1_e2_e2T = dt1_e2m.product(dt1_e2m.transpose());
        dt1_e3_e3T = dt1_e3m.product(dt1_e3m.transpose());

        dt2_e1_e1T = dt2_e1m.product(dt2_e1m.transpose());
        dt2_e2_e2T = dt2_e2m.product(dt2_e2m.transpose());
        dt2_e3_e3T = dt2_e3m.product(dt2_e3m.transpose());

    }

    // This can be changed, but probably won't be
    public static double trace = 2.1E-9;


    /**
     * Output manager
     */
    private static OutputManager om;


    public static void main(String[] args) {


        CL_Initializer.CL_init(args);
        CL_Initializer.initImagingScheme();
        OutputManager.outputDataType = "float";

        // 5E-4 step gives 1700 tensors in range FA = 0.1 - 0.9
        
        double oneDT_FAMin = 0.05;
        double oneDT_FAMax = 0.9;

        double oneDT_FAStep = 2.5E-4; 

        // 1E-2 step gives 540 two-dt voxels in FA range 0.3 - 0.9 with angles between 0 and PI / 4 
        // in angular steps of PI / 36 (5 degrees)

        // mix step defaults to 0.05, giving total of 3240 two-tensor trials in range 0.5 - 0.8
        
        double twoDT_FAMin = 0.3;
        double twoDT_FAMax = 0.9;

        double twoDT_FAStep = 0.01; 

        // dt2RotAngle means the same as in SyntheticData
        double dt2RotAngleMin = 0.0;
        double dt2RotAngleMax = Math.PI / 4.0;
        double dt2RotAngleStep = Math.PI / 36.0;
        
        double dt2MixMin = 0.2;
        double dt2MixMax = 0.8;

        double dt2MixStep = 0.05;
        
        // args from CL

        double snr = CL_Initializer.SNR;
        DW_Scheme imPars = CL_Initializer.imPars;

        // seed for the rotation of each test function
        int rotationIndex = CL_Initializer.rotationIndex;
   
        Random rotationRan = null;

        if (rotationIndex >= 0) {
            rotationRan = new Random(rotationIndex);
        }


        // CL_Intitializer has fixed default seed
        int noiseSeed = 0;//(int)System.currentTimeMillis();

       
        String infoOutputFile = "pico_calibration_info.txt";
            
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-seed")) {
                noiseSeed = Integer.parseInt(args[i+1]); 
            }
            if (args[i].equals("-onedtfarange")) {
                oneDT_FAMin = Double.parseDouble(args[i+1]);
                oneDT_FAMax = Double.parseDouble(args[i+2]);
                CL_Initializer.markAsParsed(i, 3);
            }
            if (args[i].equals("-onedtfastep")) {
                oneDT_FAStep = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-twodtfarange")) {
                twoDT_FAMin = Double.parseDouble(args[i+1]);
                twoDT_FAMax = Double.parseDouble(args[i+2]);
                CL_Initializer.markAsParsed(i, 3);
            }
            if (args[i].equals("-twodtfastep")) {
                twoDT_FAStep = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-twodtanglerange")) {
                dt2RotAngleMin = Double.parseDouble(args[i+1]);
                dt2RotAngleMax = Double.parseDouble(args[i+2]);
                CL_Initializer.markAsParsed(i, 3);
            }
            if (args[i].equals("-twodtanglestep")) {
                dt2RotAngleStep = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-twodtmixmax")) {
                dt2MixMax = Double.parseDouble(args[i+1]);

                dt2MixMin = 1.0 - dt2MixMax;

                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-twodtmixstep")) {
                dt2MixStep = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-trace")) {
                trace = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-infooutputfile")) {
                infoOutputFile = args[i+1];
                CL_Initializer.markAsParsed(i, 2);                
            }

        }

        Random noiseRan = new Random(noiseSeed);

        CL_Initializer.checkParsing(args);

        om = new OutputManager();

        // some quick sanity checks

        if (oneDT_FAMax < oneDT_FAMin) {
            throw new LoggedException("Specifed one-fibre FA range is negative: " + oneDT_FAMin 
                                      + " to " + oneDT_FAMax);
        }

        if (dt2RotAngleMax < dt2RotAngleMin || dt2RotAngleMax > Math.PI / 2.0) {
            throw new 
                LoggedException("Two-fibre rotation range is wrong (did you specify degrees instead of radians?)");
        }

        if (twoDT_FAMax < twoDT_FAMin) {
            throw new LoggedException("Specifed two-fibre FA range is negative: " + twoDT_FAMin 
                                      + " to " + twoDT_FAMax);
        }

        if (dt2MixMax < 0.0 || dt2MixMax > 1.0) {
            throw new LoggedException("Mixing parameter must range between 0.0 and 1.0, specified maximum was" + dt2MixMax);
        }


        int oneDT_Trials = (int)(1 + (oneDT_FAMax - oneDT_FAMin) / oneDT_FAStep);

        double[] oneDT_FA = new double[oneDT_Trials];

        for (int i = 0; i < oneDT_Trials; i++) {
            oneDT_FA[i] = oneDT_FAMin + i * oneDT_FAStep;
        }

        DT[] dts = getOneTensorBlock(oneDT_FA, rotationRan);
        
        for (int i = 0; i < oneDT_Trials; i++) {
            GaussianMixture g = new GaussianMixture(new DT[] {dts[i]}, new double[] {1.0});

            DataSynthesizer synth = new DataSynthesizer(g, imPars, snr, 1, noiseRan);

            om.output(synth.nextVoxel());
        }


        // now do two-tensors

        int twoDT_FASize = (int)(1 + (twoDT_FAMax - twoDT_FAMin) / twoDT_FAStep);

        double[] twoDT_FA = new double[twoDT_FASize];

        for (int i = 0; i < twoDT_FASize; i++) {
            twoDT_FA[i] = twoDT_FAMin + i * twoDT_FAStep;
        }

        // get angle range
        int twoDT_AngleSize = (int)(1 + (dt2RotAngleMax - dt2RotAngleMin) / dt2RotAngleStep);
        
        
        double[] twoDT_Angles = new double[twoDT_AngleSize];

        for (int i = 0; i < twoDT_AngleSize; i++) {
            twoDT_Angles[i] = dt2RotAngleMin + i * dt2RotAngleStep;
        }


        // finally mix
        int dt2MixSize = (int)(1 + (dt2MixMax - dt2MixMin) / dt2MixStep);
                
        double[] dt2Mix = new double[dt2MixSize];

        for (int i = 0; i < dt2MixSize; i++) {
            dt2Mix[i] = dt2MixMin + i * dt2MixStep;
        }


        int twoDT_BlockSize = (twoDT_FASize * (twoDT_FASize + 1) / 2) * dt2MixSize;
       
        for (int a = 0; a < twoDT_AngleSize; a++) {
            
            DT[][] block = getTwoTensorBlock(twoDT_FA, dt2Mix, twoDT_Angles[a], rotationRan);

            for (int i = 0; i < twoDT_FASize; i++) {
                for (int j = 0; j <= i; j++) {
                    for (int k = 0; k < dt2MixSize; k++) { 

                        int index = dt2MixSize * (i * (i+1) / 2) + dt2MixSize * j + k;

                        GaussianMixture g = new GaussianMixture(block[index], 
                                                                new double[] {1.0 - dt2Mix[k], dt2Mix[k]});

                        DataSynthesizer synth = new DataSynthesizer(g, imPars, snr, 1, noiseRan);

                        om.output(synth.nextVoxel());
                    }
                }
            }
        }
        
        
        // write the info we need for later stages
        FileOutput out = new FileOutput(infoOutputFile);

	out.writeString("VERSION\t1.0\n");
        out.writeString("ROTATION_SEED\t" + rotationIndex + "\n");
        out.writeString("ONE_DT_BLOCK_SIZE\t" + oneDT_Trials + "\n");
        out.writeString("MIN_DT2_ROT_ANGLE\t" + dt2RotAngleMin + "\n");
        out.writeString("MAX_DT2_ROT_ANGLE\t" + dt2RotAngleMax + "\n");
        out.writeString("DT2_ROT_ANGLE_STEP\t" + dt2RotAngleStep + "\n");
        out.writeString("TWO_DT_BLOCK_SIZE\t" + twoDT_BlockSize + "\n");
        out.writeString("NUM_TWO_DT_BLOCKS\t" + twoDT_AngleSize + "\n");
        out.writeString("DT1_E1_THETA\t" + DT1_E1_THETA + "\n");
        out.writeString("DT1_E1_PHI\t" + DT1_E1_PHI + "\n");
        out.writeString("DT2_E1_THETA\t" + DT2_E1_THETA + "\n");
        out.writeString("DT2_E1_PHI\t" + DT2_E1_PHI + "\n");
        out.writeString("DT2_ROT_AXIS_THETA\t" + DT1_E2_THETA + "\n");
        out.writeString("DT2_ROT_AXIS_PHI\t" + DT1_E2_PHI + "\n");
    
        out.close();

        om.close();
    }



    /**
     * @param rotationRan Random object to use to construct a random rotation matrix for each tensor.
     * If this is <code>null</code>, then no rotation is applied.
     * @return an array of tensors <code>tensors</code>, where <code>tensors[i]</code> has anisotropy
     * <code>fa[i]</code>. Each tensor is reoriented to a random rotation.
     */
    public static DT[] getOneTensorBlock(double[] fa, 
                                         Random rotationRan) {


        int numTensors = fa.length;

        DT[] tensors = new DT[numTensors];

        for (int i = 0; i < numTensors; i++) {

            double[] ev = eigenvalues(fa[i], trace);

            RealMatrix d1 = dt1_e1_e1T.scalarMult(ev[0]);
            RealMatrix d2 = dt1_e2_e2T.scalarMult(ev[1]);
            RealMatrix d3 = dt1_e3_e3T.scalarMult(ev[2]);
	
            double[][] a = (d1.add(d2).add(d3)).entries;
            
            tensors[i] = new DT(a[0][0], a[0][1], a[0][2], a[1][1], a[1][2], a[2][2]);

            if (rotationRan != null) {
                RealMatrix rotMat = Rotations.randomRotMat(rotationRan);
                
                tensors[i] = tensors[i].transform(rotMat);
            }
 
        }

        return tensors;

    }



    /**
     * Compute two tensors over a range of anisotropy, with a fixed crossing angle. 
     * Computes all combinations of <code>fa</code> and <code>mix</code>.
     * 
     *
     * @param rotationRan Random object to use to construct a random rotation matrix for each pair of
     * tensors. If this is <code>null</code>, then no rotation is applied.
     *
     * @return array of diffusion tensors, with dimensions <code>[number of combinations][2]</code>.
     * 
     */
    public static DT[][] getTwoTensorBlock(double[] fa, double[] mix,
                                           double dt2RotAngle, Random rotationRan) {

        RealMatrix dt2Rot = Rotations.getRotMat(DT1_E2, dt2RotAngle);
        RealMatrix dt2RotT = dt2Rot.transpose();

        //     rotate dt2 by dt2rotangle
        
        RealMatrix dt2_e1_e1T_rot = dt2RotT.product(dt2_e1_e1T).product(dt2Rot);
        RealMatrix dt2_e2_e2T_rot = dt2RotT.product(dt2_e2_e2T).product(dt2Rot);
        RealMatrix dt2_e3_e3T_rot = dt2RotT.product(dt2_e3_e3T).product(dt2Rot);
  
        DT[][] tensors = new DT[(fa.length * (fa.length + 1) / 2) * mix.length][2];

        double[][] evs = new double[fa.length][];
        
        for (int i = 0; i < fa.length; i++) {
            evs[i] = eigenvalues(fa[i], trace);
        }
        

        for (int i = 0; i < fa.length; i++) {

            RealMatrix d1_1 = dt1_e1_e1T.scalarMult(evs[i][0]);
            RealMatrix d2_1 = dt1_e2_e2T.scalarMult(evs[i][1]);
            RealMatrix d3_1 = dt1_e3_e3T.scalarMult(evs[i][2]);
            
            double[][] a1 = d1_1.add(d2_1).add(d3_1).entries;

            DT dt1 = new DT(a1[0][0], a1[0][1], a1[0][2], a1[1][1], a1[1][2], a1[2][2]);

            for (int j = 0; j <= i; j++) {
            
                RealMatrix d1_2 = dt2_e1_e1T_rot.scalarMult(evs[j][0]);
                RealMatrix d2_2 = dt2_e2_e2T_rot.scalarMult(evs[j][1]);
                RealMatrix d3_2 = dt2_e3_e3T_rot.scalarMult(evs[j][2]);

                double[][] a2 = d1_2.add(d2_2).add(d3_2).entries;
                 
                DT dt2 = new DT(a2[0][0], a2[0][1], a2[0][2], a2[1][1], a2[1][2], a2[2][2]);

                for (int k = 0; k < mix.length; k++) {
                    
                    int index = mix.length * (i * (i+1) / 2) + mix.length * j + k;

                    tensors[index][0] = dt1;                    
                    tensors[index][1] = dt2;

                    if (rotationRan != null ) {
                        RealMatrix rotMat = Rotations.randomRotMat(rotationRan);
                        tensors[index][0] = tensors[index][0].transform(rotMat);
                        tensors[index][1] = tensors[index][1].transform(rotMat);
                        
                    }
                }
            }
        }

        return tensors;
    }



    /**
     * Calculates eigenvalues of prolate DT with specified FA and trace.
     *
     */
    public static double[] eigenvalues(double fa, double t) {

        if (fa == 0.0) {
            return new double[] {t / 3.0, t / 3.0, t / 3.0};
        }
        else if (fa > 0.0 && fa < 1.0) {
            
            double l1 = (-3.0 * t + 2 * fa * fa * t - 2.0
                         * Math.sqrt(3.0 * fa * fa * t * t - 2 * fa * fa * fa * fa * t * t))
                /(3.0 * (-3.0 + 2.0 * fa * fa));
            
            double l2 = (t - l1) / 2.0;
            
            return new double[] {l1, l2, l2};
            
        }
        else { // catches NaN
            throw new IllegalArgumentException("Require 0 <= FA < 1");
        }
    }

}

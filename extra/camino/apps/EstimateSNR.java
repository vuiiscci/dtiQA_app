package apps;

import data.*;
import imaging.*;
import misc.*;
import tools.*;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;



/**
 *
 * Estimates SNR using different methods explained in the paper "Measurement of signal-to-noise ratios in MR
 * images: Influence of multichannel coils, parallel imaging, and reconstruction filters" by Olaf Dietrich
 * et al, JMRI 26:376-385 (2007).
 * <p>
 * The standard method is to divide the mean of a tissue ROI by the standard deviation of a background ROI, but
 * this can be unreliable with modern imaging techniques. We therefore examine the signal variance across multiple
 * b=0 images to estimate the SNR.
 * </p>
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class EstimateSNR extends Executable{
    
    //public static void main(String[] args) {
    public EstimateSNR(String[] args){
        super(args);
    }
    
    DataSource bgROI_Source;
    int[] b0Indices;
    DW_Scheme scheme;
	int numMeas;
	
    public void initDefaultVals(){
        bgROI_Source = null;
        b0Indices = null;
        
    }
    
    public void initOptions(String[] args){
              
        CL_Initializer.inputDataType = "float";
        CL_Initializer.CL_init(args);
        
        CL_Initializer.initImagingScheme();
        CL_Initializer.initDataSynthesizer();
		
		for (int i = 0; i < args.length; i++) {
            
            if(args[i].equals("-noiseroi")) {
                try {
                    bgROI_Source = ImageHeader.readHeader(args[i+1]).getImageDataSource();
                }
                catch (IOException e) {
                    throw new LoggedException("Cannot read background ROI " + args[i+1]);
                }
                
                CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-approxb0meas")) {
                String[] tmp = new String[1000];
                int t = 0;
                
                while (t + i + 1 < args.length && !args[i+t+1].startsWith("-")) {
                    tmp[t] = args[i+t+1];
                    t++;
                }
                
                b0Indices = new int[t];
                
                for (int index = 0; index < t; index++) {
                    b0Indices[index] = Integer.parseInt(tmp[index]);
                }
                
                CL_Initializer.markAsParsed(i, t+1);
            }
            
        }
        
        CL_Initializer.checkParsing(args);
        
        if (CL_Initializer.bgMask == null) {
            throw new LoggedException("This program requires a foreground ROI with -bgmask");
        }
        
    }
        
	public void initVariables(){	
        scheme = CL_Initializer.imPars;
		numMeas = scheme.numMeasurements();
    }   
	
    public void execute(OutputManager om){
        // flat list of foreground data in voxel order, ie all measurements from the
        // same voxel are together
        ArrayList<Double> fg = new ArrayList<Double>();
        
        ArrayList<Double> bg = new ArrayList<Double>();
        
//        int numMeas = scheme.numMeasurements();
        
        boolean[] isB0 = new boolean[numMeas];
        
        int numB0 = 0;
        
        if (b0Indices != null) {
            
            // override scheme file and treat specified indices as b=0
            
            numB0 = b0Indices.length;
            
            for (int i = 0; i < numB0; i++) {
                isB0[b0Indices[i]] = true;
            }
        }
        else {
            for (int i = 0; i < numMeas; i++) {
                isB0[i] = scheme.zero(i);
                
                if (isB0[i]) {
                    numB0++;
                }
            }
            
        }
        if (numB0 == 0) {
            throw new LoggedException("Can't calculate SNR without any data at b=0");
        }
        
        
        // number of voxels in the background ROI
        int bgVoxels = 0;
        
        int fgVoxels = 0;
        
        while (CL_Initializer.data.more()) {
            double[] voxel = CL_Initializer.data.nextVoxel();
            
            // check if voxel is in foreground ROI
            if (CL_Initializer.bgMask.nextVoxel()[0] > 0.0) {
                
                double[] b0 = new double[numB0];
                
                fgVoxels++;
                
                for (int i = 0; i < numMeas; i++) {
                    if (isB0[i]) {
                        fg.add(new Double(voxel[i]));
                    }
                }
                
            }
            
            // check if voxel is in the background ROI
            if (bgROI_Source != null && bgROI_Source.nextVoxel()[0] > 0.0) {
                // use first measurement for bg
                bg.add(new Double(voxel[0]));
                bgVoxels++;
            }
            
        }
        
        
        // now got data, convert to arrays
        double[][] data = new double[fgVoxels][numB0];
        
        Double[] tmp = new Double[fgVoxels*numB0];
        
        fg.toArray(tmp);
        
        int nextFG = 0;
        
        for (int i = 0; i < fgVoxels; i++) {
            for (int j = 0; j < numB0; j++) {
                data[i][j] = tmp[nextFG++].doubleValue();
            }
        }
        
        double[] bgData = null;
        
        if (bgROI_Source != null) {
            
            tmp = new Double[bgVoxels];
            
            bg.toArray(tmp);
            
            bgData = new double[bgVoxels];
            
            for (int i = 0; i < bgVoxels; i++) {
                bgData[i] = tmp[i].doubleValue();
            }
        }
        
        
        System.out.println("\nNumber of b=0 images:\t" + numB0 + "\n");
        
        System.out.println("\nNumber of foreground voxels:\t" + fgVoxels);
        System.out.println("\nNumber of background voxels:\t" + bgVoxels + "\n");
        
        System.out.println("\nSignal to noise in foreground ROI\n");
        
        double[] snrStdv = null;
        double[] snrDiff = null;
        double[] snrMult = null;
        
        if (bgData != null) {
            
            snrStdv = snrStdv(data, bgData);
        }
        if (numB0 > 1) {
            snrDiff = snrDiff(data);
            snrMult = snrMult(data);
        }
        
        DecimalFormat df = new DecimalFormat("0.00");
        
        if (snrStdv != null) {
            System.out.println("SNR stdv:\t\t" + df.format(snrStdv[0] / snrStdv[1]));
        }
        
        if (numB0 == 1) {
            System.out.println("SNR mult:\t\t undefined");
            System.out.println("SNR diff:\t\t undefined");
        }
        else {
            
            System.out.println("SNR diff:\t\t" + df.format(snrDiff[0] / snrDiff[1]));
            System.out.println("SNR mult:\t\t" + df.format(snrMult[0] / snrMult[1]));
        }
        
        System.out.println("\nNoise standard deviation\n");
        
        if (snrStdv != null) {
            System.out.println("sigma stdv:\t\t" + df.format(snrStdv[1]));
        }
        
        if (numB0 == 1) {
            System.out.println("sigma mult:\t\t undefined");
            System.out.println("sigma diff:\t\t undefined");
        }
        else {
            System.out.println("sigma diff:\t\t" + df.format(snrDiff[1]));
            System.out.println("sigma mult:\t\t" + df.format(snrMult[1]));
        }
        
        System.out.println();
    }
    
    
    
    /**
     * Computes SNR from a single ROI, imaged twice.
     *
     * @param data the data, dimensions [voxels][2].
     *
     * @return {S_{diff}, \sigma_{diff}}.
     */
    public static double[] snrDiff(double[][] data) {
        
        int N = data.length;
        
        // sMean = (1 / 2N) (data[i][0] + data[i][1]), 0 <= i < N, which is S_diff in the Dietrich paper
        double sMean = 0.0;
        
        // sDiff[i] = (data[i][0] - data[i][1])
        double[] sDiff = new double[N];
        
        double sigmaDiff = 0.0;
        
        for (int i = 0; i < N; i++) {
            // actually 2 * N * sMean here
            sMean += data[i][0] + data[i][1];
            
            sDiff[i] = data[i][0] - data[i][1];
        }
        
        sMean = sMean / (2.0 * N);
        
        sigmaDiff = ( 1.0 / Math.sqrt(2.0) ) * Math.sqrt(ArrayOps.var(sDiff, ArrayOps.mean(sDiff)));
        
        return new double[] {sMean, sigmaDiff};
    }
    
    
    /**
     * Computes SNR from a single ROI, imaged many times.
     *
     * @param data the data, dimensions [voxels][K] where there are K measurements.
     *
     * @return {S_{mult}, \sigma_{mult}}.
     */
    public static double[] snrMult(double[][] data) {
        
        // number of voxels
        int N = data.length;
        
        // number of measurements
        int K = data[0].length;
        
        // mean signal over all measurements and all voxels
        double sMult = 0.0;
        
        double[] std = new double[N];
        
        for (int i = 0; i < N; i++) {
            
            // normalize later
            sMult += ArrayOps.sum(data[i]);
            
            std[i] = Math.sqrt(ArrayOps.var(data[i], ArrayOps.mean(data[i])));
        }
        
        sMult = sMult / (N * K);
        
        return new double[] {sMult, ArrayOps.mean(std)};
        
    }
    
    
    /**
     * Computes SNR from a single ROI with noise variance estimated from background.
     *
     * @param fg the data for the foreground ROI, with dimension [voxels][1].
     * @param bg the data for the background ROI.
     *
     * @return {S, \sigma}.
     */
    public static double[] snrStdv(double[][] fg, double[] bg) {
        
        double[] fg1D = new double[fg.length];
        
        for (int i = 0; i < fg.length; i++) {
            fg1D[i] = fg[i][0];
        }
        
        double fgMean = ArrayOps.mean(fg1D);
        
        double bgStd = Math.sqrt(2.0 / (4.0 - Math.PI)) * Math.sqrt(ArrayOps.var(bg, ArrayOps.mean(bg)));
        
        return new double[] {fgMean, bgStd};
    }
    
    
    
    
}
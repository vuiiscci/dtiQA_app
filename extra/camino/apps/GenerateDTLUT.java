package apps;

import data.*;
import imaging.*;
import inverters.ModelIndex;
import numerics.*;
import tools.*;

import tractography.*;

import java.io.*;
import java.util.Random;

import misc.LoggedException;

/**
 * Makes lookup tables for PICo from diffusion tensor data.
 *
 * @author Philip Cook
 * @version $Id$
 */
public final class GenerateDTLUT extends Executable{
    
    public GenerateDTLUT(String[] args){
        super(args);
    }
    
    Random ran;
    
    private double xMin;
    private double xMax;
    private double step;
    private int samples; // sample this many tensors
    private double prop;
    private double cross;
    private int crossDeg;
    private boolean watson;
    private boolean bingham;
    private boolean acg;
    private boolean twoTensor;
    // random seed for adding noise and defining tensors
    private long seed;
    private boolean twoD;
	private double trace;
	
	private String outputRoot;
	private double snr;
	private	ModelIndex inversionIndex;
    private DW_Scheme imPars;
	
    public void initDefaultVals(){
        xMin = 0.0;
        xMax = 0.0;
        step = 0.0;
        samples = 2000;
        
        seed = System.currentTimeMillis();
        prop = 0.5;
        cross = 0.0;
        crossDeg = 0;
        
        watson = false;
        bingham = false;
        acg = false;
        
        twoTensor = false;
		
		twoD = true;		
		
		outputRoot = "";
		
		snr = 0;
		inversionIndex = null;
        imPars = null;

	trace = 2.1E-9;
    }
    
    // public static void main(String[] args) {
    
/*	Random ran;
 
    double xMin = 0.0;
    double xMax = 0.0;
 
    double step = 0.0;
 
    // sample this many tensors
    int samples = 2000;
 
    // random seed for adding noise and defining tensors
    long seed = System.currentTimeMillis();
 
    double prop = 0.5;
    double cross = 0.0;
    int crossDeg = 0;
 
    boolean watson = false;
    boolean bingham = false;
    boolean acg = false;
 
    boolean twoTensor = false;
 */
    public void initOptions(String[] args){
        if (args.length == 0) {
            System.exit(0);
        }
        CL_Initializer.CL_init(args);
        CL_Initializer.initImagingScheme();    
		
        for (int i = 0; i < args.length; i++) {
            
            if (args[i].equals("-lrange")) {
                xMin = Double.parseDouble(args[i+1]);
                xMax = Double.parseDouble(args[i+2]);
                CL_Initializer.markAsParsed(i,3);
            }
            else if (args[i].equals("-frange")) {
                xMin = Double.parseDouble(args[i+1]);
                xMax = Double.parseDouble(args[i+2]);
                CL_Initializer.markAsParsed(i,3);
            }
            else if (args[i].equals("-step")) {
                step = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-trace")) {
                trace = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-seed")) {
                seed = Long.parseLong(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-samples")) {
                samples = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-watson")) {
                watson = true;
                CL_Initializer.markAsParsed(i);
            }
            else if (args[i].equals("-bingham")) {
                bingham = true;
                CL_Initializer.markAsParsed(i);
            }
            else if (args[i].equals("-acg")) {
                acg = true;
                CL_Initializer.markAsParsed(i);
            }
            else if (args[i].equals("-prop")) {
                prop = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-cross")) {
                cross = (Math.PI / 180.0) * Double.parseDouble(args[i+1]);
                crossDeg = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-outputroot")) {
                outputRoot = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
                        
        }
		CL_Initializer.checkParsing(args);
	}
		
		
    public void initVariables() {
            
		inversionIndex = CL_Initializer.inversionIndices[0];
        imPars = CL_Initializer.imPars;
		snr = CL_Initializer.SNR;
		/*ModelIndex inversionIndex = CL_Initializer.inversionIndices[0];
        DW_Scheme imPars = CL_Initializer.imPars;
		double snr = CL_Initializer.SNR;*/
		
    }
       
       
	public void execute(OutputManager om){
        if (snr == -1) {
            // no SNR specified
            throw new LoggedException("Cannot generate LUT without valid SNR");
        }
        
        if (cross > 0.0 || inversionIndex.numDTs == 2) {
            twoTensor = true;
            
            if (cross == 0.0) {
                throw new LoggedException("Attempting to create two-tensor LUT with zero degree crossing angle");
            }
            else if (inversionIndex.numDTs < 2) {
                throw new LoggedException("Attempting to create two-tensor LUT with single tensor inversion.");
            }
        }
        
        
        // compute default range and step if none of these have been given
        if (xMin == 0.0 && xMax == 0.0) {
            if (twoTensor) {
                xMin = 0.3;
                xMax = 0.94;
            }
            else {
                xMin = 1.0;
                xMax = 10.0;
            }
        }
        if (step == 0.0) {
            if (xMax < 1.0) {
                // LUT indexed by FA
                step = 0.02;
            }
            else {
                step = 0.2;
            }
        }
        
        if (xMax < 1.0 && inversionIndex.numDTs == 1) {
            // fRange used and one Tensor inversion == 1D LUT with cyl. symm tensors
            twoD = false;
        }
        
        // Now do some checks for user mistakes
        
        // no lut selected
        if (!(watson || bingham || acg)) {
            bingham = true;
        }
        
        // enforce that range and step are sensible
        
        if ( (xMin - xMax) / step - (int)((xMin - xMax) / step) > 1E-8) {
            throw new LoggedException("LUT range must divide into a whole number of steps");
        }
        if (xMin > xMax) {
            throw new LoggedException("Illegal range specified: max < min");
        }        
        if (twoTensor) {
            if (xMin < 0.0 || xMax > 1.0) {
                throw new LoggedException("FA range must be between 0 and 1");
            }
            
        }
        else {
            
            if (!twoD) {
                if (xMin < 0.0 || xMax > 1.0) {
                    throw new LoggedException("FA range must be between 0 and 1");
                }
            }
            else {
                if (xMin < 1.0) {
                    throw new LoggedException("min(L1 / L3) must not be less than 1");
                }
            }
        }
        
        
        ran = new MTRandom(seed);
        
        DT_LookupTableGenerator generator = null;
        
        if (twoTensor) {
            generator = new TwoTensorLUTGenerator(imPars, snr, trace, prop, cross, ran);
        }
        else {
            generator = new OneTensorLUTGenerator(imPars, snr, trace, ran);
        }
        
        double[][][][] lut = generator.generateLUT(xMin, xMax, step, samples, inversionIndex,
        watson, bingham, acg);
        
        String outputSuffix = "_snr" + (int)snr + "_" + inversionIndex;
        
        if (crossDeg > 0) {
            outputSuffix += "_" + crossDeg;
        }
        
        // if only one LUT is selected, output to CL_Initializer.out
        // else to files

        boolean onePDFSelected = !(watson && bingham) && ( (watson || bingham) ^ acg );
        
        try {
            double[] dims = new double[] {xMin, xMax, xMin, xMax};
            double[] steps = new double[] {step, step};
            
            
            if (!twoD) {
                dims = new double[] {0.0, 0.0, xMin, xMax};
                steps = new double[] {0.0, step};
            }
            
            
            if (onePDFSelected) {
                
                // set up output
                //OutputManager om = new OutputManager();
                
                
                DataOutputStream dot = om.getOutputStream();
                
                if (watson) {
                    
                    writeBinaryOutput(lut[DT_LookupTableGenerator.WATSON], dims, steps, dot);
                    
                }
                if (bingham) {
                    
                    writeBinaryOutput(lut[DT_LookupTableGenerator.BINGHAM], dims, steps, dot);
                    
                }
                if (acg) {
                    
                    writeBinaryOutput(lut[DT_LookupTableGenerator.ACG], dims, steps, dot);
                    
                }
                
                //om.close();
                
            }
            else {
                
                if (watson) {
                    
                    FileOutputStream fout = new FileOutputStream(outputRoot + "watson" + outputSuffix);
                    DataOutputStream dot = new DataOutputStream(new BufferedOutputStream(fout, 1024 * 1024));
                    
                    writeBinaryOutput(lut[DT_LookupTableGenerator.WATSON], dims, steps, dot);
                    dot.close();
                }
                if (bingham) {
                    
                    FileOutputStream fout = new FileOutputStream(outputRoot + "bingham" + outputSuffix);
                    DataOutputStream dot = new DataOutputStream(new BufferedOutputStream(fout, 1024 * 1024));
                    
                    writeBinaryOutput(lut[DT_LookupTableGenerator.BINGHAM], dims, steps, dot);
                    dot.close();
                }
                if (acg) {
                    
                    FileOutputStream fout = new FileOutputStream(outputRoot + "acg" + outputSuffix);
                    DataOutputStream dot = new DataOutputStream(new BufferedOutputStream(fout, 1024 * 1024));
                    
                    writeBinaryOutput(lut[DT_LookupTableGenerator.ACG], dims, steps, dot);
                    dot.close();
                }
                
            }
            
        }
        catch(IOException e) {
            throw new LoggedException("I/O error while writing LUT");
        }
        

	}
	    
    // writes scanner order
    private void writeBinaryOutput(double[][][] data, double[] range, double[] step,
    DataOutputStream dot) {
        try {
            
            int values = data[0][0].length;
            
            // header == xrange yrange zrange xstep ystep zstep valuesPerPosition
            
            dot.writeDouble(0.0);
            dot.writeDouble(0.0);
            dot.writeDouble(range[0]);
            dot.writeDouble(range[1]);
            dot.writeDouble(range[2]);
            dot.writeDouble(range[3]);
            
            dot.writeDouble(0.0);
            dot.writeDouble(step[0]);
            dot.writeDouble(step[1]);
            
            dot.writeDouble(values);
            
            for (int v = 0; v < values; v++) {
                for (int j = 0; j < data[0].length; j++) {
                    for (int i = 0; i < data.length; i++) {
                        dot.writeDouble(data[i][j][v]);
                    }
                }
            }
            
            dot.close();
            
        }
        catch(IOException e) {
            throw new LoggedException("I/O error while writing LUT");
        }
        
    }
    
}











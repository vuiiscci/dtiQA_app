package apps;

import data.*;
import misc.*;
import tools.*;

import java.io.*;
import java.util.Random;


/**
 * Orders electrostatic points according to the method in Cook et al (ISMRM 2006, p. 1035).
 * The application will work with any point set; its energy function is based on electrostatics.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class OrderElectrostaticPoints extends Executable{
    
    private int numPoints;// = 0;
    
    // this should only be changed by the arg -temperature or by the file if resuming a run
    private double temperature;//= -1.0;
    
    private double coolingFactor;// = 1E-7;
    private int iterations;// = 1000;
    
    private boolean initFromFile;// = false;
    private boolean randomInit;// = false;
    
    private double[][] points;// = null;
    
    private String saveStateRoot;// = "OrderElectrostaticPoints";
    
	private long seed;// = System.currentTimeMillis();
	
    public OrderElectrostaticPoints(String[] args){
        super(args);
    }
    
    public void initDefaultVals(){
        numPoints = 0;
        
        // this should only be changed by the arg -temperature or by the file if resuming a run
        temperature = -1.0;
        
        coolingFactor = 1E-7;
        iterations = 1000;
        
        initFromFile = false;
        randomInit = false;
        
        double[][] points = null;
        
        saveStateRoot = "OrderElectrostaticPoints";
		
		seed = System.currentTimeMillis();
    }
    
    public void initOptions(String[] args){
        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            
            if (args[i].equals("-numpoints")) {
                numPoints = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-seed")) {
                seed = Long.parseLong(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-coolingfactor")) {
                coolingFactor = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-trialsbetweencooling")) {
                iterations = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-temperature")) {
                temperature = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-savestate")) {
                saveStateRoot = args[i+1];
                CL_Initializer.markAsParsed(i, 2);
            }
            if(args[i].equals("-resume")) {
                initFromFile = true;
                
                // read in format saved by OrderedAcqWeightedMinimizer
                
                FileInput in = new FileInput(args[i+1]);
                String[] line = in.readString().split("\\s");
                
                numPoints = Integer.parseInt(line[0]);
                
                // reset temp from file, unless it is specified as an argument
                if (temperature == -1.0) {
                    temperature = Double.parseDouble(line[1]);
                }
                
                points = new double[numPoints][3];
                
                for (int j = 0; j < numPoints; j++) {
                    String[] point = in.readString().split("\\s");
                    
                    points[j][0] = Double.parseDouble(point[0]);
                    points[j][1] = Double.parseDouble(point[1]);
                    points[j][2] = Double.parseDouble(point[2]);
                }
                
                in.close();
                
                CL_Initializer.markAsParsed(i,2);
                
            }
            if (args[i].equals("-randominit")) {
                randomInit = true;
                CL_Initializer.markAsParsed(i);
            }
            
        }
        CL_Initializer.checkParsing(args);
        
    }
    
    public void initVariables(){
    }
    
    public void execute(OutputManager om){
        
        
//    public static void main(String[] args) {
        
        // read from stdin or inputFile if no points specified by other means
        if (numPoints == 0) {
            initFromFile = true;
            
            try {
                points =
                PointSetToScheme.normalizePoints(PointSetToScheme.readPoints(CL_Initializer.inputFile));
                
                numPoints = points.length;
            }
            catch (IOException e) {
                throw new LoggedException(e);
            }
            
        }
        
        Random ran = new Random(seed);
        
        OrderedAcqWeightedMinimizer c = null;
        
        if (initFromFile) {
            c = new OrderedAcqWeightedMinimizer(points, coolingFactor, iterations, saveStateRoot, ran);
            
            if (randomInit) {
                c.initializeRandom();
            }
            
        }
        else {
            c = new OrderedAcqWeightedMinimizer(numPoints, coolingFactor, iterations, saveStateRoot, ran);
        }
        
        if (temperature < 0.0) {
            // need to either calibrate or set a temperature
            c.setTemperature(c.calibratedTemp(0.7));
        }
        else {
            c.setTemperature(temperature);
        }
        
        c.minimize();
        
        if (OutputManager.outputFile == null) {
            System.out.print(c.lowestEnergyState());
        }
        else {
            
            FileOutput out = new FileOutput(OutputManager.outputFile);
            
            out.writeString(c.lowestEnergyState());
            
            out.close();
        }
    }
    
}



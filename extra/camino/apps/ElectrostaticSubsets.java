package apps;

import data.*;
import misc.*;
import tools.*;

import java.io.*;
import java.util.Random;
import java.util.logging.Logger;

/**
 * Orders electrostatic points into subsets following the method in Cook et al (ISMRM 2005, p. 1305).
 * The application will work with any point set; its energy function is based on electrostatics.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class ElectrostaticSubsets extends Executable{
    
    private static Logger logger = Logger.getLogger("camino.apps.ElectrostaticSubsets");
    
    public ElectrostaticSubsets(String[] args){
        super(args);
    }
    
	private int numPoints;
	private int numSubsets;
	private int[] pairsPerSubset;
	
	// this should only be changed by the arg -temperature or by the file if resuming a run
    private double temperature;
    
    private double coolingFactor;
    private int iterations;    
    private boolean initFromFile;  
    private boolean randomInit;    
    private double[][] points;    
    private String saveStateRoot;   
    // if true, we optimize a single subset and ignore the other points
    private boolean singleSubset;
    private int singleSubsetSize;
	private long seed;
	private OrderedAcqMinimizer minimizer;
	private Random ran;
	
    public void initDefaultVals(){
		numPoints = 0;
		numSubsets = 0;
    
		pairsPerSubset = null;
	// this should only be changed by the arg -temperature or by the file if resuming a run
		temperature = -1.0;
    
		coolingFactor = 1E-7;
		iterations = 1000;
    
		initFromFile = false;
   		randomInit = false;
		points = null;
		saveStateRoot = "ElectrostaticSubsets";
    
    // if true, we optimize a single subset and ignore the other points
		singleSubset = false;
		singleSubsetSize = 0;

 		seed = System.currentTimeMillis();
		
		minimizer = null;
	
    }
    
    public void initOptions(String[] args){
        //public static void main(String[] args) {
        
        if (args.length == 0) {
            // this has multiple wrappers so no usage
            // should attempt to read from inputfile if numpoints is missing
            System.exit(0);
        }
        
        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            
            if (args[i].equals("-numpoints")) {
                numPoints = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-singlesubset")) {
                singleSubset = true;
                singleSubsetSize = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-numsubsets")) {
                numSubsets = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-pointspersubset")) {
                
                int[] tmpPairsPerSubset = new int[1000];
                
                int c = 0;
                
                while (i+c+1 < args.length && !args[i+1+c].startsWith("-")) {
                    tmpPairsPerSubset[c] = Integer.parseInt(args[i+c+1]);
                    c++;
                }
                
                pairsPerSubset = new int[c];
                numSubsets = c;
                
                for (int ppc = 0; ppc < pairsPerSubset.length; ppc++) {
                    pairsPerSubset[ppc] = tmpPairsPerSubset[ppc];
                }
                
                CL_Initializer.markAsParsed(i, numSubsets + 1);
                
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
                
                // read in format saved by OrderedAcqSubsetMinimizer
                
                FileInput in = new FileInput(args[i+1]);
                String[] line = in.readString().split("\\s");
                
                numPoints = Integer.parseInt(line[0]);
                
                numSubsets = Integer.parseInt(line[1]);
                
                pairsPerSubset = new int[numSubsets];
                
                for (int c = 0; c < numSubsets; c++) {
                    pairsPerSubset[c] = Integer.parseInt(line[c+2]);
                }
                
                // reset temp from file, unless it is specified as an argument
                if (temperature == -1.0) {
                    temperature = Double.parseDouble(line[2+numSubsets]);
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
        
		if (singleSubset && numSubsets > 0) {
			logger.warning("Single subset specified, but user also passed -numsubsets or -pointspersubset args.");
		}
    
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
    
    	ran = new Random(seed);
    }
    
	public void execute(OutputManager om){	
	System.err.println(seed);
    if (singleSubset) {       
        if (initFromFile) {			
            minimizer = new OrderedAcqSingleSubsetMinimizer(points, singleSubsetSize, coolingFactor,
            iterations, saveStateRoot, ran);
            if (randomInit) {
                minimizer.initializeRandom();
            }
        }
        else {
            minimizer = new OrderedAcqSingleSubsetMinimizer(numPoints, singleSubsetSize,
            coolingFactor, iterations,
            saveStateRoot, ran);
        }
        
    }
    else {
        if (pairsPerSubset == null) {
            if (numPoints % numSubsets != 0) {
                throw new IllegalArgumentException("Can't divide " + numPoints +
                " pairs equally into " + numSubsets + "subsets");
            }
            
            pairsPerSubset = new int[numSubsets];
            
            for (int c = 0; c < numSubsets; c++) {
                pairsPerSubset[c] = numPoints / numSubsets;
            }            
        }      
        
        // make sure we have the right number of points
        
        int expectedPoints = 0;
        
        for (int c = 0; c < numSubsets; c++) {
            expectedPoints += pairsPerSubset[c];
        }
        
        if (expectedPoints != numPoints) {
            throw new LoggedException("Specified subsets require " + expectedPoints +
            " points but there are " + numPoints + " points");
        }
        
        if (initFromFile) {
            minimizer = new OrderedAcqSubsetMinimizer(points, pairsPerSubset, coolingFactor,
            iterations, saveStateRoot, ran);
            
            if (randomInit) {
                minimizer.initializeRandom();
            }
        }
        else {
            minimizer = new OrderedAcqSubsetMinimizer(numPoints, pairsPerSubset, coolingFactor,
            iterations, saveStateRoot, ran);
        }
        
    }
    if (temperature < 0.0) {
        // need to either calibrate or set a temperature
        minimizer.setTemperature(minimizer.calibratedTemp(0.7));
    }
    else {
        minimizer.setTemperature(temperature);
    }
    
    minimizer.minimize();
    
    if (OutputManager.outputFile == null) {
        System.out.print(minimizer.lowestEnergyState());
    }
    else {        
        FileOutput out = new FileOutput(OutputManager.outputFile);        
        out.writeString(minimizer.lowestEnergyState());        
        out.close();
    }   
}
}
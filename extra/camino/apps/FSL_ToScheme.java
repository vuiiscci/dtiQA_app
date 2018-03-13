package apps;

import data.OutputManager;
import misc.*;
import tools.*;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;
import java.text.*;


 /**
  * Converts FSL b-vectors and b-values to Camino scheme file format.
  * 
  * @author Philip Cook
  * @version $Id$
  */
public class FSL_ToScheme extends Executable {


    /**
     * logging object
     */
    private static final Logger logger = Logger.getLogger("apps.FSL_ToScheme");

    private String bvecfile;
    private String bvalfile;
    
    private String outputFile;
    
    private double bScale;
    
    private boolean flipX;
    private boolean flipY;
    private boolean flipZ;
    
    private int numScans;

    private boolean interleave;
    
    /** 
     * Some scanners define b-value as |g|^2 * \beta where \beta is what 
     * they CLAIM the b-value is. If you use normalized gradient directions, 
     * you need to increase b accordingly to make the same ADC calculation.
     */
    private boolean useGradMod;


    // Set anything <= this to zero
    // Used to ensure that scans with minimal weighting get used as zero measurements for normalization purposes
    //
    private double zeroB_Value;
    

    public FSL_ToScheme(String[] args) {
        super(args);
        Locale.setDefault(Locale.UK);
    }


    public void initDefaultVals() {

        bvecfile = null;
        bvalfile = null;
        
        outputFile = null;
        
        bScale = 1.0;
        
        flipX = false;
        flipY = false;
        flipZ = false;
        
        numScans = 1;
        
        interleave = false;

        useGradMod = false;

        zeroB_Value = 0.0;

    }

    
    public void initOptions(String[] args) {

	CL_Initializer.CL_init(args);

	for (int i = 0; i < args.length; i++) {
	    if (args[i].equals("-bvecfile") || args[i].equals("-bvecs")) {
		bvecfile = args[i+1];
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-bvalfile") || args[i].equals("-bvals")) {
		bvalfile = args[i+1];
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-outputfile")) {
		outputFile = args[i+1];
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-bscale")) {
		bScale = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-zerobval")) {
	        zeroB_Value = Double.parseDouble(args[i+1]);
                logger.info("B-values <= " + zeroB_Value + " will be set to zero");	
		CL_Initializer.markAsParsed(i,2);
	    }
            else if (args[i].equals("-numscans")) {
                numScans = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-interleave")) {
                interleave = true;
                CL_Initializer.markAsParsed(i);
            }
	    else if (args[i].equals("-usegradmod")) {
		useGradMod = true;
		logger.info("Gradient direction magnitude will be incorporated into b-value");	
		CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-flipx")) {
		flipX = true;
		CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-flipy")) {
		flipY = true;
		CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-flipz")) {
		flipZ = true;
		CL_Initializer.markAsParsed(i);
	    }
	}

        CL_Initializer.checkParsing(args);
        
    }

    public void initVariables() {
        // Still not seeing the point of this method
    }
    

    public void execute(OutputManager om) {

	String scheme = "";

        
        double[][] bvecs = readB_Vectors(bvecfile);
       

	int numMeasSingleScan = bvecs[0].length;
	int numMeas = numScans * numMeasSingleScan;

        double[] bvals = readB_Values(bvalfile);
       
	if (bvals.length != numMeasSingleScan) {
	    throw new misc.LoggedException("Number of b-vectors and b-values do not match");
	}

	DecimalFormat df = new DecimalFormat("   0.000000;  -0.000000");

	DecimalFormat v0DF = new DecimalFormat("0.000000");

	DecimalFormat bValDF = new DecimalFormat("   0.000E00");

	scheme += "# Scheme file created by fsl2scheme\nVERSION: BVECTOR\n";

	boolean warnAboutGradDirMod = false;


        double[] vx = new double[numMeas];
        double[] vy = new double[numMeas];
        double[] vz = new double[numMeas];
        double[] bv = new double[numMeas];
        
	// add repeated scans if necessary
	if (numScans > 1) {
            

	    if (interleave) {
		for (int i = 0; i < numMeasSingleScan; i++) {
		    for (int s = 0; s < numScans; s++) {
			vx[i*numScans + s] = bvecs[0][i];
			vy[i*numScans + s] = bvecs[1][i];
			vz[i*numScans + s] = bvecs[2][i];
			bv[i*numScans + s] = bvals[i];
		    }
		}
	    }
	    
	    else {
		for (int s = 0; s < numScans; s++) {
		    System.arraycopy(bvecs[0], 0, vx, numMeasSingleScan * s, numMeasSingleScan);
		    System.arraycopy(bvecs[1], 0, vy, numMeasSingleScan * s, numMeasSingleScan);
		    System.arraycopy(bvecs[2], 0, vz, numMeasSingleScan * s, numMeasSingleScan);
		    System.arraycopy(bvals, 0, bv, numMeasSingleScan * s, numMeasSingleScan);
		}
	    }

	}
        else {
            vx = bvecs[0];
            vy = bvecs[1];
            vz = bvecs[2];
            bv = bvals;
        }

        for (int i = 0; i < numMeas; i++) {
	
            double modV = Math.sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

            if (modV == 0.0 || Math.abs(1.0 - modV) < 1E-4)
                modV = 1.0;

            if (!useGradMod && modV != 1.0) {
                warnAboutGradDirMod = true;
            }

            vx[i] = (vx[i] / modV);
		
            if (flipX && vx[i] != 0.0) {
                vx[i] = -vx[i];
            }
		
            vy[i] = (vy[i] / modV);
		
            if (flipY && vy[i] != 0.0) {
                vy[i] = -vy[i];
            }
		
            vz[i] = (vz[i] / modV);
		
            if (flipZ && vz[i] != 0.0) {
                vz[i] = -vz[i];
            }

            // get rid of -0.0
            if (vx[i] == 0.0) { 
                vx[i] = Math.abs(vx[i]);
            }
            if (vy[i] == 0.0) { 
                vy[i] = Math.abs(vy[i]);
            }
            if (vz[i] == 0.0) { 
                vz[i] = Math.abs(vz[i]);
            }

            double outputB = 0.0;
            
            if (useGradMod) {
                outputB = bv[i] * bScale * modV * modV;
            }
            else {
                outputB = bv[i] * bScale;
            }
            
            if (outputB <= zeroB_Value) {
                outputB = 0.0;
            }

            if (outputB == 0.0) {
                logger.info("Measurement " + i + " has b = 0, setting gradient to (0,0,0)");
                vx[i] = 0.0;
                vy[i] = 0.0;
                vz[i] = 0.0;
            }
				
            scheme += df.format(vx[i]);
            scheme += df.format(vy[i]);
            scheme += df.format(vz[i]);
		
            
            scheme += bValDF.format(outputB) + "\n";
            
        }

    
	
        if (warnAboutGradDirMod) {
            logger.warning("Some measurements have non-unit gradient directions. Directions have been " + 
                           "normalized to unit length");
        }
	
        
        if (outputFile != null) {
            FileOutput out = new FileOutput(outputFile);
	    
            out.writeString(scheme);
            
            out.close();
        }
        else {
            System.out.println(scheme);
        }


    } 


    private double[][] readB_Vectors(String vecFile) {
        
        // 3xN array, initialized once we know how many vectors there are
        double[][] vectors; 
        

        ArrayList<String> lines = new ArrayList<String>();

        try {

            File vecIn = new File(vecFile);

            Scanner scanner = new Scanner(vecIn);

            scanner.useDelimiter("\r\n|\n");
       
            while (scanner.hasNext()) {
                String s = scanner.next();
                
                // ignore empty lines
                if (s.length() > 0) {
                    lines.add(s);            
                }
            }
        }
        catch(IOException e) {
            throw new LoggedException("Couldn't read file " + vecFile + " : " + e);
        }

        // Assuming there are at least 7 measurements, file is in old FSL format if there are
        // only three lines

        int totalLines = lines.size();
    

        if (totalLines == 3) {
            
            // old FSL format

            String[] x = lines.get(0).trim().split("\\s+");

            int numMeas = x.length;

            vectors = new double[3][numMeas];

            for (int comp = 0; comp < 3; comp++) {
                String[] elements = lines.get(comp).trim().split("\\s+");
                
                for (int i = 0; i < numMeas; i++) {
                    vectors[comp][i] = Double.parseDouble(elements[i]);
                }
            }
            
        }
        else {
            
            int numMeas = totalLines;

            vectors = new double[3][numMeas];

            for (int i = 0; i < numMeas; i++) {
                String[] elements = lines.get(i).trim().split("\\s+");
                
                for (int comp = 0; comp < 3; comp++) {
                    vectors[comp][i] = Double.parseDouble(elements[comp]);
                }
            }
            

        }
    
        return vectors;

    }

    
    private double[] readB_Values(String valFile) {

        double[] bVals;

        // File is either nx1 or 1xn values

        ArrayList<String> lines = new ArrayList<String>();

        try {

            File valIn = new File(valFile);

            Scanner scanner = new Scanner(valIn);

            scanner.useDelimiter("\r\n|\n");
       
            while (scanner.hasNext()) {

                String s = scanner.next();

                // ignore empty lines
                if (s.length() > 0) {
                    lines.add(s);            
                }
            }
        }
        catch(IOException e) {
            throw new LoggedException("Couldn't read file " + valFile + " : " + e);
        }

        int totalLines = lines.size();        

       
        if (totalLines == 1) {
            String[] x = lines.get(0).trim().split("\\s+");
            
            int numMeas = x.length;
            
            bVals = new double[numMeas];
            
            for (int i = 0; i < numMeas; i++) {
                bVals[i] = Double.parseDouble(x[i]);
            }
        }
        else {
            int numMeas = totalLines;

            bVals = new double[numMeas];
            
            for (int i = 0; i < numMeas; i++) {
                bVals[i] = Double.parseDouble(lines.get(i));
            }
        }

        return bVals;

    }
    
}

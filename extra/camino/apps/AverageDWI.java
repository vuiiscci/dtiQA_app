package apps;

import data.*;
import imaging.*;
import misc.*;
import tools.*;
import tractography.*;

import java.io.*;
import java.util.*;
import java.util.logging.*;


/**
 * Computes an average over some or all of the DWI data
 *
 * 
 * @author Philip Cook
 *
 */
public class AverageDWI extends Executable {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.AverageDWI");


    
    private boolean median;

    private double minB;

    private double maxB;
   
    public AverageDWI(String[] args) {
        super(args);
    }


    public void initDefaultVals() {
        median = false;
        minB = 0.0;
        maxB = Double.MAX_VALUE;
    }
    

    public void initOptions(String[] args) {

	CL_Initializer.CL_init(args);

	for (int i = 0; i < args.length; i++) {
	    if (args[i].equals("-median")) {
		median = true;
		CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-minbval")) {
		minB = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-maxbval")) {
		maxB = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
	}

        CL_Initializer.checkParsing(args);
       
        CL_Initializer.initImagingScheme();

        CL_Initializer.initDataSynthesizer();
    }

    public void initVariables() {
 
    }


    public void execute(OutputManager om) {

        DW_Scheme scheme = CL_Initializer.imPars;

        int numMeas = scheme.numMeasurements();

        int numIncludedMeas = 0;

        List<Integer> includedMeasurements = new ArrayList<>();
        
        for (int i = 0; i < numMeas; i++) {

            double b = scheme.getB_Value(i);

            if (b >= minB && b <= maxB) {
                includedMeasurements.add(i);
                numIncludedMeas++;
            }
        }
        
        
        int[] indices = new int[numIncludedMeas];

        for (int i = 0; i < numIncludedMeas; i++) {
            indices[i] = includedMeasurements.get(i).intValue();
        }
       
        
        while( CL_Initializer.data.more() ) {
            
            double[] voxel = CL_Initializer.data.nextVoxel();
            
            double[] dataToAverage = new double[numIncludedMeas];
            
            for (int n = 0; n < numIncludedMeas; n++) {
                dataToAverage[n] = voxel[indices[n]];
            }
            
            double value = 0.0;
            
            if (median) {
                value = ArrayOps.median(dataToAverage);
            }
            else {
                value = ArrayOps.mean(dataToAverage);
            }
            
            om.output(new double[] {value});
        }
        

        om.close();
   
        
    }



}

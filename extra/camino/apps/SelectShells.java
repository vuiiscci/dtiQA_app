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
 * Extracts data and scheme files corresponding to selected shells. Also outputs any b=0 images
 *
 * Different from subsetscheme, which is used to subset data by individual directions.
 *
 *
 * 
 * @author Philip Cook
 *
 */
public class SelectShells extends Executable {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.SelectShells");

    private double minB;

    private double maxB;
    
    private double unweightedB;

    private String outputRoot;
    
    private boolean removeZeroMeas;


    public SelectShells(String[] args) {
        super(args);
    }


    public void initDefaultVals() {
        minB = 0.0;
        maxB = Double.MAX_VALUE;
	unweightedB = 0.0;
        outputRoot = "";
        removeZeroMeas = false;
    }
    

    public void initOptions(String[] args) {

        OutputManager.outputDataType = "float";

	CL_Initializer.CL_init(args);

	for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-minbval")) {
		minB = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-maxbval")) {
		maxB = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-unweightedb")) {
		unweightedB = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-outputroot")) {
		outputRoot = args[i+1];
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-removezeromeas")) {
                removeZeroMeas = true;
		CL_Initializer.markAsParsed(i);
	    }
	}

        CL_Initializer.checkParsing(args);
       
        CL_Initializer.initImagingScheme();

        CL_Initializer.initDataSynthesizer();

        if (ImageHeader.imageExists(CL_Initializer.inputFile)) {
            OutputManager.outputFile = outputRoot + ImageHeader.getFileExtension(CL_Initializer.inputFile);
        }
        else {
            OutputManager.outputFile = outputRoot + ".B" + OutputManager.outputDataType;
            
        }
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
            
            if (b <= unweightedB) {
                if (!removeZeroMeas) {
                    includedMeasurements.add(i);
                    numIncludedMeas++;
                }
            }
            else if ( b >= minB && b <= maxB ) {
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
            
            double[] selectedData = new double[numIncludedMeas];
            
            for (int n = 0; n < numIncludedMeas; n++) {
                selectedData[n] = voxel[indices[n]];
            }
          
            om.output(selectedData);
        }
        

        om.close();


        // Now get the scheme
        DW_Scheme subsetScheme = CL_Initializer.imPars.getSubsetScheme(indices);
        FileOutput out = new FileOutput(outputRoot + ".scheme");
        out.writeString(subsetScheme.toString());
        out.close();

        
   
        
    }



}

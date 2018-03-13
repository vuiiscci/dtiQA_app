package apps;

import imaging.*;
import misc.*;
import tools.*;

import java.io.*;
import java.text.*;
import java.util.Locale;

/**
 * Converts Camino scheme files to FSL b-vectors and b-values.
 * 
 * @author Philip Cook
 * @version $Id$
 */
public class SchemeToFSL {


    public static void main(String[] args) {
	
	Locale.setDefault(Locale.UK);

	String bvecfile = null;
	String bvalfile = null;

	String outputRoot = "scheme";

	double diffusionTime = 1.0;

	boolean flipX = false;
	boolean flipY = false;
	boolean flipZ = false;

	double bScale = 1E-6;

	CL_Initializer.CL_init(args);

	for (int i = 0; i < args.length; i++) {
	    if (args[i].equals("-outputroot")) {
		outputRoot = args[i+1];
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-bscale")) {
		bScale = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
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
	    
	DecimalFormat vecFormat = new DecimalFormat("0.000000");
	DecimalFormat valFormat = new DecimalFormat("0.000");
       
        // should use inputfile but allow schemefile
        String inputFile = CL_Initializer.inputFile != null ? CL_Initializer.inputFile : CL_Initializer.schemeFile;
        
        DW_Scheme input = DW_Scheme.readScheme(inputFile);
            
        if (flipX) {
            input = input.flipX();
        }
        if (flipY) {
            input = input.flipY();
        }
        if (flipZ) {
            input = input.flipZ();
        }
            
        String[] bvecs = new String[] {"", "", ""};

        String bvals = "";

        for (int i = 0; i < input.numMeasurements(); i++) {
                
            bvals += valFormat.format(input.getB_Value(i) * bScale) + " ";
	
            double[] gHat = input.getG_Dir(i);

            bvecs[0] += vecFormat.format(gHat[0]) + " ";
            bvecs[1] += vecFormat.format(gHat[1]) + " ";
            bvecs[2] += vecFormat.format(gHat[2]) + " ";

        }

        FileOutput out = new FileOutput(outputRoot + ".bvals");
	    
        out.writeString(bvals.trim() + "\n");
        out.close();

        out = new FileOutput(outputRoot + ".bvecs");
        out.writeString(bvecs[0] + "\n");
        out.writeString(bvecs[1] + "\n");
        out.writeString(bvecs[2] + "\n");
        out.close();


    } 




}

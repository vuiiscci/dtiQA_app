package apps;

import imaging.*;
import misc.LoggedException;
import tools.*;

import java.io.*;
import java.text.*;

import data.*;

/**
 * Splits a 4D Nifti image into 3D component volumes.
 *
 * @author Philip Cook
 */
public class Split4D_NiiImage extends Executable {

    public Split4D_NiiImage(String[] args) {
        super(args);
    }
    
    private String outputRoot;
     
    public void initDefaultVals() {
        outputRoot = "";
    }
    
    
    public void initOptions(String[] args) {
        
        CL_Initializer.CL_init(args);
	
		for (int i = 0; i < args.length; i++) {
                    if (args[i].equals("-outputroot")) {
                        outputRoot = args[i+1];
                        CL_Initializer.markAsParsed(i, 2);
                    }
                }
                
                CL_Initializer.checkParsing(args);

                if (CL_Initializer.schemeFile != null) {
                    CL_Initializer.initImagingScheme();
                }
    }
    

    public void initVariables() {

    }


    public void execute(OutputManager om) { 	
        
	Nifti1Dataset input = new Nifti1Dataset(CL_Initializer.inputFile);
        
	if (input.components() == 1) {
	    throw new LoggedException("This image is already three-dimensional");
	}
        
	try {
	    Nifti1Dataset.convertTo3D(CL_Initializer.inputFile, outputRoot);
	}
	catch (IOException e) {
	    throw new LoggedException(e);
	}

        if (CL_Initializer.imPars != null) {
            
            // produce CSV file with scheme info and imaging volumes

            DW_Scheme scheme = CL_Initializer.imPars;

            int numMeas = scheme.numMeasurements();

            if (numMeas != input.components()) {
                throw new LoggedException("Number of measurements in scheme file does not match number of volumes");
            }

            FileOutput out = new FileOutput(outputRoot + "components.csv");

            DecimalFormat dfCounter = new DecimalFormat("0000");


            DecimalFormat dfGrad = new DecimalFormat("0.000000");
            DecimalFormat dfB = new DecimalFormat("0.000E00");

            String ext = ImageHeader.getFileExtension(CL_Initializer.inputFile);

            String header = "GradX,GradY,GradZ,B_Value,ImageFile";

            out.writeString(header + "\n");
            
            for (int i = 0; i < numMeas; i++) {

                double[] grad = scheme.getG_Dir(i);
            
		//    String line = String.join(",", dfGrad.format(grad[0]), dfGrad.format(grad[1]), dfGrad.format(grad[2]),
		//                              dfB.format(scheme.getB_Value(i)), outputRoot + dfCounter.format(i+1) + ext);
		
                // do it inefficiently to avoid introducing a dependence on Java 8
                String line = dfGrad.format(grad[0]) + "," + dfGrad.format(grad[1]) + "," + dfGrad.format(grad[2]) + "," + 
		    dfB.format(scheme.getB_Value(i)) + "," + outputRoot + dfCounter.format(i+1) + ext;
		
                out.writeString(line + "\n");
                
            }

            out.close();
            
        }

    }



}

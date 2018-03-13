package apps;

import data.*;
import imaging.*;
import misc.*;
import tools.*;
import tractography.*;

import java.io.*;
import java.util.logging.*;


/**
 * Computes voxelwise statistics on an arbitrary set of input data. All images must be in the same physical space.
 *
 * 
 * @author Philip Cook
 * @version $Id: CP_Stats.java 521 2008-04-01 20:23:20Z ucacpco $
 *
 */
public class VoxelwiseImageStats extends Executable {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.VoxelwiseImageStats");


    
    private String[] imageFiles;
    
    private String outputRoot;
    
    private String operation;

    public VoxelwiseImageStats(String[] args) {
        super(args);
    }


    public void initDefaultVals() {
        imageFiles = null;
        outputRoot = null;
        operation = null;
    }
    

    public void initOptions(String[] args) {

        OutputManager.outputDataType = "float";

	CL_Initializer.CL_init(args);

	for (int i = 0; i < args.length; i++) {
	    if (args[i].equals("-images")) {
		String[] tmp = new String[100000];
                
                int t = 0;
                
                while (t + i + 1 < args.length && !args[i+t+1].startsWith("-")) {
		    tmp[t] = args[i+t+1];
		    t++;
                }
                
		imageFiles = new String[t];

		System.arraycopy(tmp, 0, imageFiles, 0, t);
                
		CL_Initializer.markAsParsed(i, t+1);
	    }
	    else if (args[i].equals("-stat")) {
		operation = args[i+1];
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-outputroot")) {
		outputRoot = args[i+1];
		CL_Initializer.markAsParsed(i,2);
	    }
	}
        
        if (outputRoot.equals("")) {
            if (OutputManager.outputFile != null) {
                outputRoot = ImageHeader.getFileRoot(OutputManager.outputFile);
            }
            else {
                throw new LoggedException("No output image file specified");
            }
        }


        CL_Initializer.checkParsing(args);

        
        CL_Initializer.headerTemplateFile = imageFiles[0];
	CL_Initializer.initInputSpaceAndHeaderOptions();


    }

    public void initVariables() {

    }


    public void execute(OutputManager om) {

        VoxelwiseStatisticalImage image;
        
        // if we are computing something that can be done incrementally, save
        // time by using a faster implementation
        
        // otherwise save memory by using a sparse vector image to store all the 
        // information we need

        // var can be done incrementally, or at least approximated, but we don't currently 
        // have that capability

        if (operation.equals("median") || operation.equals("var") || operation.equals("std") ) {
            image = new SparseVectorImage( CL_Initializer.headerTemplate.getDataDims(),  CL_Initializer.headerTemplate.getVoxelDims(), 
                                           imageFiles.length);
        }
        else {
            image = new DynamicScalarImage( CL_Initializer.headerTemplate.getDataDims(),  CL_Initializer.headerTemplate.getVoxelDims());
        }

        try {
            
            for (int n = 0; n < imageFiles.length; n++) {
                
                ImageHeader hdr = ImageHeader.readHeader(imageFiles[n]);
                
                if (! CL_Initializer.headerTemplate.sameSpace(hdr) ) {
                    throw new LoggedException("Image " + imageFiles[n] + " is not in the same physical space as " +  
                                              CL_Initializer.headerTemplate);
                } 
                
                logger.info("Processing image " + imageFiles[n]);
                
                double[][][] volume = hdr.readSingleVolumeData();
                
                for (int k = 0; k < CL_Initializer.dataDims[2]; k++) {
                    for (int j = 0; j < CL_Initializer.dataDims[1]; j++) {
                        for (int i = 0; i < CL_Initializer.dataDims[0]; i++) {

                            image.addValue(i,j,k, volume[i][j][k]);   
                        }
                    }
                }
                
            }

            double[][][] result = image.getVoxelStatistic(operation);

            ImageHeader output =  CL_Initializer.headerTemplate.writeScalarImage(result, outputRoot);

            logger.info("Output written to " + output.getHeaderFilename());


        }
	catch (IOException e) {
	    throw new LoggedException(e);
	}
    }



}

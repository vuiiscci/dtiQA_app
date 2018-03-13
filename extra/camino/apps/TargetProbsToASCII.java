package apps;

import imaging.*;
import misc.*;
import tools.*;
import data.*;
import tractography.*;

import java.io.*;
import java.text.*;
import java.util.*;



/**
 * Converts the results of tracking with target volumes to a text file.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class TargetProbsToASCII {


    private static int minTargetIndex = 0;
    private static int maxTargetIndex = 0;

    // targetIndices[i] gives order of target. Ordering is zero-indexed and ascending.
    // If there are targets labelled 1 3 5 7 9 11 in the volume, then targetIndices[1] == 0 and 
    // targetIndices[7] == 3.
    private static int[] targetIndices = null;

    // targetExists == true if there is a target with label i
    private static boolean[] targetExists = null;


    public static void main(String[] args) {

	Locale.setDefault(Locale.UK);
        
        String inputRoot = null;
        String targetFile = null;
	String seedFile = null;

	int regionIndex = 1;

	int pd = 1;

        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-inputroot")) {
                inputRoot = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            else if  (args[i].equals("-targetfile")) {
                targetFile = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
	    else if (args[i].equals("-seedfile")) {
		seedFile = args[i + 1];
		CL_Initializer.markAsParsed(i, 2);
	    } 
	    else if (args[i].equals("-pd")) {
		pd = Integer.parseInt(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
 	    else if (args[i].equals("-regionindex")) {
                regionIndex = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
	    } 

        }

        CL_Initializer.checkParsing(args);

        CL_Initializer.headerTemplateFile = seedFile;
        
        CL_Initializer.initInputSpaceAndHeaderOptions();
	
	ImageHeader ah = null;

        String ext = null;

        if (seedFile.endsWith(".gz")) {
            ext = seedFile.substring(seedFile.length() - 7);
        }
        else {
            ext = seedFile.substring(seedFile.length() - 4);
        }

	try {
	    ah = ImageHeader.readHeader(seedFile);
	}
	catch (IOException e) {
	    throw new LoggedException("Cannot read header of seed file");
	}
		    
	int xDataDim = ah.xDataDim();
	int yDataDim = ah.yDataDim();
	int zDataDim = ah.zDataDim();

	double xVoxelDim = Math.abs(ah.xVoxelDim());
	double yVoxelDim = Math.abs(ah.yVoxelDim());
	double zVoxelDim = Math.abs(ah.zVoxelDim());

	int[][][] seeds = ProcessStreamlines.readIntVolume(seedFile);

	inputRoot = inputRoot + Integer.toString(regionIndex) + "_";
	
	double[][][] targets = null;	
	
	try {
	    targets = ImageHeader.readHeader(targetFile).readSingleVolumeData();
	}
	catch (IOException e) {
	    throw new LoggedException("Cannot read header of target file");
	}
 
        setTargetIndices(targets);

        VoxelROI allROIs = new VoxelROI(seedFile, ah);

	Voxel[] seedVoxels = allROIs.getRegion(regionIndex).getSeedVoxels(xVoxelDim, yVoxelDim, zVoxelDim);
	

        // loop over all CP images, get connection probs to each target volume
        // output as
        // (seed)         (probs)
        // x \t y \t z \t cp_1 \t cp_2 ... cp_n \n


        // Set up the scientific notation format.
        String zeros = "";

        int decimalPlaces = 6;

        for (int i = 0; i < decimalPlaces; i++) {
            zeros += "0";
        }
        DecimalFormat sciForm = new DecimalFormat("0." + zeros + "E00;-0." + zeros + "E00");


	StringBuffer buffer = new StringBuffer();

        int volume = 1;
	
        String imageRoot = inputRoot + volume + "_" + pd;
	
        while (new File(imageRoot + ext).exists()) {
            
            double[] targetProbs = new double[maxTargetIndex - minTargetIndex + 1];
            
            double[][][] probs = null;
            
            try {
                probs = ImageHeader.readHeader(imageRoot + ext).readSingleVolumeData();
            }
            catch (IOException e) {
                throw new LoggedException("Cannot read header " + imageRoot + ext);
            }
            
            for (int k = 0; k < zDataDim; k++) {
                for (int j = 0; j < yDataDim; j++) {
			for (int i = 0; i < xDataDim; i++) {
                            if (probs[i][j][k] > 0.0) {
                                targetProbs[targetIndices[(int)targets[i][j][k]]] = probs[i][j][k];
                            }
                        }
                }
            } 
            
            String delim = "\t";
            
            buffer.append(seedVoxels[volume-1].x + delim + seedVoxels[volume-1].y + delim + 
                          seedVoxels[volume-1].z + delim);
		
            // output value for targets in target file only
            
            for (int i = minTargetIndex; i <= maxTargetIndex; i++) {
                if (targetExists[i]) { 
                    buffer.append(sciForm.format(targetProbs[targetIndices[i]]));
                    
                    if (i < maxTargetIndex) {
                        buffer.append(delim);
			
                    }
                }
            }
            
            buffer.append("\n");
            
            volume++;
            
            imageRoot = inputRoot + volume + "_" + pd;
        }
	    
	


	if (OutputManager.outputFile != null) {
	    FileOutput out = new FileOutput(OutputManager.outputFile);
	    out.writeString(buffer.toString());
	    out.close();
	}
	else {
	    System.out.print(buffer.toString());
	}

        
    }



    /**
     * @return the minimum (but > 0) and maximum values in the volume. Assumes volumes can be 
     * represented by shorts (this caps maximum memory consumption).
     */
    private static void setTargetIndices(double[][][] vol) {
        
        int minIndex = Short.MAX_VALUE;
        int maxIndex = 0;

    	int xDataDim = vol.length;
	int yDataDim = vol[0].length;
	int zDataDim = vol[0][0].length;

        targetExists = new boolean[Short.MAX_VALUE];
               
	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
		    if (vol[i][j][k] < minIndex && vol[i][j][k] > 0) {
                        minIndex = (short)vol[i][j][k];
		    }
                    else if (vol[i][j][k] > maxIndex) {
                        maxIndex = (short)vol[i][j][k];                        
                    }

                    if (vol[i][j][k] > 0) {
                        targetExists[(short)vol[i][j][k]] = true;
                    }
		}
	    }
	}

        minTargetIndex = minIndex;
        maxTargetIndex = maxIndex;
        

        int targetCounter = 0;

        targetIndices = new int[maxTargetIndex + 1];

        for (int i = minTargetIndex; i <= maxTargetIndex; i++) {
            if (targetExists[i]) {
                targetIndices[i] = targetCounter++;
            }
        }
        
    }


}

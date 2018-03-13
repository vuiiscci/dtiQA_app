package apps;


import java.io.*;

import imaging.*;

import misc.LoggedException;

/**
 * Counts the number of seed voxels (voxels with values > 0).
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class CountSeeds {

    
    public static void main(String[] args) {


	int regionIndex = 0;
	
	if (args.length == 2) {
	    regionIndex = Integer.parseInt(args[1]);
	}


        ImageHeader ih = null;

        try {
            ih = ImageHeader.readHeader(args[0]);
        }
        catch(IOException e) {
            throw new LoggedException(e);
        }

        
        double[][][] seeds = ih.readSingleVolumeData();
        
        int xDataDim = ih.xDataDim();
        int yDataDim = ih.yDataDim();
        int zDataDim = ih.zDataDim();
                
        int seedCounter = 0;
        
        for (int k = 0; k < zDataDim; k++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
		    
		    if (regionIndex > 0) {
			if (seeds[i][j][k] == regionIndex) {
			    seedCounter++;
			}
		    }
		    else {
			if (seeds[i][j][k] > 0.0) {
			    seedCounter++;
			}
		    }
                }
            }
        }

        System.out.println(seedCounter);

    }





}

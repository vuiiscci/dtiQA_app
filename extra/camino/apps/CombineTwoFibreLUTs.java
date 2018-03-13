package apps;

import java.io.*;


/**
 *
 * Combines LUTs with different crossing angles into one.
 * 
 * @see apps.GenerateDTLUT
 * @author Philip Cook
 * @version $Id$
 */
public final class CombineTwoFibreLUTs {

  
    public static void main(String[] args) {


	String[] luts = null;

	double[] angleRange = null;
	double angleStep = 0.0;

	if (args.length == 0) {
	    System.exit(0);
	}


	for (int i = 0; i < args.length; i++) {

            if (args[i].equals("-luts")) {
		
		int numLUTs = 0;

		while (i + numLUTs + 1 < args.length && !args[i + numLUTs + 1].startsWith("-")) {
		    numLUTs++;
		}
		
		luts = new String[numLUTs];

		for (int j = 0; j < numLUTs; j++) {
		    luts[j] = args[i+j+1];
		}
		
	    }
	    else if (args[i].equals("-anglerange")) {
		angleRange = new double[2];
		angleRange[0] = Double.parseDouble(args[i+1]);
		angleRange[1] = Double.parseDouble(args[i+2]);
	    }

	    
	}

	// degrees
	angleStep = (angleRange[1] - angleRange[0]) / (luts.length - 1.0);

	try {
	    combineTwoFibreLUTs(luts, angleRange, angleStep);
	}
	catch(IOException e) {
	    System.err.println(e);
	}

    }




    /**
     * Combine two tensor LUTs. Each LUT should be two-dimensional, giving k1, k2 as a function of FA(D_1) 
     * and FA(D_2) for a fixed crossing angle. This method combines the LUTs into 3D, where the third 
     * dimension is the crossing angle. The output is suitable for reading by PICoLUTs.
     *
     * @param luts should be in order from min crossing angle to max.
     * @param angleRange in degrees. Will be converted to radians.
     * @param angleStep in degrees. Will be converted to radians.
     * @param dataDims {xMin xMax yMin yMax zMin zMax}.
     */
    public static void combineTwoFibreLUTs(String[] luts, double[] angleRange, double angleStep) throws IOException {
	
	DataInputStream[] streams = new DataInputStream[luts.length];

	double[] dataDims = new double[6];
	double[] steps = new double[3];

	dataDims[4] = Math.PI * angleRange[0] / 180.0;
	dataDims[5] = Math.PI * angleRange[1] / 180.0;

	steps[2] = Math.PI * angleStep / 180.0 ;

	int values = 0;

	for (int s = 0; s < streams.length; s++) {

	    FileInputStream fit = new FileInputStream(luts[s]);
	    streams[s] = new DataInputStream(new BufferedInputStream(fit, 1024*1024*2));
	    
	    if (s == 0) { // set up data dims and steps
		streams[s].readDouble();
		streams[s].readDouble();

		for (int i = 0; i < 4; i++) {
		    dataDims[i] = streams[s].readDouble();
		}
		streams[s].readDouble();
		for (int i = 0; i < 2; i++) {
		    steps[i] = streams[s].readDouble();
		}
	
		values = (int)streams[s].readDouble();
	    }
	    else { // read through header info
		for (int i = 0; i < 10; i++) {
		    streams[s].readDouble();
		}
	    }

	}

	// actual dims of array
	int xDataDim = (int)(1 + Math.round((dataDims[1] - dataDims[0]) / steps[0]));
	int yDataDim = (int)(1 + Math.round((dataDims[3] - dataDims[2]) / steps[1]));
	int zDataDim = (int)(1 + Math.round((dataDims[5] - dataDims[4]) / steps[2]));

	double[][][] data = new double[xDataDim][yDataDim][values * zDataDim];


	for (int v = 0; v < values; v ++) {
	    for (int k = 0; k < zDataDim; k++) {
		for (int j = 0; j < yDataDim; j++) {
		    for (int i = 0; i < xDataDim; i++) {
			data[i][j][values * k + v] = streams[k].readDouble();
		    }
		}
	    }
	}
	
	for (int s = 0; s < streams.length; s++) {
	    streams[s].close();
	}

	DataOutputStream dout = new DataOutputStream(System.out);
	
	for (int i = 0; i < 6; i++) {
	    dout.writeDouble(dataDims[i]);
	}
	for (int i = 0; i < 3; i++) {
	    dout.writeDouble(steps[i]);
	}

	System.err.println("Data dims: " + dataDims[0] + " " + dataDims[1] + " " + dataDims[2] + " " +
			   dataDims[3] + " " + dataDims[4] + " " + dataDims[5]);

	System.err.println("Steps: " + steps[0] + " " + steps[1] + " " + steps[2]);

	dout.writeDouble(values); // values per position
	
	for (int v = 0; v < values; v ++) {
	    for (int k = 0; k < zDataDim; k++) {
		for (int j = 0; j < yDataDim; j++) {
		    for (int i = 0; i < xDataDim; i++) {
			dout.writeDouble(data[i][j][values * k + v]);
		    }
		}
	    }
	}


	dout.close();

    } 

   
}

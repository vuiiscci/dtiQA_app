package tractography;

import data.*;
import imaging.*;
import inverters.*;

import Jama.*;
import misc.DT;
import numerics.*;

import tools.*;

import Jama.*;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Random;
import java.util.Arrays;


/**
 * Fits some distribution to samples of tensors with added noise. 
 *
 * @author Philip Cook
 * @version $Id$
 */
public class OneTensorLUTGenerator extends DT_LookupTableGenerator {

 

    public OneTensorLUTGenerator(DW_Scheme ip, double snr, double trace, Random r) {
	super(ip, snr, trace, r);
    }



    /**
     * Generates the LUT. 
     * 
     * @param xMin minimum L1 / L3, the x index of the LUT. Minimum value 1.0.
     * @param xMax maximum L1 / L3, the x index of the LUT. 
     *
     * @param step must divide evenly into (yMax - yMin) and (xMax - xMin).
     */
    public double[][][][] generateLUT(double xMin, double xMax, double step, int samples, 
                                      ModelIndex inversionIndex,
                                      boolean watson, boolean bingham, boolean acg) {


	if (xMax < 1.0) {
	    return generateFA_LUT(xMin, xMax, step, samples, inversionIndex, watson, bingham, acg);
	}
	else {
	    return generateEV_LUT(xMin, xMax, step, samples, inversionIndex, watson, bingham, acg);
	}

    }


    /**
     * Generates the LUT. 
     * 
     * @param xMin minimum FA, must be >= 0.0.
     * @param xMax maximum FA, must be < 1.0.
     *
     * @param step must divide evenly into (yMax - yMin) and (xMax - xMin).
     */
    public double[][][][] generateFA_LUT(double xMin, double xMax, double step, int samples, 
					 ModelIndex inversionIndex, boolean watson, boolean bingham, boolean acg) {
	
	
	int smallEigenvalueCounter = 0;
	
	// x == FA

	double x = 0.0;

	
	double[][][][] kappaLUT = 
	    new double[3][1][1 + (int)Math.round(( (xMax - xMin ) / step ))][];

	int xCounter = 0;

	x = xMin;

        double totalIterations = kappaLUT[0][0].length;

        DecimalFormat df = new DecimalFormat("00.00");

	while (xCounter < kappaLUT[0][0].length) {
	    
            System.err.print("\r ... " + 
                             df.format( 100.0 * xCounter / (totalIterations) )
                             + "% completed");

		
		DT[] tensors = getNoisyTensors(new DT[] {getTensor(x, trace, ran)}, 
					       new double[] {1.0}, inversionIndex, imPars, 
					       snr, samples, ran)[0];


		// e1s of noisy tensors
		Vector3D[] sampleVecs = new Vector3D[samples];
		
		int fittingErrors = 0;

		boolean[] includeSample = new boolean[samples];

		for (int i = 0; i < samples; i++) {

		    if (tensors[i] == null) {
			fittingErrors++;
			continue;
		    }

		    double[][] eig = tensors[i].sortedEigenSystem();

                    includeSample[i] = true;
                    sampleVecs[i] = new Vector3D(eig[1][0], eig[2][0], eig[3][0]);

		    
		}

                
		if (fittingErrors > 5.0 * samples / 100.0) {
		    double percentageFailed = fittingErrors * 100.0 / (double)(samples); 
		    System.err.println("Warning: discarding " + percentageFailed + 
				       " % of samples at point : (" + x + ")");
		}

		// filter bad tensor fits from results
		Vector3D[] tmpSampleVecs = new Vector3D[samples - fittingErrors];
		int sampleCounter = 0;
		for (int i = 0; i < samples; i++) {
		    if (includeSample[i]) {
			tmpSampleVecs[sampleCounter++] = sampleVecs[i];
		    }
		}

		sampleVecs = tmpSampleVecs;

		if (watson) {
		    kappaLUT[WATSON][0][xCounter] = getWatsonConcentrationParams(sampleVecs);
		}
		else {
		    kappaLUT[WATSON][0][xCounter] = new double[1];
		}

		if (bingham) {
		    kappaLUT[BINGHAM][0][xCounter] = getBinghamConcentrationParams(sampleVecs);
		}
		else {
		    kappaLUT[BINGHAM][0][xCounter] = new double[2];
		}

		if (acg) {
		    kappaLUT[ACG][0][xCounter] = getACGConcentrationParams(sampleVecs);
		}
		else {
		    kappaLUT[ACG][0][xCounter] = new double[3];
		}


	    xCounter++;
	    x += step;
	}    
	
        System.err.println();

	return kappaLUT;
	
    }


    /**
     * Generates the LUT. 
     * 
     * @param xMin minimum L1 / L3, the x index of the LUT. Minimum value 1.0.
     * @param xMax maximum L1 / L3, the x index of the LUT. 
     *
     * @param step must divide evenly into (yMax - yMin) and (xMax - xMin).
     */
    public double[][][][] generateEV_LUT(double xMin, double xMax, double step, int samples, 
                                      ModelIndex inversionIndex,
                                      boolean watson, boolean bingham, boolean acg) {

	int smallEigenvalueCounter = 0;
	
	// x == L1 / L3
	// y == L2 / L3
	double x, y;

        double yMin = xMin;        
        double yMax = xMax;
	
	// kappaLUT[x][y] = params for x,y
	double[][][][] kappaLUT = 
	    new double[3][1 + (int)Math.round( (xMax - xMin ) / step )][1 + (int)Math.round( (yMax - yMin ) / step )][];

	


	double l1, l2, l3;

	int xCounter = 0;
	int yCounter = 0;
	

	x = xMin;
	y = yMin;

	// increment x and y by step until they exceed the max

        double totalIterations = kappaLUT[0].length * (kappaLUT[0].length + 1) / 2.0;

        DecimalFormat df = new DecimalFormat("00.00");

	while (xCounter < kappaLUT[0].length) {
	    
            System.err.print("\r ... " + 
                             df.format( 100.0 * xCounter * (xCounter+1) / (2.0 * totalIterations) )
                             + "% completed");

	    while (yCounter < kappaLUT[0][0].length) {

		
		if (xCounter < yCounter) {
		    kappaLUT[WATSON][xCounter][yCounter] = new double[1];
		    kappaLUT[BINGHAM][xCounter][yCounter] = new double[2];
		    kappaLUT[ACG][xCounter][yCounter] = new double[3];
	    
		    yCounter++;
		    y += step;
		    continue;
		}
		
		DT[] tensors = getNoisyTensors(new DT[] {getTensor(x, y, trace, ran)}, 
					       new double[] {1.0}, inversionIndex, imPars, 
					       snr, samples, ran)[0];


		// e1s of noisy tensors
		Vector3D[] sampleVecs = new Vector3D[samples];
		
		int fittingErrors = 0;

		boolean[] includeSample = new boolean[samples];

		for (int i = 0; i < samples; i++) {

		    if (tensors[i] == null) {
			fittingErrors++;
			continue;
		    }

		    double[][] eig = tensors[i].sortedEigenSystem();

                    includeSample[i] = true;
                    sampleVecs[i] = new Vector3D(eig[1][0], eig[2][0], eig[3][0]);

		    
		}

                
		if (fittingErrors > 5.0 * samples / 100.0) {
		    double percentageFailed = fittingErrors * 100.0 / (double)(samples); 
		    System.err.println("Warning: discarding " + percentageFailed + 
				       " % of samples at point : (" + x + ", " + y + ")");
		}

		// filter bad tensor fits from results
		Vector3D[] tmpSampleVecs = new Vector3D[samples - fittingErrors];
		int sampleCounter = 0;
		for (int i = 0; i < samples; i++) {
		    if (includeSample[i]) {
			tmpSampleVecs[sampleCounter++] = sampleVecs[i];
		    }
		}

		sampleVecs = tmpSampleVecs;

		if (watson) {
		    kappaLUT[WATSON][xCounter][yCounter] = getWatsonConcentrationParams(sampleVecs);
		}
		else {
		    kappaLUT[WATSON][xCounter][yCounter] = new double[1];
		}

		if (bingham) {
		    kappaLUT[BINGHAM][xCounter][yCounter] = getBinghamConcentrationParams(sampleVecs);
		}
		else {
		    kappaLUT[BINGHAM][xCounter][yCounter] = new double[2];
		}

		if (acg) {
		    kappaLUT[ACG][xCounter][yCounter] = getACGConcentrationParams(sampleVecs);
		}
		else {
		    kappaLUT[ACG][xCounter][yCounter] = new double[3];
		}



		
		yCounter++;
		y += step;
		
	    }
	
	    xCounter++;
	    x += step;
	    
	    yCounter = 0;
	    y = yMin;
	
	    
	}    
	
        System.err.println();


	return kappaLUT;
	
    }





   protected double[] getWatsonConcentrationParams(Vector3D[] sampleVecs) {
	
	EigenSystem3D eig = WatsonFitter.tBarEigenSystem(sampleVecs);
	
	return new double[] {WatsonFitter.fitKappa(eig, sampleVecs)};

    }

   
    protected double[] getBinghamConcentrationParams(Vector3D[] sampleVecs) {
	
	EigenSystem3D eig = BinghamFitter.tBarEigenSystem(sampleVecs);
	
	try {
	    double[] bngpar = BinghamFitter.bngpar(eig.eigenvalues[1], eig.eigenvalues[2]);

	    return new double[] {bngpar[0], bngpar[1]};

	}
	catch (ConvergenceException e) {
	    return new double[] {0.0, 0.0};
	}

    }


   protected double[] getACGConcentrationParams(Vector3D[] sampleVecs) {

	RealMatrix A = ACG_Fitter.findA(sampleVecs);
	
	Jama.Matrix aj = new Jama.Matrix(3,3);
	
	aj.set(0,0, A.entries[0][0]);
	aj.set(0,1, A.entries[0][1]);
	aj.set(0,2, A.entries[0][2]);
	aj.set(1,0, A.entries[1][0]);
	aj.set(1,1, A.entries[1][1]);
	aj.set(1,2, A.entries[1][2]);
	aj.set(2,0, A.entries[2][0]);
	aj.set(2,1, A.entries[2][1]);
	aj.set(2,2, A.entries[2][2]);
	
	EigenSystem3D eigA = EigenSystem3D.sort(aj.eig());
	
	return new double[] {eigA.eigenvalues[0], eigA.eigenvalues[1], eigA.eigenvalues[2]};

    }

}				   
	
				   







    

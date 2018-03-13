package tractography;

import Jama.*;

import imaging.*;
import inverters.ModelIndex;
import misc.DT;
import numerics.*;
import optimizers.*;

import java.text.DecimalFormat;
import java.util.Random;
import java.util.Arrays;


/**
 * Fits some distribution to samples of tensors with added noise. 
 *
 * @author Philip Cook
 * @version $Id$
 */
public class TwoTensorLUTGenerator extends DT_LookupTableGenerator {


    protected double prop;
    protected double rotAngle;
    
    /**
     * @param prop the mixing proportion of the two tensors.
     * @param cross the crossing angle (in radians) of the tensors.
     */
     public TwoTensorLUTGenerator(DW_Scheme ip, double snr, double trace, 
				     double prop, double cross, Random r) {

	super(ip, snr, trace, r);

	this.prop = prop;
	rotAngle = cross;
    }


    /**
     * Generates the LUT. 
     * 
     * @param xMin minimum fa of dt1, the x index of the LUT.
     * @param xMax maximum fa of dt1, the x index of the LUT. 
     *
     * @param step must divide evenly into (yMax - yMin) and (xMax - xMin).
     */
    public double[][][][] generateLUT(double xMin, double xMax, double step, int samples, 
                                      ModelIndex inversionIndex,
				    boolean watson, boolean bingham, boolean acg) {

        double yMin = xMin;        
        double yMax = xMax;

	
	// kappaLUT[x][y] = params for x,y
	double[][][][] kappaLUT = 
	    new double[3][1 + (int)( Math.round((xMax - xMin ) / step) )]
	    [1 + (int)( Math.round((yMax - yMin ) / step) )][];

	
	double[] propArray = new double[] {prop, 1.0 - prop};
	
	RealMatrix dt2Rot = Rotations.getRotMat(e2Vec, rotAngle);

	Vector3D dt2E1Vec = Rotations.rotateVector(e1Vec, dt2Rot);

	double x = xMin;
	double y = yMin;

	int xCounter = 0;
	int yCounter = 0;

        double totalIterations = kappaLUT[0].length * (kappaLUT[0].length + 1) / 2.0;

        if (prop != 0.5) {
            totalIterations = kappaLUT[0].length * kappaLUT[0][0].length;
        }

        DecimalFormat df = new DecimalFormat("00.00");
	
	while (xCounter < kappaLUT[0].length) {

            System.err.print("\r ... " + 
                             df.format( 100.0 * xCounter * (xCounter+1) / (2.0 * totalIterations) )
                             + "% completed");
	    
	    while (yCounter < kappaLUT[0][0].length) {

		if (xCounter < yCounter && prop == 0.5) {
		    // if prop is equal then LUT is symmetric. Save time by not doing both halves
		    yCounter++;
		    y += step;
		    continue;
		}
		
		Vector3D[] tmpSampleVecs = new Vector3D[2 * samples];

		DT dt1 = getTensor(x, trace, ran);
		DT dt2 = getTensor(y, trace, ran);
		
		dt2 = dt2.transform(dt2Rot);

		DT[][] tensors = getNoisyTensors(new DT[] {dt1, dt2}, propArray, 
						 inversionIndex, imPars, snr, samples, ran);
		
		
		int successfulTensorFits = 0;
		
		for (int i = 0; i < samples; i++) {
		    
		    if (tensors[0][i] == null) {
			// fitting failed
			continue;
		    }
		    
		    double[][] eig1 = tensors[0][i].sortedEigenSystem();
		    double[][] eig2 = tensors[1][i].sortedEigenSystem();

		    // check for fitting errors
		    // specifically negative eigenvalues
		    if (eig1[0][0] < 0.0 || eig1[0][1] < 0.0 || eig1[0][2] < 0.0 || 
			eig2[0][0] < 0.0 || eig2[0][1] < 0.0 || eig2[0][2] < 0.0 ) {
			continue;
		    }

		    tmpSampleVecs[2 * successfulTensorFits] = new Vector3D(eig1[1][0],
									   eig1[2][0],
									   eig1[3][0]);
		    
		    tmpSampleVecs[2 * successfulTensorFits + 1] = new Vector3D(eig2[1][0],
									       eig2[2][0],
									       eig2[3][0]);
		    successfulTensorFits++;
		}

		if ( ((double)successfulTensorFits) / ((double)samples) < 0.95) {
		    double percentageFailed = ( 1.0 - ((double)successfulTensorFits) / 
						((double)samples) ) * 100.0;

		    System.err.println("Warning: discarding " + percentageFailed + 
				       "% of samples at point : (" + x + ", " + y + ")");
		}


		Vector3D[] sampleVecs = new Vector3D[2 * successfulTensorFits];


		

		for (int i = 0; i < 2 * successfulTensorFits; i++) {
			sampleVecs[i] = tmpSampleVecs[i];
                }
		
		if (watson) {
		    kappaLUT[WATSON][xCounter][yCounter] = 
			getWatsonConcentrationParams(sampleVecs, e1Vec, dt2E1Vec);
		}
		else {
		    kappaLUT[WATSON][xCounter][yCounter] = new double[2];
		}

		if (bingham) {
		    kappaLUT[BINGHAM][xCounter][yCounter] = 
			getBinghamConcentrationParams(sampleVecs, e1Vec, dt2E1Vec);
		}
		else {
		     kappaLUT[BINGHAM][xCounter][yCounter] = new double[4];
		}

		if (acg) {
		    kappaLUT[ACG][xCounter][yCounter] = 
			getACGConcentrationParams(sampleVecs, e1Vec, dt2E1Vec);
		}
		else {
		    kappaLUT[ACG][xCounter][yCounter] = new double[6];
		}

	
		if (prop == 0.5 && xCounter != yCounter) {

		    double[] tmp = new double[2];

		    tmp[0] = kappaLUT[WATSON][xCounter][yCounter][1];
		    tmp[1] = kappaLUT[WATSON][xCounter][yCounter][0];

		    kappaLUT[WATSON][yCounter][xCounter] = tmp;

		    tmp = new double[4];

		    tmp[0] = kappaLUT[BINGHAM][xCounter][yCounter][2];
		    tmp[1] = kappaLUT[BINGHAM][xCounter][yCounter][3];
	    	    tmp[2] = kappaLUT[BINGHAM][xCounter][yCounter][0];
		    tmp[3] = kappaLUT[BINGHAM][xCounter][yCounter][1];
		    
		    kappaLUT[BINGHAM][yCounter][xCounter] = tmp;

		    tmp = new double[6];

		    tmp[0] = kappaLUT[ACG][xCounter][yCounter][3];
		    tmp[1] = kappaLUT[ACG][xCounter][yCounter][4];
	    	    tmp[2] = kappaLUT[ACG][xCounter][yCounter][5];
		    tmp[3] = kappaLUT[ACG][xCounter][yCounter][0];
		    tmp[4] = kappaLUT[ACG][xCounter][yCounter][1];
	    	    tmp[5] = kappaLUT[ACG][xCounter][yCounter][2];

		    kappaLUT[ACG][yCounter][xCounter] = tmp;

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


    /**
     * @return <code>{k1, k2}</code>, where <code>k1</code> comes from the distribution 
     * aligned closest to <code>e1Vec</code>.
     *
     */
    protected double[] getWatsonConcentrationParams(Vector3D[] sampleVecs, Vector3D e1Vec, 
					      Vector3D dt2E1Vec) {
	

	Vector3D[][] sortedSamples = sortAxes(sampleVecs, e1Vec, dt2E1Vec);

	EigenSystem3D[] tBar = new EigenSystem3D[2];

	double[] k = new double[2];
	Vector3D[] mu = new Vector3D[2];

	for (int i = 0; i < 2; i++) {

	    tBar[i] = WatsonFitter.tBarEigenSystem(sortedSamples[i]);

	    k[i] = WatsonFitter.fitKappa(tBar[i], sortedSamples[i]);
	    
	    if (k[i] > 0.0) {
		mu[i] = tBar[i].eigenvectors[0];
	    }
	    else {
		mu[i] = tBar[i].eigenvectors[2];
	    }
	}

        int mostConcentratedMu = k[0] > k[1] ? 0 : 1;

	double poss1 = Math.abs(mu[mostConcentratedMu].dot(e1Vec));
	double poss2 = Math.abs(mu[mostConcentratedMu].dot(dt2E1Vec));

	if (poss1 > poss2) {
            return new double[] {k[mostConcentratedMu], k[1 - mostConcentratedMu]};
	}
	else {
            return new double[] {k[1 - mostConcentratedMu], k[mostConcentratedMu]};
        }
	
    }
    

    /**
     * @return <code>{k11, k12, k21, k22}</code>, where <code>k1</code> and <code>k21</code> 
     * come from the distribution aligned closest to <code>e1Vec</code>.
     *
     */
   protected double[] getBinghamConcentrationParams(Vector3D[] sampleVecs, Vector3D e1Vec, 
					      Vector3D dt2E1Vec) {

	Vector3D[][] sortedSamples = sortAxes(sampleVecs, e1Vec, dt2E1Vec);

	EigenSystem3D[] tBar = new EigenSystem3D[2];

	double[] ks = new double[4];
	Vector3D[] mus = new Vector3D[6];

	for (int i = 0; i < 2; i++) {

	    tBar[i] = BinghamFitter.tBarEigenSystem(sortedSamples[i]);

	    mus[3*i] = tBar[i].eigenvectors[0];
	    mus[3*i + 1] = tBar[i].eigenvectors[1];
	    mus[3*i + 2] = tBar[i].eigenvectors[2];

	    try {
		double[] bngpar = BinghamFitter.bngpar(tBar[i].eigenvalues[1], tBar[i].eigenvalues[2]);
		
		ks[2*i] = bngpar[0];
		ks[2*i+1] = bngpar[1];
	    }
	    catch (ConvergenceException e) {
		// do nothing, concentration is 0
	    }
	}

        int mostConcentratedMu = Math.abs(ks[0] + ks[1]) > Math.abs(ks[2] + ks[3]) ? 0 : 1;


        // m11 (of the Bingham distribution) should be normal to e1Vec (of the tensor)
        // according to eigenvector numbering scheme of Bingham
	double poss1 = Math.abs(mus[3 * mostConcentratedMu].dot(e1Vec));
	double poss2 = Math.abs(mus[3 * mostConcentratedMu].dot(dt2E1Vec));

	if (poss1 > poss2) {
	    return new double[] {ks[2 * mostConcentratedMu], ks[2 * mostConcentratedMu + 1], 
                                 ks[2 * (1 - mostConcentratedMu)], ks[2 * (1 - mostConcentratedMu) + 1]};
	}
	else {
	    return new double[] {ks[2 * (1 - mostConcentratedMu)], ks[2 * (1 - mostConcentratedMu) + 1],
                                 ks[2 * mostConcentratedMu], ks[2 * mostConcentratedMu + 1]};

	}

    }

    /**
     * @return <code>{sigma11, sigma21, sigma31, sigma12, sigma22, sigma32}</code>, 
     * where the first three parameters come from the distribution aligned closest to 
     * <code>e1Vec</code>.
     *
     */
    protected double[] getACGConcentrationParams(Vector3D[] sampleVecs, Vector3D e1Vec, 
					      Vector3D dt2E1Vec) {


	double[] sigmas = new double[6];
	Vector3D[] e1s = new Vector3D[2];

	Vector3D[][] sortedSamples = sortAxes(sampleVecs, e1Vec, dt2E1Vec);

	
	for (int i = 0; i < 2; i++) {	
	    
	    RealMatrix A = ACG_Fitter.findA(sortedSamples[i]);
	    
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

	    sigmas[3*i] = eigA.eigenvalues[0];
	    sigmas[3*i+1] = eigA.eigenvalues[1];
	    sigmas[3*i+2] = eigA.eigenvalues[2];
	    
	    e1s[i] = eigA.eigenvectors[0];
	}
	    
	
        int mostConcentratedMu = (sigmas[0] + sigmas[1] + sigmas[2]) < (sigmas[3] + sigmas[4] + sigmas[5]) ? 0 : 1;
        
	double poss1 = Math.abs(e1s[mostConcentratedMu].dot(e1Vec));
	double poss2 = Math.abs(e1s[mostConcentratedMu].dot(dt2E1Vec));
        
	if (poss1 > poss2) {
            return new double[] {sigmas[3*mostConcentratedMu], sigmas[3*mostConcentratedMu+1],
                                 sigmas[3*mostConcentratedMu+2],
                                 sigmas[3*(1-mostConcentratedMu)], sigmas[3*(1-mostConcentratedMu)+1],
                                 sigmas[3*(1-mostConcentratedMu)+2]};
        }
	else {
            return new double[] {sigmas[3*(1-mostConcentratedMu)], sigmas[3*(1-mostConcentratedMu)+1],
                                 sigmas[3*(1-mostConcentratedMu)+2],
                                 sigmas[3*mostConcentratedMu], sigmas[3*mostConcentratedMu+1],
                                 sigmas[3*mostConcentratedMu+2]};
	}
            
        
    }


    private static Vector3D[][] sortAxes(Vector3D[] samples, Vector3D e1Vec, Vector3D dt2E1Vec) {
	
	int numSamples = samples.length / 2;

	Vector3D[] dt1Samples = new Vector3D[numSamples];
	Vector3D[] dt2Samples = new Vector3D[numSamples];

	for (int i = 0; i < numSamples; i++) {
	    Vector3D sample1 = samples[2*i];
	    Vector3D sample2 = samples[2*i+1];

	    // maximise mu.dot(sample)

	    double poss1 = Math.abs(sample1.dot(e1Vec)) + Math.abs(sample2.dot(dt2E1Vec));
	    double poss2 = Math.abs(sample2.dot(e1Vec)) + Math.abs(sample1.dot(dt2E1Vec));

	    if (poss1 > poss2) {
		dt1Samples[i] = sample1;
		dt2Samples[i] = sample2;
	    }
	    else {
		dt1Samples[i] = sample2;
		dt2Samples[i] = sample1;
	    }

	}

	EigenSystem3D dt1Eig = WatsonFitter.tBarEigenSystem(dt1Samples);
	EigenSystem3D dt2Eig = WatsonFitter.tBarEigenSystem(dt2Samples);


	Vector3D estimatedMu1 = dt1Eig.eigenvectors[0];
	Vector3D estimatedMu2 = dt2Eig.eigenvectors[0];

	boolean swapped = true;
	int counter = 0;
	while (swapped && counter < 100) {
	    swapped = false;
	    for (int i = 0; i < numSamples; i++) {
		if (Math.abs(dt1Samples[i].dot(estimatedMu1)) <
		    Math.abs(dt1Samples[i].dot(estimatedMu2))) {
		    Vector3D tmp = dt1Samples[i];
		    dt1Samples[i] = dt2Samples[i];
		    dt2Samples[i] = tmp;
		    swapped = true;
		}
	    }

	    counter++;

	    dt1Eig = WatsonFitter.tBarEigenSystem(dt1Samples);
	    dt2Eig = WatsonFitter.tBarEigenSystem(dt2Samples);

	    estimatedMu1 = dt1Eig.eigenvectors[0];
	    estimatedMu2 = dt2Eig.eigenvectors[0];
	}


	return new Vector3D[][] {dt1Samples, dt2Samples};
    }



}				   
	
				   







    

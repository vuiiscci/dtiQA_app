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
import java.util.Random;
import java.util.Arrays;


/**
 * Fits some distribution to samples of tensors with added noise. 
 *
 * @author Philip Cook
 * @version $Id$
 */
public abstract class DT_LookupTableGenerator {
   
    protected Vector3D e1Vec;
    protected Vector3D e2Vec;
    protected Vector3D e3Vec;

    // eigenvectors
    private Matrix e1;
    private Matrix e2;
    private Matrix e3;
    
    private Matrix e1T;
    private Matrix e2T;
    private Matrix e3T;
 
    private Matrix minusE1;
    private Matrix minusE2;
    private Matrix minusE3;
    
    private Matrix minusE1T;
    private Matrix minusE2T;
    private Matrix minusE3T;


    protected double snr;

    protected double trace;

    protected DW_Scheme imPars;
    protected Random ran;
    
    // indices to LUTs returned by generateLUT
    public static final int WATSON = 0;
    public static final int BINGHAM = 1;
    public static final int ACG = 2;


    /**
     * 
     *
     * @param ip imaging parameters for the synthesis of data.
     * @param snr signal to noise ratio in the q=0 image.
     * @param trace trace of the diffusion tensors used as test functions.
     * @param r a random number generator.
     *
     */
    protected DT_LookupTableGenerator(DW_Scheme ip, double snr, double trace, 
				       Random r) {
	imPars = ip;
	this.snr = snr;
	this.trace = trace;
	ran = r;
	initializeEigenvectors(ran);
    }



    /**
     * Get a prolate axi-symmetric tensor.
     *
     * @param fa the fractional anisotropy.
     * @param trace the total trace of the tensor.
     */    
    protected DT getTensor(double fa, double trace, Random ran) {

	if (fa < 0.0 || fa > 1.0) {
	    throw new IllegalArgumentException("Require 0.0 <= FA <= 1.0");
	}


	// from Mathematica, solving the FA equation for L1, using L2 = L3 = (Trace - L1) / 2

	double t = trace;

	double l1 = (-3.0 * t + 2 * fa * fa * t - 2.0 * 
		     Math.sqrt(3.0 * fa * fa * t * t - 2 * fa * fa * fa * fa * t * t))
	    /(3.0 * (-3.0 + 2.0 * fa * fa));

	double l2 = (t - l1) / 2.0;
	double l3 = l2;

	return getFullySpecifiedTensor(l1, l2, l3, ran);
    }


    /**
     * Get a diffusion tensor, oriented according to the e{1-3}Vecs, with anisotropy described by x and y.
     * @param x L1 / L3
     * @param y L2 / L3
     * @param trace the total trace of the tensor
     * @param ran used to decide sign of eigenvectors.
     */    
    protected DT getTensor(double x, double y, double trace, Random ran) {
	
	// L1 + L2 + L3 == trace
	// x * L3 + y * L3 + L3 = trace
	// (x + y + 1)L3 = trace
	
	double l3 = trace / (x + y + 1.0);
	
	double l2 = l3 * y;  
	double l1 = l3 * x;  

	return getFullySpecifiedTensor(l1, l2, l3, ran);
    }


    
    /**
     * Get a diffusion tensor, with specified eigenvalues and orientation according to the e{1-3}Vecs.
     */    
    protected DT getFullySpecifiedTensor(double l1, double l2, double l3, Random ran) {

	// D == l1 * (e1 * e1T) + l2 * (e2 * e2T) + l3 * (e3 * e3T)
	
	Matrix e1e1T = null;
	Matrix e2e2T = null;
	Matrix e3e3T = null;
	
	Matrix d;
	
	// Each tensor has random signed eigenvectors
	
	if (ran.nextDouble() < 0.5) {
	    
	    e1e1T = minusE1.times(minusE1T);
	    
	}
	else {
	    
	    e1e1T = e1.times(e1T);
	    
	}	    
	if (ran.nextDouble() < 0.5) {
	    
	    e2e2T = minusE2.times(minusE2T);
	    
	}
	else {
	    
	    e2e2T = e2.times(e2T);
	    
	}	    
	if (ran.nextDouble() < 0.5) {
	    
	    e3e3T = minusE3.times(minusE3T);
	    
	}
	else {
	    
	    e3e3T = e3.times(e3T);
	    
	}
	
	// note use of timesEquals, avoids object creation but changes e1e1T
	d = e1e1T.timesEquals(l1).plusEquals(e2e2T.timesEquals(l2).plusEquals(e3e3T.timesEquals(l3)) );
	
	double[][] a = d.getArray();
	
	DT tensor = new DT(a[0][0], a[0][1], a[0][2], a[1][1], a[1][2], a[2][2]);
	
	return tensor;
    }


    /**
     * Get a diffusion tensor, with specified FA and orientation according to the e{1-3} vectors.
     */    
    protected DT getTensor(Vector3D e1, Vector3D e2, Vector3D e3, double fa, double trace, Random ran) {

        if (fa == 0.0) {
            return getTensor(e1, e2, e3, 1.0, 1.0, trace, ran);
        }
        else if (fa > 0.0 && fa < 1.0) {
            double t = trace;
            
            double l1 = (-3.0 * t + 2 * fa * fa * t - 2.0
                         * Math.sqrt(3.0 * fa * fa * t * t - 2 * fa * fa * fa * fa * t * t))
                /(3.0 * (-3.0 + 2.0 * fa * fa));
            
            double l2 = (t - l1) / 2.0;
            
            return getTensor(e1, e2, e3, l1 / l2, 1.0, trace, ran);

        }
        else {
            throw new IllegalArgumentException("Require 0 <= FA < 1");
        }

	
    }

    /**
     * Get a diffusion tensor, with specified eigenvalues and orientation 
     * according to the e{1-3} vectors.
     */    
    protected DT getTensor(Vector3D e1, Vector3D e2, Vector3D e3, double x, double y, 
			 double trace, Random ran) {

	double l3 = trace / (x + y + 1.0);
	
	double l2 = l3 * y;  
	double l1 = l3 * x;  

	// D == l1 * (e1 * e1T) + l2 * (e2 * e2T) + l3 * (e3 * e3T)
	
	Matrix e1e1Tm = null;
	Matrix e2e2Tm = null;
	Matrix e3e3Tm = null;

	Matrix e1m = e1.toJamaMatrix();
	Matrix e2m = e2.toJamaMatrix();
	Matrix e3m = e3.toJamaMatrix();
	
	Matrix d;
	
	// Each tensor has random signed eigenvectors
	
	if (ran.nextDouble() < 0.5) {
	    e1e1Tm = e1m.times(e1m.transpose());
	}
	else {
	    e1e1Tm = e1m.timesEquals(-1.0).times(e1m.transpose());
	}	    

	if (ran.nextDouble() < 0.5) {
	    e2e2Tm = e2m.times(e2m.transpose());
	}
	else {
	    e2e2Tm = e2m.timesEquals(-1.0).times(e2m.transpose());
	}	    

	if (ran.nextDouble() < 0.5) {
	    e3e3Tm = e3m.times(e3m.transpose());
	}
	else {
	    e3e3Tm = e3m.timesEquals(-1.0).times(e3m.transpose());
	}	    

	
	// note use of timesEquals, avoids object creation but changes e1e1T
	d = e1e1Tm.timesEquals(l1).plusEquals(e2e2Tm.timesEquals(l2).plusEquals(e3e3Tm.timesEquals(l3)) );
	
	double[][] a = d.getArray();
	
	DT tensor = new DT(a[0][0], a[0][1], a[0][2], a[1][1], a[1][2], a[2][2]);
	
	return tensor;
    }


  

    /**
     * Get noisy tensors from a one-tensor or two-tensor GaussianMixture. If fitting fails, 
     * the array may contain nulls.
     * 
     * @param noiseFreeTensors should be one or two tensors from a Gaussian mixture.
     * @return an array noisyTensors which is of type and dimensions 
     * <code>DT[noiseFreeTensors.length][samples]</code>. 
     */    
    protected static DT[][] getNoisyTensors(DT[] noiseFreeTensors, double[] mix, 
					    ModelIndex inversionIndex, DW_Scheme imPars, double snr, 
					    int samples, Random ran) {
	
	
	GaussianMixture g = new GaussianMixture(noiseFreeTensors, mix);

	DataSynthesizer source = new DataSynthesizer(g, imPars, snr, samples, ran);
		
	DT[][] dts = new DT[noiseFreeTensors.length][samples];

	if (noiseFreeTensors.length == 1) {

	    DT_Inversion inv = DT_Inversion.getIndexedDT_Inversion(inversionIndex, imPars);
	    
	    for (int i = 0; i < samples; i++) {
		try {
		    double[] invertedData = inv.invert(source.nextVoxel());

		    if (invertedData[0] == 0.0) {
			dts[0][i] = new DT(invertedData[2], invertedData[3], invertedData[4], 
					   invertedData[5], invertedData[6], invertedData[7]);
		    }

		}
		catch(DataSourceException e) {
		    System.out.println(e);
		}
	    }
	    
	}
	else { // two tensors
	    
	    // Choose the inversion to run.
	    TwoTensorInversion inv = new TwoTensorInversion(imPars, inversionIndex, 
							    ModelIndex.NLDT_POS);
		
	    for (int i = 0; i < samples; i++) {
		    
		double[] fittedData = null;
		    
		try {
		    fittedData = inv.invert(source.nextVoxel());
		} catch(Exception e) {
		    System.err.println(e);
		    continue;
		}
		    
		if (fittedData[0] != 0.0) {
		    // fitting failed
		    System.err.println("fittedData[0]: " + fittedData[0]);
		    continue;
		}
		    
		dts[0][i] = new DT(fittedData[4], fittedData[5], fittedData[6], fittedData[7], fittedData[8], fittedData[9]);
		
		dts[1][i] = new DT(fittedData[11], fittedData[12], fittedData[13], fittedData[14], fittedData[15], fittedData[16]);
		
	    }
	
	}


	return dts;

    }



	
    // called when e1, e2, e3 change
    protected void initializeEigenvectors(RealMatrix r) {
	
	e1Vec = new Vector3D(1.0, 0.0, 0.0);
	e2Vec = new Vector3D(0.0, 1.0, 0.0);
	e3Vec = new Vector3D(0.0, 0.0, 1.0);

	e1Vec = Rotations.rotateVector(e1Vec, r);
	e2Vec = Rotations.rotateVector(e2Vec, r);
	e3Vec = Rotations.rotateVector(e3Vec, r);

	// eigenvectors
	e1 = e1Vec.toJamaMatrix();
	e2 = e2Vec.toJamaMatrix();
	e3 = e3Vec.toJamaMatrix();
	
	e1T = e1.transpose();
	e2T = e2.transpose();
	e3T = e3.transpose();
	
	minusE1 = e1.times(-1.0);
	minusE2 = e2.times(-1.0);
	minusE3 = e3.times(-1.0);
	
	minusE1T = minusE1.transpose();
	minusE2T = minusE2.transpose();
	minusE3T = minusE3.transpose();
	
    }

    // called when e1, e2, e3 change
    // more random 
    protected void initializeEigenvectors(Random ran) {
	
	double theta1 = Math.acos(2.0 * ran.nextDouble() - 1.0);
	double phi1 = Math.PI * 2.0 * ran.nextDouble();
	
	e1Vec = Vector3D.vectorFromSPC(1.0, theta1, phi1);

	// get axis orthogonal to e1
	e2Vec = Rotations.rotateVector(Rotations.X_AXIS, Rotations.Y_AXIS, theta1);
	e2Vec = Rotations.rotateVector(e2Vec, Rotations.Z_AXIS, phi1);

	// now rotate about e1
	e2Vec = Rotations.rotateVector(e2Vec, e1Vec, Math.PI * 2.0 * ran.nextDouble());

	e3Vec = e1Vec.cross(e2Vec);

	// eigenvectors
	e1 = e1Vec.toJamaMatrix();
	e2 = e2Vec.toJamaMatrix();
	e3 = e3Vec.toJamaMatrix();
	
	e1T = e1.transpose();
	e2T = e2.transpose();
	e3T = e3.transpose();
	
	minusE1 = e1.times(-1.0);
	minusE2 = e2.times(-1.0);
	minusE3 = e3.times(-1.0);
	
	minusE1T = minusE1.transpose();
	minusE2T = minusE2.transpose();
	minusE3T = minusE3.transpose();
	
    }




    /**
     * Generates the LUT. See subclasses for definition of indices of LUT.
     * 
     */
    public abstract double[][][][] generateLUT(double xMin, double xMax, double step, int samples, 
					     ModelIndex inversionIndex,
					     boolean watson, boolean bingham, boolean acg);


}				   
	
				   







    

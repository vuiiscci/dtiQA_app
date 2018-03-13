package data;

import java.io.*;

import imaging.*;
import misc.*;
import numerics.*;
import tools.*;

import java.util.Random;


/**
 * Provides wild bootstrap data based on diffusion tensor reconstruction. 
 *
 * For information on the technique, see B Whitcher et al, 
 * "Using the wild bootstrap to quantify uncertainty in diffusion tensor imaging", Human Brain Mapping
 * (in press)
 * 
 *
 * 
 * @author Philip Cook
 * @version $Id$
 * 
 */
public class WildBS_DT_DataSynth extends WildBootstrapDataSynth {


    /**
     * 
     * @param input
     *            A data source that should provide voxel-order data.
     *
     * @param scheme
     *            Imaging scheme for the data.
     * 
     *
     * @param samples
     *            The number of bootstrap samples to generate from each voxel.
     * 
     */
    public WildBS_DT_DataSynth(DataSource input, DW_Scheme scheme, int samples) {

	this(input, scheme, samples, new numerics.MTRandom(2350));
    }



    /**
     * 
     * @param input
     *            A data source that should provide voxel-order data.
     *
     * @param scheme
     *            Imaging scheme for the data.
     *
     * @param samples
     *            The number of bootstrap samples to generate from each voxel.
     *
     * @param seed 
     *            Seed for the random number generator.
     * 
     */
    public WildBS_DT_DataSynth(DataSource input, DW_Scheme scheme, int samples, int seed) {

	this(input, scheme, samples, new numerics.MTRandom(seed));
    }



    /**
     * 
     * @param input
     *            A data source that should provide voxel-order data.
     * 
     * @param scheme
     *            Imaging scheme for the data.
     *
     * @param samples
     *            The number of bootstrap samples to generate from each voxel.
     * 
     * @param ran 
     *          A random number generator.
     */
    public WildBS_DT_DataSynth(DataSource input, DW_Scheme scheme, int samples, Random ran) {
	super(input, scheme, samples, ran);

	// tell superclass to bootstrap log data since DT is linear on log measurements
	useLogData();

	init();
    }




    /**
     * Initialize X, H, linearInv.
     *
     */
    protected void initializeReconSpecificParams() {

	// Flachaire, Chung et al [NeuroImage 33:531-541 (2006)] do things a little differently
	// than Camino inverters.
	// Given y = X \beta, they calculate \beta = (X^T X)^{-1} X^T y
	
	// since  X^T y = X^T X \beta, and then (X^T X)^{-1} X^T y = (X^T X)^{-1} X^T X \beta = \beta

	// (X^T X)^{-1} is inverted directly.


	// Camino computes pseudoinverse X^* and then takes \beta = X^* y.

	// We follow the procedure that Whitcher appears to follow, which is not
	// spelled out in his paper but inferred from Flachaire, ie we use H = X (X^T X)^{-1} X^T and 
	// set linearInv = (X^T X)^{-1} X^T, so H = X * linearInv


	// Set up the matrix X.
        X = ip.getB_Matrix();

	// now get H = X (X^TX)^{-1} X^T

	RealMatrix XTX = X.transpose().product(X);

	RealMatrix XTX_INV = null;

	try {
	    XTX_INV = XTX.inverse();
        }
        catch (Exception e) {
            throw new LoggedException(e);
        }


	linearInv = XTX_INV.product(X.transpose());
	
	H = X.product(linearInv);

    }


   
}

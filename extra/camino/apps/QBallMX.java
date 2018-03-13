package apps;

import java.util.logging.Logger;
import tools.*;
import data.*;
import misc.*;
import imaging.*;
import sphfunc.*;
import inverters.*;
import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Generates the transition matrix for QBall
 * 
 * <dt>Description:  The default option uses the method 
 * descibed by Tuch (MRM, 2004) to generate
 * the transition matrix for QBall.  Alternatively, this method
 * can generate the transition matrix using the method described
 * by Alexander (Annals of the New York Academy of Sciences, 2005). 
 * 
 * <dd> QBall.
 * 
 * </dl>
 * 
 * @author Kiran Seunarine
 * @version $Id$
 *  
 */
public class QBallMX extends Executable{

	public QBallMX(String[] args){
		super(args);
	}
	
    /**
     * Logger object
     */
    private static Logger logger = Logger.getLogger("apps/QBallMX");


    /**
     * The width of the radial basis functions used to smooth the data
     */
    private static double dataSmoothSigma;
	//basis functions
	private LinearBasisSum lbs;
	private LinearBasisFunction lbf;
	private	DW_Scheme imp;
	private int basis;
	private double[][] normQs;
	private int numQPoints;
	
	private double [] [] odfBasisCentres;
	private int numOdfBasisPoints;
	private RealMatrix yMat;


    /**
     * Output manager
     */
    //private static OutputManager om;

	public void initDefaultVals(){
	    dataSmoothSigma = 0.1308996938995747; // 7.5 degrees

	}
	
	public void initOptions(String[] args){	
//    public static void main(String[] args) {	
        // Parse the command line arguments
        CL_Initializer.CL_init(args);

	// parse technique specific args
		for(int i=0; i<args.length; i++) {
			if(args[i].equals("-smoothingsigma")) {
				dataSmoothSigma = Double.parseDouble(args[i + 1]);
				CL_Initializer.markAsParsed(i, 2);
			}
		}
		
        CL_Initializer.checkParsing(args);
        CL_Initializer.initImagingScheme();
	}

	// om = new OutputManager();

	public void initVariables(){
		// Find out what type of basis function to use
		basis = CL_Initializer.basisType;

		// Create the basis functions
		lbs = BasisSumFactory.getBasisSum(basis);
//		LinearBasisFunction lbf = null;
		logger.info("BASIS FUNCTION SETTINGS:");
		logger.info(lbs.getSettings());

		// Get gradient points from Scheme
		imp = CL_Initializer.imPars;

		// use gradient dirs in place of normed Qs. The former was not 
		// consistent across scheme formats
		normQs = imp.getNonZeroG_Dirs();
		numQPoints = normQs.length;

		// Get the ODF sample points
		odfBasisCentres = SphericalPoints.getElecPointSet(lbs.numBasisFunctions());
		numOdfBasisPoints = odfBasisCentres.length;

		// Calculate the matrix Y (basis functions, centred on odfBasisFunctionCentres)
		yMat = new RealMatrix(numQPoints, lbs.numBasisFunctions());
	}

	public void execute(OutputManager om){

		for(int i=0; i < lbs.numBasisFunctions(); i++) {
			lbf = lbs.basisFunction(i);

			for(int j=0; j< numQPoints; j++) {
				yMat.setEntry(j, i, lbf.getRadius(normQs[j][0],normQs[j][1],normQs[j][2]));
			}
		}

		// Calculate the matrix Phi (the great circle integrals)
		RealMatrix phiMat = new RealMatrix(numOdfBasisPoints, lbs.numBasisFunctions());
		for(int i=0; i< lbs.numBasisFunctions(); i++) {
			lbf =lbs.basisFunction(i);
	    
			for(int j=0; j< numOdfBasisPoints ; j++) 
			{
				phiMat.setEntry(j, i, lbf.greatCircleIntegral(odfBasisCentres[j]));
			}
		}

		// the output matrix
		RealMatrix odfMat;

	/*
	 * The matrix Theta and the calculation of the ODF depend on whether
	 * Tuch's method or Alexander's method is used.
	 * Here, Tuch's method is used when rbf's are specified as the basis,
	 * otherwise Alexander's method is used.
	 */
	if(basis != BasisSumFactory.TUCH_RBF) {
	    // Calculate the matrix Theta
	    RealMatrix thetaMat = new RealMatrix(numOdfBasisPoints, lbs.numBasisFunctions());
	    for(int i=0; i< lbs.numBasisFunctions(); i++) {
		lbf = lbs.basisFunction(i);

		for(int j=0; j< numOdfBasisPoints ; j++) {
		    thetaMat.setEntry(j, i, lbf.getRadius(odfBasisCentres[j][0],odfBasisCentres[j][1],odfBasisCentres[j][2]));
		}
	    }

	    // Get the linear transformation matrix
	    RealMatrix thetaMat_Inv = thetaMat.pseudoInv();

	    if(basis == BasisSumFactory.TUCH_RBF) // never true!
	    {
		thetaMat_Inv = tuchsNormalization(thetaMat_Inv, numOdfBasisPoints);
	    }

	    RealMatrix yMat_Inv = yMat.pseudoInv();
	    odfMat = thetaMat_Inv.product(phiMat);
	    odfMat = odfMat.product(yMat_Inv);
	}
	else {
	    // Get new Basis Sum so that theta works
	    RBF_Sum.setPoints(normQs);
	    TuchRBF_Sum.setSigma(dataSmoothSigma);
	    lbs = BasisSumFactory.getBasisSum(basis);

	    logger.info("SMOOTHING FUNCTION SETTINGS:");
	    logger.info(lbs.getSettings());

	    // Calculate the matrix Theta
	    RealMatrix thetaMat = new RealMatrix(numQPoints, lbs.numBasisFunctions());
	    for(int i=0; i< lbs.numBasisFunctions(); i++) {
		lbf = lbs.basisFunction(i);

	   	for(int j=0; j< numQPoints ; j++) {
		    thetaMat.setEntry(j, i, lbf.getRadius(normQs[j][0],normQs[j][1],normQs[j][2]));
		}
	    }
	    	
	    // Get the matrix
	    yMat = yMat.pseudoInv();
	    thetaMat = tuchsNormalization(thetaMat, numQPoints);
	    odfMat = phiMat.product(yMat);
	    odfMat = odfMat.product(thetaMat);
	}

       	// Output the matrix
       	double[] row = new double[odfMat.columns()];
        for (int j = 0; j < odfMat.rows(); j++) {
	    for (int i = 0; i < odfMat.columns(); i++) {
		row[i] = odfMat.entries[j][i];
	    }
	    om.output(row);
	}

        om.close();
    }

    

    /** 
     * Normalization used by Tuch
     * @param M the matrix to normalize
     * @param matLength the length of the matrix (it is assumed the matrix is square)
     */
    private static RealMatrix tuchsNormalization(RealMatrix M, int matLength) {

	double [] thetaSum = new double [matLength]; 
	for(int j=0; j< matLength ; j++) {
	    thetaSum [j] = 0;
	}
	//adjust thetaMat
	for(int i=0; i< matLength; i++) {
	    for(int j=0; j< matLength ; j++) {
		thetaSum [i] += M.entry(i, j);
	    }
	}
		
	RealMatrix tempMat = new RealMatrix(matLength, matLength);
	double tempVal = 0;
	for(int i=0; i< matLength; i++) {
	    for(int j=0; j< matLength ; j++) {
		tempVal = M.entry(i,j)/thetaSum[j];
		tempMat.setEntry(j, i,tempVal);
	    }
	}

	return tempMat;
    }

}

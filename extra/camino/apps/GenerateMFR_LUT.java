package apps;
import java.util.logging.Logger;
import misc.*;
import numerics.*;
import tractography.Histogram;
import tractography.HistBin;
import tools.*;
import data.*;
import apps.*;

import java.io.*;
import java.util.Random;
import java.util.Set;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Iterator;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd> Generates a LUT that is indexed using the Hessian of the peaks of
 * multi-fibre methods.  Calibration for PICo.
 * 
 * <dt>Description:  Generates a lookup table.
 *  
 * 
 * <dd> blah.
 * 
 * </dl>
 * 
 * @author Kiran Seunarine
 * @version $Id$
 *  
 */
public class GenerateMFR_LUT {
    
    /**
     * Logger object
     */
    private static Logger logger = Logger.getLogger("camino.GenerateMFR_LUT");

    /**
     * Distribution types
     */
    private static final String BINGHAM = "bingham";
    private static final String WATSON = "watson";
    private static final String ACG = "acg";
    
    /**
     * test data flag
     */
    private static boolean testDataOutput = false;

    /**
     * Flag data!
    */
    private static int order = 2;
    private static int numFOEsPerVoxel = 3;
    private static double binSize = 1;
    private static int numVectsPerBin = 50;
    private static boolean useLog = true;
    private static Histogram hist_1Fibre, hist_2Fibre;

    /**
     * Calibration Data Settings (version 1.0)
     */
    private static  Random rand;
    private static int oneDT_BlockSize = 1;
    private static double minDT2RotAngle = 0.0;
    private static double maxDT2RotAngle = 0.5;
    private static double dt2RotAngleStep = 0.1;
    private static int twoDT_BlockSize = 1;
    private static int numTwoDT_Blocks = 1;
    private static double dt1E1Theta = 0.0;
    private static double dt1E1Phi = 0.0;
    private static double dt2E1Theta = 0.0;
    private static double dt2E1Phi  = 0.0;
    private static double dt2RotAxisTheta = 0.0;
    private static double dt2RotAxisPhi = 0.0;

    private static String outputStem = "PICoCalibOut";
    private static String distributionType = BINGHAM;

    /*
     * Reads the info file created by SphFuncPICoCalibrationData.java
     * The data in this file must be in the form FLAG \t [value]
     */
    private static void readInfoFile(String infoFilename){
	
	// try and open the file
	FileInput fInData = new FileInput(infoFilename);

	// read first line of data, which should be the version information
	String line =  fInData.readString();
	String flag = line.substring(0,line.indexOf("\t")).trim();
	String value = line.substring(line.indexOf("\t")+1,line.length());
	
	double version = 0.0;

	/* if the first line contains the version information, process the info file using the flags for that version,
	 * else assume that the info file is version one
	 */
	if(flag.equals("VERSION")) {
	    version = Double.parseDouble(value);

	    // read the next line ready for processing (assuming that there is more than the version information in the file!)
	    line = fInData.readString();
	    flag = line.substring(0,line.indexOf("\t")).trim();
	    value = line.substring(line.indexOf("\t")+1,line.length());
        
	}
	else {
		// if no version information found, assume version 1.0
	    logger.warning("no version information found in info file.  Assuming version 1.0. flag: "+flag);
		version = 1.0;
	}

	// parse the rest of the arguments
	boolean done = false;
	while(!done) {

	    // process the data using version 1 flags
	    processVersion1flags(flag, value);

	    // if not the end of the file, read the next line of data
	    if(!fInData.eof()) {
		line = fInData.readString();
			
		if(line != null){
		    flag = line.substring(0,line.indexOf("\t")).trim();
		    value = line.substring(line.indexOf("\t")+1,line.length());
		}
            
	    }
	    else {
		done = true;
		fInData.close();
	    }
	}
    }
   
    /*
     * Processes the flags for info files that are either VERSION 1.0 or have no version information
     */
    private static void processVersion1flags(String flag, String value)
    {
	if(flag.equals("ROTATION_SEED")) {
	    //create a new Random object, seeded using the specified integer
	    rand = new Random(Integer.parseInt(value));
	    if(testDataOutput)
		System.err.println("rotation seed: " + value);
	}
       	else if(flag.equals("ONE_DT_BLOCK_SIZE")) {
       	    oneDT_BlockSize = Integer.parseInt(value);
	    if(testDataOutput)
		System.err.println("one block size: " + value);
	}
		else if(flag.equals("MIN_DT2_ROT_ANGLE")) {
			minDT2RotAngle = Double.parseDouble(value);
			if(testDataOutput)
			    System.err.println("min dt2 rot angle: " + value);
		}
		else if(flag.equals("MAX_DT2_ROT_ANGLE")) {
			maxDT2RotAngle = Double.parseDouble(value);
			if(testDataOutput)
			    System.err.println("max dt2 rot angle: " + value);
		}
		else if(flag.equals("DT2_ROT_ANGLE_STEP")) {
			dt2RotAngleStep = Double.parseDouble(value);
			if(testDataOutput)
			    System.err.println("dt2 rot angle step: " + value);
		}
		else if(flag.equals("TWO_DT_BLOCK_SIZE")) {
			twoDT_BlockSize = Integer.parseInt(value);
			if(testDataOutput)
			    System.err.println("dt2 block size: " + value);
		}
		else if(flag.equals("NUM_TWO_DT_BLOCKS")) {
			numTwoDT_Blocks = Integer.parseInt(value);
			if(testDataOutput)
			    System.err.println("num dt2 blocks: " + value);
		}
		else if(flag.equals("DT1_E1_THETA")) {
			dt1E1Theta = Double.parseDouble(value);
			if(testDataOutput)
			    System.err.println("dt1 e1 theta: " + value);
		}
		else if(flag.equals("DT1_E1_PHI")) {
			dt1E1Phi = Double.parseDouble(value);
			if(testDataOutput)
			    System.err.println("dt1 e1 phi: " + value);
		}
		else if(flag.equals("DT2_E1_THETA")) {
			dt2E1Theta = Double.parseDouble(value);
			if(testDataOutput)
			    System.err.println("dt2 e1 theta: " + value);
		}
		else if(flag.equals("DT2_E1_PHI")) {
			dt2E1Phi = Double.parseDouble(value);
			if(testDataOutput)
			    System.err.println("dt2 e1 phi: " + value);
		}
		else if(flag.equals("DT2_ROT_AXIS_THETA")) {
			dt2RotAxisTheta = Double.parseDouble(value);
			if(testDataOutput)
			    System.err.println("dt2 rot axis theta: " + value);
		}
		else if(flag.equals("DT2_ROT_AXIS_PHI")) {
			dt2RotAxisPhi = Double.parseDouble(value);
			if(testDataOutput)
			    System.err.println("dt2 rot axis phi: " + value);
		}
		else {
			// if non of the above work, somethings messed up
			logger.warning("could not parse info file entry: " + flag + "\t" + value);
		}
    }

    /*
     * Loads data voxel by voxel and bins the data as required. 
     *
     */
    private static void processData(String filename, int numComponents){
	    
	    // prepare the data source
	    VoxelOrderDataSource dataIn = new VoxelOrderDataSource(filename, numComponents, "double");

	    // variables to keep track of the test function parameters
	    int voxelNum = 0;
	    double rotAngle=0;
	    double dt2rotangle=0;
	     
	    /* create vectors that hold the true-fibre orientations for
	     * the one and two fibre cases respectively
	     */
	    Vector3D [] trueFibreOneDT = new Vector3D[1];
	    Vector3D [] trueFibresTwoDTs = new Vector3D[2];

	    int totalDataSize = oneDT_BlockSize + (twoDT_BlockSize * numTwoDT_Blocks);
		
	    if(testDataOutput)
		System.err.println("total data size: " + totalDataSize);

	    // create Vectors for the two test function directions (prior to any rotations)
	    Vector3D dir1 = Vector3D.vectorFromSPC(1.0, dt1E1Theta, dt1E1Phi);
	    Vector3D dir2 = Vector3D.vectorFromSPC(1.0, dt2E1Theta, dt2E1Phi);
	    
	    // define the rotation axis for the two fibre cases
	    Vector3D dt2RotAxisVect = Vector3D.vectorFromSPC(1.0, dt2RotAxisTheta, dt2RotAxisPhi);
	    double dt2angle;
	    Vector3D dt2Vect = null;
	    RealMatrix transformed = null;
	    RealMatrix rotationMat = null;
	    double [] data = null;
	    // for each set of fibre-orientation estimates
	    while(voxelNum < totalDataSize){
		data = dataIn.nextVoxel();
		
		if(data[2]>0){
		    //first deal with the single fibre cases
		    if(voxelNum < oneDT_BlockSize) {
			// get the true fibre orientation
			rotationMat = Rotations.randomRotMat(rand);
			trueFibreOneDT[0] = Rotations.rotateVector(dir1, rotationMat);

			//bin the foe into histogram
			binFOEs(data, trueFibreOneDT, 1);
		    }
		    // the two fibre cases
		    else if(voxelNum >= oneDT_BlockSize) {

			//need to determine what block you're in...
			int blocknum = (int)(voxelNum - oneDT_BlockSize)/twoDT_BlockSize;

			//calculate the rotation of the 2nd diffusion tensor
			dt2angle = -(minDT2RotAngle + (dt2RotAngleStep * blocknum));

			if(dt2angle > maxDT2RotAngle)
			    logger.warning("dt2angle exceeded maximum angle specified! angle: " + dt2angle );
				
			dt2Vect = Rotations.rotateVector(dir2, dt2RotAxisVect, dt2angle);

			// calculate the true-fibre orientations, including random rotations
			rotationMat = Rotations.randomRotMat(rand);
			trueFibresTwoDTs[0] = Rotations.rotateVector(dir1, rotationMat);
			trueFibresTwoDTs[1] = Rotations.rotateVector(dt2Vect, rotationMat);
				
			//bin the foe into histogram
			binFOEs(data, trueFibresTwoDTs, 2);
		    }
		}
		else {
		    // iterate random rotation anyway
		    rotationMat = Rotations.randomRotMat(rand);
		}
		//increment the voxel counter
		voxelNum++;
	    }
	}

    /*
     * Bins the fibre-orientation estimates into a histogram 
     *
     */
   private static void binFOEs(double [] data, Vector3D [] trueFibres, int numFibreDirs) {

	    /* 	for each combination
	     *      if sum of angles is less than current minimum
	     *	    set fibreOrder to current order
	     *	    for each true fibre
	     *	    calculate stats against fibre[fibreOrder[i]]
	     *	    save output
	     */

	     double dotProdSum = 0;
	     int [] fibreOrder = null;
	     int [] [] combination = getCombinations(numFOEsPerVoxel);

	     Vector3D dataVec = null;
	     // find correct allocations
	     for(int j=0;j<combination.length;j++) {
		 
		 // check this!!!!!
		 double sum = 0.0;
		     
		 for(int k=0;k<numFibreDirs;k++) {
		     dataVec = new Vector3D(data[6 + combination[j][k]*8], data[7 + combination[j][k]*8], data[8 + combination[j][k]*8]);
		     // if fibre is (0,0,0), i,e, it wasn't recovered, dot product is zero, so does not contribute to sum
		     sum+=Math.abs(trueFibres[k].dot(dataVec));
		 }

		 // if current sum is greater than the previous best, use the current combination instead.		
		 if(sum > dotProdSum)
		 {
		     dotProdSum = sum;
		     fibreOrder = combination[j];
		 }
	     }

	     Vector3D fibre = null;
	     RealMatrix hessianMat = new RealMatrix(2,2);
	     double meanF = 1.0;

	     // calculate stats
	     for(int j=0;j<numFibreDirs;j++) {
		 //extract the mean of the function mean(f)
		 meanF = data[4];

		 // extract the fibre from the data
		 fibre = new Vector3D(data[6 + fibreOrder[j]*8], data[7 + fibreOrder[j]*8], data[8 + fibreOrder[j]*8]);
		 

		 //generate hessian matrix for the peak
		 hessianMat.setEntry(0,0, data[10 + fibreOrder[j]*8]);
		 hessianMat.setEntry(0,1, data[11 + fibreOrder[j]*8]);
		 hessianMat.setEntry(1,0, data[12 + fibreOrder[j]*8]);
		 hessianMat.setEntry(1,1, data[13 + fibreOrder[j]*8]);
					
		 // if the magnitude of the fibre orientation estimate is > 0
		 if((fibre.x*fibre.x + fibre.y*fibre.y + fibre.z*fibre.z) > 0) {
		     // transform the fibre and extract the rotated true fibre orientation and eigenvalues of the hessian
		     transformFibre(fibre, hessianMat, trueFibres[j], meanF, numFibreDirs);    
		 }
		 else {
		     //logger.warning("null fibre-orientation estimate.  Program will continue");
		 }
	     }
   }

private static void transformFibre(Vector3D fibreVect, RealMatrix hessianMat, Vector3D trueVect, double meanF, int numFibreDirs)
	{
	    Vector3D zVect = new Vector3D(0.0, 0.0, 1.0);

		// rotate the fibre to the z-axis
		RealMatrix fibreRotationMat = Rotations.getRotMat(fibreVect, zVect);

		//calculate the eigenvects (index 1) and eigenvals (index 0)
		RealMatrix [] eigHessian = hessianMat.jacobi();

		//rotate the fibre - don't assume first eigenVector corresponds to largest eigenValue!
		Vector3D eigenVect;
		double eig1 = 0.0;
		double eig2 = 0.0;
		
		if(Math.abs(eigHessian[0].entry(0,0))>Math.abs(eigHessian[0].entry(1,1)))
		    {
			eigenVect = new Vector3D(eigHessian[1].entry(0,0), eigHessian[1].entry(1,0), 0.0);
			eig1 = Math.abs(eigHessian[0].entry(0,0));
			eig2 = Math.abs(eigHessian[0].entry(1,1));
			
		    }
		else
		    {
			eigenVect = new Vector3D(eigHessian[1].entry(1,0), eigHessian[1].entry(1,1), 0.0);
			eig1 = Math.abs(eigHessian[0].entry(1,1));
			eig2 = Math.abs(eigHessian[0].entry(0,0));
		    }

		Vector3D xVect = new Vector3D(1.0, 0.0, 0.0);

		// find the rotation matrix that aligns the largest eigenVector with the x axis
		RealMatrix hessianRotationMat = Rotations.getRotMat(eigenVect, xVect);

		// combine the two rotation matrices
		RealMatrix rotationMatrix = hessianRotationMat.product(fibreRotationMat);

		//rotate the true fibre orientations using the combined rotation tranform
		Vector3D transformedFibre = Rotations.rotateVector(trueVect, rotationMatrix);

		if(trueVect.dot(fibreVect) < 0)
		    transformedFibre = transformedFibre.negated();

		// meanF is used for scaling the eigenvalues of the Hessian, so we need to make sure that they are not zero!
		if(meanF == 0.0)
		    {
			meanF=1.0;
			logger.warning("mean of function is zero!  Setting to " + meanF + " and continuing");
		    }
		

		// add to the processed list
		if(numFibreDirs==1){
		    hist_1Fibre.add(new Vector3D(transformedFibre.x, transformedFibre.y, transformedFibre.z), eig1/meanF, eig2/meanF);
		}
		else if(numFibreDirs==2){
		    hist_2Fibre.add(new Vector3D(transformedFibre.x, transformedFibre.y, transformedFibre.z), eig1/meanF, eig2/meanF);
		}

	}

   /*
    * Returns the possible organization mappings between the fibre-orientation
    * estimates and the true fibres
    */
    private static int [] [] getCombinations(int numFibres) {
	if(numFibres==1) {
	    int [][] result = {{0}};
	    return result;
	}
	else if(numFibres==2) {
	    int [][] result = {{0,1},{1,0}};
	    return result;
	}
	else if(numFibres==3) {
	    int [][] result = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
	    return result;
	}
	else {
	    //warning!
	    logger.warning("numpds greater than 3.  Currently not supported!  Program will attempt to continue");
	}
	return null;
   }

    /*
    * Fits a distribution to the vectors in a bin
    *
    */
   private static double [] [] fitDistribution(String type, Histogram histogram){
       
       // output list
       double list [] [];

       Set histogramSet = histogram.getEntrySet();
       Iterator itor = histogramSet.iterator();

       // find out how many bins there are
       int numBins = histogramSet.size();

       // choose distribution and set up array
	if(type.equals(BINGHAM) || type.equals(ACG)){

	    list = new double [numBins][4];
	}
	else if (type.equals(WATSON)){
	    list = new double [numBins][2];
	}
	else {
	    throw new LoggedException("distribution type " + type + " unknown");
	}

	// get iterator
	int listCount = 0;

	HistBin bin;
	double [] tmp, eigs;
	Vector3D [] dirs;
	while(itor.hasNext()){
	    // get next bin
	    bin = (HistBin)(((Entry)itor.next()).getValue());
	    dirs = bin.getDirsList();

	    // if the bin has at least numVectsPerBin then calculate the distribution parameters
	    if(dirs.length >= numVectsPerBin) { 
		if(type.equals(BINGHAM)){
		    
		    //store eigs and bingham params
		    tmp = SphPDF_Fit.getBinghamConcentrationParams(dirs);
		    if(useLog)
		    	eigs = bin.getLogEigs();
		    else
			eigs = bin.getEigs();

		    list[listCount][0] = eigs[0];
		    list[listCount][1] = eigs[1];
		    list[listCount][2] = Math.log(-tmp[0]);
		    list[listCount][3] = Math.log(-tmp[1]);
		    //System.err.println(list[listCount][0] + " " + list[listCount][1] + " " + list[listCount][2] + " "  + list[listCount][3]);
		}
		else if (type.equals(WATSON)){
		    // store trace and watson param
		    tmp = SphPDF_Fit.getWatsonConcentrationParams(dirs);
		    if(useLog)
			list[listCount][0] = bin.getLogTrace();
		    else
			list[listCount][0] = bin.getTrace();
		    
		    list[listCount][1] = tmp[0];

		}
		else if (type.equals(ACG)){ //three params!
		    //store eigs and ACG params
		    tmp = SphPDF_Fit.getACGConcentrationParams(dirs);
		    //System.err.println("ACG temp size: " + tmp.length);
		    System.err.println("ACG eigs: " + tmp[0] + " " + tmp[1] + " " + tmp[2]); 
		    if(useLog)
			eigs = bin.getLogEigs();
		    else
			eigs = bin.getEigs();

		    list[listCount][0] = eigs[0];
		    list[listCount][1] = eigs[1];
		    list[listCount][2] = Math.log(tmp[0]/tmp[2]);
		    list[listCount][3] = Math.log(tmp[1]/tmp[2]);
		}
		listCount++;
	    }
	}

	//resize array to account for the bins that have less than numVectsPerBin fibre orientations per bin and return
	return trimArray(list, listCount);
   }

    /*
     * trims a 2D array to a specified size
     */
    private static double [][] trimArray(double [][] a, int trimLength) {
	int numComponents = a[0].length;

	double [][] b = new double [trimLength][numComponents];

	for(int i=0; i<trimLength;i++)
	    for(int j=0;j<numComponents;j++)
		b[i][j] = a[i][j];

	return b;
    }

    /*
    * Fits surfaces to the distribution parameters.  For use with
    * generalized PICo
    */
    private static double [] fitSurface(double [] [] orderedData, int zInd){
	//System.err.println("OrderedData length = " + orderedData.length);

	// calculate the number of coeffs required
	int tmpInd = ((order+1)*(order+2))/2;

	        if(orderedData.length<1)
	    		throw new LoggedException("no fibre-orientation estimates in histogram!  Increase calibration data set size (best) or reduce -minvectsperbin for testing.");
		else if(orderedData.length < tmpInd)
		    throw new LoggedException("not enough populated bins to enable surface fit.  Number of bins: " + orderedData.length);
	    	if(zInd >= orderedData[0].length)
	    		throw new LoggedException("out of bounds! zInd = " + zInd + " orderedDataSize = " + orderedData[0].length);

		// convert to 1D arrays or Matrices as appropriate!!
		double [] x = new double[orderedData.length];
		double [] y = new double[orderedData.length];

		RealMatrix zMat = new RealMatrix(1,orderedData.length);

		for(int i=0; i<orderedData.length; i++) {
			// x coord  (eigenvalue 1)
			x[i] = orderedData[i][0];

			// y coord (eigenvalue 2)
			y[i] = orderedData[i][1];
            
			// distribution coeff
			zMat.setEntry(0,i,orderedData[i][zInd]);
		}

		//f(x,y)=sum_i(sum_j(a_ij x^i y^j)) where (i+j)<=3
		RealMatrix xyMat = new RealMatrix(tmpInd, x.length);
		int indexForXyMat=0;
		for(int i=0; i <= order; i++) {
			for(int j=0;j <= order; j++) {
				if((i+j)<=order) {
					for(int k=0;k < x.length; k++) {
					    xyMat.setEntry(indexForXyMat,k,Math.pow(x[k],i) * Math.pow(y[k],j));
					}
					indexForXyMat++;
				}
			}
		}
		RealMatrix xyInvMat = xyMat.pseudoInv();
		RealMatrix coeffsMat = zMat.product(xyInvMat);

		// convert result to array form and return
		double [] coeffs = new double [coeffsMat.columns()];

		for(int i=0;i<coeffs.length;i++) {
		    coeffs[i]=coeffsMat.entry(0,i);
		}

		if(testDataOutput) {
		    System.err.println("Coeffs: ");
		    for(int i=0;i<coeffs.length;i++) {
			System.err.println(coeffs[i]);
		    }
		    System.err.println();
		}
		return coeffs;
    }

    /*
    * Fits line to the distribution parameter.  For use with Parker and Alexanders method
    */
    private static double [] fitLine(double [] [] orderedData){
	
	if(orderedData.length==0)
	    		throw new LoggedException("no fibre-orientation estimates in histogram!");

	// distribution coeff
	double [] x = new double[orderedData.length];
	RealMatrix yMat = new RealMatrix(1,orderedData.length);

	for(int i=0;i<orderedData.length;i++) {
	    x[i] = orderedData[i][0];
	    yMat.setEntry(0,i,orderedData[i][1]);
	}

	RealMatrix xMat = new RealMatrix(order+1, x.length);

	for(int i=0;i < order+1;i++) {
	    for(int j=0; j< x.length; j++) {
		xMat.setEntry(i,j,Math.pow(x[j],i));
	    }
	}

	RealMatrix xInvMat = xMat.pseudoInv();
	RealMatrix coeffsMat = yMat.product(xInvMat);

	// convert result to array form and return
	double [] coeffs = new double [coeffsMat.columns()];

	for(int i=0;i<coeffs.length;i++) {
	    coeffs[i]=coeffsMat.entry(0,i);
	}

	if(testDataOutput) {
	    System.err.println("Coeffs: ");
	    for(int i=0;i<coeffs.length;i++) {
		System.err.println(coeffs[i]);
	    }
	    System.err.println();
	}

	return coeffs;
    }

    /*
    * outputs coeffs to specified file
    */
    private static void saveCoeffs(double [] coeffs, String filename)
    {
	OutputManager.outputFile = filename;
	OutputManager om = new OutputManager();
	
	// output data as doubles
	om.output(coeffs);

        om.close();
    }

    private static void outputCoeffs(double [] coeffs, String filename)
    {

	double [] finalCoeffs = new double [coeffs.length+2];

	if(distributionType.equals(BINGHAM) || distributionType.equals(ACG))
	    finalCoeffs[0]=2;
	else if(distributionType.equals(WATSON))
	    finalCoeffs[0]=1;
	    
	finalCoeffs[1]=order;
	for(int i=0; i < coeffs.length; i++) {
	    finalCoeffs[i+2] = coeffs[i];
	}

	saveCoeffs(finalCoeffs, filename);
    }

    private static void outputCoeffs(double [] coeffsA, double [] coeffsB, String filename)
    {
	// coeffs arrays will be the same size!
	if(coeffsA.length != coeffsB.length) {
	    throw new LoggedException("surface coefficient arrays not the same length for 2 fibre case");
	}

	// concatenate the lists and output
	double [] finalCoeffs = new double [(coeffsA.length * 2)+2];
	if(distributionType.equals(BINGHAM) || distributionType.equals(ACG))
	    finalCoeffs[0]=2;
	else if(distributionType.equals(WATSON))
	    finalCoeffs[0]=1;
	    
	finalCoeffs[1]=order;
	for(int i=0; i < coeffsA.length; i++) {
	    finalCoeffs[i+2] = coeffsA[i];
	    finalCoeffs[i+coeffsA.length+2] = coeffsB[i];
	}

	saveCoeffs(finalCoeffs, filename);
    }

    public static void main(String[] args) {

	// Parse the command line arguments
	CL_Initializer.CL_init(args);
	String filename = CL_Initializer.inputFile;
	numFOEsPerVoxel = CL_Initializer.numPDsIO;
		//need to know the number of PDs per voxel!
	String infoFile = "";

	for(int i=0; i<args.length; i++) {
	    if(args[i].equals("-order")) {
			order = Integer.parseInt(args[i + 1]);
			CL_Initializer.markAsParsed(i, 2);
	    }
	    else if(args[i].equals("-binincsize")) {
			binSize = Double.parseDouble(args[i + 1]);
			CL_Initializer.markAsParsed(i, 2);
	    }
	    else if(args[i].equals("-minvectsperbin")) {
			numVectsPerBin = Integer.parseInt(args[i + 1]);
			CL_Initializer.markAsParsed(i, 2);
	    }
	    else if(args[i].equals("-directmap")) {
		useLog = false;    
		CL_Initializer.markAsParsed(i, 1);
	    }
	    else if(args[i].equals("-infofile")) {
			infoFile = args[i + 1];
			CL_Initializer.markAsParsed(i, 2);
	    }
	    else if(args[i].equals("-outputstem")) {
			outputStem = args[i + 1];
			CL_Initializer.markAsParsed(i, 2);
	    }
	    else if(args[i].equals("-pdf")) {
			distributionType = args[i + 1];
			CL_Initializer.markAsParsed(i, 2);
	    }
	    else if(args[i].equals("-test")) {
		testDataOutput = true;
			CL_Initializer.markAsParsed(i, 1);
	    }
	}
	
	CL_Initializer.checkParsing(args);

	// Now do some checks for user mistakes
	if(!distributionType.equals(BINGHAM) && !distributionType.equals(WATSON) && !distributionType.equals(ACG))
	    throw new LoggedException("distribution not recognized.  requested: " + distributionType);
	    
	// need to account for situations where more than one distribution selected
	if(order < 0)
	    throw new LoggedException("order must be >= 0.  order=" + order);
	if(binSize <= 0.0 )
	    throw new LoggedException("bin size must be > 0.0.  binincsize="+binSize);
	    
	if(distributionType.equals(BINGHAM) || distributionType.equals(ACG)) {
	    hist_1Fibre = new Histogram(2, binSize, useLog);
	    hist_2Fibre = new Histogram(2, binSize, useLog);
	}
	else {
	    hist_1Fibre = new Histogram(1, binSize, useLog);
	    hist_2Fibre = new Histogram(1, binSize, useLog);
	}
	int numElementsPerVoxel = 6 + (numFOEsPerVoxel*8);
	readInfoFile(infoFile);
	processData(filename, numElementsPerVoxel);


	// calculate the distribution parameters for each bin in the histograms
	double [] [] dist1 = fitDistribution(distributionType, hist_1Fibre);
	
	if(testDataOutput){
	    System.err.println("Dist 1");
	    for(int i=0;i< dist1.length;i++){
		for(int j=0;j<dist1[0].length;j++){
			System.err.print(dist1[i][j] + " ");

		}
		System.err.println();
	    }
	}
	double [] [] dist2 = fitDistribution(distributionType, hist_2Fibre);
	if(testDataOutput) {
	    System.err.println("Dist 2");
	    for(int i=0;i< dist2.length;i++){
		for(int j=0;j<dist2[0].length;j++){
			System.err.print(dist2[i][j] + " ");
		}
		System.err.println();
	    }
	}
	/* fit a surface to each set of distribution parameters and save the results into two
	 * separate files (for the one fibre case and the two fibre case respectively.
	 */
	if(distributionType.equals(BINGHAM) || distributionType.equals(ACG)) {
	    
	    double [] oneFibreSurface1 = fitSurface(dist1 , 2);
	    double [] oneFibreSurface2 = fitSurface(dist1 , 3);
	    outputCoeffs(oneFibreSurface1, oneFibreSurface2, outputStem + "_oneFibreSurfaceCoeffs.Bdouble");

	    double [] twoFibreSurface1 = fitSurface(dist2 , 2);
	    double [] twoFibreSurface2 = fitSurface(dist2 , 3);
	    outputCoeffs(twoFibreSurface1, twoFibreSurface2, outputStem + "_twoFibreSurfaceCoeffs.Bdouble");

	}
	else if(distributionType.equals(WATSON)){
	    double [] oneFibreLine = fitLine(dist1);
	    outputCoeffs(oneFibreLine, outputStem + "_oneFibreLineCoeffs.Bdouble");


	    double [] twoFibreLine = fitLine(dist2);
	    outputCoeffs(twoFibreLine, outputStem + "_twoFibreLineCoeffs.Bdouble");
	}
    }
}


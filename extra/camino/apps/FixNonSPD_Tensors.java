package apps;

import data.*;
import imaging.*;
import inverters.*;
import misc.*;
import tools.*;




/**
 * 
 * Fix negative eigenvalues and output diagnostics
 *
 * @author Philip Cook
 */
public class FixNonSPD_Tensors extends Executable {


    // dimensions of DT space
    private int xDataDim;
    private int yDataDim;
    private int zDataDim;
    
    // DWI data
    private DataSource dataSource;

    private boolean haveDWI;
    
    // If true, replace the tensor with one constrained to be SPD 
    private boolean fitCholTensor;

    // If true, winsorize bad DWI data before re-fitting 
    private boolean winsorizeDWI;

    private double unweightedB;

    // true if b[i] <= unweightedB
    boolean[] unweightedMeas;

    // Don't refit, just edit the tensor
    //
    // Replace negative eigenvalue l_i with min[abs(l_i), l_{i-1} - delta]
    //
    // Preserve ordering by not allowing the replacement value to be greater than the first positive eigenvalue
    //
    // If L2 and L3 are both < 0, use mean abs value for both
    // 
    // If L1 is also < 0, something seriously wrong, mask out the voxel
    private boolean editEigenvalues;


    // Diagnostic stuff


    // Number of negative eigenvalues, 1-3
    private double[][][] numNegEigenvalues;

    // Number of measurements implying negative diffusion, ie A > mean(A_0)
    private double[][][] numNegDiffCoeff;

    // is there an outlier b0 at this voxel
    private double[][][] outlierB0;

    // Output the corrected data, mask, and some diagnostics
    private String outputRoot;

    // Brain mask, might be edited compared to input
    private double[][][] brainMask;

    // Imaging scheme
    private DW_Scheme scheme;

    private int numMeas;

    private NonLinearDT_Inversion cholInv;

    
    public FixNonSPD_Tensors(String[] args) {
        super(args);
    }

    
    public void initDefaultVals() {
	xDataDim = 0;
        yDataDim = 0;
        zDataDim = 0;
        dataSource = null;
        unweightedB = 0.0;
        unweightedMeas = null;
        haveDWI = false;
	fitCholTensor = false;
	winsorizeDWI = false;
	editEigenvalues = false;
	numNegDiffCoeff = null;
        numNegEigenvalues = null;
	outlierB0 = null;
	outputRoot = "";
	scheme = null;
	numMeas = 0;
	cholInv = null;
    }

    
    public void initOptions(String[] args) {
              
        CL_Initializer.inputDataType = "double";
        CL_Initializer.CL_init(args);

	// -inputfile is the tensors (required), this is optional
	String dwiFile = "";

	for (int i = 0; i < args.length; i++) {
            
            if(args[i].equals("-dwifile")) {
		dwiFile = args[i+1];
                haveDWI = true;
                CL_Initializer.markAsParsed(i,2);
            }
	    if(args[i].equals("-outputroot")) {
		outputRoot = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
	    if(args[i].equals("-fitchol")) {
		fitCholTensor = true;
                CL_Initializer.markAsParsed(i);
            }
	    if(args[i].equals("-winsorizedwi")) {
		winsorizeDWI = true;
                CL_Initializer.markAsParsed(i);
            }
	    if(args[i].equals("-editeigenvalues")) {
		editEigenvalues = true;
                CL_Initializer.markAsParsed(i);
            }
	    if (args[i].equals("-unweightedb")) {
		unweightedB = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i,2);
	    }
	    
        }
        
        CL_Initializer.checkParsing(args);
        
	CL_Initializer.inputModel = "dt";

	CL_Initializer.initTensorDataSource();

	if (haveDWI) {
            
            CL_Initializer.initImagingScheme();
            
            scheme = CL_Initializer.imPars;
            
            numMeas = scheme.numMeasurements();

            // This relies on DWI data being the same type as tensors, or a NIFTI image (type determined by header)
	    dataSource = ExternalDataSource.getDataSource(dwiFile, scheme.numMeasurements(), CL_Initializer.inputDataType);

            cholInv = new NonLinearDT_Inversion(scheme);

	}
        else {
            editEigenvalues = true;
        }

        if (!(editEigenvalues || fitCholTensor)) {
            throw new LoggedException("Choice of replacement strategy required (edit eigenvalues or refit tensor)");
        }
        
        if (CL_Initializer.headerTemplate == null) {

	    if (ImageHeader.imageExists(dwiFile)) {
		CL_Initializer.headerTemplateFile = dwiFile;
	    } 
	    else if (ImageHeader.imageExists(CL_Initializer.bgMaskFile)) {
		CL_Initializer.headerTemplateFile = CL_Initializer.bgMaskFile;
	    } 
	    else {
		throw new LoggedException("This program requires a header template with -header or a NIFTI brain mask");
	    }
	}

	CL_Initializer.initInputSpaceAndHeaderOptions();

	xDataDim = CL_Initializer.dataDims[0];
	yDataDim = CL_Initializer.dataDims[1];
	zDataDim = CL_Initializer.dataDims[2];

        OutputManager.outputFile = outputRoot + "dtSPD.B" + OutputManager.outputDataType;

    }

    
    public void initVariables() {
	
	numNegDiffCoeff = new double[xDataDim][yDataDim][zDataDim];

        numNegEigenvalues = new double[xDataDim][yDataDim][zDataDim];

	brainMask = new double[xDataDim][yDataDim][zDataDim];

        if (haveDWI) {

            unweightedMeas = new boolean[numMeas];
            
            for (int m = 0; m < numMeas; m++) {
                if (scheme.getB_Value(m) <= unweightedB) {
                    unweightedMeas[m] = true;
                }
            }
        }
        
    }   

    
    public void execute(OutputManager om) {

        
	for (int k = 0; k < zDataDim; k++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
		    
		    double[] voxelTensorData = CL_Initializer.data.nextVoxel();

		    boolean background = CL_Initializer.bgMask.nextVoxel()[0] == 0.0;
            
		    double[] voxelData = null;

		    if (haveDWI) {
			voxelData = dataSource.nextVoxel();
		    }

		    if ( background != (voxelTensorData[0] == -1.0) ) {
			throw new LoggedException("Brain mask does not match tensor exit code background classification"); 
		    }

		    numNegEigenvalues[i][j][k] = 0.0;
		    numNegDiffCoeff[i][j][k] = 0.0;
		    
		    if (background) {
			brainMask[i][j][k] = 0.0;
			om.output(voxelTensorData);
			continue;
		    }
                    // Code what was done to the data 
                    //
                    // -100    Bad data flagged by original fit, removed from brain mask
                    // -2    Attempted Chol fit failed, edited eigenvalues
                    // -1    Background, nothing done
                    // 0     Good data, nothing done
                    // 2     Original nonlinear fit failed, but eigenvalues OK
                    // 6     Bad data changed in original fit, but eigenvalues OK
                    // 10    Modified by changing eigenvalues
                    // 12    Modified with Cholesky fit
                    
		    if (voxelTensorData[0] == -100.0) {
                        brainMask[i][j][k] = 0.0;
                        om.output(voxelTensorData);
			continue;
		    }
		    
		    DT inputTensor = dtFromModelFit(voxelTensorData);
		    
		    double[][] seig = inputTensor.sortedEigenSystem();

                    double meanB0 = 0.0;
                    
                    if (haveDWI) {
                        meanB0 = geoMean(voxelData, unweightedMeas);
                        numNegDiffCoeff[i][j][k] = numNegDiffCoeffs(voxelData, meanB0);
                    }
                                                
                    // By here, we know this is a voxel to keep in the mask
		    brainMask[i][j][k] = 1.0;

		    if (seig[0][2] > 0.0) {
			// Data all good, output it
			om.output(voxelTensorData);
		    }
		    else {
			
			double[] fixedTensorData = null;

			if (seig[0][0] < 0.0) {
                            numNegEigenvalues[i][j][k] = 3.0;
                        }
			else if (seig[0][1] < 0.0) {
			    numNegEigenvalues[i][j][k] = 2.0;
			}
			else {
			    numNegEigenvalues[i][j][k] = 1.0;
			}
                       
			if (editEigenvalues) {
			    fixedTensorData = editEigenvalues(seig);
                            // insert original estimate of ln(S0)
                            fixedTensorData[1] = voxelTensorData[1];
			}
			else {

			    // Attempt to re-fit tensor 
			    
			    if (winsorizeDWI) {
                                fixedTensorData = fitCholTensor(winsorizeDWI(voxelData, meanB0), seig);
			    }
			    else {
				fixedTensorData = fitCholTensor(voxelData, seig);
			    }
			    
			}

			om.output(fixedTensorData);			    
			    
		    }
   
		}
		
	    }    
	    
	}

        om.close();

        

	// now output image data

        try {
            ImageHeader shortHeader = ImageHeader.readHeader(CL_Initializer.headerTemplateFile);
            
            shortHeader.setDataType("short");
            
            if (haveDWI) {
                shortHeader.writeScalarImage(numNegDiffCoeff, outputRoot + "numNegDiffCoeff");
            }
            
            shortHeader.writeScalarImage(numNegEigenvalues, outputRoot + "numNegEV");

            shortHeader.writeScalarImage(brainMask, outputRoot + "brainMask");
            
        }
        catch(java.io.IOException e) {
            throw new LoggedException("Could not write output images, caught exception " + e);
        }

    }

    
    private DT dtFromModelFit(double[] voxelTensorData) {
	return new DT(voxelTensorData[2], voxelTensorData[3], voxelTensorData[4], voxelTensorData[5], 
		      voxelTensorData[6], voxelTensorData[7]);
    }


    private double geoMean(double[] voxelData, boolean[] includeMeasurements) {
        double gm = 1.0;

        double numMeas = voxelData.length;

        int numIncludedMeas = 0;
        
        for (int i = 0; i < numMeas; i++) {
            if (includeMeasurements[i]) {
                gm = gm * voxelData[i];
                numIncludedMeas += 1;
            }
        }
        
        return Math.pow(gm, 1.0 / numIncludedMeas);
        
    }

    
    // Expects a fixable tensor (ie at least one positive eigenvalue)
    private double[] editEigenvalues(double[][] seig) {

	double code = 10.0;

	if (seig[0][1] > 0.0) {
	    // Just one eigenvalue to change

	    seig[0][2] = Math.min(Math.abs(seig[0][2]), seig[0][1] * 0.99);
	    
	}
	else {
	    
	    if (seig[0][0] > 0.0) {
                double fixEig = ( Math.abs(seig[0][1]) +  Math.abs(seig[0][2]) ) / 2.0;

                if (fixEig == 0.0) {
                    // Chol fit will sometimes converge to zero eigenvalues, which cannot be represented in the
                    // log space
                    fixEig = seig[0][0] / 5.0;
                }
                
                seig[0][1] = fixEig;
                seig[0][2] = fixEig;
            }
            else {
                double md = Math.abs(seig[0][0]);
                
                seig[0][0] = md;
                seig[0][2] = md;
                seig[0][1] = md;
            }

	}
	
	double[] dtComps = DT.dtFromEig(seig).getComponents();

	double[] fixedTensor = new double[8];

	fixedTensor[0] = code;
	fixedTensor[1] = 1.0;
	
	for (int i = 0; i < 6; i++) {
	    fixedTensor[2+i] = dtComps[i];
	}
	
	return fixedTensor;
	
    }

    
    private double[] winsorizeDWI(double[] inputData, double meanB0) {
	double[] fixedData = new double[numMeas];
	
	System.arraycopy(inputData, 0, fixedData, 0, numMeas);

	for (int i = 0; i < numMeas; i++) {
	    if (!(unweightedMeas[i]) && fixedData[i] > meanB0) {
		fixedData[i] = meanB0 * 0.9999;
	    } 
	}

        return fixedData;
    }

    
    // Attempt chol fit, edit eigenvalues directly if that fails
    private double[] fitCholTensor(double[] dwiData, double[][] seig) {
	
	double[] fixedTensorComponents = cholInv.invert(dwiData);

        double code = 12.0;
        
	if (fixedTensorComponents[0] == 2.0) {
	    // fit failed
	    fixedTensorComponents = editEigenvalues(seig);
	    code = -2.0;
	}

        fixedTensorComponents[0] = code;
        
	
	// Chol fit can return 0 eigenvalues
        DT cholTensor = dtFromModelFit(fixedTensorComponents);
        double[][] cholSeig = cholTensor.sortedEigenSystem();

        if (cholSeig[0][2] <= 0.0) {
            double[] cholFix = editEigenvalues(cholSeig);
            cholFix[0] = code;
            return(cholFix);
        }

        return(fixedTensorComponents);
    }

 
    private int numNegDiffCoeffs(double[] data, double meanB0) {

	int numNeg = 0;

	for (int i = 0; i < numMeas; i++) {
	    if (!(unweightedMeas[i]) && data[i] > meanB0) {
		numNeg += 1;
	    } 
	}

	return numNeg;
    }

    
    
    
}

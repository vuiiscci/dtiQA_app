package apps;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;

import java.io.*;
import java.util.ArrayList;

/**
 * Fits a spherical PDF to a collection of axes, and outputs a result in the same format
 * as PICoPDFs (number of PDs is assumed to be 1).
 *
 * @author Philip Cook
 * @version $Id$
 *
 * @see apps.PICoPDFs
 */
public class SphPDF_Fit extends Executable {
    
    
    /**
     * Output manager
     */
    //private static OutputManager om;
    
    public SphPDF_Fit(String[] args){
        super(args);
    }
    
    // dyads default to upper triangular, NIFTI uses lower
    boolean upperTriangular;
    boolean outputParamsOnly;
    
    String pdf;
    Vector3D e1;
    Vector3D e2;
    Vector3D e3;
    int vectorsPerVoxel;
    
    // dyadic tensors can come with headers
    String inputImageType;
	
    boolean lastRun;
    // vector or dyad data
    DataSource dataSource;
    DataSource bgSource;
    EigenSystem3D nextTBarEig;
    Vector3D[] vectors;
    int counter;	
    
    public void initDefaultVals(){
        upperTriangular = true;
        outputParamsOnly = true;
        pdf = "bingham";
        e1 = null;
        e2 = null;
        e3 = null;
        vectorsPerVoxel = 0;

        inputImageType = "raw";

        lastRun = false;
        dataSource = null;
        bgSource = null;
        nextTBarEig = null;
        vectors = null;      
        counter = 0;
    }
    
    public void initOptions(String[] args){

        CL_Initializer.inputDataType = "double";
        CL_Initializer.inputModel = "vector";
        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            
            if (args[i].equals("-pdf")) {
                pdf = args[i+1];
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-evecs")) {
                
                e1 = Vector3D.vectorFromSPC(1.0, Double.parseDouble(args[i+1]),
                Double.parseDouble(args[i+2]));
                
                e2 = Vector3D.vectorFromSPC(1.0,  Double.parseDouble(args[i+3]),
                Double.parseDouble(args[i+4]));
                
                e3 = Vector3D.vectorFromSPC(1.0,  Double.parseDouble(args[i+5]),
                Double.parseDouble(args[i+6]));
                
                CL_Initializer.markAsParsed(i, 7);
            }
            if (args[i].equals("-outputpicoformat")) {
                outputParamsOnly = false;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-vectorspervoxel")) {
                vectorsPerVoxel = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            
            
        }
        
        CL_Initializer.checkParsing(args);
        
        if (CL_Initializer.inputFile != null) {
            
        }
    }
    
    public void initVariables(){
       
        if (pdf.equals("acg") && !CL_Initializer.inputModel.equals("vector")) {
            throw new LoggedException("Need vectors to fit ACG");
        }
        
       
        if (CL_Initializer.bgMaskFile != null) {
            CL_Initializer.initMaskSource();
            bgSource = CL_Initializer.bgMask;
            
        }
                
        if (CL_Initializer.inputModel.equals("vector")) {
            
            // vector data must be raw and voxel ordered, too much memory use otherwise

            dataSource = new VoxelOrderDataSource(CL_Initializer.inputFile, 3, CL_Initializer.inputDataType);
            
            if (vectorsPerVoxel == 0) {
                
                ArrayList<Vector3D> vectorList = new ArrayList<Vector3D>(2000);
                
                //int counter = 0;
                
                while (dataSource.more()) {
                    vectorList.add(new Vector3D(dataSource.nextVoxel()));
                    counter++;
                }
                
                vectorsPerVoxel = counter;
                
                vectors = new Vector3D[vectorsPerVoxel];
                
                vectorList.toArray(vectors);
                
                nextTBarEig = SphericalDistributionFitter.tBarEigenSystem(vectors);
                
            }
            else {
                
                vectors = new Vector3D[vectorsPerVoxel];
                
                for (int i = 0; i < vectorsPerVoxel; i++) {
                    vectors[i] = new Vector3D(dataSource.nextVoxel());
                }
                
                
                nextTBarEig = SphericalDistributionFitter.tBarEigenSystem(vectors);
                
                
            }
            
        }
        else if (CL_Initializer.inputModel.equals("dyad")) {

 
            if (CL_Initializer.inputFile != null) {
                if (CL_Initializer.inputFile.endsWith(".mha") || CL_Initializer.inputFile.endsWith(".mhd")) {
                    inputImageType = "meta";
                }
                else {
                    if (CL_Initializer.inputFile.endsWith(".nii") || CL_Initializer.inputFile.endsWith(".nii.gz")) {
                        inputImageType = "nifti";
                        upperTriangular = false;
                    }
                }
                
            }
                        
            if (inputImageType.equals("raw")) {
                dataSource = new VoxelOrderDataSource(CL_Initializer.inputFile, 6, CL_Initializer.inputDataType);
            }
            else if (inputImageType.equals("meta")) {
                try {
                    MetaImageHeader mh = MetaImageHeader.readHeader(CL_Initializer.inputFile);
                    dataSource = mh.getImageDataSource();
                }
                catch (IOException e) {
                    throw new LoggedException(e);
                }
            }
            else if (inputImageType.equals("nifti")) {
                Nifti1Dataset nd = new Nifti1Dataset(CL_Initializer.inputFile);
                dataSource = nd.getImageDataSource();
            }
            
            
            nextTBarEig = SphericalDistributionFitter.tBarEigenSystem(dataSource.nextVoxel(), upperTriangular);
        }
        
        
        int counter = 0;
        
//        boolean lastRun = false;
        
    }
    
    public void execute(OutputManager om){
        while (!lastRun) {
            
            lastRun = !dataSource.more();
            
            double[] concentrationParams = null;
            
            boolean background = false;
            
            if ((bgSource != null && bgSource.nextVoxel()[0] == 0.0) || nextTBarEig.eigenvalues[0] == 0.0) {
                
                background = true;
                
                if (pdf.equals("bingham")) {
                    concentrationParams = new double[2];
                }
                else if (pdf.equals("watson")) {
                    concentrationParams = new double[1];
                }
                else if (pdf.equals("acg")) {
                    concentrationParams = new double[3];
                }
                else {
                    throw new LoggedException("Unrecognized pdf: " + pdf);
                }
            }
            else {
                // fit concentration
                try {
                    
                    if (pdf.equals("bingham")) {
                        concentrationParams = getBinghamConcentrationParams(nextTBarEig);
                    }
                    else if (pdf.equals("watson")) {
                        concentrationParams = getWatsonConcentrationParams(nextTBarEig, vectors);
                    }
                    else if (pdf.equals("acg")) {
                        concentrationParams = getACGConcentrationParams(vectors);
                    }
                    else {
                        throw new LoggedException("Unrecognized pdf: " + pdf);
                    }
                    
                }
                catch(LoggedException e) {
                    // what to do about errors? Should probably warn
                    for (int i = 0; i < concentrationParams.length; i++) {
                        concentrationParams[i] = 0.0;
                    }
                    System.err.println("Couldn't fit params in voxel " + counter +
                    ", scatter matrix eig\n" + nextTBarEig);
                }
                
                // screen out insanely large params, get these in background regions
                if (!(Math.abs(concentrationParams[0]) < 1E6)) {
                    for (int i = 0; i < concentrationParams.length; i++) {
                        concentrationParams[i] = 0.0;
                    }
                    System.err.println("Couldn't fit params in voxel " + counter +
                    ", scatter matrix eig\n" + nextTBarEig);
                }
                
            }
            
            
            if (outputParamsOnly) {
                // params
                om.output(concentrationParams);
                
            }
            else {
                // output eigenvectors specified on command line, but if these aren't given, output
                // vectors from scatter matrix of samples
                double[] evecs = new double[9];
                
                if (e1 != null) {
                    evecs[0] = e1.x;
                    evecs[1] = e1.y;
                    evecs[2] = e1.z;
                    
                    evecs[3] = e2.x;
                    evecs[4] = e2.y;
                    evecs[5] = e2.z;
                    
                    evecs[6] = e3.x;
                    evecs[7] = e3.y;
                    evecs[8] = e3.z;
                    
                }
                else {
                    
                    
                    for (int i = 0; i < 3; i++) {
                        evecs[3*i] = nextTBarEig.eigenvectors[i].x;
                        evecs[3*i+1] = nextTBarEig.eigenvectors[i].y;
                        evecs[3*i+2] = nextTBarEig.eigenvectors[i].z;
                    }
                }
                
                // numPDs, mix
                om.output(new double[] {background ? 0.0 : 1.0, 1.0});
                
                // evecs
                om.output(evecs);
                
                // params
                om.output(concentrationParams);
            }
            
            counter++;
            
            if (!lastRun) {
                
                if (CL_Initializer.inputModel.equals("vector")) {
                    vectors = new Vector3D[vectorsPerVoxel];
                    
                    for (int i = 0; i < vectorsPerVoxel; i++) {
                        vectors[i] = new Vector3D(dataSource.nextVoxel());
                    }
                    
                    nextTBarEig = SphericalDistributionFitter.tBarEigenSystem(vectors);
                }
                else if (CL_Initializer.inputModel.equals("dyad")) {
                    nextTBarEig = SphericalDistributionFitter.tBarEigenSystem(dataSource.nextVoxel(),
                    upperTriangular);
                }
                
            }
            
        }
        
        om.close();
    }
    
    
    
    protected static double[] getWatsonConcentrationParams(Vector3D[] sampleVecs) {
        
        EigenSystem3D eig = SphericalDistributionFitter.tBarEigenSystem(sampleVecs);
        
        if (sampleVecs != null) {
            return new double[] {WatsonFitter.fitKappa(eig, sampleVecs)};
        }
        else {
            return new double[] {WatsonFitter.fitKappa(eig)};
        }
        
    }
    
    
    
    protected static double[] getWatsonConcentrationParams(EigenSystem3D eig, Vector3D[] sampleVecs) {
        
        if (sampleVecs != null) {
            return new double[] {WatsonFitter.fitKappa(eig, sampleVecs)};
        }
        else {
            return new double[] {WatsonFitter.fitKappa(eig)};
        }
        
    }
    
    
    protected static double[] getBinghamConcentrationParams(EigenSystem3D eig) {
        
        try {
            double[] bngpar = BinghamFitter.bngpar(eig.eigenvalues[1], eig.eigenvalues[2]);
            
            return new double[] {bngpar[0], bngpar[1]};
            
        }
        catch (ConvergenceException e) {
            return new double[] {0.0, 0.0};
        }
        
    }
    
    protected static double[] getBinghamConcentrationParams(Vector3D[] sampleVecs) {
        
        EigenSystem3D eig = SphericalDistributionFitter.tBarEigenSystem(sampleVecs);
        
        try {
            double[] bngpar = BinghamFitter.bngpar(eig.eigenvalues[1], eig.eigenvalues[2]);
            
            return new double[] {bngpar[0], bngpar[1]};
            
        }
        catch (ConvergenceException e) {
            return new double[] {0.0, 0.0};
        }
        
    }
    
    
    
    protected static double[] getACGConcentrationParams(Vector3D[] sampleVecs) {
        
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

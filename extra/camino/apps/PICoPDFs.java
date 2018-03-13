package apps;

import java.util.logging.Logger;
import data.*;
import misc.*;
import numerics.*;
import tools.*;
import tractography.*;

import java.io.*;



/**
 * Calculate PICo PDFs in each voxel.
 *
 * @version $Id$
 *
 * @author Philip Cook
 */
public class PICoPDFs extends Executable{
    
    public PICoPDFs(String[] args){
        super(args);
    }
    
    /**
     * Concentration parameters per PD.
     */
    private static int paramsPerPD;
    
    /**
     * Logger object
     */
    private static Logger logger = Logger.getLogger("apps/PICoPDFs");
    
    public static boolean flipX;
    public static boolean flipY;
    public static boolean flipZ;
    
	private String[] lutFiles;        
    private String pdf;        
    private boolean log;
        
    private int valuesPerVoxel;
    private int numLUTs;
    
    /**
     * Output manager
     */
    //private static OutputManager om;
    
    public void initDefaultVals(){
		paramsPerPD = 0;
		flipX = false;
		flipY = false;
		flipZ = false;
		lutFiles = new String[10];
		pdf = "bingham";
		log = true;
    }
    
    /**
     * Get PICo PDFs for SF data.
     */
    public static void mapSFData(DataSource data, String pdf, boolean log,
    String[] funcFiles, OutputManager om) {
        
        boolean watson = pdf.equals("watson");
        boolean bingham = pdf.equals("bingham");
        boolean acg = pdf.equals("acg");
        
        if (!(watson || bingham || acg)) {
            throw new LoggedException("Unrecognized PDF: " + pdf);
        }
        
        if (watson) {
            paramsPerPD = 1;
        }
        else if (bingham) {
            paramsPerPD = 2;
        }
        else {
            paramsPerPD = 3;
        }
        
        // functions[number of pds][function for each param of PDF]
        PolynomialFunction[][] functions = new PolynomialFunction[CL_Initializer.numPDsIO][];
        
        for (int i = 0; i < CL_Initializer.numPDsIO; i++) {
            
            if (funcFiles[i] == null) {
                throw new LoggedException("Number of LUTs must be " + CL_Initializer.numPDsIO +
                ", the maximum number of PDs in each voxel");
            }
            
            functions[i] = PolynomialFunction.readFunctions(funcFiles[i]);
        }
        
        while (data.more()) {
            
            double[] voxel = data.nextVoxel();
            
            double[][] params = new double[CL_Initializer.numPDsIO][paramsPerPD];
            
            // all three vectors in one array
            double[][] vectors = new double[CL_Initializer.numPDsIO][9];
            
            double[] mix = new double[CL_Initializer.numPDsIO];
            
            int numPDs = 0; // number of PDs in this voxel
            
            numPDs = (int)voxel[2];
            
            if (numPDs > CL_Initializer.numPDsIO) {
                numPDs = CL_Initializer.numPDsIO;
            }
            
            if (voxel[0] < 0.0) {
                numPDs = 0;
            }
            
            // starts at 6 (x, y, z, f, H00, H01, H10, H11)
            pdloop:
                for (int i = 0; i < numPDs; i++) {
                    
                    Vector3D peak =
                    new Vector3D(voxel[6 + i * 8], voxel[7 + i * 8], voxel[8 + i * 8]);
                    
                    if (flipX) {
                        peak = new Vector3D(-peak.x, peak.y, peak.z);
                    }
                    if (flipY) {
                        peak = new Vector3D(peak.x, -peak.y, peak.z);
                    }
                    if (flipZ) {
                        peak = new Vector3D(peak.x, peak.y, -peak.z);
                    }
                    
                    double h11 = voxel[10 + i * 8];
                    double h12 = voxel[11 + i * 8];
                    double h21 = voxel[12 + i * 8];
                    double h22 = voxel[13 + i * 8];
                    
                    double sqrt = Math.sqrt( 4.0 * h12 * h21 + (h11 - h22) * (h11 - h22) );
                    
                    double ev1 = 0.5 * ( (h11 + h22) + sqrt );
                    double ev2 = 0.5 * ( (h11 + h22) - sqrt );
                    
                    double fmean = voxel[4];
                    
                    // if abs(ev1) or abs(ev2) or f are < 0, or NaN, or Infinity, don't fit anything
                    if ( !(Math.abs(ev1) > 0.0 && Math.abs(ev2) > 0.0 && fmean > 0.0) ) {
                        
                        // output zero concentration
                        
                        // peak axes are unit axes
                        
                        vectors[i][0] = 1.0;
                        vectors[i][1] = 0.0;
                        vectors[i][2] = 0.0;
                        vectors[i][3] = 0.0;
                        vectors[i][4] = 1.0;
                        vectors[i][5] = 0.0;
                        vectors[i][6] = 0.0;
                        vectors[i][7] = 0.0;
                        vectors[i][8] = 1.0;
                        
                        continue pdloop;
                        
                    }
                    
                    
                    /**
                     *
                     * hessian = [d^2f/ds^2 d^2f/dsdt]
                     * [d^2f/dtds d^2f/dt^2]
                     *
                     * = [H00 H01]
                     * [H10 H11]
                     *
                     * where s and t are orthogonal coordinates local to the peak.
                     *
                     * s and t are actually distances along the coordinate system with axes parallel
                     * to the vectors t_1 = (y, -x, 0) and t_2 = (-xz, -yz, x^2+y^2), where (x, y, z)
                     * is the peak direction so that x^2+y^2+z^2 = 1. Suppose \hat{t}_i = t_i/|t_i|
                     * are unit vectors along the axes parallel to t_1 and t_2. If an eigenvector of
                     * your Hessian matrix is (a, b), then the direction along that axis of the peak
                     * in the usual 3D coordinates is parallel to e = a \hat{t}_1 + b \hat{t}_2; ie
                     * the direction you need is \hat{e} = e/|e|.
                     *
                     * (PAC) - the axes t_1 and t_2 are the rotation of (0, -1, 0) and (-1, 0, 0) to the peak
                     * direction (x, y, z). When the peak is (0, 0, 1), these rotations are undefined so we use
                     * the unrotated axes.
                     *
                     */
                    
                    // axes of hessian coordinate system
                    Vector3D t1;
                    Vector3D t2;
                    
                    if((peak.y*peak.y + peak.x*peak.x)==0) {
                        t1 = new Vector3D(0.0,-1.0,0.0);
                        t2 = new Vector3D(-1.0,0.0,0.0);
                    }
                    else{
                        t1 = new Vector3D(peak.y, -1.0 * peak.x, 0.0).normalized();
                        
                        t2 = new Vector3D(peak.x * peak.z * (-1.0),
                        peak.y * peak.z * (-1.0),
                        peak.x * peak.x +
                        peak.y * peak.y).normalized();
                    }
                    // eigenvectors of hessian in 3D
                    double x1 = 1.0;
                    double y1 = 0.0;
                    if(h12!=0.0){
                        x1 = 1.0;
                        y1 = (ev1 - h11) / h12;
                    }
                    
                    double x2 = 0.0;
                    double y2 = 1.0;
                    if((ev2 - h22) != 0.0){
                        x2 = 1.0;
                        y2 = h21 / (ev2 - h22);
                    }
                    
                    Vector3D e1 =
                    t1.scaled(x1).plus(t2.scaled(y1));
                    
                    Vector3D e2 =
                    t1.scaled(x2).plus(t2.scaled(y2));
                    
                    e1 = e1.normalized();
                    e2 = e2.normalized();
                    
                    // smallest eigenvalue first
                    // largest Hessian eigenvalue has highest concentration, and thus is last
                    // this might be confusing because the largest hessian eigenvector corresponds
                    // to the smallest eigenvector of the scatter matrix of samples from the PDF
                    // this is because Hessian measures curvature and a large Hessian means
                    // low uncertainty
                    if (Math.abs(ev1) > Math.abs(ev2)) {
                        Vector3D tmp = e1;
                        e1 = e2;
                        e2 = tmp;
                        
                        double tmp2 = ev1;
                        ev1 = ev2;
                        ev2 = tmp2;
                    }
                    
                    // PDF axes are {peak e1 e2}
                    
                    vectors[i][0] = peak.x;
                    vectors[i][1] = peak.y;
                    vectors[i][2] = peak.z;
                    vectors[i][3] = e1.x;
                    vectors[i][4] = e1.y;
                    vectors[i][5] = e1.z;
                    vectors[i][6] = e2.x;
                    vectors[i][7] = e2.y;
                    vectors[i][8] = e2.z;
                    
                    
                    if (watson) {
                        // normalise vals by mean of function
                        double fx = (Math.abs(ev1) + Math.abs(ev2)) / fmean;
                        
                        if (log) {
                            fx = Math.log(fx);
                            params[i][0] = functions[i][0].evaluate(new double[] {fx});
                        }
                        else {
                            params[i][0] = functions[i][0].evaluate(new double[] {fx});
                        }
                        
                        // if kappa out of range, impose limit
                        if(params[i][0]>1E5){
                            logger.warning("kappa out of range, value: "+params[i][0]+ ".  Setting to 1E5");
                            params[i][0]=1E5;
                        }
                        
                    }
                    else if (bingham) {
                        double fx = Math.abs(ev2) / fmean; // corresponds to ev2 (biggest)
                        double fy = Math.abs(ev1) / fmean; // corresponds to ev1
                        
                        if (log) {
                            fx = Math.log(fx);
                            fy = Math.log(fy);
                        }
                        
                        params[i][0] =
                        -1.0 * Math.exp(functions[i][0].evaluate(new double[] {fx, fy}));
                        
                        params[i][1] =
                        -1.0 * Math.exp(functions[i][1].evaluate(new double[] {fx, fy}));
                        
                        if(Math.abs(params[i][0])<Math.abs(params[i][1])) {
                            params[i][0]=params[i][1];
                        }
                        
                        
                        // check for errors
                        if (Double.isNaN(params[i][0]) || Double.isNaN(params[i][1])) {
                            logger.warning("NaN kappa detected, setting concentration parameters to zero");
                            
                            logger.info("NaN kappa detected. PD index == " + i + "\nfx == " + fx + "\nfy == " +
                            fy + "\n log == " + log + "\nfunctions[i][0] == " +
                            functions[i][0].evaluate(new double[] {fx, fy}) + "\nfunctions[i][1] == " +
                            functions[i][1].evaluate(new double[] {fx, fy}));
                            
                            params[i][0] = 0.0;
                            params[i][1] = 0.0;
                        }
                        
                        
                        if(params[i][0]>0)
                            params[i][0]=0;
                        else if(params[i][0]<-1E5){
                            logger.warning("kappa1 out of range, value: "+params[i][0]+ ".  Setting to -1E5");
                            params[i][0]=-1E5;
                        }
                        
                        if(params[i][1]>0)
                            params[i][1]=0;
                        else if(params[i][1]<-1E5){
                            logger.warning("kappa2 out of range, value: "+params[i][0]+ ".  Setting to -1E5");
                            params[i][1]=-1E5;
                        }
                        
                        if(Math.abs(params[i][0])<Math.abs(params[i][1]))
                            params[i][0]=params[i][1];
                        
                    }
                    else { // ACG
                        
                        double fx = Math.abs(ev2) / fmean; // s1, corresponds to ev2
                        double fy = Math.abs(ev1) / fmean; // s2, corresponds to ev1
                        
                        if (log) {
                            fx = Math.log(fx);
                            fy = Math.log(fy);
                            
                            params[i][0] =
                            Math.exp(functions[i][0].evaluate(new double[] {fx, fy}));
                            
                            params[i][1] =
                            Math.exp(functions[i][1].evaluate(new double[] {fx, fy}));
                            
                        }
                        else {
                            
                            params[i][0] =
                            functions[i][0].evaluate(new double[] {fx, fy});
                            
                            params[i][1] =
                            functions[i][1].evaluate(new double[] {fx, fy});
                            
                            params[i][2] = 1.0;
                        }
                        
                    }
                    
                }
                
                
                om.output(new double[] {numPDs});
                
                for (int i = 0; i < CL_Initializer.numPDsIO; i++) {
                    if (numPDs == 0) {
                        om.output(new double[] {0.0});
                    }
                    else {
                        om.output(new double[] {1.0 / (double)numPDs});
                    }
                    
                    om.output(vectors[i]);
                    om.output(params[i]);
                }
                
        }
        
    }
    
    
    /**
     * Get PICo PDFs for DT data.
     */
    public static void mapDTData(DataSource data, LookupTable[] luts, OutputManager om) {
        
        while (data.more()) {
            
            double[] voxel = data.nextVoxel();
            
            double[][] params = new double[CL_Initializer.numPDsIO][paramsPerPD];
            
            double[][] vectors = new double[CL_Initializer.numPDsIO][9];
            // all three vectors in one array
            
            double[] mix = new double[CL_Initializer.numPDsIO];
            
            double numPDs = 0.0;
            
            if (voxel[0] < 0.0) {
                // background or bad data, do not fit parameters
                
            }
            else if (voxel.length == CL_Initializer.ONETENVALSPERVOXEL) {
                
                numPDs = 1.0;
                
                mix[0] = 1.0;
                
                boolean gotNaN = false;
                
                // should really check this for all images, but occurs most often in
                // dt images that have been transformed
                
                for (int i = 0; i < 6; i++) {
                    if ( Double.isNaN(voxel[2+i]) ) {
                        gotNaN = true;
                    }
                }
                
                if (gotNaN) {
                    numPDs = 0.0;
                }
                else {
                    
                    DT tensor = new DT(voxel[2], voxel[3], voxel[4],
                    voxel[5], voxel[6], voxel[7]);
                    
                    
                    // default LUT is based on L1 / L3 and L2 / L3
                    boolean fa = false;
                    
                    if (luts[0].yMin() == luts[0].yMax()) {
                        fa = true;
                    }
                    
                    double[][] seig = tensor.sortedEigenSystem();
                    
                    for (int i = 0; i < 3; i++) {
                        vectors[0][3*i+0] = seig[1][i];
                        vectors[0][3*i+1] = seig[2][i];
                        vectors[0][3*i+2] = seig[3][i];
                    }
                    
                    double y = seig[0][0] / seig[0][2];
                    double z = seig[0][1] / seig[0][2];
                    
                    if ( Double.isNaN(y) || Double.isInfinite(y) || Double.isNaN(z) ||
                    Double.isInfinite(z) ) {
                        
                        // do nothing, output is zeros
                    }
                    else {
                        
                        if (fa) {
                            
                            double f = tensor.fa();
                            
                            if (f >= 0.0 && f < 1.0) {
                                try {
                                    params[0] = luts[0].getValues(f, true, true);
                                }
                                catch (OutsideLUTRangeException e) {
                                    throw new LoggedException("Got LUT range exception with clamping on: " + e);
                                }
                                
                            }
                            
                        }
                        else {
                            try {
                                params[0] = luts[0].getValues(y, z, true, true);
                            }
                            catch (OutsideLUTRangeException e) {
                                throw new LoggedException("Got LUT range exception with clamping on: " + e);
                            }
                        }
                        
                    } // else (eigenvalues not infinite, zero, or NaN)
                    
                }  // matches else (no values are NaN)
                
            }
            else if (voxel.length == CL_Initializer.TWOTENVALSPERVOXEL) {
                
                numPDs = voxel[2];
                
                if (numPDs == 0.0) {
                    // do nothing, output is zeros
                    
                }
                else if (numPDs == 1.0) {
                    
                    DT dt1 = new DT(voxel[4], voxel[5], voxel[6], voxel[7], voxel[8], voxel[9]);
                    
                    
                    double[][] seig1 = dt1.sortedEigenSystem();
                    
                    for (int i = 0; i < 3; i++) {
                        vectors[0][3*i+0] = seig1[1][i];
                        vectors[0][3*i+1] = seig1[2][i];
                        vectors[0][3*i+2] = seig1[3][i];
                    }
                    
                    mix[0] = voxel[3];
                    mix[1] = voxel[10];
                    
                    double y = seig1[0][0] / seig1[0][2];
                    double z = seig1[0][1] / seig1[0][2];
                    
                    if ( Double.isNaN(y) || Double.isInfinite(y) || Double.isNaN(z) ||
                    Double.isInfinite(z) ) {
                        
                        // do nothing, output is zeros
                    }
                    else {
                        try {
                            params[0] = luts[0].getValues(y, z, true, true);
                        }
                        catch (OutsideLUTRangeException e) {
                            throw new LoggedException("Got LUT range exception with clamping on: "
                            + e);
                        }
                    }
                    
                }
                else if (numPDs == 2.0) {
                    
                    DT dt1 = new DT(voxel[4], voxel[5], voxel[6], voxel[7], voxel[8], voxel[9]);
                    
                    DT dt2 = new DT(voxel[11], voxel[12], voxel[13],
                    voxel[14], voxel[15], voxel[16]);
                    
                    double fa1 = dt1.fa();
                    double fa2 = dt2.fa();
                    
                    double[][] seig1 = dt1.sortedEigenSystem();
                    double[][] seig2 = dt2.sortedEigenSystem();
                    
                    mix[0] = voxel[3];
                    mix[1] = voxel[10];
                    
                    Vector3D[] evecs1 = new Vector3D[3];
                    Vector3D[] evecs2 = new Vector3D[3];
                    
                    // first eigenvector is the e1 of each tensor
                    evecs1[0] = new Vector3D(seig1[1][0], seig1[2][0], seig1[3][0]);
                    evecs2[0] = new Vector3D(seig2[1][0], seig2[2][0], seig2[3][0]);
                    
                    double dot = Math.abs(evecs1[0].dot(evecs2[0]));
                    double cross = Math.acos(dot);
                    
                    // check for fitting failure, zero crossing angle
                    if(dot > 0.9999) {
                        // fitting failure, output identical PDFs
                        evecs1[1] = new Vector3D(seig1[1][1], seig1[2][1], seig1[3][1]);
                        evecs2[1] = new Vector3D(seig2[1][1], seig2[2][1], seig2[3][1]);
                        
                        evecs1[2] = new Vector3D(seig1[1][2], seig1[2][2], seig1[3][2]);
                        evecs2[2] = new Vector3D(seig2[1][2], seig2[2][2], seig2[3][2]);
                        
                    }
                    else {
                        // third eigenvector is cross product of the two firsts,
                        // ie normal to plane of crossing
                        evecs1[2] = evecs1[0].cross(evecs2[0]).normalized();
                        evecs2[2] = evecs1[2];
                        
                        // second is the first rotated about the third by 90 degrees
                        evecs1[1] = Rotations.rotateVector(evecs1[0], evecs1[2], Math.PI / 2.0);
                        evecs2[1] = Rotations.rotateVector(evecs2[0], evecs2[2], Math.PI / 2.0);
                    }
                    
                    if (fa1 < 0.0 || fa2 < 0.0 || fa1 > 1.0 || fa2 > 1.0) {
                        // output zeros
                    }
                    else if ( Double.isNaN(fa1) || Double.isInfinite(fa1) || Double.isNaN(fa2) ||
                    Double.isInfinite(fa2) ) {
                        // output zeros
                    }
                    else if ( Double.isNaN(cross) || Double.isInfinite(cross) ) {
                        // output zeros
                    }
                    else try {
                        double[] tmp = luts[1].getValues(fa1, fa2, cross, true, true);
                        
                        for (int p = 0; p < 2; p++) {
                            for (int i = 0; i < paramsPerPD; i++) {
                                params[p][i] = tmp[p * paramsPerPD + i];
                            }
                        }
                    }
                    catch (OutsideLUTRangeException e) {
                        throw new LoggedException("Got LUT range exception with clamping on: " + e);
                    }
                    
                    // write evecs into array
                    for (int i = 0; i < 3; i++) {
                        vectors[0][3*i+0] = evecs1[i].x;
                        vectors[0][3*i+1] = evecs1[i].y;
                        vectors[0][3*i+2] = evecs1[i].z;
                        
                        vectors[1][3*i+0] = evecs2[i].x;
                        vectors[1][3*i+1] = evecs2[i].y;
                        vectors[1][3*i+2] = evecs2[i].z;
                    }
                    
                    
                }
                else {
                    throw new LoggedException("Invalid number of PDs, " + numPDs + ", in voxel");
                }
                
                
                
            }
            else if (voxel.length == CL_Initializer.THREETENVALSPERVOXEL) {
                throw new LoggedException("Three tensor PICo is not available.");
            }
            else {
                throw new LoggedException("Got unknown data with " + voxel.length +
                " values per voxel");
            }
            
            om.output(new double[] {numPDs});
            
            for (int i = 0; i < CL_Initializer.numPDsIO; i++) {
                om.output(new double[] {mix[i]});
                
                om.output(vectors[i]);
                om.output(params[i]);
            }
            
        }
        
    }
    
    
    public void initOptions(String[] args) {
        //public static void main(String[] args) {
        
        // numPDsIO set by CL_Initializer according to inputmodel
        
        CL_Initializer.inputDataType = "double";
        CL_Initializer.CL_init(args);
        
        // set up output
        //   om = new OutputManager();      
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-luts")) {
                int t = 0;
                while (t + i + 1 < args.length && !args[t+i+1].startsWith("-")) {
                    lutFiles[t] = args[i+t+1];
                    t++;
                }
                CL_Initializer.markAsParsed(i, t+1);
            }
            else if (args[i].equals("-flipx")) {
                flipX = true;
                CL_Initializer.markAsParsed(i);
            }
            else if (args[i].equals("-flipy")) {
                flipY = true;
                CL_Initializer.markAsParsed(i);
                
            }
            else if (args[i].equals("-flipz")) {
                flipZ = true;
                CL_Initializer.markAsParsed(i);
            }
            else if (args[i].equals("-pdf")) {
                pdf = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-directmap")) {
                log = false;
                CL_Initializer.markAsParsed(i);
            }
            
        }
        
        CL_Initializer.checkParsing(args);
    }
    
    public void initVariables(){
        
        valuesPerVoxel = 0;
        
        numLUTs = 1;
        
    }
    
    
    public void execute(OutputManager om){
        
        if (CL_Initializer.inputModel == null) {
            throw new LoggedException("This program requires an input model. Specify one with the -inputmodel option");
        }
        
        if (CL_Initializer.inputModel.equals("pds")) {
            
            valuesPerVoxel = 6 + 8 * CL_Initializer.numPDsIO;
            
            DataSource data =
                ExternalDataSource.getDataSource(CL_Initializer.inputFile, valuesPerVoxel, CL_Initializer.inputDataType);
            
            mapSFData(data, pdf, log, lutFiles, om);
        }
        else {
            
            
            if (CL_Initializer.inputModel.equals("dt")) {
                valuesPerVoxel = CL_Initializer.ONETENVALSPERVOXEL;
            }
            else if (CL_Initializer.inputModel.equals("multitensor")) {
                valuesPerVoxel = CL_Initializer.TWOTENVALSPERVOXEL;
                numLUTs = 2;
                
                if (CL_Initializer.maxTensorComponents > 2) {
                    throw new LoggedException("Attempted to use unsupported tensor input model with " +
                    CL_Initializer.maxTensorComponents + " components");
                }
            }
            else {
                throw new LoggedException("Unsupported input model " + CL_Initializer.inputModel);
            }
            
            DataSource data =
                ExternalDataSource.getDataSource(CL_Initializer.inputFile, valuesPerVoxel, CL_Initializer.inputDataType);
            
            
            LookupTable[] luts = new LookupTable[numLUTs];
            
            luts[0] = LookupTable.readLUT(lutFiles[0]);
            
            paramsPerPD = luts[0].valuesPerPosition();
            
            for (int i = 1; i < luts.length; i++) {
                luts[i] = LookupTable.readLUT(lutFiles[i]);
                
                if (luts[i].valuesPerPosition() != (i+1) * paramsPerPD) {
                    throw new LoggedException("LUTs must be the same PDF type for all data");
                }
            }
            
            mapDTData(data, luts, om);
            
        }
        
        om.close();
        
        
        
    }
    
}
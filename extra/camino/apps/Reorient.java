package apps;

import java.util.logging.Logger;

import tools.*;
import misc.*;
import imaging.*;
import data.*;
import numerics.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Does diffusion tensor reorientation to accompany an image
 * transformation.
 * 
 * <dt>Description:
 * 
 * <dd>Will do PPD or FS reorientation, specified by the
 * -reorientation option, for affine transformations specified with
 * the -trans option with an affine transformation in the format
 * output by FSL's flirt.
 * 
 * 
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class Reorient {

    private static Logger logger = Logger.getLogger("camino.apps.Reorient");


    /**
     * Linear transformation matrix
     */
    private static RealMatrix linear = null;

    /**
     * Rotation matrix for global FS.
     */
    private static RealMatrix commonRot = null;

    /**
     * X-component of warp.
     */
    private static ScalarImage transX = null;

    /**
     * Y-component of warp.
     */
    private static ScalarImage transY = null;

    /**
     * Y-component of warp.
     */
    private static ScalarImage transZ = null;

    /**
     * Flag for linear transformation.
     */
    private static boolean linearTrans = false;

    /**
     * Flag for non-linear transformation.
     */
    private static boolean nonlinearTrans = false;


    public static void main(String[] args) {

        // Read and process command line options.
        initOptions(args);

	// Set up output 
        OutputManager om = new OutputManager();

        // Initialize all the variables that control what this
        // program does.
        initVariables();

        // Do the processing.
        execute(om);

        // Tidy up.
        om.close();
    }

    public static void initOptions(String[] args) {

        // Input type defaults to single DT. Input datatype defaults to
        // double.
        CL_Initializer.inputModel = "dt";
        CL_Initializer.inputDataType = "double";

        CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
        // This program only deals with tensor data.
        CL_Initializer.initTensorDataSource();

    }


    public static void initVariables() {

        // Get the image transformations
        if(CL_Initializer.transformFile == null && CL_Initializer.transformFileX == null) {
            throw new LoggedException("No transformation specified.  Use the -trans option.");
        }

        // Load linear part if specified.
        linear = new RealMatrix(3, 3);

        // If it is an image file, it is a non-linear warp.
        if(CL_Initializer.transformFile != null && (!ImageHeader.imageExists(CL_Initializer.transformFile))) {
            RealMatrix affine = readFSL_Affine(CL_Initializer.transformFile);
            for(int i=0; i<3; i++) {
                for(int j=0; j<3; j++) {
                    linear.entries[i][j] = affine.entries[i][j];
                }
            }
            linearTrans = true;
        }

        // Compute the finite strain rotation for the linear
        // transformation if required.
        commonRot = new RealMatrix(3,3);
        if(linearTrans && CL_Initializer.reorientation.equals("fs")) {
            // We need to transpose the result here as DT.iTransform
            // computes R D R^T rather than R^T D R
            commonRot = getFS_Reorientation(linear).transpose();
        }
       

        // Load non-linear part if specified.
        if(CL_Initializer.transformFileX != null) try {
            transX = new ScalarImage(CL_Initializer.transformFileX);
            transY = new ScalarImage(CL_Initializer.transformFileY);
            transZ = new ScalarImage(CL_Initializer.transformFileZ);
            nonlinearTrans = true;
        }
        catch(Exception e) {
            LoggedException.logExceptionSevere(e, Thread.currentThread().getName());
            System.exit(1);
        }
        else if (CL_Initializer.transformFile != null && ImageHeader.imageExists(CL_Initializer.transformFile)) try {
            // Assume the transformation is in a single multi-component file
            // as output by FSLs fnirt.
            transX = new ScalarImage(CL_Initializer.transformFile, 0);
            transY = new ScalarImage(CL_Initializer.transformFile, 1);
            transZ = new ScalarImage(CL_Initializer.transformFile, 2);
            nonlinearTrans = true;
        }
        catch(Exception e) {
            LoggedException.logExceptionSevere(e, Thread.currentThread().getName());
            System.exit(1);
        }
            
            

        // If a non-linear transformation is specified, we need
        // the image dimensions, so check they have been provided.
        if(nonlinearTrans && (CL_Initializer.dataDims[0] == 0 || CL_Initializer.voxelDims[0] == 0)) {
            throw new LoggedException("Must specify -datadims and -voxeldims for non-linear transformations.");
        }

    }


    public static void execute(OutputManager om) {

        // Loop over the data
        int voxX = 0;
        int voxY = 0;
        int voxZ = 0;
        while (CL_Initializer.data.more())
            try {

                double[] vox = CL_Initializer.data.nextVoxel();
                DT[] dts = FracAnis.getTensorList(vox, CL_Initializer.inputModel);

                // Check if the voxel is background.  
     		boolean notBG = !ModelFit.isBG(Math.exp(vox[1]));
		
		if (vox[0]>=0 && notBG) {

                    // Reorient according to the linear
                    // transformation.
                    if(linearTrans) {
                        for (int i = 0; i < dts.length; i++) {
                            if(CL_Initializer.reorientation.equals("ppd")) {
                                dts[i] = dts[i].ppd(linear);
                            }
                            else if(CL_Initializer.reorientation.equals("fs")) {
                                dts[i] = dts[i].iTransform(commonRot);
                            }
                            else {
                                throw new LoggedException("Unrecognized reorientation strategy: " + CL_Initializer.reorientation);
                            }
                        }
                    }

                    // Reorient according to the non-linear
                    // transformation.
                    if(nonlinearTrans) {

                        for (int i = 0; i < dts.length; i++) {

                            // Get the local Jacobian
                            double posX = (voxX+0.5)*CL_Initializer.voxelDims[0];
                            double posY = (voxY+0.5)*CL_Initializer.voxelDims[1];
                            double posZ = (voxZ+0.5)*CL_Initializer.voxelDims[2];
                            RealMatrix localLinear = getJacobian(transX, transY, transZ, posX, posY, posZ);

                            // We need the inverse transform, since
                            // the warp maps the target to the source.
                            RealMatrix iLin = localLinear.inverse();

                            //System.out.println(localLinear);
//                             Point3D mmPoint = new Point3D(posX, posY, posZ);
//                             double dispX = transX.valueAt(mmPoint);
//                             double dispY = transY.valueAt(mmPoint);
//                             double dispZ = transZ.valueAt(mmPoint);
//                             System.out.println(voxX + " " + voxY + " " + voxZ + " " + dispX + " " + dispY + " " + dispZ);
                            
                            if(CL_Initializer.reorientation.equals("ppd")) {
                                dts[i] = dts[i].ppd(iLin);
                            }
                            else if(CL_Initializer.reorientation.equals("fs")) {
                                RealMatrix rot = getFS_Reorientation(iLin).transpose();
                                dts[i] = dts[i].iTransform(rot);
                            }
                            else {
                                throw new LoggedException("Unrecognized reorientation strategy: " + CL_Initializer.reorientation);
                            }
                        }
                    }

                    replaceTensorList(vox, dts, CL_Initializer.inputModel);
                }

                om.output(vox);

                // Need to keep track of the current
                // position in the image.
                voxX += 1;
                if(voxX == CL_Initializer.dataDims[0]) {
                    voxX = 0;
                    voxY += 1;
                    if(voxY == CL_Initializer.dataDims[1]) {
                        voxY = 0;
                        voxZ += 1;
                    }
                }
                else if(nonlinearTrans && voxZ == CL_Initializer.dataDims[2]) {
                            logger.warning("More voxels on input than -datadims specifies.  Probable error.");
                }

            }
            catch (Exception e) {
                LoggedException.logExceptionSevere(e, Thread.currentThread().getName()); 
            }
    }


    /**
     * Reads an affine transformation from a file in the format output
     * by FSL's flirt program into a matrix object.
     *
     * @param filename
     *            The name of the file containing the transformation.
     *
     * @return
     *            A 4x4 matrix containing the transformation.
     */
    public static RealMatrix readFSL_Affine(String filename) {
        FileInput f = new FileInput(filename);
        RealMatrix aff = new RealMatrix(4,4);
        for(int i=0; i<4; i++) try {
            String[] nums = f.readString().split("\\s+");
            for(int j=0; j<4; j++) {
                aff.entries[i][j] = Double.parseDouble(nums[j]);
            }
        }
        catch(Exception e) {
            throw new LoggedException("Could not parse transformation file: " + filename + e);
        }

        return aff;
    }
            

    /**
     * Computes the finite strain strategy reorientation matrix.
     *
     * @param linear  The linear part F of the affine transformation.
     *
     * @return The rotation matrix R = [F F^T]^-0.5 F
     */
    public static RealMatrix getFS_Reorientation(RealMatrix F) {
        RealMatrix fft = F.product(F.transpose());
        RealMatrix[] esys = fft.jacobi();
        
        // Raise to the power -0.5.
        for(int i=0; i<esys[0].rows(); i++) {
            if(esys[0].entries[i][i] == 0.0) {
                throw new LoggedException("Transformation is singular: " + F);
            }
            esys[0].entries[i][i] = 1.0/Math.sqrt(esys[0].entries[i][i]);
        }

        // Reconstruct the matrix raised to power -0.5.
        RealMatrix fftp = esys[1].product(esys[0].product(esys[1].transpose()));
        return fftp.product(F);
    }


    /**
     * Computes the Jacobian of the image transformation at point
     * (posX, posY, posZ).
     *
     * @param transX X-Component of warp.
     *
     * @param transY Y-Component of warp.
     *
     * @param transZ Z-Component of warp.
     *
     * @param posX X-coord of point.
     *
     * @param posY Y-coord of point.
     *
     * @param posZ Z-coord of point.
     *
     * @return The jacobian matrix [[dTx/dx dTx/dy dTx/dz]; [dTy/dx
     * dTy/dy dTy/dz]; [dTz/dx dTz/dy dTz/dz]].
     */
    public static RealMatrix getJacobian(ScalarImage transX, ScalarImage transY, ScalarImage transZ, double posX, double posY, double posZ) {
        Point3D p = new Point3D(posX, posY, posZ);
        double[] dTxdp = transX.derivAt(p);
        double[] dTydp = transY.derivAt(p);
        double[] dTzdp = transZ.derivAt(p);

        RealMatrix J = new RealMatrix(3, 3);
        for(int i=0; i<3; i++) {
            J.entries[0][i] = dTxdp[i];
            J.entries[1][i] = dTydp[i];
            J.entries[2][i] = dTzdp[i];
            
            // Need to add the identity for the jacobian.
            J.entries[i][i] += 1.0;

        }

        return J;
    }


    /**
     * Puts a list of diffusion tensors in a single voxel into a data
     * array for output.
     * 
     * @param vox
     *            The data.
     * 
     * @param dts 
     *            The array of diffusion tensors.
     *
     * @param model
     *            The model.
     */
    public static void replaceTensorList(double[] vox, DT[] dts, String model) {

        if (model.equals("dt")) {
            double[] comp1 = dts[0].getComponents();
            for(int i=0; i<6; i++) {
                vox[i+2] = comp1[i];
            }
        }
        else if (model.equals("twotensor")) {
            double[] comp1 = dts[0].getComponents();
            for(int i=0; i<6; i++) {
                vox[i+4] = comp1[i];
            }
            double[] comp2 = dts[1].getComponents();
            for(int i=0; i<6; i++) {
                vox[i+11] = comp2[i];
            }
        }
        else if (model.equals("threetensor")) {
            double[] comp1 = dts[0].getComponents();
            for(int i=0; i<6; i++) {
                vox[i+4] = comp1[i];
            }
            double[] comp2 = dts[1].getComponents();
            for(int i=0; i<6; i++) {
                vox[i+11] = comp2[i];
            }
            double[] comp3 = dts[2].getComponents();
            for(int i=0; i<6; i++) {
                vox[i+18] = comp3[i];
            }
        }
        else if (model.equals("multitensor")) {
            for(int i=0; i<dts.length; i++) {
                double[] comp = dts[i].getComponents();
                for(int j=0; j<6; j++) {
                    vox[4+7*i+j] = comp[j];
                }
            }
        }
        else {
            throw new LoggedException("Cannot extract tensors from model: " + model);
        }

    }

}

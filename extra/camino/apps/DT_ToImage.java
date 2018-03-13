package apps;

import data.*;
import imaging.*;
import misc.*;
import tools.*;

import java.io.IOException;

/**
 * Converts Camino dts to image format.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class DT_ToImage extends Executable {
    

    private ImageHeader header;

    private int xDataDim;
    private int yDataDim;
    private int zDataDim;

    private int numDTs;
    

    private static String outputRoot = "camino_";
    
    /** Header for diffusion data, used to initialize tensor header.*/
    private static String headerFile = null;

    public DT_ToImage(String[] args) {
        super(args);
    }
    
    public void initDefaultVals() {	
        xDataDim = 0;
        yDataDim = 0;
        zDataDim = 0;

        header = null;

        numDTs = 1;
    }
    
    public void initOptions(String[] args) {
        
        CL_Initializer.inputDataType = "double";
        OutputManager.outputDataType = "float";

        CL_Initializer.inputModel = "dt";

        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-outputroot")) {
                outputRoot = args[i+1];
                CL_Initializer.markAsParsed(i, 2);
            }
        }
        
        CL_Initializer.checkParsing(args);
        
        CL_Initializer.initTensorDataSource();

        numDTs = CL_Initializer.maxTensorComponents;

        CL_Initializer.initInputSpaceAndHeaderOptions();

        header = CL_Initializer.headerTemplate;
        
        xDataDim = header.xDataDim();
        yDataDim = header.yDataDim();
        zDataDim = header.zDataDim();

    }
    
    public void initVariables() {
    }
	
    public void execute(OutputManager om) {	
      
           
        // data array for image, in Nifti format
        double[][][][][] tensorData = new double[numDTs][xDataDim][yDataDim][zDataDim][6];
            
        double[][][] exit = new double[xDataDim][yDataDim][zDataDim];
            
        double[][][] lnS0 = new double[xDataDim][yDataDim][zDataDim];
            
        double[][][][] mix = new double[numDTs][xDataDim][yDataDim][zDataDim];

            
        for (int k = 0; k < zDataDim; k++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
                    double[] voxel = CL_Initializer.data.nextVoxel();
                        
                    exit[i][j][k] = voxel[0];
                    lnS0[i][j][k] = voxel[1];

                    DT[] dts = FracAnis.getTensorList(voxel, CL_Initializer.inputModel);
                    
                    for (int c = 0; c < numDTs; c++) {
                        double[] comps = dts[c].getComponents();
                        
                        System.arraycopy(comps, 0, tensorData[c][i][j][k], 0, 6);

                        if (numDTs > 1) {
                            mix[c][i][j][k] = voxel[3+c*7];
                        }

                    }
                }
            }
            
        }
                
            
        String tensorFileRoot = outputRoot + "dt";

        
        String mixFileRoot = outputRoot + "mix";
        
        if (numDTs > 1) {
            for (int c = 0; c < numDTs; c++) {
                header.writeTensorImage(tensorData[c], tensorFileRoot + (c+1));

                header.writeScalarImage(mix[c], mixFileRoot + (c+1));
            }
        }
        else {
            header.writeTensorImage(tensorData[0], tensorFileRoot);
        }
                                    
        String exitFileRoot = outputRoot + "exitcode";
        
        header.writeScalarImage(exit, exitFileRoot);
    
        String lnS0FileRoot = outputRoot + "lns0";
            
        header.writeScalarImage(lnS0, lnS0FileRoot);
    }  
 
}

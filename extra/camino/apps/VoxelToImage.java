package apps;

import data.*;
import imaging.*;
import misc.*;
import tools.*;

import java.io.*;

import java.text.*;

import java.util.logging.Logger;


/**
 * Converts Camino voxel-order files to 3D images for visualization.
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class VoxelToImage extends Executable {
    
    public VoxelToImage(String[] args) {
        super(args);
    }
    
    private String outputRoot;
    private int components;

    private ImageHeader header;

    private int xDataDim;
    private int yDataDim;
    private int zDataDim;
    

    private boolean outputVector;
    

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.VoxelToImage");
    
    
    public void initDefaultVals(){
	outputRoot = "voxel2image";
	components = 1;

	header = null;

        outputVector = false;
	
	xDataDim = 0;
	yDataDim = 0;
	zDataDim = 0;
    }
    
    public void initOptions(String[] args){
        CL_Initializer.inputDataType = "double";
        OutputManager.outputDataType = "float";
        CL_Initializer.CL_init(args);


        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-outputroot")) {
		outputRoot = args[i+1];
		CL_Initializer.markAsParsed(i,2);
            }
	    if (args[i].equals("-components")) {
                components = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i,2);
            }
	    if (args[i].equals("-header")) {
		CL_Initializer.headerTemplateFile = args[i+1];
		CL_Initializer.markAsParsed(i,2);
            }
            if (args[i].equals("-outputvector")) {
                outputVector = true;
		CL_Initializer.markAsParsed(i);
            }
	    
	    
        }
        
        CL_Initializer.checkParsing(args);
	
	CL_Initializer.initInputSpaceAndHeaderOptions();
    
	header = CL_Initializer.headerTemplate;

	xDataDim = header.xDataDim();
	yDataDim = header.yDataDim();
	zDataDim = header.zDataDim();

    }
    
    
    public void initVariables() {
    
    }
    
    public void execute(OutputManager om) {
        
	VoxelOrderDataSource input = 
	    new VoxelOrderDataSource(CL_Initializer.inputFile, components, CL_Initializer.inputDataType);
	

        double[][][][] data;

        if (outputVector) {
            data = new double[xDataDim][yDataDim][zDataDim][components];
        }
        else {

            data = new double[components][xDataDim][yDataDim][zDataDim];

        }

        for (int k = 0; k < zDataDim; k++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
                    
                    double[] voxel = input.nextVoxel();
                    
                    if (outputVector) {
                        data[i][j][k] = voxel;
                    }
                    else {
                        
                        for (int n = 0; n < components; n++) {
                            data[n][i][j][k] = voxel[n];
                        }

                    }
                    
                }
            }
        }
            

        if (outputVector) {
            header.writeVectorImage(data, outputRoot);
        }
        else {
            DecimalFormat df = new DecimalFormat("0000");
            
            for (int n = 0; n < components; n++) {
                
                String fileRoot = outputRoot;
                
                if (components > 1) {
                    fileRoot = outputRoot + df.format(n+1);
                }

                header.writeScalarImage(data[n], fileRoot);
            }
            
        }
	
        
    }
    
    
}



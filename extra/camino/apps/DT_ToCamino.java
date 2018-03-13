package apps;

import data.*;
import imaging.*;
import tools.*;

/**
 * Converts a NIFTI, ITK or Analyze tensor volume to Camino format. Just does a straight conversion,
 * does not alter the coordinate system of the data. Scaling is applied according to the header,
 * though additional scaling may be applied.
 *  
 * @author Philip Cook
 * @version $Id$
 */
public class DT_ToCamino extends Executable {

    public DT_ToCamino(String[] args) {
        super(args);
    }
    
    private DataSource tensorSource;
    private DataSource s0Source;
    
    private String s0File;
    private boolean fileIsLnS0;
    private double scaleSlope;
    private double scaleInter;
    
    // true for ITK, false for NIFTI, undefined for Analyze
    private boolean upperTriangular;
    private ImageHeader dtHeader;
	
    public void initDefaultVals() {
	tensorSource = null;
	s0Source = null;
        s0File = null;
	fileIsLnS0 = false;
        scaleSlope = 1.0;
	scaleInter = 0.0;
	// true for ITK, false for NIFTI, undefined for Analyze
	upperTriangular = false;
    }

    public void initOptions(String[] args) {
	
	CL_Initializer.inputDataType = "double";
	CL_Initializer.CL_init(args);
	CL_Initializer.initMaskSource();	
	
	for (int i = 0; i < args.length; i++) {
	    if (args[i].equals("-s0")) {
		s0File = args[i+1];
		CL_Initializer.markAsParsed(i,2);
	    }
	    if (args[i].equals("-lns0")) {
		s0File = args[i+1];
		fileIsLnS0 = true;
		CL_Initializer.markAsParsed(i,2);
	    }
	    if (args[i].equals("-uppertriangular")) {
		upperTriangular = true;
		CL_Initializer.markAsParsed(i);
	    } 
	    if (args[i].equals("-lowertriangular")) {
		upperTriangular = false;
		CL_Initializer.markAsParsed(i);
	    } 
	    
	}

	CL_Initializer.checkParsing(args);

        CL_Initializer.initMaskSource();
    }

    public void initVariables() {
	// get data
	try {
	    dtHeader = ImageHeader.readHeader(CL_Initializer.inputFile);
	    tensorSource = dtHeader.getImageDataSource();
	    
            if (dtHeader.components() != 6) {
              throw new misc.LoggedException("Input is not a tensor image");
            }
	}
	catch (java.io.IOException e) {
	    throw new misc.LoggedException(e);
	}

	if (s0File != null) {
            s0Source = ExternalDataSource.getDataSource(s0File, 1, CL_Initializer.inputDataType);
	}


    }	

    public void execute(OutputManager om) { 

	while (tensorSource.more()) {

	    double[] input = tensorSource.nextVoxel();

	    double[] caminoFormat = new double[8];

            if (CL_Initializer.bgMask != null) {
		caminoFormat[0] = CL_Initializer.bgMask.nextVoxel()[0] == 0.0 ? -1.0 : 0.0;
	    }
	    

	    if (caminoFormat[0] >= 0.0) {
		if (upperTriangular) {
		    caminoFormat[2] = input[0] * scaleSlope + scaleInter;
		    caminoFormat[3] = input[1] * scaleSlope + scaleInter;
		    caminoFormat[4] = input[2] * scaleSlope + scaleInter;
		    caminoFormat[5] = input[3] * scaleSlope + scaleInter;
		    caminoFormat[6] = input[4] * scaleSlope + scaleInter;
		    caminoFormat[7] = input[5] * scaleSlope + scaleInter;
		}	
		else {
		    caminoFormat[2] = input[0] * scaleSlope + scaleInter;
		    caminoFormat[3] = input[1] * scaleSlope + scaleInter;
		    caminoFormat[4] = input[3] * scaleSlope + scaleInter;
		    caminoFormat[5] = input[2] * scaleSlope + scaleInter;
		    caminoFormat[6] = input[4] * scaleSlope + scaleInter;
		    caminoFormat[7] = input[5] * scaleSlope + scaleInter;
		}

	    }

	    if (s0Source != null) {

		double lnS0 = 0.0;

		if (fileIsLnS0) {
		    lnS0 = s0Source.nextVoxel()[0];
		}
		else {
		    double s0 = s0Source.nextVoxel()[0];
		    
		    if (s0 > 0.0) {
			lnS0 = Math.log(s0);
		    }
		}

		if (caminoFormat[0] >= 0.0) {
		    caminoFormat[1] = lnS0;
		}
	    }
		
	    om.output(caminoFormat);
	}

	om.close();
	
    }




}

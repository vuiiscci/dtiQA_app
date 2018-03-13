package apps;

import data.*;
import imaging.*;
import misc.LoggedException;
import tools.*;

import java.io.*;
import java.text.*;
import java.util.logging.Logger;

/**
 * Unary or binary math ops on images
 *
 *
 * @author Philip Cook
 */
public class ImageMath extends Executable {

    // I have avoided writing this for a long time because I don't want to waste time
    // maintaining a competitor to c3d / ANTS ImageMath, but I need something internal to let me diff
    // images for tests


    // Works like this
    // imagemath -inputfile image1.nii.gz -operation image2.nii.gz -outputfile result.nii.gz

    // For summary metrics (eg ssd), output number to stdout

  

   
    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.ImageMath");


    // first image is from inputfile
    private String secondImageFile;
    
    private String operation;

    public ImageMath(String[] args) {
        super(args);
    }


    public void initDefaultVals() {
        secondImageFile = null;
        operation = null;
    }
    

    public void initOptions(String[] args) {

	CL_Initializer.CL_init(args);

	for (int i = 0; i < args.length; i++) {
	    if (args[i].equals("-ssd")) {

                operation = "ssd";

                secondImageFile = args[i+1];

		CL_Initializer.markAsParsed(i, 2);
	    }
	 
	}
        
        CL_Initializer.checkParsing(args);

        
        CL_Initializer.headerTemplateFile = CL_Initializer.inputFile;
	CL_Initializer.initInputSpaceAndHeaderOptions();


        

    }

    public void initVariables() {

    }


    public void execute(OutputManager om) {

        String firstImageFile = CL_Initializer.inputFile;

        DecimalFormat df = new DecimalFormat("0.000E00");

        try {
            
            ImageHeader hdr = ImageHeader.readHeader(firstImageFile);
            
            // Assumes two images per operation and checks image space consistency
            // Voxelwise operations only

            if (! hdr.sameSpace(secondImageFile) ) {
                throw new LoggedException("Images must be in the same physical and voxel space");
            }

            if (operation.equals("ssd")) {
                System.out.println(df.format(ssd(firstImageFile, secondImageFile)));
            }
            else {
                throw new LoggedException("Unrecognized operation " + operation);
            }
            

        }
        catch (IOException e) {
            throw new LoggedException(e);
        }
         

        
        
    }
    

    public double ssd(String image1, String image2) {

        
        DataSource ds1 = null;
        DataSource ds2 = null;
        
        double ssd = 0.0;

        
        try {
            ds1 = ImageHeader.readHeader(image1).getImageDataSource();
        
            ds2 = ImageHeader.readHeader(image2).getImageDataSource();
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }
         
        
        while (ds1.more()) {

            double diff = ds1.nextVoxel()[0] - ds2.nextVoxel()[0];

            ssd += diff * diff;
        }


        return ssd;

    }


}
package apps;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;
import tractography.*;


/**
 *
 * Generates statistics of streamline tracts, either of the tracts themselves or of scalar values
 * along the tracts. The program can either output the raw statistics or combine them into an Analyze
 * image.
 *
 *
 * @author Philip Cook
 * @version $Id$
 *
 * @see tractography.TractStatisticFilter
 * @see tractography.TractStatisticImage
 * 
 *
 */
public class TractStatistics {


    public static void main(String[] args) {
	
	boolean outputImage = false;

	String scalarFile = null;

	CL_Initializer.inputModel = "raw";
	
	CL_Initializer.CL_init(args);

	String outputRoot = "tractstats";

	// If false, only compute the statistic in the voxel in which the fibre was seeded,
	// else statistic is placed in all voxels the fibre intersects
	boolean countIntersect = false;

	// if true, interpolate scalar data set
	boolean interpolated = false;

	// The statistic to compute from the tract: none (output raw scalar values), length (of tract), or 
	// [mean | max | min | median | std] (from image)
	String tractStat = "mean";

	// The statistic to compute in the image, ie what scalar to derive in each voxel :
	// mean, max, min, median, std.
	String imageStat = "mean";

	// Example: we compute mean FA of streamlines, tractStat == mean, countIntersect == true.
	// We get a 4D image where each voxel contains the mean FA of all tracts that intersect the voxel.
	// The imageStat defines what to do with those values in order to produce a 3D scalar output
	// image. If the imageStat is "mean" or "var", these are always weighted, so tracts that intersect
	// more of the voxel have higher weight.


	for (int i = 0; i < args.length; i++) {

	    // image args
            if(args[i].equals("-scalarfile")) {
		scalarFile = args[i + 1];
		CL_Initializer.markAsParsed(i, 2);
	    } 
	    else if (args[i].equals("-countintersect")) {
		countIntersect = true;
	    	CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-interpolate")) {
		interpolated = true;
		CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-tractstat")) {
		tractStat = args[i+1];
	    	CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-imagestat")) {
		imageStat = args[i+1];
	    	CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-outputroot")) {
		outputRoot = args[i+1];
	    	CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-outputimage")) {
		outputImage = true;
	    	CL_Initializer.markAsParsed(i);
	    }
	}	    

	CL_Initializer.checkParsing(args);

	if (outputImage && tractStat.equals("none")) {
	    throw new LoggedException("Tract statistic must be specified when an image is to be generated.");
	}

	if (tractStat.equals("length") || tractStat.equals("endpointsep")) {
          if (scalarFile != null) {
            throw new LoggedException("Scalar image defined, but non-scalar statistic requested. Pass scalar file with -header <file> to define output space");
          }
        } 

	ImageHeader scalarHeader = null;

	try {

	    ScalarImage scalars = null;

	    if (scalarFile != null) {
                
                CL_Initializer.headerTemplateFile = scalarFile;
                
		// derive stats from the tracts and the scalars
		scalarHeader = ImageHeader.readHeader(scalarFile);
	
		double[][][] scalarVol = scalarHeader.readSingleVolumeData();
		
		scalars = new ScalarImage(scalarVol, scalarHeader.getVoxelDims());

	    }

            CL_Initializer.initInputSpaceAndHeaderOptions();


            TractSource source = new TractSource(CL_Initializer.inputFile, CL_Initializer.headerTemplate);

	    if (outputImage) {

		TractStatisticImage statImage = null;

		if (scalars != null) {
		    		    
		    statImage = new TractStatisticImage(scalars);

		}
		else {
		    statImage = new TractStatisticImage(CL_Initializer.dataDims, CL_Initializer.voxelDims);
		}

		statImage.setTractStatistic(tractStat);
		statImage.setImageStatistic(imageStat);
		statImage.setInterpolate(interpolated);
		statImage.setCountIntersect(countIntersect);
		
		while (source.more()) {
		    Tract t = source.nextTract();
		    
		    statImage.processTract(t);
		}
	    
		double[][][] result = statImage.getImageStatistic();
		
                CL_Initializer.headerTemplate.writeScalarImage(result, outputRoot);

	    }
	    else {

		// output raw values for each tract

		OutputManager om = new OutputManager();

		TractStatisticFilter statFilter = null;
		
		if (scalars != null) {
		    
		    statFilter = new TractStatisticFilter(scalars);
		}
		else {
		    statFilter = new TractStatisticFilter(CL_Initializer.dataDims, CL_Initializer.voxelDims);
		}

		statFilter.setTractStatistic(tractStat);
		statFilter.setInterpolate(interpolated);

		while (source.more()) {
		    Tract t = source.nextTract();
		    
		    om.output(statFilter.processTract(t));
		}

		om.close();
	    }
	    

	}
	catch (java.io.IOException e) {
	    throw new LoggedException(e);
	}

	
    }







}

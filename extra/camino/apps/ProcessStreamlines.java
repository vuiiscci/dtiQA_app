package apps;

import imaging.*;

import misc.LoggedException;

import numerics.*;
import tools.*;
import data.*;
import tractography.*;

import java.util.logging.*;
import java.util.Random;
import java.util.zip.*;


import java.io.*;


/**
 * Processes streamline output from StreamlineTractography or from a previous invocation of this
 * program. The program can output streamlines or connection probability images. It can also be
 * used to change the format of streamlines stored on disk.  
 *
 * @author Philip Cook
 * @version $Id$
 * @see apps.StreamlineTractography
 */
public class ProcessStreamlines extends Executable {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.ProcessStreamlines");

    private static int MAX_PDS = 1000;
	
    private boolean outputCP; // if true, output connection probability maps
	
    // if true, implies outputCP
    private boolean outputCBS;

    // if true, implies outputCP 
    // "anatomical connectivity map", global streamline count
    private boolean outputACM;
	
    // if true, implies outputCP
    private boolean normalizeCP;

    // dimensions of seed space
    private int xDataDim;
    private int yDataDim;
    private int zDataDim;

    private double xVoxelDim;
    private double yVoxelDim;
    private double zVoxelDim;

    private String seedFile;
    private String targetFile;
    private String waypointFile;
    private String exclusionFile;
    private String endFile;

    private RealMatrix voxelToPhysical;


    private boolean countFirstEntry; // by default, break on entry to target region


    // discard streamlines on entry to exclusion ROI
    // default behaviour
    private boolean discardOnExclusionEntry; 

    // truncate at second waypoint entry point
    private boolean truncateLoops;

    // discard if fibers loop
    private boolean discardLoops;

    // seed ROIs
    private RegionOfInterest[] allROIs;


    // regionIndex counts from zero, but for I/O purposes we
    // index from 1. If it's -1 then we process all regions in the seed file in sequence
    private int regionIndex;


    private int minTractPoints;	

    private int maxTractPoints;

    private double minTractLength;	

    private double maxTractLength;	

    private String outputRoot;

    private double resampleStepSize;

    // gzip output if true
    private boolean gzip;


    // number of iterations for tracking
    private int iterations;
	
    // do not print to stderr if true
    private boolean silent;


    public ProcessStreamlines(String[] args){
        super(args);
    }


    /**
     * initialise the default values of class-level
     * fields, as they would be at declaration.
     * 
     * It is unclear if this is really necessary, and if it is, why.
     *
     */
    public void initDefaultVals() {
        
	
        outputCP = false; // if true, output connection probability maps
	
	// if true, implies outputCP
        outputCBS = false;

 	// if true, implies outputCP
	// "anatomical connectivity map", global streamline count
        outputACM = false;
	
	// if true, implies outputCP
        normalizeCP = false;

	// dimensions of seed space
        xDataDim = 0;
        yDataDim = 0;
        zDataDim = 0;

        xVoxelDim = 0.0;
        yVoxelDim = 0.0;
        zVoxelDim = 0.0;

        seedFile = null;
        targetFile = null;
        waypointFile = null;
        exclusionFile = null;
        endFile = null;
        
        voxelToPhysical = null;
        
        countFirstEntry = true; // by default, break on entry to target region


	// discard streamlines on entry to exclusion ROI
	// default behaviour
        discardOnExclusionEntry = true; 

	// truncate at second waypoint entry point
        truncateLoops = false;

	// discard if fibers loop
        discardLoops = false;

        allROIs = null;


        // regionIndex counts from zero, but for I/O purposes we
        // index from 1. If it's -1 then we process all regions in the seed file in sequence
        regionIndex = -1;


        minTractPoints = 0;	

        maxTractPoints = 0;

        minTractLength = 0.0;	

        maxTractLength = 0;	

        outputRoot = null;

        resampleStepSize = 0.0;

	// number of iterations for tracking
        iterations = 1;
	
	// do not print to stderr if true
        silent = false;

    }
    
    
    /**
     *
     * Command line parsing
     *
     */
    public void initOptions(String[] args) {

        // your CL_Initializer defaults here

	CL_Initializer.inputModel = "raw";

	// default number of PDs is 1
	CL_Initializer.numPDsIO = 1;

        // Affects image output
        OutputManager.outputDataType = "float";

        // Parse the command line arguments
        CL_Initializer.CL_init(args);

        // gzip output if true
        gzip = OutputManager.gzipOut;

        // your args parsed here	

        for (int i = 0; i < args.length; i++) {
	    // image args
            if (args[i].equals("-resamplestepsize")) { 
		resampleStepSize = Double.parseDouble(args[i + 1]);
		CL_Initializer.markAsParsed(i, 2);
	    }
	    else if (args[i].equals("-noresample")) { 
		resampleStepSize = -1.0;
		CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-seedfile")) {
		seedFile = args[i + 1];
		CL_Initializer.markAsParsed(i, 2);
	    } 
 	    else if (args[i].equals("-regionindex")) {
                regionIndex = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
	    } 
            else if (args[i].equals("-waypointfile")) {
                waypointFile = args[i + 1];
                CL_Initializer.markAsParsed(i,2);
	    }
            else if (args[i].equals("-exclusionfile")) {
                exclusionFile = args[i + 1];
                CL_Initializer.markAsParsed(i,2);
	    }
            else if (args[i].equals("-targetfile")) {
                targetFile = args[i + 1];
                CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-endpointfile")) {
		endFile = args[i + 1];
		CL_Initializer.markAsParsed(i, 2);
	    } 
            else if (args[i].equals("-mintractpoints")) {
		minTractPoints = Integer.parseInt(args[i + 1]);
		CL_Initializer.markAsParsed(i, 2);
	    }
            else if (args[i].equals("-mintractlength")) {
		minTractLength = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i, 2);
	    }
            else if (args[i].equals("-maxtractpoints")) {
		maxTractPoints = Integer.parseInt(args[i + 1]);
		CL_Initializer.markAsParsed(i, 2);
	    }
            else if (args[i].equals("-maxtractlength")) {
		maxTractLength = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i, 2);
	    }
	    else if (args[i].equals("-truncateinexclusion")) {
		discardOnExclusionEntry = false;
		CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-truncateloops")) {
		truncateLoops = true;
                CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-discardloops")) {
		discardLoops = true;
                CL_Initializer.markAsParsed(i);
	    }
            else if (args[i].equals("-allowmultitargets")) {
                countFirstEntry = false;
                CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-outputtracts")) {
                outputCP = false;
                CL_Initializer.markAsParsed(i, 1);
            }
	    else if  (args[i].equals("-outputcp")) {
		outputCP = true;
		normalizeCP = true;
		CL_Initializer.markAsParsed(i);
	    }
	    else if  (args[i].equals("-outputcbs")) {
		outputCBS = true;
		outputCP = true;
		CL_Initializer.markAsParsed(i);
	    }
	    else if  (args[i].equals("-outputsc")) {
		outputCP = true;
		normalizeCP = false;
		CL_Initializer.markAsParsed(i);
	    }
	    else if  (args[i].equals("-outputacm")) {
		outputCP = true;
		outputACM = true;
		CL_Initializer.markAsParsed(i);
	    }
	    else if (args[i].equals("-iterations")) {
		iterations = Integer.parseInt(args[i + 1]);
		CL_Initializer.markAsParsed(i, 2);
	    }
	    else if (args[i].equals("-outputroot")) {
		outputRoot = args[i + 1];
		CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-silent")) {
		silent = true;
	    	CL_Initializer.markAsParsed(i);
	    }
	}

        CL_Initializer.checkParsing(args);

    }
    

    /**
     * Initialize class variables 
     */
    public void initVariables() {
        
    }
    

    /**
     * abstract method execution for commands
     */
    public void execute(OutputManager om) {

        // IMPORTANT - input tracts will be in PHYSICAL SPACE
        
        // the filter brings them into Camino space for processing

        // if outputting tracts, put them back into physical space


	if (truncateLoops && discardLoops) {
	    logger.warning("Discard and truncate loops both selected. Program will truncate loops " + 
			   "(run again with discard option to remove them if desired)");
	    
	    discardLoops = false;
	}

       
        if (CL_Initializer.headerTemplateFile != null || CL_Initializer.dataDims[0] > 0) { 
            CL_Initializer.initInputSpaceAndHeaderOptions();
        }
        
   
	if (seedFile != null) {
            
            if (CL_Initializer.headerTemplate == null) {
                logger.info("Determining input physical space from seed file");
                CL_Initializer.headerTemplateFile = seedFile;
                CL_Initializer.initInputSpaceAndHeaderOptions();
            }
            else {
                if ( !CL_Initializer.headerTemplate.sameSpace(seedFile) ) {
                    throw new LoggedException(" Inconsistent image space : " + seedFile);
                }
            }
 
            
	    allROIs = new VoxelROI(seedFile, CL_Initializer.headerTemplate).getAllRegions();

            // seed file but no output root
	    if (outputRoot == null) {
                if (OutputManager.outputFile != null) {
                    logger.info("Have seed file but no output root, output will be to output file " + 
                                OutputManager.outputFile);
                }
                else {
                    logger.info("Have seed file but no output root, output will be to stdout");
                }
            }
            
            
        }
	else {
            if (outputCP && !outputACM) {
                throw new LoggedException("Seed points must be specified for connection " + 
                                          "probability images");
            }
	    else if (!outputCP && outputRoot != null) {
		throw new LoggedException("Without seed ROI information, output is to a single stream " + 
					  "(default stdout). Specify a file with -outputfile");
	    }
            
	}

	if (outputRoot == null) {
	    om = new OutputManager();
	}

        
        int[][][] waypointVol = null;
        
        if (waypointFile != null) {
	   
             if (CL_Initializer.headerTemplate == null) {
                 logger.info("Determining input physical space from waypoint file");
                 CL_Initializer.headerTemplateFile = waypointFile;
                 CL_Initializer.initInputSpaceAndHeaderOptions();
             }
             else {
                 if ( !CL_Initializer.headerTemplate.sameSpace(waypointFile) ) {
                     throw new LoggedException(" Inconsistent image space : " + waypointFile);
                 }
             }
            
            waypointVol = readIntVolume(waypointFile);
        }

        int[][][] exclusionVol = null;

        if (exclusionFile != null) {

            if (CL_Initializer.headerTemplate == null) {
                logger.info("Determining input physical space from exclusion file");
                CL_Initializer.headerTemplateFile = exclusionFile;
                CL_Initializer.initInputSpaceAndHeaderOptions();
            }
            else {
                if ( !CL_Initializer.headerTemplate.sameSpace(exclusionFile) ) {
                    throw new LoggedException(" Inconsistent image space : " + exclusionFile);
                }
            }
            
            exclusionVol = readIntVolume(exclusionFile);
        }

        int[][][] targetVol = null;
        
        if (targetFile != null) {

            if (CL_Initializer.headerTemplate == null) {
                logger.info("Determining input physical space from target file");
                CL_Initializer.headerTemplateFile = targetFile;
                CL_Initializer.initInputSpaceAndHeaderOptions();
            }
            else {
                if ( !CL_Initializer.headerTemplate.sameSpace(targetFile) ) {
                    throw new LoggedException(" Inconsistent image space : " + targetFile);
                }
            }

            targetVol = readIntVolume(targetFile);

        }



        int[][][] endVol = null;

        if (endFile != null) {
            
            if (CL_Initializer.headerTemplate == null) {
                logger.info("Determining input physical space from endpoint file");
                CL_Initializer.headerTemplateFile = endFile;
                CL_Initializer.initInputSpaceAndHeaderOptions();
            }
            else {
                if ( !CL_Initializer.headerTemplate.sameSpace(endFile) ) {
                    throw new LoggedException(" Inconsistent image space : " + endFile);
                }
            }
            
            endVol = readIntVolume(endFile);
        }
        

        // If no header template at this point, we have no definition of physical space
        if (CL_Initializer.headerTemplate == null) {

            throw new LoggedException("Definition of physical space required, use -header <image> with any image in " + 
                                      "the correct space");
        }

        voxelToPhysical = CL_Initializer.headerTemplate.getVoxelToPhysicalTransform();

	xDataDim = CL_Initializer.dataDims[0];
	yDataDim = CL_Initializer.dataDims[1];
	zDataDim = CL_Initializer.dataDims[2];

	xVoxelDim = CL_Initializer.voxelDims[0];
	yVoxelDim = CL_Initializer.voxelDims[1];
	zVoxelDim = CL_Initializer.voxelDims[2];


	TractSource tractSource = new TractSource(CL_Initializer.inputFile, CL_Initializer.headerTemplate);
	

	if (outputCBS && targetVol == null) {
	    throw new LoggedException("Cannot do connectivity segmentation without targets");
	}
	
	if (xDataDim == 0.0 || yDataDim == 0.0 || zDataDim == 0.0 || xVoxelDim == 0.0 || yVoxelDim == 0.0 || zVoxelDim == 0.0) {
        
            throw new LoggedException("This program requires a definition of the image space");
	}
        
	StreamlineROI_Filter filter = 
	    new StreamlineROI_Filter(CL_Initializer.dataDims, CL_Initializer.voxelDims);

	if (waypointFile != null) {
	    filter.setWaypoints(waypointVol);
	    filter.setTruncateLoops(truncateLoops);
	    filter.setDiscardLoops(discardLoops);
        }

	if (exclusionFile != null) {
	    filter.setExclusionROIs(exclusionVol);
	    filter.setDiscardOnExclusionEntry(discardOnExclusionEntry);
	    
	}        

	if (endFile != null) {
	    filter.setEndZones(endVol);
	}        


	// if resampleStepSize == 0.0 use default resampling scheme
	// if < 0 then do not resample
	// else use specified step size
	if (resampleStepSize > 0.0) {
	    filter.setResampleStepSize(resampleStepSize);
	}
	else { 
            filter.setResampleTracts(false);
	}
	
	if (minTractPoints > 0) {
	    filter.setMinTractPoints(minTractPoints);
	}
	if (minTractLength > 0.0) {
	    filter.setMinTractLength(minTractLength);
	}
	if (maxTractPoints > 0) {
	    filter.setMaxTractPoints(maxTractPoints);

	    // makes no sense to resample tracts if we are restricting the maximum number of points
	    filter.setResampleTracts(false);

            if (resampleStepSize > 0.0) {
                logger.info("Max tract points set, resampling disabled");
            }
            
	}
	if (maxTractLength > 0.0) {
	    filter.setMaxTractLength(maxTractLength);
	}

	
	try {

	    int numRegions = 1;

	    if (allROIs != null) {
		numRegions = allROIs.length;
	    }
	    
	    
	    // output numbering for CBS images is
	    // outputRoot%d_%d, outputRegionID, pd
	    // output numbering for targetCP images is
	    // outputRoot%d_%d_%d, outputRegionID, sp+1, pd,

	    // output numbering for raw CP images is
	    // outputRoot%d_%d_%d, regionIndex, sp+1, pd
	    
	    if (outputACM) {
		
		int tractCounter = 0;

		if (!silent) {
	    	    System.err.println();
		}
                
                // untested feature - allow acm with targets
                // something of a hack for now.
                if (targetVol != null) {

                    
                    TargetCP_Image tacm = 
                        new TargetCP_Image(targetVol, xVoxelDim, yVoxelDim, zVoxelDim);


                    tacm.setCountFirstEntry(countFirstEntry);


                    while (tractSource.more()) {
                        
                        if (!silent) {		    
                            System.err.print("\rProcessing streamline " + (tractCounter+1));
                        }
                        
                        TractCollection nextTC = new TractCollection(2, 100.0);
                        
                        nextTC.addTract(tractSource.nextTract());
                        
                        nextTC = filter.processTracts(nextTC);
                        
                        tacm.processTracts(nextTC);
			
                        tractCounter++;
                    }                                   
                  
                    if (normalizeCP) {
                        double[][][] targetCPImage = tacm.getConnectionProbabilities();
                        CL_Initializer.headerTemplate.writeScalarImage(targetCPImage, outputRoot + "acm_target_cp");
                    }
                    else {
                        double[][][] targetCPImage = tacm.getStreamlineCounts();
                        CL_Initializer.headerTemplate.writeScalarImage(targetCPImage, outputRoot + "acm_target_sc");
                    }
	
                    
                }
                else {
                    ConnectionProbabilityImage acm = 
                        new ConnectionProbabilityImage(xDataDim, yDataDim, zDataDim, 
                                                       xVoxelDim, yVoxelDim, zVoxelDim);
                    
                    while (tractSource.more()) {
		    
                        if (!silent) {		    
                            System.err.print("\rProcessing streamline " + (tractCounter+1));
                        }
                        
                        TractCollection nextTC = new TractCollection(2, 100.0);
                        
                        nextTC.addTract(tractSource.nextTract());
                        
                        nextTC = filter.processTracts(nextTC);
                        
                        acm.processTracts(nextTC);
			
                        tractCounter++;
                    }
                    
                    // now output
                    if (normalizeCP) {
                        double[][][] cp = acm.getConnectionProbabilities();
                        
                        CL_Initializer.headerTemplate.writeScalarImage(cp, outputRoot + "acm_cp");
                    }
                    else {
                        double[][][] sc = acm.getStreamlineCounts();
                        
                       CL_Initializer.headerTemplate.writeScalarImage(sc, outputRoot + "acm_sc");
                    }

                }
                    
                if (!silent) {
                    System.err.println();
                }
                
		logger.info("Processed " + tractCounter + " streamlines for anatomical connectivity map");
		
	    } // end if ACM
	    else if (!outputCP && outputRoot == null) {

		DataOutputStream dout = om.getOutputStream();
		
		// if no outputRoot is given, output to stdout or single file
		
		// just read tracts, filter and and output
		// no need to worry about tracts per voxel, etc
		int tractCounter = 0;

		if (!silent) {
	    	    System.err.println();
		}
    
		while (tractSource.more()) {
		
		    if (!silent) {		    
		        System.err.print("\rProcessing streamline " + (tractCounter+1));
		    }
	
		    TractCollection nextTC = new TractCollection(2, 100.0);
		    
		    nextTC.addTract(tractSource.nextTract());
		    
		    nextTC = filter.processTracts(nextTC);
		    
		    for (int i = 0; i < nextTC.numberOfTracts(); i++) {
			
			Tract t = nextTC.getTract(i);

                        t.transformToPhysicalSpace(voxelToPhysical, xVoxelDim, yVoxelDim, zVoxelDim);
			    
                        t.writeRaw(dout);
                        
		    }
			
			
		    tractCounter++;
		}
		
		dout.close();
		
	    }
	    else { // output images, but not ACM

		int seedCounter = 0;

		// read first Tract outside of loop, that way we can tell if a chunk of streamlines
		// come from the same seed and different PD or a different seed
		
		// don't filter here because we need this tract regardless of its waypoint
		// or exclusion status
		
		TractCollection nextTC = new TractCollection(2, 100.0);
		    
		nextTC.addTract(tractSource.nextTract());
		
                		
		regions : for (int region = 0; region < numRegions; region++) {
                    
		    RegionOfInterest roi = allROIs[region];

                    int seedPointsThisROI = roi.getSeedPoints().length;

		    // number ROIs from 1 upwards to maintain correspondence to 
		    // values in seed file.
		    int outputRegionID = roi.getRegionLabel();
		
                    
                    if (regionIndex > -1) {
                        if (outputRegionID != regionIndex) {
                            continue regions;
                        }
                    }
		   
                    if (!silent) { 
                        System.err.println("\nProcessing ROI " + (region+1) + " of " + numRegions);
                    }
		
		    if (outputCBS) {

			ConnectivitySegmentedImage[] cbsImages = new ConnectivitySegmentedImage[MAX_PDS];
			
			for (int sp = 0; sp < seedPointsThisROI; sp++) {
			 
			    if (!silent) {	   
			        System.err.print("\rProcessing seed " + (sp + 1) + " of " + seedPointsThisROI); 
			    }

			    // index of current PD
			    int pd = 0;
			    
			    Point3D seedPoint = nextTC.getTract(0).getSeedPoint();
                            
                            Point3D ntSeedPoint = nextTC.getTract(0).getSeedPoint();
			    
			    // a bunch of (iterations) tracts sharing the same seed are
			    // interpreted as belonging to a different PD in the same seed
			    // this means we can only have one CP map per voxel.
			    
			    while (nextTC != null && seedPoint.equals(ntSeedPoint)) {

				// initialize output for this PD 
				if (cbsImages[pd] == null) {
				    cbsImages[pd] = new ConnectivitySegmentedImage(roi, targetVol, 
                                                                                   xVoxelDim, yVoxelDim, zVoxelDim);
				    
				    cbsImages[pd].setCountFirstEntry(countFirstEntry);
				}
				
				cbsImages[pd].processTracts(sp, filter.processTracts(nextTC));
				
				for (int i = 1; i < iterations; i++) {
				    TractCollection tc = filter.processTract(tractSource.nextTract());
				    cbsImages[pd].processTracts(sp, tc);
				}
				
				if (tractSource.more()) {
				    pd++;
				    
				    nextTC = new TractCollection(2, 100.0);
				
				    nextTC.addTract(tractSource.nextTract());
				
				    ntSeedPoint = nextTC.getTract(0).getSeedPoint();
				}
				else {
				    nextTC = null;
                                    seedPoint = null;
				    ntSeedPoint = null;
                                    
				    if (sp < seedPointsThisROI - 1) {
					throw new 
					    LoggedException("No more tracts in input after processing seed " + 
							    (sp+1) + " of " + seedPointsThisROI);
				    }
				}
				
			    }
			    
			}
			
			
			// now output CBS
			for (int p = 0; p < MAX_PDS; p++) {
			    if (cbsImages[p] != null) {
				
				double[][][] labelledSeeds = cbsImages[p].getSegmentedSeeds();
			    
				double[][][] targetCPImage = null;


                                CL_Initializer.headerTemplate.writeScalarImage(labelledSeeds, outputRoot + "labels_" + 
                                                                               outputRegionID + "_" + (p+1));


				if (normalizeCP) {
				    targetCPImage = cbsImages[p].getMaxTargetCP();
                                    
                                    CL_Initializer.headerTemplate.writeScalarImage(targetCPImage, outputRoot + "labelcp_" + 
                                                                                   outputRegionID + "_" + (p+1));
				}
				else {
				    targetCPImage = cbsImages[p].getMaxTargetSC();

                                    CL_Initializer.headerTemplate.writeScalarImage(targetCPImage, outputRoot + "labelsc_" + 
                                                                                   outputRegionID + "_" + (p+1));
				}
				
			    }

			}
			
		    }
		    else if (outputCP && targetVol != null) {

			TargetCP_Image cpImage = 
			    new TargetCP_Image(targetVol, xVoxelDim, yVoxelDim, zVoxelDim);
			
			cpImage.setCountFirstEntry(countFirstEntry);
			

			for (int sp = 0; sp < seedPointsThisROI; sp++) {

			    if (!silent) {	
			    	System.err.print("\rProcessing seed " + (sp + 1) + " of " + seedPointsThisROI); 
			    }		
			   
			    // index of current PD
			    int pd = 0;
		
			    Point3D seedPoint = nextTC.getTract(0).getSeedPoint();
                            
                            Point3D ntSeedPoint = nextTC.getTract(0).getSeedPoint();
			
			    while (nextTC != null && seedPoint.equals(ntSeedPoint)) {
				
				cpImage.reset();
			
				cpImage.processTracts(filter.processTracts(nextTC));

				for (int i = 1; i < iterations; i++) {
				    TractCollection tc = filter.processTract(tractSource.nextTract());
				    cpImage.processTracts(tc);
				}


		
				if (normalizeCP) {
				    double[][][] targetCPImage = cpImage.getConnectionProbabilities();

                                    CL_Initializer.headerTemplate.writeScalarImage(targetCPImage, outputRoot + outputRegionID + 
                                                                                   "_" + (sp+1) +  "_" + (pd+1));

				}
				else {
				    double[][][] targetCPImage = cpImage.getStreamlineCounts();

                                    CL_Initializer.headerTemplate.writeScalarImage(targetCPImage, outputRoot + outputRegionID + 
                                                                                   "_" + (sp+1) + "_" + (pd+1));
				}
                                
				pd++;
			
				if (tractSource.more()) {
				    nextTC = new TractCollection(2, 100.0);
			
				    nextTC.addTract(tractSource.nextTract());
			    			    
				    ntSeedPoint = nextTC.getTract(0).getSeedPoint();
				}
				else {
				    nextTC = null;
				    seedPoint = null;
				    ntSeedPoint = null;

				    if (sp < seedPointsThisROI - 1) {
					throw new 
					    LoggedException("No more tracts in input after processing seed " + 
							    (sp+1) + " of " + seedPointsThisROI);
				    }

				}
			
			    }

		    
			}
		    }
		    else if (outputCP) { // no targets
			
			for (int sp = 0; sp < seedPointsThisROI; sp++) {
			    
				
			    if (!silent) {
			        System.err.print("\rProcessing seed " + (sp + 1) + " of " + seedPointsThisROI); 
			    }

			    // index of current PD
			    int pd = 0;
		
                            Point3D seedPoint = nextTC.getTract(0).getSeedPoint();
                            
                            Point3D ntSeedPoint = nextTC.getTract(0).getSeedPoint();

			    while (nextTC != null && seedPoint.equals(ntSeedPoint)) {

				ConnectionProbabilityImage cpImage = 
				    new ConnectionProbabilityImage(xDataDim, yDataDim, zDataDim, 
								   xVoxelDim, yVoxelDim, zVoxelDim);
			
				cpImage.processTracts(filter.processTracts(nextTC));

				for (int i = 1; i < iterations; i++) {
				    TractCollection tc = filter.processTract(tractSource.nextTract());
				    cpImage.processTracts(tc);
				}

				// now output
				if (normalizeCP) {
				    double[][][] cp = cpImage.getConnectionProbabilities();

				    CL_Initializer.headerTemplate.writeScalarImage(cp, outputRoot + outputRegionID
							     + "_" + (sp+1) + "_" + (pd+1));
				}
				else {
				    double[][][] sc = cpImage.getStreamlineCounts();
			    
				     CL_Initializer.headerTemplate.writeScalarImage(sc, outputRoot + outputRegionID
							     + "_" + (sp+1) + "_" + (pd+1));
				}
			
				
				if (tractSource.more()) {
				    nextTC = new TractCollection(2, 100.0);
			    
                                    nextTC.addTract(tractSource.nextTract());

                                    ntSeedPoint = nextTC.getTract(0).getSeedPoint();
                                    
				}
				else {
				    nextTC = null;
                                    seedPoint = null;
                                    ntSeedPoint = null;
                                    
				    if (sp < seedPointsThisROI - 1) {
					throw new 
					    LoggedException("No more tracts in input after processing seed " + 
							    (sp+1) + " of " + seedPointsThisROI);
				    }

				}
				
				pd++;

			    }
			} // for sp
                        
		    } 
		    else { // output tracts, have an outputRoot
                        
			// if outputroot is specified it implies that the output is to be structured
			// by ROI, as in the old track
                        
			FileOutputStream fout = null;
			DataOutputStream dout = null;

                        
                        if (gzip) {
                            fout = new FileOutputStream(outputRoot + outputRegionID + ".Bfloat.gz");
                            dout = new DataOutputStream
				    (new GZIPOutputStream(fout, 1024*1024*16));
                        }
                        else {
                            fout = new FileOutputStream(outputRoot + outputRegionID + ".Bfloat");
                            dout = new DataOutputStream
                                (new BufferedOutputStream(fout, 1024*1024*16));
                        }
			
		 	if (!silent) {		
			    System.err.println();
			}

			for (int sp = 0; sp < seedPointsThisROI; sp++) {
			    
		    
			    if (!silent) {				    
			    	System.err.print("\rProcessing seed " + (sp + 1) + " of " + seedPointsThisROI); 
		    	    }
                            

                            Point3D seedPoint = nextTC.getTract(0).getSeedPoint();
                            
                            Point3D ntSeedPoint = nextTC.getTract(0).getSeedPoint();
                            
			    while (nextTC != null && seedPoint.equals(ntSeedPoint)) {
                                
				TractCollection filtered = filter.processTracts(nextTC);

				for (int tr = 0; tr < filtered.numberOfTracts(); tr++) {
				    
                                    Tract t = filtered.getTract(tr);
                                    t.transformToPhysicalSpace(voxelToPhysical, xVoxelDim, yVoxelDim, zVoxelDim);
                                    t.writeRaw(dout);
				}
				
				
				if (tractSource.more()) {
                                    nextTC = new TractCollection(2, 100.0);
                                    
                                    nextTC.addTract(tractSource.nextTract());
                                    
                                    ntSeedPoint = nextTC.getTract(0).getSeedPoint();
                                    
				}
				else {
                                    nextTC = null;
                                    seedPoint = null;
                                    ntSeedPoint = null;

				    if (sp < seedPointsThisROI - 1) {
					throw new 
					    LoggedException("No more tracts in input after processing seed " + 
							    (sp+1) + " of " + seedPointsThisROI);
				    }

				}
				
			    }
			}
		    	
			dout.close();
			
		    } // else (outputting tracts with outputroot specified)
		
		} // end for regions
		
	    } // end else (ie if outputCP || tracts by region)

	    if (tractSource.more()) {
		logger.warning("Tract processing finished but there are more Tracts in input. Check -seedfile and -iterations options.");
	    }
	  
	    if (!silent) {	
	    	System.err.println("\n");
	    }
	}
	catch(IOException e) {
	    throw new LoggedException(e);
	}
    
    }


    protected static int[][][] readIntVolume(String file) {
        try {
         
            ImageHeader ih = ImageHeader.readHeader(file);
            
            int xDataDim = ih.xDataDim();
            int yDataDim = ih.yDataDim();
            int zDataDim = ih.zDataDim();
            
            int[][][] vol = new int[xDataDim][yDataDim][zDataDim];
            
            DataSource vin = ih.getImageDataSource();
            
            for (int k = 0; k < zDataDim; k++) {
                for (int j = 0; j < yDataDim; j++) {
                    for (int i = 0; i < xDataDim; i++) {
                        
                        double value = vin.nextVoxel()[0];
                        if (value > 0.0) {
                            if (value > Integer.MAX_VALUE) {
                                throw new LoggedException("Maximum value allowed for seeds / waypoints / targets is " + 
                                                          Integer.MAX_VALUE);
                            }
                            vol[i][j][k] = (int)Math.round(value);
                        }
                    }
                }
            }
            return vol;
            
        } catch (IOException e) {
            throw new LoggedException(e);
        }
        
    }



}

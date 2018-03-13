package apps;

import data.*;

import imaging.*;
import inverters.ModelIndex;
import misc.LoggedException;

import numerics.*;
import tools.*;
import tractography.*;

import java.util.Random;
import java.util.zip.*;
import java.util.logging.*;

import java.io.*;

/**
 *
 * Does streamline tractography.
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class StreamlineTractography extends Executable {
    
    public StreamlineTractography(String[] args){
        super(args);
    }
    
    // Not the input image type, but the image type that will be used to track
    public enum ImageType {
        
	DT("dt"),
	MULTITENSOR("multitensor"),
	SF_PEAK("sfpeak"),
	WILDBS_DT("wildbs_dt"),
	BAYESDIRAC("bayesdirac"),
	BAYESDIRAC_DT("bayesdirac_dt"),
	REPBS_DT("repbs_dt"),
	REPBS_MULTITENSOR("repbs_multitensor"),
	PICO("pico"),
	DWI_DT("dwi_dt"),
	DWI_MULTITENSOR("dwi_multitensor"),
	VECTOR("vector"),
	BEDPOSTX("bedpostx"),
	BEDPOSTX_DYAD("bedpostx_dyad");

        ImageType(String s) {
	    modelString = s;
        }

	public static ImageType getImageType(String s) {
            
	    for (ImageType model : ImageType.values()) {
		if (model.modelString.equals(s)) {
		    return model;
		}
		
	    }  
            
	    throw new IllegalArgumentException("Unrecognized input model " + s);
	}
        
	public final String modelString;
        
    }
    

    public enum TrackingAlgorithm {
      
	FACT("fact"),
	EULER("euler"),
	RK4("rk4");
        
        TrackingAlgorithm(String s) {
	    algString = s;
        }

	public static TrackingAlgorithm getTrackingAlgorithm(String s) {

	    for (TrackingAlgorithm alg : TrackingAlgorithm.values()) {
		if (alg.algString.equals(s)) {
		    return alg;
		}
		
	    }  
	    
	    throw new IllegalArgumentException("Unrecognized tracking algorithm " + s);
	
	}

	public final String algString;

    }


    public enum DataInterpolation {
	NEAREST_NEIGHBOUR("nn"),
	NEIGHBOUR_CHOICE("prob_nn"),
	DWI_TRILINEAR("dwi_linear"),
	VECTOR_TRILINEAR("linear"),
	TEND_NN("tend"),
	TEND_NC("tend_prob_nn");
        
        DataInterpolation(String s) {
	    interpModel = s;
        }
        
	public static DataInterpolation getDataInterpolation(String s) {
            
	    for (DataInterpolation inter : DataInterpolation.values()) {
		if (inter.interpModel.equals(s)) {
		    return inter;
		}
		
	    }  
	    
	    throw new IllegalArgumentException("Unrecognized interpolation algorithm " + s);
	
	}
        
        
	private String interpModel;
        
    }
    
	
    private double stepSize;

    private int regionIndex;
                     
    private int xDataDim;
    private int yDataDim;
    private int zDataDim;
       
    private double xVoxelDim;
    private double yVoxelDim;
    private double zVoxelDim;

    // image containing seed ROI
    private String seedFile;
	
    // text file containing physical space seed points
    private String seedList;
  
    // Bayesian options
    private double curvePriorK;
    private double curvePriorG;

    // random seed, defaults to system time
    private long seed;
    
    private double anisThresh;

    // curvature threshold
    private double ipThresh;

    // distance to track before checking curvature threshold
    private double checkCurveLength;


    // if no one-DT data is provided, then an ANIS map is needed
    // otherwise no ANIS threshold can be applied.
    private String anisMapFile;
	
    // contains all ROIs in one volume or seed list
    private RegionOfInterest[] allROIs;

    // path to image containing f for TEND
    private String tendF_File;
        
    // constant f used if file is not given
    private double tendF;
        
    // use a constant g, since low anisotropy DTs automatically deflect the direction less
    private double tendG;

    // if true, don't print anything to stderr
    private boolean silent;

    private int iterations;
  
    private PICoPDF picoPDF;

    private Random ran;
    private String outputRoot;

    // type of image and tracking algorithm
    private ImageType imageType;

    private TrackingAlgorithm trackingAlgorithm;

    private DataInterpolation dataInterpolation;


    // gzip output if true
    private boolean gzip;
	
    // path to a PICo tractography image
    private String externalPriorImageFile;
    private String externalPriorDataType;
  	
    // voxel space of diffusion image to physical space
    private RealMatrix voxelToPhysicalTrans;


    // List of files containing vectors for deterministic tracking. At least one image of
    // vectors but optionally more with extra vectors in some or all voxels
    private String[] vectorFiles;

    // Path to bedpostx output
    private String bedpostxDir;

    // Minimum compartment size for bedpost sticks to be used for tracking
    private double bedpostxMinF;
	

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.StreamlineTractography");
    
    
    public void initDefaultVals(){
		
        stepSize = 1.0;
		        
        xDataDim = 0;
        yDataDim = 0;
        zDataDim = 0;
        
        xVoxelDim = 0.0;
        yVoxelDim = 0.0;
        zVoxelDim = 0.0;

        regionIndex = -1;
        
        seedFile = null;
        
        // text file containing mm seed points
        seedList = null;
  
        curvePriorK = 0.0;
        curvePriorG = 0.0;

        seed = System.currentTimeMillis();
    
        anisThresh = 0.0;
        ipThresh = 0.0;

        checkCurveLength = 5.0;
		
        anisMapFile = null;

        allROIs = null;

        tendF_File = null;
        tendF = 0.0;
        tendG = 1.0;
		
        silent = false;

        iterations = 1;
  
        picoPDF = PICoPDF.BINGHAM;
		
        ran = null;
        outputRoot = null;
		
        imageType = null;
        trackingAlgorithm = TrackingAlgorithm.FACT;
        dataInterpolation = DataInterpolation.NEAREST_NEIGHBOUR;

        gzip = OutputManager.gzipOut;
        
        externalPriorImageFile = null;
        externalPriorDataType = "double";
  
        voxelToPhysicalTrans = RealMatrix.identity(4);

	vectorFiles = null;
	
	bedpostxDir = null;

	bedpostxMinF = 0.01;

    }
    
    public void initOptions(String[] args){
        
        if (args.length == 0) {
            System.exit(0);
        }
        
        CL_Initializer.inputDataType = "double";
        
        // for Bayesian tracking
        CL_Initializer.pointSetInd = 1;
        
        CL_Initializer.CL_init(args);
     	

	if (CL_Initializer.inputModel == null) {
            throw new LoggedException("An input model is required.");
        }

	imageType = ImageType.getImageType(CL_Initializer.inputModel);

	// Set default input type to float for anything that uses raw DWI data
	
	switch (imageType) {

	case DWI_DT : case DWI_MULTITENSOR : case WILDBS_DT : case REPBS_DT : case REPBS_MULTITENSOR : case BAYESDIRAC : case BAYESDIRAC_DT :
	    CL_Initializer.inputDataType = "float";
	}

        for (int i = 0; i < args.length; i++) {
            // input data
            if (args[i].equals("-inputdatatype")) {
                // needed because we change default for bootstrap data
                CL_Initializer.inputDataType = args[i+1];
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-regionindex")) {
                regionIndex = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-seedfile")) {
                seedFile = args[i + 1];
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-seedlist")) {
                seedList = args[i + 1];
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-randomseed") || args[i].equals("-seed")) {
                // don't want CL default of 0, so parse here
                // also named differently to avoid confusion with tracking seed
                // but allow CL default
                seed = Long.parseLong(args[i + 1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-anisthresh")) {
                anisThresh = Double.parseDouble(args[i + 1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-anisfile")) {
                anisMapFile = args[i + 1];
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-ipthresh")) {
                ipThresh = Double.parseDouble(args[i + 1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-curvethresh")) {
                ipThresh = Math.cos(Math.PI * Double.parseDouble(args[i + 1]) / 180.0);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-curveinterval")) {
                checkCurveLength = Double.parseDouble(args[i + 1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-stepsize") || args[i].equals("-steplength")) {
                stepSize = Double.parseDouble(args[i + 1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-iterations")) {
                iterations = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-pdf")) {
                CL_Initializer.markAsParsed(i);
                if (i == args.length - 1) {
                    throw new LoggedException("Missing PDF type");
                }
                picoPDF = PICoPDF.getPDF(args[i+1]);
                
                CL_Initializer.markAsParsed(i+1);
            }
            else if (args[i].equals("-curvepriork")) {
                // Bayesian tracking option
                curvePriorK = Double.parseDouble(args[i+1]);
                curvePriorG = 0.0;
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-curvepriorg")) {
                // Bayesian tracking option
                curvePriorG = Double.parseDouble(args[i+1]);
                curvePriorK = 0.0;
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-extpriorfile")) {
                // Bayesian tracking option
                externalPriorImageFile = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-extpriordatatype")) {
                // Bayesian tracking option
                externalPriorDataType = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-outputroot")) {
                outputRoot = args[i + 1];
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-interpolator")) {
                
                dataInterpolation = DataInterpolation.getDataInterpolation(args[i+1]);
                
                CL_Initializer.markAsParsed(i, 2);
                
            }
            else if (args[i].equals("-tracker")) {
                
                trackingAlgorithm = TrackingAlgorithm.getTrackingAlgorithm(args[i+1]);
                
                CL_Initializer.markAsParsed(i, 2);
                
            }
            else if (args[i].equals("-tendffile")) {
                tendF_File = args[i + 1];
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-tendf")) {
                tendF = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-tendg")) {
                tendG = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-silent")) {
                silent = true;
                CL_Initializer.markAsParsed(i);
            }
	    if (args[i].equals("-vectorfiles")) {

		String[] tmp = new String[100];

		int t = 0;

		while (t + i + 1 < args.length && !args[t + i + 1].startsWith("-")) {
		    tmp[t] = args[i + t + 1];
		    t++;
		}

		vectorFiles = new String[t];
		
		System.arraycopy(tmp, 0, vectorFiles, 0, t);

		CL_Initializer.markAsParsed(i, t + 1);
		
	    }
	    if (args[i].equals("-bedpostxdir")) {
                bedpostxDir = args[i+1];
                CL_Initializer.markAsParsed(i, 2);
            }
	    if (args[i].equals("-bedpostxminf")) {
		bedpostxMinF = Double.parseDouble(args[i+1]);
		CL_Initializer.markAsParsed(i, 2);	
	    }
	}
        CL_Initializer.checkParsing(args);
	
    }
    
    
    public void initVariables() {

        // do some default setting and sanity checks on args

	// If the user specifies -header, use that
	// Else look for a header to use
	if (CL_Initializer.headerTemplateFile == null) {

	    if (ImageHeader.imageExists(CL_Initializer.inputFile)) {
		CL_Initializer.headerTemplateFile = CL_Initializer.inputFile;
	    }
	    else if (imageType == imageType.BEDPOSTX || imageType == imageType.BEDPOSTX_DYAD) {

		String bedpostxRoot = BedpostxTractographyImage.getBedpostxInputRoot(bedpostxDir);
		String ext = BedpostxTractographyImage.getBedpostxImageExtension(bedpostxRoot);
		
		CL_Initializer.headerTemplateFile = bedpostxRoot + "dyads1" + ext;
	    
	    }
	    
	}

	
        if (CL_Initializer.headerTemplateFile != null) { 
            logger.info("Defining input physical space from " + CL_Initializer.headerTemplateFile);
            CL_Initializer.initInputSpaceAndHeaderOptions();
        }
        else if (CL_Initializer.voxelDims[0] > 0.0) { 
            CL_Initializer.initInputSpaceAndHeaderOptions();
        }
        else if (seedFile != null) {
            logger.info("Defining input physical space from seed file");
            CL_Initializer.headerTemplateFile = seedFile;
            CL_Initializer.initInputSpaceAndHeaderOptions();
        }
     


        xDataDim = CL_Initializer.dataDims[0];
        yDataDim = CL_Initializer.dataDims[1];
        zDataDim = CL_Initializer.dataDims[2];
        
        xVoxelDim = CL_Initializer.voxelDims[0];
        yVoxelDim = CL_Initializer.voxelDims[1];
        zVoxelDim = CL_Initializer.voxelDims[2];

        if (xVoxelDim == 0.0) {
            // failed to get a definition of input space from anywhere
            throw new LoggedException("Definition of input space required, use -header");
        }


        voxelToPhysicalTrans = CL_Initializer.headerTemplate.getVoxelToPhysicalTransform();


	if (anisThresh > 0.0)  {
	    switch(imageType) {
		
	    case DT : case MULTITENSOR : case BEDPOSTX : case BEDPOSTX_DYAD :
		// No problem since these input formats have anisotropy information built in
		break;
		
	    default: 
		
		if (anisMapFile == null) {
		    throw new LoggedException("Input data does not contain anisotropy, anisotropy map (-anisfile) must be " +
					      "supplied when -anisthresh is used");
		}
	    }
	}


        // no interpolation with FACT, by definition
        if (trackingAlgorithm == TrackingAlgorithm.FACT) {

            if (dataInterpolation != DataInterpolation.NEAREST_NEIGHBOUR) { 
                logger.warning("Interpolation is not compatible with FACT tracking, using Euler tracker with step size " + stepSize + " mm");
            }

            trackingAlgorithm = TrackingAlgorithm.EULER;
        }

        
        ran = new MTRandom(seed);
        
                
        // get seeds
        if (seedFile != null) {
            
            VoxelROI imageROIs = new VoxelROI(seedFile, CL_Initializer.headerTemplate);

            allROIs = imageROIs.getAllRegions();
            
        }
        else if (seedList != null) {
            allROIs = new RegionOfInterest[] {PointListROI.readPoints(seedList, CL_Initializer.headerTemplate)};
        }
	else {
	    throw new LoggedException("No seed points specified");
	}
       
       
    }
    
    
    public void execute(OutputManager om) {
        

        // All internal operations are in Camino space; we warp tracts to physical space before
        // writing them out

        int numRegions = allROIs.length;
        
        
        // get the voxel classification
        // either from a voxel class map or from a brain / background segmentation
        int[][][] vc = new int[xDataDim][yDataDim][zDataDim];
        
        // number of PDs in each voxel (for Bayesian images)
        int[][][] voxelNumPDs = new int[xDataDim][yDataDim][zDataDim];
        
        for (int k = 0; k < zDataDim; k++) {
            for (int j = 0; j < yDataDim; j++) {
                for (int i = 0; i < xDataDim; i++) {
                    vc[i][j][k] = 2;
                    voxelNumPDs[i][j][k] = 1;
                }
            }
        }
        
        if (CL_Initializer.bgMaskFile != null) {
            
            if (ImageHeader.imageExists(CL_Initializer.bgMaskFile)) {
                if ( !CL_Initializer.headerTemplate.sameSpace(CL_Initializer.bgMaskFile) ) {
                    throw new LoggedException("Brain mask must be in the same voxel space as the input data");
                }
            }

            CL_Initializer.initMaskSource();
            
            for (int k = 0; k < zDataDim; k++) {
                for (int j = 0; j < yDataDim; j++) {
                    for (int i = 0; i < xDataDim; i++) {
                        
                        double maskValue = CL_Initializer.bgMask.nextVoxel()[0];
                        
                        vc[i][j][k] = maskValue > 0.0 ? 2 : -1;
                        voxelNumPDs[i][j][k] = maskValue > 0.0 ? 1 : 0;
                    }
                }
            }
            
        }
        
        
        // use VC from file if we have it - note this overrides the bgmask

        // This needs generalizing to support NIfTI files / other formats. Ideally it should be decoupled from the
        // SH order (which is ambiguous over order 4 anyhow), it should be simply the number of PDs in the voxel
        //
        // Actually doing this will require work to the TractographyImage tree and other classes
        //
        if (CL_Initializer.voxelClassMap != null) {
            vc = new int[xDataDim][yDataDim][zDataDim];
            
            DataSource vcSource = ExternalDataSource.getDataSource(CL_Initializer.voxelClassMap, 1, "int");
            
            for (int k = 0; k < zDataDim; k++) {
                for (int j = 0; j < yDataDim; j++) {
                    for (int i = 0; i < xDataDim; i++) {
                        vc[i][j][k] = (int)vcSource.nextVoxel()[0];
                    }
                }
            }
        }
        
        // get the anisotropy map, if any
        //
        // Could be FA or any scalar image where we cease tracking below a threshold
        double[][][] anisMap = null;
        
        if (anisMapFile != null) {
            if (ImageHeader.imageExists(anisMapFile)) {
                try {
                    
                    ImageHeader ih = ImageHeader.readHeader(anisMapFile);

                    
                    if (!CL_Initializer.headerTemplate.sameSpace(ih)) {
                        throw new LoggedException("Anisotropy image must be in the same voxel space as the input data");
                    }


                    anisMap = ih.readSingleVolumeData();
                }
                catch (IOException e) {
                    throw new LoggedException(e);
                }
            }
            else {
                anisMap = new double[xDataDim][yDataDim][zDataDim];
                
                DataSource anisSource =
                    ExternalDataSource.getDataSource(anisMapFile, 1, CL_Initializer.inputDataType);
                
                for (int k = 0; k < zDataDim; k++) {
                    for (int j = 0; j < yDataDim; j++) {
                        for (int i = 0; i < xDataDim; i++) {
                            anisMap[i][j][k] = anisSource.nextVoxel()[0];
                        }
                    }
                }
                
            }
        }


        // set up the image
        
        TractographyImage image = null;
        
        switch (imageType) {

	case BEDPOSTX : case BEDPOSTX_DYAD :

	    boolean probabilistic = ( imageType == ImageType.BEDPOSTX ); 
	    
	    image = BedpostxTractographyImage.getTractographyImage(bedpostxDir, probabilistic, bedpostxMinF, anisMap, anisThresh,
								   new int[] {xDataDim, yDataDim, zDataDim},
								   new double[] {xVoxelDim, yVoxelDim, zVoxelDim}, ran);
	    
	    break;

                
        case REPBS_DT: case REPBS_MULTITENSOR: 
            
            CL_Initializer.initImagingScheme();
            
                
            image =
                RepBS_DWI_TractographyImage.getTractographyImage(CL_Initializer.bsDataFiles,
                                                                 CL_Initializer.inputDataType,
                                                                 CL_Initializer.imPars,
                                                                 CL_Initializer.inversionIndices,
                                                                 vc, anisMap, anisThresh,
                                                                 new int[] {xDataDim, yDataDim, zDataDim},
                                                                 new double[] {xVoxelDim, yVoxelDim, zVoxelDim}, ran);
            
            break;
            
            
        case WILDBS_DT:
	    CL_Initializer.initImagingScheme();
                
	    image =
		DT_WildBS_DWI_TractographyImage.getTractographyImage(CL_Initializer.inputFile,
								     CL_Initializer.inputDataType,
								     CL_Initializer.imPars,
								     vc, anisMap, anisThresh,   
								     new int[] {xDataDim, yDataDim, zDataDim},
								     new double[] {xVoxelDim, yVoxelDim, zVoxelDim}, ran);
                
	    break;
                
                
                
        case PICO :
                    
            image = PICoTractographyImage.getTractographyImage
                (CL_Initializer.inputFile, CL_Initializer.inputDataType, CL_Initializer.numPDsIO,
                 picoPDF, anisMap, anisThresh, new int[] {xDataDim, yDataDim, zDataDim},
                 new double[] {xVoxelDim, yVoxelDim, zVoxelDim}, ran);
            
            break;
            
        case BAYESDIRAC: case BAYESDIRAC_DT :
            
            CL_Initializer.initImagingScheme();
            
            // Friman Bayesian method
            // add anisotropy map option

            BayesDataModel dataModel = imageType == ImageType.BAYESDIRAC ? BayesDataModel.BALL_STICK : BayesDataModel.CYL_SYMM_DT;

            BayesDiracTractographyImage bi = BayesDiracTractographyImage.getTractographyImage
                (CL_Initializer.inputFile, CL_Initializer.inputDataType,
                 CL_Initializer.imPars, dataModel, voxelNumPDs, anisMap, anisThresh,
                 new int[] {xDataDim, yDataDim, zDataDim},
                 new double[] {xVoxelDim, yVoxelDim, zVoxelDim}, CL_Initializer.pointSetInd, ran);
            
            if (curvePriorK > 0.0) {
                bi.setCurvePriorKappa(curvePriorK);
            }
            if (curvePriorG > 0.0) {
                bi.setCurvePriorGamma(curvePriorG);
            }
            
            if (externalPriorImageFile != null) {
                PICoTractographyImage ePrior = PICoTractographyImage.getTractographyImage
                    (externalPriorImageFile, externalPriorDataType, CL_Initializer.numPDsIO,
                     picoPDF, anisMap, anisThresh, new int[] {xDataDim, yDataDim, zDataDim},
                     new double[] {xVoxelDim, yVoxelDim, zVoxelDim}, ran);
                
                bi.setExternalPriors(ePrior);
            }
            
            image = bi;
            
            break;
            
        case SF_PEAK :
            image =
                SF_TractographyImage.getTractographyImage(CL_Initializer.inputFile, CL_Initializer.inputDataType,
                                                          CL_Initializer.numPDsIO, anisMap, anisThresh,
                                                          new int[] {xDataDim, yDataDim, zDataDim},
                                                          new double[] {xVoxelDim, yVoxelDim, zVoxelDim});
            break;
            
        case DT: case MULTITENSOR :
            image =
                DT_TractographyImage.getTractographyImage(CL_Initializer.inputFile, CL_Initializer.inputDataType,
                                                          CL_Initializer.maxTensorComponents, anisMap, anisThresh,
                                                          new int[] {xDataDim, yDataDim, zDataDim},
                                                          new double[] {xVoxelDim, yVoxelDim, zVoxelDim});
            break;

        case DWI_DT : case DWI_MULTITENSOR:

            CL_Initializer.initImagingScheme();

            image =
                DWI_TractographyImage.getTractographyImage(CL_Initializer.inputFile, CL_Initializer.inputDataType,
                                                           CL_Initializer.imPars, CL_Initializer.inversionIndices,
                                                           vc, anisMap, anisThresh,                                                         
                                                           new int[] {xDataDim, yDataDim, zDataDim},
                                                           new double[] {xVoxelDim, yVoxelDim, zVoxelDim});
            break;

	case VECTOR : 

	    image = PD_TractographyImage.getTractographyImage(vectorFiles, CL_Initializer.inputDataType,
							      anisMap, anisThresh,   
							      new int[] {xDataDim, yDataDim, zDataDim},
							      new double[] {xVoxelDim, yVoxelDim, zVoxelDim});
	    break;
            
        default : throw new LoggedException("Unsupported image type : " + imageType);
            
        }
        
        
        // Set up the interpolation

        ImageInterpolator interp = null;

       
        switch (dataInterpolation) {


        case NEAREST_NEIGHBOUR: 
            
            interp = new NearestNeighbourInterpolator(image);

            break;

        case NEIGHBOUR_CHOICE: 

            interp = new NeighbourChoiceInterpolator(image, ran);

            break;

        case DWI_TRILINEAR: 

            interp = new DWI_LinearInterpolator((DWI_TractographyImage)image);

            break;

        case VECTOR_TRILINEAR: 

            interp = new VectorLinearInterpolator(image);

            break;

        case TEND_NN: case TEND_NC:

            // how to interpolate the tensor data itself, in order to provide
            // the tensor for the TEND term
            TensorInterpolator dataInterp = null;

            if (dataInterpolation == DataInterpolation.TEND_NC) {
                dataInterp = new DT_NC_Interpolator((TensorTractographyImage)image, ran);
            }
            else {
                dataInterp = new DT_NN_Interpolator((TensorTractographyImage)image);
            }

            if (tendF_File != null) {
                try {                                
                    ImageHeader ih = ImageHeader.readHeader(tendF_File);
                    
                    interp = new TendInterpolator((TensorTractographyImage)image, dataInterp, ih.readSingleVolumeData(), tendG);
                }
                catch (IOException e) {
                    throw new LoggedException(e);

                }
            
            }
            else {
                interp = new TendInterpolator((TensorTractographyImage)image, dataInterp, tendF, tendG);   
            }
            
            break;

        default: throw new LoggedException("Unsupported interpolation : " + dataInterpolation);


        }

        
        // set up the tracker

        FibreTracker tracker = null;


        switch (trackingAlgorithm) {

        case FACT:
 
            tracker = new FACT_FibreTracker(image);
            
            break;
            
        case EULER: 
            
            tracker = new EulerFibreTracker(interp, stepSize); 

            break;

        case RK4: 

            tracker = new RK4FibreTracker(interp, stepSize);

            break;
        }

        tracker.setCurveCheckInterval(checkCurveLength);
        tracker.setIP_Threshold(ipThresh);
        
        
        // And finally, do the tracking
        
        regions : for (int region = 0; region < allROIs.length; region++) {

            RegionOfInterest roi = allROIs[region];
            
            if (regionIndex > -1) {
                if (roi.getRegionLabel() != regionIndex) {
                    continue;
                }
            }
            
            
            int outputRegionID = roi.getRegionLabel();
            
            
            // points defined in Camino space
            Point3D[] seeds = roi.getSeedPoints();
            
            if (!silent) {
                System.err.println("Processing ROI " + (region + 1) + " of " + numRegions);
            }
          
            FileOutputStream fout = null;
            DataOutputStream dout = null;
            
            try {
                
                if (outputRoot == null) {
                    dout = om.getOutputStream();
                }
                else {
                    
                    if (gzip) {
                        fout = new FileOutputStream(outputRoot + outputRegionID + ".Bfloat.gz");
                        dout = new DataOutputStream(new GZIPOutputStream(fout, 1024*1024*16));
                        
                    }
                    else {
                        fout = new FileOutputStream(outputRoot + outputRegionID + ".Bfloat");
                        dout = new DataOutputStream
                            (new BufferedOutputStream(fout, 1024*1024*16));
                    }
                }
                
                
                
                
                seed: for (int sp = 0; sp < seeds.length; sp++) {
                    
                    if (!silent) {
                        System.err.print("\rProcessing seed " + (sp + 1) + " of " + seeds.length);
                    }
                    
                    Point3D seedPoint = seeds[sp];
                    
                    if (!tracker.inBounds(seedPoint)) {
                        logger.warning("Seed point \n\t" + seedPoint + "\n is outside the diffusion image space, ignoring");
                        continue seed;
                    }

                    int xVox = (int)(seedPoint.x / xVoxelDim);
                    int yVox = (int)(seedPoint.y / yVoxelDim);
                    int zVox = (int)(seedPoint.z / zVoxelDim);
                    
                    int numPDs = image.numberOfPDs(xVox, yVox, zVox);
                    
                    // if number of PDs is zero, track once
                    // tracker will return the seed point
                    numPDs = numPDs == 0 ? 1 : numPDs;
                    
                    for (int p = 0; p < numPDs; p++) {
                        
                        for (int i = 0; i < iterations; i++) {
                            
                            Tract t = tracker.trackFromSeed(seedPoint, p);
                            
                            // warp tracts to physical space
                            t.transformToPhysicalSpace(voxelToPhysicalTrans, xVoxelDim, yVoxelDim, zVoxelDim);
                            
                            t.writeRaw(dout);
                            
                        }
                        
                    }
                    
                } // end for seedpoints
                
                if (!silent) {
                    System.err.println("\n");
                }
                
                if (outputRoot != null) {
                    // want to close file
                    dout.close();
                }
                
            }
            catch(IOException e) {
                throw new LoggedException(e);
            }
            
        } // end for all ROIs
        
        // close om stream
        if(om != null)
            om.close();
        
    }
    
   
}

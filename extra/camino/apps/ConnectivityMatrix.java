package apps;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;
import tractography.*;

import java.io.*;
import java.text.*;
import java.util.logging.Logger;

/**
 * Matrices of connectivity derived from streamlines
 *
 * @author Philip Cook
 */
public class ConnectivityMatrix extends Executable {
   
    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.ConnectivityMatrix");


    // for defining nodes of connectivity matrix - defines input space
    private String targetFile;
    
    // for stats on tracts
    private String statFile;

    // Names the labels present in the target file
    private String labelDefinitionFile;

    // Optional seed file for connectivity segmentation
    private String seedFile;
    
    private String tractStat;

    private String outputRoot;
    
    private DecimalFormat df;

    // true if we should do CBS, not pairwise connectivity
    private boolean doCBS;
    
    // for pairwise conmat
    private StreamlineROI_Filter filter;

    public ConnectivityMatrix(String[] args) {
        super(args);
    }


    public void initDefaultVals() {
        targetFile = null;
	seedFile = null;
        statFile = null;
        labelDefinitionFile = null;
        tractStat = null;
        outputRoot = "conmat_";
	doCBS = false;
        filter = null;
	df = new DecimalFormat("0.00000E00");
    }
    

    public void initOptions(String[] args) {
       	
	// Affects image output
        OutputManager.outputDataType = "float";

        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-targetfile")) {
                targetFile = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
	    if (args[i].equals("-seedfile")) {
                seedFile = args[i+1];
                CL_Initializer.markAsParsed(i,2);

            }
            else if (args[i].equals("-scalarfile")) {
		statFile = args[i+1];
		CL_Initializer.markAsParsed(i, 2);

	    } 
	    else if (args[i].equals("-tractstat")) {
		tractStat = args[i+1];
	    	CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-targetnamefile")) {
		labelDefinitionFile = args[i+1];
	    	CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-outputroot")) {
		outputRoot = args[i+1];
	    	CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-outputcbs")) {
		doCBS = true;
		CL_Initializer.markAsParsed(i);
	    }
	    
        }
	
	if (doCBS && seedFile == null) {
	    throw new LoggedException("Seed file required for connectivity segmentation");
	}

        if (statFile != null && tractStat == null) {
            // default to mean
            tractStat = "mean";
        } 

        CL_Initializer.headerTemplateFile = targetFile;

        CL_Initializer.initInputSpaceAndHeaderOptions();
        
        CL_Initializer.checkParsing(args);
    }

    public void initVariables() {

    }


    public void execute(OutputManager om) {
        
       
        if (doCBS) {
            writeCBS_Matrix();
        }
        else {
            writeConnectivityGraph();
        }




    }

    public void writeCBS_Matrix() {
        
        ImageHeader targetHdr;

	ImageHeader seedHdr;

        try {
            targetHdr = ImageHeader.readHeader(targetFile);
	    seedHdr = ImageHeader.readHeader(seedFile);

	    if (!targetHdr.sameSpace(seedHdr)) {
		throw new LoggedException("Target and seed images must be in the same space");
	    } 
	    
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }

        int[][][] targets = ProcessStreamlines.readIntVolume(targetFile);

	int[][][] seeds = ProcessStreamlines.readIntVolume(seedFile);

	
	double[] voxelDims = targetHdr.getVoxelDims();
	int[] dataDims = targetHdr.getDataDims();


        TractSource source = new TractSource(CL_Initializer.inputFile, targetHdr);
        
        TargetCP_Image image;

        if (labelDefinitionFile != null) {
            image = new TargetCP_Image(targets, voxelDims, labelDefinitionFile);
        }
        else {
            image = new TargetCP_Image(targets, voxelDims);
        }

	int numNodes = image.numNodes();

	String[] targetNames = image.getTargetNames();
	int[] targetLabels = image.getTargetLabels();

	boolean done = false;

	Tract t = source.nextTract();

	FileOutput out = new FileOutput(outputRoot + "sc.csv");        

	double[][][] seedConnLabels = new double[dataDims[0]][dataDims[1]][dataDims[2]];
	double[][][] seedMaxSC = new double[dataDims[0]][dataDims[1]][dataDims[2]];

	String[] header = new String[numNodes + 5];

	header[0] = "SeedX";
	header[1] = "SeedY";
	header[2] = "SeedZ";
	header[3] = "SeedLabel";
	header[4] = "InputStreamlines";
	
	System.arraycopy(targetNames, 0, header, 5, numNodes);

	out.writeString(join(header, ",") + "\n");
          
	while (!done) {
  
	    TractCollection tc = new TractCollection(1000, 100.0);
	    
	    Point3D seed = t.getSeedPoint();    
	    
	    Point3D nextSeed = seed;
	    
	    while (!done && seed.equals(nextSeed)) {

		tc.addTract(t);

		if (source.more()) {
		    t = source.nextTract();
		    nextSeed = t.getSeedPoint();
		}
		else {
		    done = true;
		}
	    }


	    // process all tracts from this seed point
	    image.processTracts(tc);

	    // output a line of the matrix - seedX, seedY, seedZ, seedLabel, <connectivity>
	    int[] sc = image.getStreamlineCountList();

	    double[] line = new double[sc.length + 5];

            // Got to put seed point back into physical space
            // It is currently in Camino space
	    // 
            Point3D seedVoxel = new Point3D( seed.x / seedHdr.xVoxelDim() - 0.5, seed.y / seedHdr.yVoxelDim() - 0.5, seed.z / seedHdr.zVoxelDim() - 0.5 );

            Point3D seedPointPhys = seedVoxel.transform(seedHdr.getVoxelToPhysicalTransform());

	    line[0] = seedPointPhys.x;
	    line[1] = seedPointPhys.y;
	    line[2] = seedPointPhys.z;

	    int[] seedVox = new int[3];

	    seedVox[0] = (int)(seed.x / voxelDims[0]);
	    seedVox[1] = (int)(seed.y / voxelDims[1]);
	    seedVox[2] = (int)(seed.z / voxelDims[2]);

	    line[3] = seeds[seedVox[0]][seedVox[1]][seedVox[2]];
	    
	    line[4] = image.totalStreamlines();


	    // different types so no arraycopy
	    for (int i = 0; i < sc.length; i++) {
		line[5+i] = sc[i];
	    }

	    String outLine = join(line, ",") + "\n";

	    out.writeString(outLine);


	    int maxSC = 0;
	    int maxIndex = -1;
	    
	    for (int i = 0; i < numNodes; i++) {
		if (sc[i] > maxSC) {
		    maxSC = sc[i];
		    maxIndex = i; 
		} 
	    }

	    seedConnLabels[seedVox[0]][seedVox[1]][seedVox[2]] = maxIndex > -1 ? targetLabels[maxIndex] : 0;

	    seedMaxSC[seedVox[0]][seedVox[1]][seedVox[2]] = maxSC;


	    image.reset();

	}

        out.close();

	// write images as ints
	CL_Initializer.headerTemplate.writeScalarImage(seedConnLabels, outputRoot + "cbs_labels");
	CL_Initializer.headerTemplate.writeScalarImage(seedMaxSC, outputRoot + "cbs_sc");
	

    }


    public void writeConnectivityGraph() {


        ImageHeader targetHdr;

        try {
            targetHdr = ImageHeader.readHeader(targetFile);
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }

        int[][][] targets = ProcessStreamlines.readIntVolume(targetFile);

        filter = new StreamlineROI_Filter(CL_Initializer.dataDims, CL_Initializer.voxelDims);

        TractSource source = new TractSource(CL_Initializer.inputFile, targetHdr);
    
        StreamlineConnectivityGraph graph;

        if (labelDefinitionFile != null) {
            graph = new StreamlineConnectivityGraph(targets, CL_Initializer.voxelDims, labelDefinitionFile);
        }
        else {
            graph = new StreamlineConnectivityGraph(targets, CL_Initializer.voxelDims);
        }


        if (tractStat != null) {

            TractStatisticFilter statFilter;

            if (statFile != null) {
                ImageHeader scalarHdr;
                
                try {
                    scalarHdr = ImageHeader.readHeader(statFile);

                    ScalarImage image = new ScalarImage(statFile);
                    
                    statFilter = new TractStatisticFilter(image);
                    
                }
                catch (IOException e) {
                    throw new LoggedException(e);
                }
                
                if (!targetHdr.sameSpace(scalarHdr)) {
                    throw new LoggedException("Target and stat file must be in the same space");
                }
                
            }
            else {
                statFilter = new TractStatisticFilter(targetHdr.getDataDims(), targetHdr.getVoxelDims());
            }
            
            statFilter.setTractStatistic(tractStat);

            statFilter.setInterpolate(true);

            graph.setTractStatFilter(statFilter);
        }       

           
        while (source.more()) {
            graph.processTracts(filter.processTract(source.nextTract()));
        }


        // output streamline counts and tract statistics

        FileOutput out = new FileOutput(outputRoot + "sc.csv");

        out.writeString(graph.getStreamlineCountMatrixCSV());

        out.close();

        if (tractStat != null) {
            out = new FileOutput(outputRoot + "ts.csv");
            
            out.writeString(graph.getTractStatisticMatrixCSV());
            
            out.close();
            
        }
        
    } 
    

    private String join(double[] values, String delim) {
	StringBuffer buff = new StringBuffer(values.length * 11);

	for (int i = 0; i < values.length - 1; i++) {
	    buff.append(df.format(values[i]));
	    buff.append(delim);
	}

	buff.append(df.format(values[values.length - 1]));

	return buff.toString();

    }


    private String join(String[] values, String delim) {
	StringBuffer buff = new StringBuffer(values.length * 11);

	for (int i = 0; i < values.length - 1; i++) {
	    buff.append(values[i]);
	    buff.append(delim);
	}

	buff.append(values[values.length - 1]);

	return buff.toString();

    }


}


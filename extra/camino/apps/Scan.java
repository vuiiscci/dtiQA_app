package apps;

import imaging.SimulableScheme;

import java.util.logging.Logger;

import misc.LoggedException;

import simulation.SimulationParams;
import simulation.measurement.AgnosticScan;
import simulation.measurement.ScanFactory;
import tools.CL_Initializer;

import data.OutputManager;


/**
 * Executable module for the scan command. For some reason this
 * particular command wasn't included with the rest of the 
 * camino commands in the entrypoint framework, and for some 
 * other reason wasn't updated with the new schemefile 
 * framework. So well done there.
 * 
 * This takes a trajfile and a chemefile and generates noise-free
 * signals.
 * 
 * @author Matt Hall (matt.hall@ucl.ac.uk)
 *
 */
public class Scan extends Executable {

	// logging object
	private static final Logger logger= Logger.getLogger("apps.Scan");

	/**
	 * constructor, passed commandline args
	 * 
	 * @param args
	 */
	public Scan(String[] args){
		super(args);
	}
	
	
	/** 
	 * nothing to do here
	 */
	public void initDefaultVals() {

	}

	/** 
	 * nothing here either
	 */
	public void initVariables() {

	}

	/** 
	 * over-ride the default initOptions. This is in-line with what's
	 * in SyntheticData. Here we don't care about input data type as
	 * this command doesn't read from stdin, but we do care about
	 * output data type as this should be in line with the rest of 
	 * the output of SyntheticData.
	 * 
	 * Beyond that we just need to parse the commandline and initialise
	 * the imaging scheme.
	 * 
	 */
	public void initOptions(String[] args) {

        // The output defaults to float type.
        OutputManager.outputDataType = "float";
        
        // Parse the command line arguments
        CL_Initializer.CL_init(args);
        
        CL_Initializer.checkParsing(args);
        
        // Now initialize the acquisition scheme and model parameters.
        CL_Initializer.initImagingScheme();
        //CL_Initializer.initDataSynthesizer();
     
    }
	
	
	/** 
	 * executes the command. Gets the scheme object, checks for compatibility. Then
	 * gets the trajfile, creates an agnostic scan object and gets the signals. The 
	 * output manager then reads out.
	 * 
	 * Main data synthesis functionality is contained in AgnosticScan (@see AgnosticScan,
	 * {@link AgnosticScan}).
	 * 
	 * @param om OutputManager object for readout.
	 */
	public void execute(OutputManager om) {
			
		try{
			SimulableScheme simScheme= (SimulableScheme)CL_Initializer.imPars;
		}
		catch(ClassCastException cce){
			logger.severe("The scheme file provided doesn't contain enough information to generate measurements from simulation");
        	logger.severe("Simulation-based data synthesis requires gradent strength, pulse durations and timings to be specified");
        	logger.severe("Compatible scheme file formats include STEJSKALTANNER, TRSE, and genwave");
        	
			throw new LoggedException("Specified schemefile is not compatible with data synthesis from simulation. Scheme must be Simulable");
		}
		
		SimulableScheme scheme= (SimulableScheme)CL_Initializer.imPars;
		
		String trajfile= SimulationParams.trajFile;
		
		logger.info("Synthesising measurements from "+trajfile+" using "+CL_Initializer.schemeFile);

		AgnosticScan scan= new AgnosticScan(scheme, trajfile);
		
		double[] signals= scan.getSignalsFromTrajectories();
		
		om.output(signals);

		om.close();
	}

}

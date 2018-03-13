/* DiffusionSimulation.java created on 25-Nov-2005
 * (simulation)
 * 
 * author: matt (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation;

import imaging.RectGradSteTanScheme;
import imaging.SimulableScheme;

import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.Random;
import java.util.logging.Logger;

import numerics.MTRandom;

import misc.LoggedException;

import simulation.dynamics.StepGenerator;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.Walker;
import simulation.dynamics.WalkerFactory;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.geometry.elements.Cylinder;
import simulation.geometry.elements.SquashyCylinder;
import simulation.geometry.substrates.CylinderSubstrate;
import simulation.geometry.substrates.SquashyInflammationSubstrate;
import simulation.geometry.substrates.Substrate;
import simulation.geometry.substrates.SubstrateFactory;
import simulation.measurement.ScanFactory;
import simulation.measurement.StatisticsModule;
import simulation.measurement.StatisticsModuleFactory;
import simulation.measurement.SyntheticScan;
import tools.CL_Initializer;


import data.DataSource;
import data.DataSourceException;
import data.OutputManager;

/**
 * top level class in diffusion simulation. This implements the DataSource
 * interface and contains the main monte carlo loop over all walkers and
 * timesteps, calling the step generator and the geometry/substrate objects
 * to amend steps.
 *
 * it also uses the synthetic scan to assemble noisy signals and multiple 
 * voxels-worth of readings.
 *
 * @author matt
 *
 */
public class DiffusionSimulation implements DataSource {

    /** logging object */
    private Logger logger= Logger.getLogger(this.getClass().getName());
    
    /** the dimensionality of the simulation. 
     * this is the central value from which all others are taken */
    public static final int D=3; 

    /** the self-diffusion constant of water at room temperature */
    public static final double DIFF_CONST = CL_Initializer.DIFF_CONST;
    
    /** the membrane transition probability */
    private final double p;

    /** random number generator */
    //private final MTRandom twister = new MTRandom((1321839371<<32)|(129375817));
    private MTRandom twister = new MTRandom(CL_Initializer.seed+189);
    
    /** diffusion simulation parameters */
    private final SimulationParams simParams;
    
    /** synthetc scan parameters */
    private SyntheticScan synthScan;
       
    /** array of walkers */
    private final Walker[] walker;
    
    /** diffusion substrate */
    private Substrate substrate;
    
    /** likely border in which to clone bits of substrate */
    public static double border;
    
    /** object for generating random steps in the walk */
    private StepGenerator stepGenerator;
    
    /** number of voxels to generate */
    private final int numVoxels;
    
    /** index of the current voxel */
    private int voxel=0;
    
    /** flag indicating if we want a separate simulation for each voxel or  not */
    private final boolean separateRuns = SimulationParams.sim_separate_runs;
    
    /** number of timesteps in simulation */
    private final int tmax;
    
    /** time increment associated with a timestep */
    public final double dt;
    
    /** duration of simulation */
    private final double duration;
    
    /** counter for the number of times the substrate initialiser function is called */
    public static int calls=-1;
    
    /**
     * are we selecting only one iteration to output?
     * -1 if not, otherwise equal to iteration number
     */
    private final int onlyRun= SimulationParams.sim_onlyRun;
    
    /** file writer for trajectories */
    private final DataOutputStream trajWriter;
    
    /** file writer for statistical measures of dynamics */
    private final DataOutputStream statsWriter;

    /** statistics module */
    private final StatisticsModule statsMod;
    
    
    /** debug file writers */
    //private BufferedWriter intraWriter=null;
    //private BufferedWriter extraWriter=null;
    
    public DiffusionSimulation(SimulationParams simParams, SimulableScheme imParams){
        
        this.simParams=simParams;

        this.stepGenerator=StepGeneratorFactory.getStepGenerator(simParams);
        
        DiffusionSimulation.border= stepGenerator.getBorder();
        
        this.substrate=SubstrateFactory.getSubstrate(simParams.getGeometryType(), simParams);
    
        this.numVoxels=CL_Initializer.numVoxels;
    
        this.p=simParams.getP();
               
        this.walker=new Walker[simParams.getN_walkers()];
        
        this.synthScan=ScanFactory.getMeasurementModule(simParams, imParams, substrate, walker);

        this.tmax=simParams.getTmax();
        
        this.dt=simParams.getDt();

        SimulationParams.duration=tmax*dt;
        
        this.duration= SimulationParams.duration;
        
        this.trajWriter= null;
        
        if(SimulationParams.sim_statsfile!=null){
            
            FileOutputStream fos=null;
            try{
               fos = new FileOutputStream(SimulationParams.sim_statsfile);
            }
            catch(Exception e){
                throw new LoggedException(e);
            }
            
            statsMod= StatisticsModuleFactory.getStatsModule(walker);
            statsWriter= new DataOutputStream(new BufferedOutputStream(fos));
        }
        else{
            statsMod= null;
            statsWriter= null;
        }
        
        logger.info("running simulation: "+walker.length+" walkers, "+simParams.getTmax()+" timesteps, p= "+p);
        logger.info("dynamics duration "+imParams.getDuration());
    }

    
    public DiffusionSimulation(SimulationParams simParams){
        
        this.simParams=simParams;

        this.stepGenerator=StepGeneratorFactory.getStepGenerator(simParams);
        
        DiffusionSimulation.border= stepGenerator.getBorder();
        
        this.substrate=SubstrateFactory.getSubstrate(simParams.getGeometryType(), simParams);
    
        this.numVoxels=CL_Initializer.numVoxels;
    
        this.p=simParams.getP();
               
        this.walker=new Walker[simParams.getN_walkers()];
        
        this.synthScan=null;

        this.tmax=simParams.getTmax();
        
        this.dt=simParams.getDt();
        
        this.duration= SimulationParams.duration;
        
        /** initalise file output -- data output writer wrapped round a
         * 							 buffered stream, wrapped round a file stream!
         */
        try{
        	FileOutputStream fstream = null;
        	try{
        		String fname= SimulationParams.trajFile;
        		fstream= new FileOutputStream(fname);
        	}
        	catch(Exception e){
        		throw new LoggedException(e);
        	}
        	
        	this.trajWriter= new DataOutputStream(new BufferedOutputStream(fstream, simParams.buffsize));
        }
        catch(Exception ioe){
        	throw new LoggedException(ioe);
        }
    
        
        
        if(SimulationParams.sim_statsfile!=null){
            
            FileOutputStream fos;
            try{
                fos = new FileOutputStream(SimulationParams.sim_statsfile);
            }
            catch(Exception e){
                throw new LoggedException(e);
            }
            
            statsWriter= new DataOutputStream(new BufferedOutputStream(fos));
            statsMod= StatisticsModuleFactory.getStatsModule(walker);
            
        }
        else{
            statsMod= null;
            statsWriter= null;
        }
        
        
        
        logger.info("running simulation: "+walker.length+" walkers, "+simParams.getTmax()+" timesteps, p= "+p);
        logger.info("no scheme used, trajectories created instead. duration= "+duration);
    }
    
    /**
     * test constructor which is handed a substrate instead of getting one from the subtrate factory
     * 
     * @param simParams
     * @param imParams
     * @param substrate
     */
    public DiffusionSimulation(SimulationParams simParams, SimulableScheme imParams, Substrate substrate){
        
        this.simParams=simParams;

        this.stepGenerator=StepGeneratorFactory.getStepGenerator(simParams);
        
        DiffusionSimulation.border= stepGenerator.getBorder();
        
        this.substrate=substrate;
    
        this.numVoxels=CL_Initializer.numVoxels;
    
        this.p=simParams.getP();
               
        this.walker=new Walker[simParams.getN_walkers()];
        
        this.synthScan=ScanFactory.getMeasurementModule(simParams, imParams, substrate, walker);

        this.tmax=simParams.getTmax();
        
        this.dt=simParams.getDt();

        SimulationParams.duration=tmax*dt;
        
        this.duration= SimulationParams.duration;
        
        this.trajWriter= null;
        
        if(SimulationParams.sim_statsfile!=null){
            
            FileOutputStream fos=null;
            try{
               fos = new FileOutputStream(SimulationParams.sim_statsfile);
            }
            catch(Exception e){
                throw new LoggedException(e);
            }
            
            statsMod= StatisticsModuleFactory.getStatsModule(walker);
            statsWriter= new DataOutputStream(new BufferedOutputStream(fos));
        }
        else{
            statsMod= null;
            statsWriter= null;
        }
        
        logger.info("running simulation: "+walker.length+" walkers, "+simParams.getTmax()+" timesteps, p= "+p);
        logger.info("dynamics duration "+imParams.getDuration());
    }

    
    /** 
     * sets all walkers in initial conditions
     */
    public void initialiseWalkers(){	

    	logger.info("initialising spins...");
    	
        if(simParams.getInitialConditions()==SimulationParams.SPIKE){
            // initially delta-peaked at centre of substrate
            double midway=substrate.getPeakCoord();
            
            double[] midPoint=new double[D];
            
            
            for(int i=0; i<D; i++){
                midPoint[i]=midway;
            }
            
            for(int i=0; i<walker.length; i++){
                //walker[i]=new Walker(midPoint, stepGenerator, substrate, synthScan, trajWriter);
                walker[i]= WalkerFactory.getWalker(midPoint, stepGenerator, substrate, synthScan, trajWriter, simParams);
            }
        }
        else if(simParams.getInitialConditions()==SimulationParams.UNIFORM){
            // initially uniformly distributed across substrate
            double[] substrateSize=substrate.getSubstrateSize();
            
            double[] bottomLeft= new double[]{Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE};
            double[] topRight= new double[]{0.0, 0.0, 0.0};
            
            for(int i=0; i<walker.length; i++){
                double[] r0= new double[D];

                do{
	                for(int j=0; j<D; j++){
	                    r0[j]= twister.nextDouble()*substrateSize[j];
	                    
	                }
                }
                while(!substrate.positionOk(r0,stepGenerator.getWalkerRadius()));
                
                for(int j=0; j<D; j++){
                	// get top right and bottom left corners of envloping cube
	                if(r0[j]<bottomLeft[j]){
	                	bottomLeft[j]=r0[j];
	                }
	                if(r0[j]>topRight[j]){
	                	topRight[j]=r0[j];
	                }
                }
                
                
                //walker[i]= new Walker(r0, stepGenerator, substrate, synthScan, trajWriter);
                walker[i]= WalkerFactory.getWalker(r0, stepGenerator, substrate, synthScan, trajWriter, simParams);
            }
        }
        else if(simParams.getInitialConditions()==SimulationParams.INTRACELLULAR){
            // initially uniformly distributed across substrate
            double[] substrateSize=substrate.getSubstrateSize();
            
            boolean extracellular=true;
            
            double[] bottomLeft= new double[]{Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE};
            double[] topRight= new double[]{0.0, 0.0, 0.0};
            
        	for(int i=0; i<walker.length; i++){
        	    extracellular=true;
        	    
        		while(extracellular){
            		double[] r0= new double[D];
            		for(int j=0; j<D; j++){
            			r0[j]= twister.nextDouble()*substrateSize[j];
                    
            			// get top right and bottom left corners of enveloping cube
            			if(r0[j]<bottomLeft[j]){
            				bottomLeft[j]=r0[j];
            			}
            			if(r0[j]>topRight[j]){
            				topRight[j]=r0[j];
            			}
            		}
            		//walker[i]= new Walker(r0, stepGenerator, substrate, synthScan, trajWriter);
                    walker[i]= WalkerFactory.getWalker(r0, stepGenerator, substrate, synthScan, trajWriter, simParams);

            		
            		extracellular=!substrate.intracellular(walker[i]);
            	}
            }
        }
        else if(simParams.getInitialConditions()==SimulationParams.EXTRACELLULAR){
            // initially uniformly distributed across substrate
            double[] substrateSize=substrate.getSubstrateSize();
            boolean intracellular=true;
            
            double[] bottomLeft= new double[]{Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE};
            double[] topRight= new double[]{0.0, 0.0, 0.0};
            
            for(int i=0; i<walker.length; i++){
                intracellular=true;
                
            	while(intracellular){
	                double[] r0= new double[D];
	                for(int j=0; j<D; j++){
	                    r0[j]= twister.nextDouble()*substrateSize[j];
	                    
	                    // get top right and bottom left corners of envloping cube
	                    if(r0[j]<bottomLeft[j]){
	                    	bottomLeft[j]=r0[j];
	                    }
	                    if(r0[j]>topRight[j]){
	                    	topRight[j]=r0[j];
	                    }
	                }
	                
	                //walker[i]= new Walker(r0, stepGenerator, substrate, synthScan, trajWriter);
	                walker[i]= WalkerFactory.getWalker(r0, stepGenerator, substrate, synthScan, trajWriter, simParams);

	                intracellular=substrate.intracellular(walker[i]);
	            }
            }
            
            FileWriter walkerPos= null;
            
            
            logger.info("writing initial walker positions");
            try{
                
                walkerPos= new FileWriter("walkerPos-extra.csv");
                
                for(int i=0; i<walker.length; i++){
                    
                    walkerPos.write(walker[i].r[0]+","+walker[i].r[1]+"\n");
                }
                
                walkerPos.flush();
                
                walkerPos.close();
            }
            catch(IOException ioe){
                
                throw new LoggedException(ioe);
            }
            logger.info("walkr positions written");
            
            
        }
        else if(simParams.getInitialConditions()==SimulationParams.SPECIAL){

        	double[] substrateSize= substrate.getSubstrateSize();
        	
        	for(int i=2; i<walker.length; i++){
                double[] r0= new double[D];
            	for(int j=0; j<D; j++){
            		r0[j]= twister.nextDouble()*substrateSize[j];
                    
            		walker[i]= new Walker(r0, stepGenerator, substrate, synthScan, trajWriter);
            	}
        	}

        }
        else{
            String errMess= new String("unrecognised initial conditions code "+
                    	         simParams.getInitialConditions());
            logger.severe(errMess);
            throw new LoggedException(errMess);
        }
        
        logger.info("done.");
                
    }
    
    
    /** 
     *  runs the simualtion. tmax updates of N_walkers walkers
     *  using the simulation parameters given.
     */
    public void runMainLoop(){
          	
    	boolean[] isIntra= new boolean[walker.length];
        
    	substrate.init();        
    	
    	// write cylinders file if specified
    	if(SimulationParams.sim_cylFile!=null){
    		
    		try{
    			CylinderSubstrate cylSubs= (CylinderSubstrate)substrate;

				Cylinder[] cylinder= cylSubs.getCylinders();

    			try{
    				FileWriter fw= new FileWriter(SimulationParams.sim_cylFile);
    				
    				String newline= System.getProperty("line.separator");
    				
    				logger.info("Writing cylinders to file");
    				
    				for(int i=0; i<cylinder.length; i++){
    					double r= cylinder[i].getRadius();
    					double[] p= cylinder[i].getPosition();
    					
    					String cylLine= new String(r+","+p[0]+","+p[1]+","+p[2]+newline);
    					
    					fw.write(cylLine);
    				}
    				
    				fw.flush();
    				fw.close();
    				
    				logger.info("done.");
    			}
    			catch(IOException ioe){
    				throw new LoggedException(ioe);
    			}
    			
    			
    		}
    		catch(ClassCastException cce){
    			logger.warning("Cylnder output file specified with non-cylinder substrate. No file will be output.");
    		}
    		
    	}

    	
    	
    	if((onlyRun==-1)||(calls==onlyRun)){
    	    for(int i=0; i<walker.length; i++){
    	        isIntra[i]=substrate.intracellular(walker[i]);
    	    }

    	    // write initial information to traj file
            if(trajWriter!=null){
                try {
                    trajWriter.writeDouble(dt*(double)tmax);
                    trajWriter.writeDouble((double)simParams.getN_walkers());
                    trajWriter.writeDouble((double)tmax);
                } catch (IOException ioe) {
                    throw new LoggedException(ioe);
                }
			}
        
            int when=0;
            int who=42;
            boolean report=false;
            
            for(int t=0; t<tmax; t++){
                if((t%100)==0){
                	System.err.print("\r"+100.0*(double)t/(double)(simParams.getTmax())+"%     ");
                }
                
                for(int i=0; i<simParams.getN_walkers(); i++){
    
                	report= false;
					/*if(t==when){
						if(i==who){
							System.err.println("catch clause reached");
							report= true;
							double[] pos= new double[D];
							substrate.getSubstrateCoords(walker[i].r, new double[]{0.0, 0.0, 0.0}, pos);
							System.err.println("** pos: "+pos[0]+","+pos[1]+","+pos[2]);
							double size[]= substrate.getSubstrateSize();
							System.err.println("** size: "+size[0]+","+size[1]+","+size[2]);
							
						}
					}*/
                    	
                	if(p==0.0){
                		if(substrate.intracellular(walker[i])!=isIntra[i]){            		
                			logger.severe("walker pos: "+walker[i].r[0]+"  "+walker[i].r[1]+"  "+walker[i].r[2]);
                			
                			throw new LoggedException("t= "+t+" i= "+i+" has crossed. isIntra="+isIntra[i]);
                		}
                	}
                	
                	//System.err.println("t="+t+", i="+i);
                	walker[i].update(t*dt, t, i, report);
                }
                
                // stats measures if we're generating them
                if(statsMod!=null){
                    double[] stats= statsMod.getRuntimeStats(t*this.dt);
                    try{
                        statsWriter.writeDouble(t*dt);
                        for(int i=0; i<stats.length; i++){
                            statsWriter.writeDouble(stats[i]);
                        }
                    }
                    catch(IOException ioe){
                        throw new LoggedException(ioe);
                    }
                }
                
                
                // get the scan to do its thing...
                if(synthScan!=null){
                	synthScan.update(t);
                }
            }
            
            // write the final walker positions to traj file
            if(trajWriter!=null){
            	for(int i=0; i<simParams.getN_walkers(); i++){
            		try{
            			trajWriter.writeDouble(tmax*dt);
    	        		trajWriter.writeInt(i);
    	        		for(int j=0; j<D; j++){
    	        			trajWriter.writeDouble(walker[i].r[j]);
    	        		}
    	        	}
    	        	catch(IOException ioe){
    	        		throw new LoggedException(ioe);
    	        	}
            	}
            
            
    	        //flush and close the file writer
    	        try{
    	        	trajWriter.flush();
    	        	trajWriter.close();
    	        }
    	        catch(IOException ioe){
    	        	throw new LoggedException(ioe);
    	        }
            }
            
            if(statsMod!=null){
                // get final runtime stats measures
                double[] stats= statsMod.getRuntimeStats(tmax*dt);
                try{
                    statsWriter.writeDouble(tmax*dt);
                    for(int i=0; i<stats.length; i++){
                        statsWriter.writeDouble(stats[i]);
                    }
                    statsWriter.flush();
                    statsWriter.close();
                }
                catch(IOException ioe){
                    throw new LoggedException(ioe);
                }
                
                // get post-processing stats
                stats= statsMod.getPostSimulationStats();
                
                if(stats!=null){
                    // report stasts into log
                    logger.info("post-processing stats:");
                    for(int j=0; j<stats.length; j++){
                        logger.info(stats[j]+" ");
                    }
                    
                    // report stats into file
                    try{
                        DataOutputStream postproWriter= new DataOutputStream(
                                new FileOutputStream(SimulationParams.sim_postproStatsFname));

                        for(int i=0; i<stats.length; i++){
                            postproWriter.writeDouble(stats[i]);
                        }
                    }
                    catch(IOException ioe){
                        throw new LoggedException(ioe);
                    }
                }
            }
        }
    }
    
    
    
    /** 
     * initialises a simulation and runs the main loop before
     * constructing a 
     * @see data.DataSource#nextVoxel()
     */
    public double[] nextVoxel() throws DataSourceException {
        
        double[] S=null;
        
        if(!separateRuns){
        	// run the simulation once. subsequent voxels
        	// are different noise realisations of the 
        	// simulation data.
        	
        	if(voxel==0){
        		initialiseWalkers();

        		runMainLoop();
        	}
        }
        else{
        	// run a separate simulation for each voxel
        	
        	stepGenerator= StepGeneratorFactory.getStepGenerator(simParams);
        	initialiseWalkers();
    	    System.err.println("onlyrun= "+onlyRun+" calls= "+calls);
    	    runMainLoop();
        	
        	CL_Initializer.seed=twister.nextInt();
        	twister=new MTRandom(CL_Initializer.seed+189);
        }
        
        
        /**
         * debug code. spin out directional mean squared disp
         */
        /*
        int N=50;
        double[] netDisp= new double[D];
        double[] n= new double[D];
        
        FileWriter msdWriter;
        try{
            msdWriter= new FileWriter("angular_msd.csv");
        }
        catch(IOException ioe){
            throw new RuntimeException(ioe);
        }
        for(int t=0; t<N; t++){
            double theta= 2.0*t*Math.PI/N;
            
            n[0]=Math.cos(theta);
            n[1]= 0.0;
            n[2]=Math.sin(theta);
            
            double msd1=0.0;
            double msd2=0.0;
            
            int msdCount1=0;
            int msdCount2=0;
            
            for(int i=0; i<walker.length; i++){
                double dp=0.0;
                for(int j=0; j<D; j++){
                    netDisp[j]= walker[i].r[j]-walker[i].r0[j];
                    dp+=netDisp[j]*n[j];
                }
                
                if(walker[i].r0[1]<substrate.getSubstrateSize()[1]/2){
                    msd1+= dp*dp;
                    msdCount1++;
                }
                else{
                    msd2+= dp*dp;
                    msdCount2++;
                }
            }
            
            double msdTot= msd1+msd2;
            
            msd1/=msdCount1;
            msd2/=msdCount2;
            msdTot/=(msdCount1+msdCount2);
            
            try{
                msdWriter.write(theta+","+msd1+","+msd2+","+msdTot+","+msdCount1+","+msdCount2+"\n");
            }
            catch(IOException ioe){
                throw new RuntimeException(ioe);
            }
            
        }
        try{
            msdWriter.flush();
            msdWriter.close();
        }
        catch(IOException ioe){
            throw new RuntimeException(ioe);
        }
        */
        
        
        
        /**
         *  TODO: synthetic scan and trajectories currently 
         *  implemented in parallel. in future remove synthscan
         *  in simulation completely. 
         *  currently main output will be {-1.0, -1.0, -1.0} if 
         *  no scheme is passed to diffusion sim. 
         */
        if(synthScan!=null){
            // check compartmental output
            if(SimulationParams.sim_compartmentSignal==SimulationParams.INTRAONLY){
                // intracellular only
                S=synthScan.getCompartmentalSignals(true);
            }
            else if(SimulationParams.sim_compartmentSignal==SimulationParams.EXTRAONLY){
                // extracellular only
                S=synthScan.getCompartmentalSignals(false);
            }
            else if(SimulationParams.sim_compartmentSignal==SimulationParams.ALLCOMPS){
                /* all compartments, so must get data from intra and extracellular as
                 * well as the total, and then concatenate them together into a single
                 * array for output. this SHOULD work, but is a bit hacky.
                 */
                
                // get signals
                double[] Sin= synthScan.getCompartmentalSignals(true);
                double[] Sout= synthScan.getCompartmentalSignals(false);
                double[] Sall= synthScan.getSignals();
                
                // make new array
                S=new double[Sin.length+Sout.length+Sall.length];
                
                // concatenate signals into a single array
                for(int i=0; i<Sin.length; i++){
                    S[i]=Sin[i];
                }
                for(int i=0; i<Sout.length; i++){
                    S[Sin.length+i]=Sout[i];
                }
                for(int i=0; i<Sall.length;i++){
                    S[Sin.length+Sout.length+i]=Sall[i];
                }
            }
            else{
                // non-compartmental.
                S=synthScan.getSignals();
            }
        }
        else{
        	S= new double[] {-1.0, -1.0, -1.0};
        }
    	voxel++;
    
    	return S;
        
    }

    /** 
     * @see data.DataSource#more()
     */
    public boolean more() {
        if(voxel<numVoxels){
            return true;
        }
        else{
            return false;
        }
    }    
  
    /** calulate the mean square displacement of walkers at their current 
     *   locations. This is designed to be used by test code after a 
     *   simulation has run.
     * 
     * @return current mean square dispalacement of walkers
     */
    public double getMeanSquareDisplacement(){
        
        double meanSquareDisp=0.0;
        
        for(int i=0; i<walker.length; i++){
            
            double disp[]=walker[i].getDisplacement();
            double squareDisp=0.0;
            
            for(int j=0; j<D; j++){
                squareDisp+=disp[j]*disp[j];
            }
            
            
            
            meanSquareDisp+=squareDisp;
            
        }
        
        meanSquareDisp/=walker.length;
        
        return meanSquareDisp;
        
        
    }
    
    /**
     * @return the array of walker in the simulation
     */
    public final Walker[] getWalkers(){
    	
    	if(walker==null){
    		throw new LoggedException("attempt to retreive un-initialised walker array");
    	}
    	
    	return walker;
    }
    
    
    /**
     * 
     * @return substrate object
     * 
     */
    public final Substrate getSubstrate(){
    	
    	if(substrate==null){
    		throw new LoggedException("attempt to retreive uninitialised substrate from simulation");
    	}
    	
    	return substrate;
    }
    
    
    /**
     * @return step generator object
     */
    public final StepGenerator getStepGenerator(){
    	
    	if(stepGenerator==null){
    		throw new LoggedException("attempt to retreive uninitialised stepo generator from simulation");
    	}
    	
    	return stepGenerator;
    }
    
    /**
     * @return synthetic scan object
     */
    public final SyntheticScan getScan(){
    	
    	return synthScan;	
    }
    
    
    public static final void testDiffusionSimulation(){

        
        SimulationParams simParams;
        DiffusionSimulation diffSim;
        
        
        double[] stepParams = new double[] {0.01};

        SimulationParams.sim_l=1.0;
        SimulationParams.sim_L=3;
        
        SimulableScheme scheme;
        
        try{
            URI url= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= url.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }
        
        
        simParams= new SimulationParams(1, 70000, 0.002, SimulationParams.SPIKE, SubstrateFactory.SubstrateType.CYL_2_FIXED, 
                	           StepType.FIXEDLENGTH,
                	           1.5, scheme);
        
        simParams.setStepParams(stepParams);
        
        diffSim= new DiffusionSimulation(simParams, scheme);
        
        diffSim.nextVoxel();
        
        simParams= new SimulationParams(1000, 100000, 0.0001, SimulationParams.SPIKE, 
                SubstrateFactory.SubstrateType.CELL_STRIPED, StepType.FIXEDLENGTH, 
                1.5, scheme);
  
        diffSim= new DiffusionSimulation(simParams, scheme);
        
        diffSim.nextVoxel();
        

        SimulationParams.sim_l=1.0;
        SimulationParams.sim_L=20.0;
        SimulationParams.sim_stripethickness=3;
        
        simParams= new SimulationParams(100, 100000, 0.0001, SimulationParams.UNIFORM, 
                SubstrateFactory.SubstrateType.CELL_STRIPED, StepType.FIXEDLENGTH, 
                10.0, scheme);

        simParams.setStepParams(stepParams);
        
        diffSim= new DiffusionSimulation(simParams, scheme);
  
        diffSim.nextVoxel();        
    }

    
    private static void testFullScale(){
        // tests diffusion statistics on a striped substrate
        SimulationParams simParams;
        DiffusionSimulation diffSim;
        
        int tmax=100000;
        
        double p=0;


        double[] stepParams = new double[] {1E-8};

        SimulationParams.sim_l=1E-6;
        SimulationParams.sim_L=20;
        SimulationParams.sim_stripethickness=1;
        
        SimulableScheme scheme;
        
        try{
            URI url= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= url.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }

               
        simParams = new SimulationParams(100, tmax, 
                p, SimulationParams.SPIKE, SubstrateFactory.SubstrateType.CELL_STRIPED,
                StepType.FIXEDLENGTH, 0.5, scheme);
        
        stepParams=StepGeneratorFactory.getStepParamsArray(StepType.FIXEDLENGTH, simParams);
        
        simParams.setStepParams(stepParams);

        diffSim= new DiffusionSimulation(simParams, scheme);
        
        double[] S= diffSim.nextVoxel();
    	
    }
    
    
   private static void testReflectionFromCylinders(){
    	
        SimulationParams simParams;
        DiffusionSimulation diffSim;
        
        int tmax=1000;
        
        double p=0.0;

        double[] p_arr= new double[1];

        CL_Initializer.gamma_k=1E-6;
        CL_Initializer.gamma_beta=2E-5;
        SimulationParams.sim_cyl_dist_size=1;
        CL_Initializer.numVoxels=1;
        SimulationParams.sim_inflamm_increments=1;
        SimulationParams.sim_L=2.0;
        
        
        double[] stepParams = new double[] {DiffusionSimulation.DIFF_CONST};

        
        SimulableScheme scheme;
        
        try{
            URI url= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= url.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }

               
        simParams = new SimulationParams(1, tmax, 
                p, SimulationParams.UNIFORM, SubstrateFactory.SubstrateType.CYL_1_INFLAM,
                StepType.FIXEDLENGTH, 1.0, scheme);
        
        stepParams=StepGeneratorFactory.getStepParamsArray(StepType.FIXEDLENGTH, simParams);
        
        simParams.setStepParams(stepParams);

        StepGenerator stepGen= StepGeneratorFactory.getStepGenerator(simParams);
        
        diffSim= new DiffusionSimulation(simParams, scheme);
    	
        diffSim.substrate.init();
        
    	// that's everything. now we'll reach inside and set up the test case.
        SquashyCylinder cylinder= new SquashyCylinder(new double[]{1.0, 1.0, 0.0}, 0.5, 0.0);
        ((SquashyInflammationSubstrate)diffSim.substrate).cylinder[0]=cylinder;
    	
        
        double[] P = cylinder.getPosition();
        double[] V = cylinder.getAxis();
        double r = cylinder.getRadius();
        
        System.err.println("cylinder position =("+P[0]+","+P[1]+","+P[2]+")");
        System.err.println("cylinder axis =("+V[0]+","+V[1]+","+V[2]+")");
        System.err.println("radius = "+r);
        
        double[] r0 = new double[D];
        
        r0[0]=P[0]+0.9*r;
        r0[1]=P[1];
        r0[2]=0.0;
        
        Walker walker = new Walker(r0);
        
        double[] step= new double[D];
        
        step[0]=0.2*r;
        step[1]=0.0;
        step[2]=0.2*r;
        
        double[] normal= new double[D];
        double[] d= new double[1];
        double origLength=0.0;
        
        for(int i=0; i<D; i++){
        	origLength+=step[i]*step[i];
        }
        origLength=Math.sqrt(origLength);
      
        
        System.err.println("walker pos = "+walker.r[0]+","+walker.r[1]+","+walker.r[2]+")");
        System.err.println("step = ("+step[0]+","+step[1]+","+step[2]+")");
        
        boolean[] in = new boolean[1];
        
        double[] newStep= new double[]{0.0, 0.0, 0.0};
        
        
        boolean crosses=false;
		try {
			crosses = diffSim.substrate.crossesMembrane(walker, newStep, step, normal, d, false, origLength, in, p_arr, false, null);
		} catch (StepRejectedException e) {
			e.printStackTrace();
		}
    	
        System.err.println("crosses membrane is "+crosses);
        
        double[] toBarrier= new double[D];
        double[] amended= new double[D];
        double[] unamended= new double[D];
        
        diffSim.substrate.testAmendment(walker, step, normal, d, origLength, toBarrier, amended, unamended);
        
        System.err.println("toBarrier =("+toBarrier[0]+","+toBarrier[1]+","+toBarrier[2]+")");
        System.err.println("amended =("+amended[0]+","+amended[1]+","+amended[2]+",");
        System.err.println("unamended =("+unamended[0]+","+unamended[1]+","+unamended[2]+",");
        
        
        
    }
   
   
   	private static void testNonIntersection(){
   	
       SimulationParams simParams;
       DiffusionSimulation diffSim;
       
       int tmax=1000;
       
       double p=0.0;

/*       try{
           DiffusionSimulation.meanSquareWriter = new BufferedWriter(new FileWriter("squareDisps_p="+p+"_stripes_1.dat"));
       }
       catch(IOException ioe){
           throw new LoggedException(ioe);
       }
*/

       double[] stepParams = new double[] {DiffusionSimulation.DIFF_CONST};

       CL_Initializer.gamma_k=1E-6;
       CL_Initializer.gamma_beta=2E-5;
       SimulationParams.sim_cyl_dist_size=1;
       CL_Initializer.numVoxels=1;
       SimulationParams.sim_inflamm_increments=1;
       SimulationParams.sim_L=2.0;
       
       
       
       
       SimulableScheme scheme;
       
       try{
           URI url= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
           
           String path= url.getPath();
           
           scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
       }
       catch(URISyntaxException urise){
           throw new LoggedException(urise);
       }

              
       simParams = new SimulationParams(1, tmax, 
               p, SimulationParams.UNIFORM, SubstrateFactory.SubstrateType.CYL_1_INFLAM,
               StepType.FIXEDLENGTH, 1.0, scheme);
       
       stepParams=StepGeneratorFactory.getStepParamsArray(StepType.FIXEDLENGTH, simParams);
       
       simParams.setStepParams(stepParams);

       StepGenerator stepGen = StepGeneratorFactory.getStepGenerator(simParams);
       
       diffSim= new DiffusionSimulation(simParams, scheme);
   	
       diffSim.substrate.init();
       
   	// that's everything. now we'll reach inside and set up the test case.
       SquashyCylinder cylinder= new SquashyCylinder(new double[]{1.0, 1.0, 0.0}, 2.0, 0.0);
       ((SquashyInflammationSubstrate)diffSim.substrate).cylinder[0]=cylinder;
   	
       
       double[] P = cylinder.getPosition();
       double[] V = cylinder.getAxis();
       double r = cylinder.getRadius();
       
       System.err.println("cylinder position =("+P[0]+","+P[1]+","+P[2]+")");
       System.err.println("cylinder axis =("+V[0]+","+V[1]+","+V[2]+")");
       System.err.println("radius = "+r);
       
       double[] r0 = new double[D];
       
       r0[0]=P[0];
       r0[1]=P[1];
       r0[2]=0.0;
       
       Walker walker = new Walker(r0);
       
       double[] step= new double[D];
       
       step[0]=0.5*r;
       step[1]=0.0;
       step[2]=0.0;//0.2*r;
       
       double[] normal= new double[D];
       double[] d= new double[1];
       double origLength=0.0;
       
       for(int i=0; i<D; i++){
       	origLength+=step[i]*step[i];
       }
       origLength=Math.sqrt(origLength);
     
       
       System.err.println("walker pos = "+walker.r[0]+","+walker.r[1]+","+walker.r[2]+")");
       System.err.println("step = ("+step[0]+","+step[1]+","+step[2]+")");
       
       boolean[] in= new boolean[1];
       double[] offset= new double[]{0.0, 0.0, 0.0};
       double[] p_arr= new double[1];
       
       boolean crosses=false;
       try {
    	   crosses = diffSim.substrate.crossesMembrane(walker, offset, step, normal, d, false, origLength, in, p_arr, false, null);
       } catch (StepRejectedException e) {
    	   e.printStackTrace();
       }

       System.err.println("crosses membrane is "+crosses);
       
       double[] toBarrier= new double[D];
       double[] amended= new double[D];
       double[] unamended= new double[D];
       
       diffSim.substrate.testAmendment(walker, step, normal, d, origLength, toBarrier, amended, unamended);
       
       System.err.println("toBarrier =("+toBarrier[0]+","+toBarrier[1]+","+toBarrier[2]+")");
       System.err.println("amended =("+amended[0]+","+amended[1]+","+amended[2]+",");
       System.err.println("unamended =("+unamended[0]+","+unamended[1]+","+unamended[2]+",");
       
       
       
   }

    private static void testSingleCylinder(){
    	
    	SimulationParams simParams;
        DiffusionSimulation diffSim;
        
        int tmax=1000;
        
        double p=0.0;

        double L=2E-6;

        double[] stepParams = new double[] {1E-7};

        CL_Initializer.gamma_k=1E-6;
        CL_Initializer.gamma_beta=2E-5;
        SimulationParams.sim_cyl_dist_size=1;
        CL_Initializer.numVoxels=1;
        SimulationParams.sim_inflamm_increments=1;
        SimulationParams.sim_L=L;
        

        
        
        
        SimulableScheme scheme;
        
        try{
            URI url= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= url.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }
        simParams = new SimulationParams(100, 100000, 
                p, SimulationParams.SPIKE, SubstrateFactory.SubstrateType.CYL_1_INFLAM,
                StepType.FIXEDLENGTH, 1.0, scheme);
        
        //stepParams=StepGeneratorFactory.getStepParamsArray(StepGeneratorFactory.FIXEDLENGTH,
        //		simParams, imParams);
        
        simParams.setStepParams(stepParams);

        diffSim= new DiffusionSimulation(simParams, scheme);
    	
        SquashyInflammationSubstrate inflam=(SquashyInflammationSubstrate)diffSim.substrate;
        
        
        inflam.P[0][0]= L/2.0;
        inflam.P[0][1]= L/2.0;
        inflam.P[0][2]= L/2.0;
        
        
        diffSim.nextVoxel();
        double[] signal=diffSim.synthScan.getSignals();
        
        for(int i=0; i<signal.length; i++){
        	System.err.println("signal["+i+"] = "+signal[i]);
        }
        
        double[] meanSq= new double[]{0.0, 0.0, 0.0};
        
        for(int i=0; i<diffSim.walker.length; i++){
        	for(int j=0; j<D; j++){
        		meanSq[j]+=(diffSim.walker[i].r[j]-diffSim.walker[i].r0[j])*(diffSim.walker[i].r[j]-diffSim.walker[i].r0[j]);
        	}
        }
        for(int j=0; j<D; j++){
        	meanSq[j]/=diffSim.walker.length;
        	System.err.println("meanSq["+j+"]= "+meanSq[j]);
        }
    }
    
    private static void testPerfectIsotropy(){
        
        CL_Initializer.brownianSimulation= true;
        SimulationParams.sim_N_walkers= 10000;
        SimulationParams.sim_tmax= 100;
        SimulationParams.sim_geomType= SubstrateFactory.SubstrateType.EMPTY;
        SimulationParams.sim_initial= SimulationParams.SPIKE;
        CL_Initializer.numVoxels= 1;
        SimulationParams.sim_l= 1e-3;
        CL_Initializer.schemeFile= "facet_scheme_1000.scheme";
        SimulationParams.sim_G= 0.022;
        SimulationParams.sim_G_set= true;
        SimulationParams.sim_delta= 0.032;
        SimulationParams.sim_delta_set= true;
        SimulationParams.sim_DELTA= 0.04;
        SimulationParams.sim_DELTA_set= true;
        OutputManager.outputFile= "perfect.bfloat";
        
        CL_Initializer.initImagingScheme();
        CL_Initializer.initDataSynthesizer();
        
        
        DiffusionSimulation diffsim= (DiffusionSimulation)CL_Initializer.data;
        SyntheticScan scan= diffsim.synthScan;
        
        double[] s= diffsim.stepGenerator.getStep(null);
        double len=0.0;
        for(int i=0; i<s.length; i++){
            len+=s[i]*s[i];
        }
        len=Math.sqrt(len);
        
        double[][] step= new double[SimulationParams.sim_N_walkers][DiffusionSimulation.D];
        
        int n = SimulationParams.sim_N_walkers;
        
        Random rng= new Random(183745634);
        
        for(int i=0; i<n; i++){
            double theta= (double)(2.0*Math.PI*((double)i/(double)n));
            
            //double theta= 2.0*Math.PI*rng.nextDouble();
            
            double x= Math.cos(theta)*len;
            double y= Math.sin(theta)*len;
            
            step[i][0]=x;
            step[i][1]=y;
            step[i][2]=0.0;
            
        }
        
        diffsim.initialiseWalkers();
        Walker[] walker= diffsim.walker;
        
        
        
        
        // first update loop
        for(int i=0; i<walker.length; i++){
            
            /*for(int j=0; j<walker[i].dPhi.length; j++){
                walker[i].dPhi[j]+=scan.getPhaseShift(walker[i].r, 0.0, j, 0.0);
            }*/
            
            double[] newStep= new double[D];
            for(int j=0; j<D; j++){
                newStep[j]= step[i][j];
            }
            
            diffsim.substrate.amend(walker[i], step[i], 0.0, i, false, null);
            
            for(int j=0; j<D; j++){
                if(newStep[j]!=step[i][j]){
                    diffsim.logger.warning("step "+i+" amended (dir "+j);
                }
            }
            
            walker[i].makeStep(step[i]);
            
        }
        
        
        
        BufferedWriter phasesWriter;
        
        try{
            phasesWriter= new BufferedWriter(new FileWriter("walkerPhases.csv"));
        }
        catch(IOException ioe){
            throw new LoggedException(ioe);
        }
        
        
        //second update loop
        for(int i=0; i<walker.length; i++){
            
            for(int j=0; j<walker[i].dPhi.length; j++){
                walker[i].dPhi[j]+=scan.getPhaseShift(walker[i], diffsim.dt, j, 0.0);
            }
            
            double[] newStep= new double[D];
            for(int j=0; j<D; j++){
                newStep[j]= step[i][j];
            }
            
            diffsim.substrate.amend(walker[i], step[i], 0.0, i, false, null);
            
            for(int j=0; j<D; j++){
                if(newStep[j]!=step[i][j]){
                    diffsim.logger.warning("step "+i+" amended (dir "+j);
                }
            }
            
            walker[i].makeStep(step[i]);
            
        }
        
        double[] signal= diffsim.synthScan.getSignals();
        
        for(int i=0; i<signal.length; i++){
            
            double theta= (double)((double)i/(double)n)*2*Math.PI;
            
            System.err.println(theta+" "+signal[i]);
        }
        
    }
    
    
    public static void main(String[] args){

    	//DiffusionSimulation.testFullScale();
    	
    	//DiffusionSimulation.testReflectionFromCylinders();
    	
    	//DiffusionSimulation.testNonIntersection();
    	
    	//DiffusionSimulation.testSingleCylinder();
        
        //DiffusionSimulation.testPerfectIsotropy();
        
        DiffusionSimulation.testDiffusionSimulation();
    }
    
}

/* Walker.java created on 25-Nov-2005
 * (simulation)
 * 
 * author: matt
 * 
 */
package simulation.dynamics;


import java.io.DataOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Logger;

import misc.LoggedException;

import simulation.DiffusionSimulation;
import simulation.geometry.substrates.StickyCylinderSubstrate;
import simulation.geometry.substrates.Substrate;
import simulation.measurement.AgnosticScan;
import simulation.measurement.SyntheticScan;


/**
 *  Camino fibre reconstruction and tracking toolkit
 * 
 * Walker (simulation.dynamics)
 * 
 * describes a single brownian walker in a diffusion simulation
 * each walker has an initial position, a current position and 
 * sums of positions during each gradient block.
 * 
 * 
 * 
 *
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public class Walker {

    /** Logging object */
    private static Logger logger=Logger.getLogger("simulation.Walker");
    
    /** dimensionality of the system */
    int D;
    
    /** initial coords */
    public double r0[];
    
    /** current coords */
    public double r[]=null;
        
    /** space to store next step */
    public double[] step;
    
    /** step generator */
    private final StepGenerator stepGen;
    
    /** substrate */
    protected final Substrate substrate;
    
    /** scan object */
    private final SyntheticScan scan;
    
    /** trajectory output stream (if used) */
    protected final DataOutputStream trajWriter;
    
    /** space to store phase shifts in each scan direction */
    public final double[] dPhi;
    
    /** space to store magnetizations */
    private final double[] M;
    
    /** last time we queried the scan */
    private double tLast=0.0;
    
    /** walker size */
    public final double R;
    
    /** public constructor. takes an array of initial coords 
     * @param r0 initial coords array
     */
    public Walker(double[] r0, StepGenerator stepGen, Substrate substrate, SyntheticScan scan, DataOutputStream trajWriter){
        
        if(r0.length!=DiffusionSimulation.D){
            logger.severe("dimension mismatch in Walker initialisation");
            throw new RuntimeException("Initial coordinates passed to Walker constructor" +
            		"is not of same length ("+r0.length+") as simualtion dimensionality "+
            		DiffusionSimulation.D);
        }
        
        this.D=r0.length;
        this.r0=new double[D];
        this.r=new double[D];
        
        this.R=stepGen.getWalkerRadius();
        
        this.stepGen=stepGen;
        this.substrate=substrate;
        this.scan=scan;
        
        this.trajWriter= trajWriter;

        if(scan!=null){
        	int numMeas= scan.getNumMeasurements();
        
        	this.dPhi= new double[numMeas];
        	this.M= new double[numMeas];
        
        	for(int i=0; i<numMeas; i++){
        		dPhi[i]=0.0;
        		M[i]=0.0;
        	}
        }
        else{
        	this.dPhi=null;
        	this.M=null;
        }
        
        for(int i=0; i<D; i++){
        	this.r0[i]=r0[i];
        	this.r[i]=r0[i];
        }
    }
    
    /**
     * initialises a walker without a step generator. this is for testing
     * purposes only and WILL NOT WORK FOR A FULL SIMULATION.
     * 
     * @param r0 inital walker position
     * 
     */
    public Walker(double[] r0){

    	this(r0, null, null, null, null);

    	logger.warning("walker intantiated without step generator, substrate or scan.\nfull simulation update NOT POSSIBLE.");
    	
    }
    
    /**
     * initialises a walker without a step generator. this is for testing
     * purposes only and WILL NOT WORK FOR A FULL SIMULATION.
     * 
     * @param r0 inital walker position
     * 
     */
    public Walker(double[] r0, Substrate substrate){

        this(r0, null, substrate, null, null);

    }
    
        /**
     * initialises a walker with a step generator but no substrate. this is for testing
     * purposes only and WILL NOT WORK FOR A FULL SIMULATION.
     * 
     * @param r0 inital walker position
     * 
     */
    public Walker(double[] r0, StepGenerator stepGen){

        this(r0, stepGen, null, null, null);

    }
    
    
    
    /**
     * update walkers position.
     *
     */
    public void update(double t, int ti, int i, boolean report){
        
    	if(trajWriter!=null){
        	try{
        		trajWriter.writeDouble(t);
        		trajWriter.writeInt(i);
        		for(int j=0; j<D; j++){
        			trajWriter.writeDouble(r[j]);
        		}
        	}
        	catch(IOException ioe){
        		throw new RuntimeException(ioe);
        	}
    	}
    	
        if(scan!=null){
        	// get phase and magnetisation from scan
        	for(int j=0; j<dPhi.length; j++){
        		dPhi[j]+=scan.getPhaseShift(this, t, j, tLast);
        		dPhi[j]=AgnosticScan.mapToCircle(dPhi[j]);
        		M[j]+=substrate.getLogMagnetisationChange(this, t, tLast);
           	}
        	        	
        	tLast=t;
        }        
        
        //logger.info("Initial position: "+r[0]+", "+r[1]+", "+r[2]);
        //int when=3039;
    	//int who=35662;
    	
        boolean stepOk=false;
        int counter=0;

        FileWriter debugWriter=null;
        
        if(report){
        	logger.info("opening debug report file");
        	try{
        		debugWriter= new FileWriter("debug report "+ti+"_"+i+".csv");
        	}
        	catch(IOException ioe){
        		throw new LoggedException(ioe);
        	}
        }
        
        while(!stepOk){

            if(counter>=10){
                throw new LoggedException("unable to correctly amend step in walker update loop. ti="+ti+", i="+i);
            }

            // get a step
    	    step = stepGen.getStep(this);
    	
    	    //if(ti==when){
		//    if(i==who){
		//	    logger.info("Initial walker position: "+r[0]+", "+r[1]+", "+r[2]);
		//	    logger.info("Step: "+step[0]+", "+step[1]+", "+step[2]);
		//	    System.err.println("catch clause reached");
		//    }
	    //}
    	    
    	    if(report){
    	    	double[] pos= new double[D];
    	 
    	    	substrate.getSubstrateCoords(r, new double[]{0.0,0.0,0.0}, pos);
    	    	
    	    	logger.info("writing walker and step information");
    	    	
	    	    try{
	    	    	double[] subsSize= substrate.getSubstrateSize();

	    	    	System.err.println("pos, ,"+pos[0]+","+pos[1]+","+pos[2]+"\n");
	    	    	System.err.println("orig step, ,"+step[0]+","+step[1]+","+step[2]+"\n\n");	    	    	
	    	    	System.err.println("subs size, ,"+subsSize[0]+","+subsSize[1]+","+subsSize[2]+"\n\n");
	    	    	
	    	    	
	    	    	debugWriter.write("pos, ,"+pos[0]+","+pos[1]+","+pos[2]+"\n");
	    	    	debugWriter.write("orig step, ,"+step[0]+","+step[1]+","+step[2]+"\n\n");	    	    	
	    	    	debugWriter.write("subs size, ,"+subsSize[0]+","+subsSize[1]+","+subsSize[2]+"\n\n");
	    	    	
	    	    	debugWriter.flush();
	    	    }
	    	    catch(IOException ioe){
	    	    	throw new LoggedException(ioe);
	    	    }
    	    }
    	    
            // interact with substrate
            stepOk=substrate.amend(this, step, t, i, report, debugWriter);
            counter++;
        }
        
        // make the step 
        makeStep(step);
        
        if(report){
        	logger.info("closing debug report");
        	try{
        		debugWriter.flush();
        		debugWriter.close();
        	}
        	catch(IOException ioe){
        		throw new LoggedException(ioe);
        	}
        }
        
        //logger.info("Final position: "+r[0]+", "+r[1]+", "+r[2]);
        
    }

    
    /**
     * returns the phase shift due to diffusion weighting in 
     * a specified direction. no distinction is made between
     * weighted and unweighted directions. the phase shift in 
     * an unweighted direction will be zero.
     * 
     * @param i index of direction
     * @return dPhi[i]
     */
    public final double getPhaseShift(int i){
    	return dPhi[i];
    }
    
    /**
     * returns log magnetisaion for a given measurement
     * 
     * @param i direction index
     * @return M[i]
     */
    public final double getLogMagnetisation(int i){
    	return M[i];
    }
    
    /** make a specified step
     * @param step the step vector
     */
    public void makeStep(double[] step){
                
        for(int i=0; i<D; i++){
            this.r[i]+=step[i];
        }
                
    }
    
    /** 
     * set the position directly. This replaces the current position
     * with the given vector.
     * 
     * THIS IS TEST CODE AND SHOULD NOT BE USED IN A SIMULATION
     * 
     * @param pos new position vector (uses same instance)
     */
    public void testReplacePositionVector(double[] pos){
    	this.r=pos;
    }
    
    /** fetches the current displacemnt vector (r - r0)
     * 
     * @return new vector containing r-r0
     */
    public double[] getDisplacement(){
        double[] disp=new double[D];
        
        for(int i=0; i<D; i++){
            disp[i]=r[i]-r0[i];
        }
        
        return disp;
    }

    /** outputs walker fields as string */
    public String toString(){
        String walkerString= new String();
        walkerString+="r0= (";
        for(int i=0; i<D; i++){
            walkerString+=""+r0[i];
            if(i<D-1){
                walkerString+=", ";
            }
        }
        walkerString+="\n";
        walkerString+="r= (";
        for(int i=0; i<D; i++){
            walkerString+=""+r[i];
            if(i<D-1){
                walkerString+=", ";
            }
        }
        
        return walkerString;
    }


    
    
    /** main is just a test method */
    public static void main(String[] args) {
        int D=DiffusionSimulation.D;
        
        double r0[]= new double[] {0.5, 0.5, 0.5};
                
        Walker johnny=new Walker(r0);
        
        System.err.println("created walker: "+johnny);
        
    }


}

package simulation.measurement;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Logger;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StickyWalker;
import simulation.dynamics.Walker;

/**
 * statistics module that calculates sticky stats as well as mean-squared 
 * displacements. This generates number of free and stuck walkers as well as
 * the usual mean-squared displacements during the simulations, and then
 * also does mean occupation time and equilibrium time in post-processing.
 *
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class StickyStatsModule extends StatisticsModule {

    /** logging object */
    private final Logger logger= Logger.getLogger(this.getClass().getName());
    
    /** dimensionality of space */
    private final int D= DiffusionSimulation.D;
    
    /** more space for stick walker values */
    private final int[] stickyStats= new int[2];
    
    /** width of window for equilibration windowed average */
    private final int windowWidth= 100;
    
    /** threshold for windowed gradient of stuck people */
    private final double equilibrationThreshold= 0.1;
    
    
    /** 
     * constructor, takes array of walkers
     * 
     * @param array of walkers from simulation
     */
    public StickyStatsModule(Walker[] walker){
        
        super(walker);
        
        Ns=5;
        
        stats= new double[Ns];
        
        SimulationParams.sim_postproStatsFname= 
            new String("postproStats_"+SimulationParams.sim_p_stick+"_"+
                    SimulationParams.sim_p_unstick);
        
    }
    
    
    
    public double[] getRuntimeStats(double t){
        /* NB -- this zeros Ns values, which happens
         * currently to be the same as D but in future this
         * may change. don't forget!
         */
        for(int j=0; j<Ns; j++){
            stats[j]=0.0;
        }
        
        // accumulate mean-squared displacements
        for(int i=0; i<walker.length; i++){
            for(int j=0; j<D; j++){
                stats[j]+= (walker[i].r[j]- walker[i].r0[j])*(walker[i].r[j]- walker[i].r0[j]);
            }
        }
        
        // divide by number of walkers and time to give diffusivity
        for(int j=0; j<D; j++){
            stats[j]/=(2*walker.length*t);
        }
        
        // get the sticky stats
        int[] stickyStats= getStickyStats();
        
        //System.err.println(stickyStats[0]+" "+stickyStats[1]);
        
        if((stickyStats[0]!=0)&&(stickyStats[1]!=0)){
            stats[3]=(double)stickyStats[0];
            stats[4]=(double)stickyStats[1];
        }
        
        return stats;
        
    }
    
    
    /**
     * if the walkers are sticky, this counts how many 
     * are stuck and how many are free, returns them in
     * an array of two ints.
     * 
     * @return array of two ints {num free, num stuck}
     */
    private int[] getStickyStats(){
        
        if(walker[0] instanceof StickyWalker){
            
            for(int i=0; i<stickyStats.length; i++){
                stickyStats[i]=0;
            }
            
            for(int i=0; i<walker.length; i++){
                if(((StickyWalker)walker[i]).free){
                    stickyStats[0]++;
                }
                else{
                    stickyStats[1]++;
                }
            }
        }
        
        return stickyStats;
        
    }


    
    
    /**
     * parses binary stats file into csv for easy excel importing
     */
    public void parseBinaryStatsFile(String fname){
        
        String fnameroot;
        
        if(fname.endsWith(".bdouble")){
            fnameroot= fname.substring(0, fname.length()-8);
        }
        else{
            fnameroot= new String(fname);
            fname+=".bdouble";
        }
        
        DataInputStream binIn;
        
        FileWriter csvOut;
                
        int cols= 6;
        
        logger.info("parsing "+fname);
        
        try{
            
            binIn= new DataInputStream(new FileInputStream(fname));
            csvOut= new FileWriter(fnameroot+".csv");
            while(binIn.available()>0){
                for(int i=0; i<cols; i++){
                    double num= binIn.readDouble();
                    char extra= ',';
                    if(i==cols-1){
                        extra='\n';
                    }
                    csvOut.write(""+num+extra);
                }
            }
            binIn.close();
            csvOut.flush();
            csvOut.close();

        }
        catch(IOException ioe){
            throw new RuntimeException(ioe);
        }

        System.err.println("done.");

    }

    /**
     * gets the time to equilibration for the the number of walkers
     * stuck to the surface. It does this by taking a windowed average
     * over the time evolution of the number of stuck walkers. The 
     * equilibration time is defined as the earliest time for which the 
     * abs value of the gradient of the average of the next n 
     * samples (n const,  @see windowWidth) is less than a given theshold.
     * 
     * The equlibration level is the mean over the window at this time.
     * 
     * @return {eq. time, eq. level}
     */
    private double[] getEquilibrationTime(){
        
        String statsFileName= SimulationParams.sim_statsfile;
        String fnameroot;
        
        if(statsFileName.endsWith(".bdouble")){
            fnameroot= statsFileName.substring(0, statsFileName.length()-8);
        }
        else{
            fnameroot= new String(statsFileName);
            statsFileName+=".bdouble";
        }
        
        logger.info("seeking equilibration time. window width is "+windowWidth);
        logger.info("parsing sticky stats file '"+statsFileName+"'");
        logger.info("file name root is '"+fnameroot+"'");
        
        
        double[] windowData= new double[windowWidth];
        double[] windowTimes= new double[windowWidth];
        
        int counter=0;
        
        DataInputStream binIn;

        int cols= Ns+1;

        double[] line= new double[cols];
        double lowestGrad=Double.MAX_VALUE;
        double lowestGradTime=0.0;
        

        
        try{
            
            double lastMean=0.0;
            double mean=0.0;
            binIn= new DataInputStream(new FileInputStream(statsFileName));
            
            while(binIn.available()>0){
                for(int i=0; i<cols; i++){
                    double num= binIn.readDouble();
                    line[i]=num;
                }
                
                // accumulate data and time into the window
                windowData[counter%windowData.length]= line[5];
                windowTimes[counter%windowTimes.length]= line[0];
                
                counter++;
                
                if(counter==windowWidth){
                    lastMean=0.0;
                    
                    
                    for(int i=0; i<windowWidth; i++){
                        lastMean+=windowData[i];
                    }
                    lastMean/=windowWidth;
                }
                
                if(counter>windowWidth){
                    
                    mean=0.0;
                    
                    for(int i=0; i<windowWidth; i++){
                        mean+=windowData[i];
                    }
                    mean/=windowWidth;
                    
                    double grad=Math.abs(mean-lastMean);
                    
                    if(grad<equilibrationThreshold){
                        // if we're in here, we've reached equilibrium
                        double time= Double.MAX_VALUE;
                        for(int i=0; i<windowWidth; i++){
                            if(windowTimes[i]<time){
                                time=windowTimes[i];
                            }
                        }
                        
                        return new double[]{time, mean};
                    }
                    
                    if(grad<lowestGrad){
                        lowestGrad=grad;
                        double time= Double.MAX_VALUE;
                        for(int i=0; i<windowWidth; i++){
                            if(windowTimes[i]<time){
                                time=windowTimes[i];
                            }
                        }
                        lowestGradTime=time;
                    }
                    
                    lastMean= mean;
                }
            }
            binIn.close();

        }
        catch(IOException ioe){
            throw new RuntimeException(ioe);
        }

        logger.warning("equilibrium not reached in sticky stats file "+statsFileName+" lowest grad= "+lowestGrad+" at "+lowestGradTime);
        
        return null;
    }
    
    
    /**
     * calculates the mean occupation time of walkers.
     * 
     * TODO: include std dev occu time, maybe do histogram
     * 
     * TODO: obtain the equilibration time from the 
     * windowed average and gradient threshold on the
     * number of walkers bound to the surface
     * 
     * @return new array containing mean occ time, equilibration time and level
     *  
     */
    public double[] getPostSimulationStats(){
        
        double[] postSimStats= new double[3];
        
        double meanOccuTime= 0.0;
        
        int freeCount=0;
        
        for(int i=0; i<walker.length; i++){
            
            StickyWalker swalker= (StickyWalker)walker[i];
            
            
            
            if(swalker.free){
                // add the mean stuck time to the population average
                meanOccuTime+=(swalker.totalTimeStuck/swalker.numSticks);
                freeCount++;
            }
            else{
                
                swalker.totalTimeStuck += (SimulationParams.duration-swalker.lastStickTime);
                meanOccuTime+=(swalker.totalTimeStuck/swalker.numSticks);
            }
            
        }
    
        meanOccuTime/=walker.length;
        double[] eqStats= getEquilibrationTime();
        
        postSimStats[0]=meanOccuTime;
        
        if(eqStats!=null){
            postSimStats[1]=eqStats[0];
            postSimStats[2]=eqStats[1];
        }
        
        return postSimStats;
    }
    
    /**
     * instead of doing the full set of stats, this entrypoint parses
     * an existing binary stats file into csv and does the post-simulation
     * stats on it. it's a quick way of checking thresholds an' ting.
     * 
     * 
     * @param args
     */
    public static void main(String[] args){
        
        // set the name of the statsfile to input
        SimulationParams.sim_statsfile=new String("allstats_0.3_0.1.bdouble");
        
        // new stats mod with null walker array -- won't work for runtume stats!
        StickyStatsModule testStatsMod= new StickyStatsModule(null);
        
        // parse binary into csv
        testStatsMod.parseBinaryStatsFile(SimulationParams.sim_statsfile);
        
        // test the equilibrium stats gathering
        double[] eqStats= testStatsMod.getEquilibrationTime();
        
        System.err.println("Output from equilibration stats:");
        for(int i=0; i<eqStats.length; i++){
            System.err.println("\t"+eqStats[i]);    
        }
    }
    
    
}

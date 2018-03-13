/* SimulationParams.java created on 28-Nov-2005
 * (simulation)
 * 
 * author: Matt Hall (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation;

import imaging.DW_Scheme;
import imaging.SimulableScheme;

import java.util.logging.Logger;

import simulation.dynamics.StepAmenderFactory;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.StepAmenderFactory.AmenderType;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.geometry.elements.CylinderFactory;
import simulation.geometry.elements.CylinderFactory.CylType;
import simulation.geometry.substrates.ParallelCylinderSubstrate;
import simulation.geometry.substrates.SubstrateFactory;
import simulation.geometry.substrates.SubstrateFactory.SubstrateType;
import simulation.measurement.ScanFactory;
import simulation.measurement.ScanFactory.ScanType;
import simulation.measurement.StatisticsModuleFactory.StatsModuleType;
import tools.CL_Initializer;

/**
 *  Camino fibre reconstruction and tracking toolkit
 * 
 * SimulationParams (simulation)
 * 
 * contains extra parameters for simulations not contained in
 * imaging params.
 * 
 * 
 *
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public class SimulationParams {

    private Logger logger = Logger.getLogger(this.getClass().getName());

    
    /**
     * indicates whether to run a separate simulation for
     * each voxel (true) or generate all voxels from the 
     * same simulation
     */
    public static boolean sim_separate_runs = false;
    
    public static int sim_inflamm_increments = 10;
    
    /**
     * number of walkers in the simulation 
     */
    public static int sim_N_walkers= 10000;
    
    /**
     * number of timesteps in the diffusion simulation
     */
    public static int sim_tmax= 100000;
    
    /**
     * membrane transition probability
     */
    public static double sim_p= 0.0;
    
    /**
     * initial conditions flag
     */
    public static int sim_initial= SimulationParams.UNIFORM;
    
    /**
     * the default geometry type
     */
    public static SubstrateType sim_geomType= SubstrateFactory.SubstrateType.CYL_1_FIXED;
    
    /**
     * the substrate size
     */
    public static double sim_L=20.0;
    
    /**
     * cell size for cellular lattice
     */
    public static double sim_l=1.0;
    
    /**
     * stripe width for striped substrate
     */
    public static int sim_stripethickness=3;
    
    /**
     * percolation prob for perc lattic
     */
    public static double sim_p_perc=0.5;
    
    /**
     * the size of the central voxel
     */
    public static double sim_voxelSize= 10.0;
    
    /**
     * the type of step distibution
     */
    public static StepType sim_stepType= StepType.FIXEDLENGTH;
    
    /**
     * has the simualtion set delta?
     */
    public static boolean sim_delta_set=false;
    
    /**
     * simulation delta value
     */
    public static double sim_delta;
    
    /**
     * has the simulation set DELTA?
     */
    public static boolean sim_DELTA_set=false;
    
    /**
     * simulation diffusion time
     */
    public static double sim_DELTA;
    
       
    /**
     * has the simulation set gradient strength?
     */
    public static boolean sim_G_set=false;
    
    /**
     * simulation G value
     */
    public static double sim_G;

    /**
     * diffusivity on surface of an object for sticky walks
     * (default is one order of magnitude lower than free 
     * diffusivity)
     */
    public static double sim_surfaceDiffusivity= 2.02E-10;
    
    
    /**
     * compartmental T2 in seconds 
     * 
     * (default vals: 2.28s for bulk CSF, 10% of that on surface)
     */
    public static double[] sim_T2= new double[]{0.5, 0.09};
    
    
    /** 
     * type of cylinders
     */
    public static CylType cylinderType= CylType.BASIC;
    
    /**
     * max recursion level for nested cylinders
     */
    public static int cyl_nest_depth_max=2;
    
    /**
     * cylinder packing 
     */
    public static int sim_cyl_pack = ParallelCylinderSubstrate.HEX;

    /**
     * number of cylinders along a side of square distributed 
     * radius cylinder substrate
     */
    public static int sim_cyl_dist_size=20;
    
    /**
     * distributed cylinder min radius (meters)
     */
    public static double sim_cyl_min_r=0.0;
    
    /**
     * distributed cylinder max radius (meters)
     */
    public static double sim_cyl_max_r=2E-5;
    
    /** 
     * fixed radius of cylinder or sphere (meters)
     */
    public static double sim_r = 1E-5;
    
    /**
     * fixed separation of cylinder or sphere (meters)
     */
    public static double sim_R = 3E-5;
    
    /**
     * separation of mesh objects in 3D.
     * if null, assume minimal bounding box
     */
    public static double[] sim_mesh_sep= null;
    
    /**
     * ratio of outer to inner ratio
     */
    public static double sim_gRatio=0.7;
    
    /** 
     * outer radius of myelinated cylinders
     */
    public static double sim_cyl_r1;
    

    
    /**
     * core diffusivity of myelinated cylinders
     */
    public static double sim_cyl_D1= CL_Initializer.DIFF_CONST;
    
    /**
     * sheath diffusivity in meylinated cyliners
     */
    public static double sim_cyl_D2= CL_Initializer.DIFF_CONST;
    
    /**
     * angle between crossing fibres on crossing substrate. 
     * 90 degrees by default.
     */
    public static double sim_cAngle= Math.PI/2;
    
    /**
     * number of cylinders on substrate
     */
    public static int sim_num_cylinders=20;
    
    /**
     * number of facets on a facetted cylinder
     */
    public static int sim_num_facets=0;
    
    /**
     * name of PLY file to read
     */
    public static String sim_plyfile=null;
    
    /**
     * name of cylinders file to write
     */
    public static String sim_cylFile=null;
    
    /**
     * fraction of outer voxel voxel on mesh substrate used 
     * to generate data
     */
    public static double sim_voxelSizeFrac= 0.75;
    
    /**
     * name of stats output file
     */
    public static String sim_statsfile= null;
    
    /**
     * type of statistics to generate
     */
    public static StatsModuleType sim_StatsModType= StatsModuleType.MS_DISP;
    
    /**
     * density of square grid used for sptial sorting
     */
    public static int sim_spatial_grid_size= 10;
    
    /**
     * flag to set drawing of cross section of cylinder substrate (TODO: squashy cyls only)
     */
    public static boolean sim_drawCrossSection= false;
    
    /**
     * Name of file to output the cross section figure
     */
    public static String crossSectionFig= "crossSec_dyn.gray";
    
    /**
     * size of cross section image in pixels
     */
    public static int sim_crossSectionImageSize= 1000;
    
    /**
     * type of step amender
     */
    public static AmenderType sim_amender_type=AmenderType.ELESTIC_REFLECTOR;
    
    /**
     * probability for sticking to a membrane 
     */
    public static double sim_p_stick= 0.75;
    
    /**
     * probability for a spin unsticking itself from a membrane
     */
    public static double sim_p_unstick= 0.0;
    
    /**
     * name of post-processing stats file
     */
    public static String sim_postproStatsFname= new String("postProcessingStats.bdouble"); 
    
    /**
     * scan type
     */
    public static ScanType scanType= ScanType.PGSE_SCAN;
    
    /** flag for delta-peaked initially posiitoned walkers */
    public static final int SPIKE=0;
    
    /** flag for uniformly distributed initial conditions */
    public static final int UNIFORM=1;
    
    /** flag for special (debug) initial conditions. one intra one extra */
    public static final int SPECIAL=2;
    
    /** flag for intracellular walkers only */
    public static final int INTRACELLULAR=3;
    
    /** flag for extracellular walkers only */
    public static final int EXTRACELLULAR=4;
    

    
    /** flag for compartmental signal: default (no compartmental segregation */
    public static final int NOCOMPS= 0;
    
    /** flag for compartmental signal: extracellular only */
    public static final int EXTRAONLY= 1;
    
    /** flag for compartmental signal: intracellular only */
    public static final int INTRAONLY= 2;
    
    /** flag for compartmental signal: output for each component */
    public static final int ALLCOMPS= 3;
    
    
    
    
    
    
    /** number of walkers */
    private final int N_walkers;
    
    /** number of timesteps (exclusing transient) */
    private final int tmax;
    
    /** the increment of time corresponding to a timestep */ 
    private final double dt;   
    
    /** the type of step distribution we'd like */
    private final StepType stepType;
    
    /** additional step parameters */
    private double[] stepParams;
       
    /** the size of the imaging voxel */
    private final double voxelSize;
    
    /** the membrane transition probability */
    private final double p;
    
    /** initial conditions flag */
    private final int initial;
    
    /** geometry type */
    private final SubstrateType geometryType;
    
    /** extra geometry parameters */
    //private final Object[] geometryParams;
    
    /** name of trajectories file */
    public static String trajFile=null;
    
    /** size of traj file output buffer (one meg by default)*/
    public int buffsize=1048576;
    
    /** duration of simulation in seconds */
    public static double duration;
    
    /** are we generating trajectories? */
    public static boolean trajectories=false;
    
    /** flag to read out cylinder info */
    public static boolean substrateInfo= false;
    
    /** only run a single iteration of the substrate? -1 if no. other number to run */
    public static int sim_onlyRun= -1;
    
    /** compartmental signal output? default no */
    public static int sim_compartmentSignal= NOCOMPS;
    
    
    public SimulationParams(int N_walkers, int tmax, double p, 
            				int initial, SubstrateType geomType, //Object[] geomParams, 
            				StepType stepType,
            				double voxelSize, SimulableScheme imPars){

        this.N_walkers=N_walkers;
        this.tmax=tmax;
        this.p=p;
        this.geometryType=geomType;
        this.stepType=stepType;
        this.voxelSize=voxelSize;
                
        /** set the echo time and dt increment properly */
        double TE=((SimulableScheme)imPars).getDuration();
        
        //double halfP90= imPars.getHalfP90(0);
        
        //this.dt=(halfP90+TE)/((double)tmax);
        
        this.dt=TE/((double)tmax);
        
        trajectories=false;
        
        this.initial=initial;
    }

    
    /** alternate constructor with no scheme object. this instantiates
     * a simulation that generates a trajectories file instead of scan
     * results.
     * 
     * @param N_walkers number of walkers
     * @param tmax number of timesteps
     * @param p permeability
     * @param initial initial conditions flag
     * @param geomType substrate type
     * @param stepType type of steps to generate
     * @param voxelSize size of voxel
     * @param duration duration of simulated dynamics in seconds
     */
    public SimulationParams(int N_walkers, int tmax, double p, 
			int initial, SubstrateType geomType, //Object[] geomParams, 
			StepType stepType,
			double voxelSize, double duration){

    	this.N_walkers=N_walkers;
		this.tmax=tmax;
		this.p=p;
		this.geometryType=geomType;
		//this.geometryParams=geomParams;
		this.stepType=stepType;
		this.voxelSize=voxelSize;
		
		//this.trajFile= new String("N="+N_walkers+"_T="+tmax+"_m"+CL_Initializer.seed+".traj");
		
		/** set the echo time and dt increment properly */
		this.dt=duration/((double)tmax);
		
		this.duration= duration;
		this.trajectories=true;
		
		this.initial=initial;
}

    
    
    
    
    
    
    /**
     * @return Returns the dt.
     */
    public double getDt() {
        return dt;
    }

    /**
     * @return Returns the n_walkers.
     */
    public int getN_walkers() {
        return N_walkers;
    }

    /**
     * @return Returns the tmax.
     */
    public int getTmax() {
        return tmax;
    }
    public static void main(String[] args) {
    }
    
    /**
     * @return Returns the membrane transition prob p.
     */
    public double getP() {
        return p;
    }
    
    /**
     * @return Returns the geometryType.
     */
    public SubstrateType getGeometryType() {
        return geometryType;
    }
    
    /**
     * @return Returns the geometryParams.
     */
    /*public Object[] getGeometryParams() {
        return geometryParams;
    }*/
    
    /**
     * @return Returns the initial.
     */
    public int getInitialConditions() {
        return initial;
    }
    /**
     * @return Returns the stepParams.
     */
    public double[] getStepParams() {
        return stepParams;
    }
    /** sets the parameters for the step generator
     * 
     * @param steParams array of step generator parameters
     */
    public void setStepParams(double[] stepParams){
        this.stepParams=stepParams;
    }
    /**
     * @return Returns the stepType.
     */
    public StepType getStepType() {
        return stepType;
    }
    /**
     * @return Returns the voxelSize.
     */
    public double getVoxelSize() {
        return voxelSize;
    }
    
}

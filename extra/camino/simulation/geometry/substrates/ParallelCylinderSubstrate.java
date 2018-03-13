package simulation.geometry.substrates;

import imaging.RectGradSteTanScheme;
import imaging.SimulableScheme;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.logging.Logger;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.Walker;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.geometry.elements.BasicCylinder;
import simulation.geometry.elements.Cylinder;
import simulation.geometry.elements.CylinderFactory;
import simulation.geometry.elements.MyelinatedCylinder;
import simulation.geometry.elements.CylinderFactory.CylType;
import simulation.geometry.elements.NestedCylinder;
import simulation.geometry.substrates.SubstrateFactory.SubstrateType;
import tools.CL_Initializer;
import misc.LoggedException;


/** 
 *  Implements a class of substrates with a single population
 * of parallel cylinders with a single constant thickness. Both
 * square and hexagonal packing are implemented. Square is easy,
 * with just a single cylinder in the cell and a square unit cell 
 * or size R (R = cylinder separation)
 * 
 *    -----
 *   |     |
 *   |  O  |
 *   |     |
 *    -----
 *   
 *   In hexagonal packing the situation is slightly (but not greatly)
 *  more complicated. here we have a rectangular unit cell and four 
 *  cylinders arranged in a diamond configuration:
 *  
 *    -o-
 *   |   |
 *   o   o
 *   |   |
 *    -o-
 *   
 *   where cell width is cylinder separation R and cell height is
 *   is 2sqrt(3/4)R. An equivalent solution has the cylinders in a
 *   domino 5 cross formation, but the diamond formulation only 
 *   requires 4 cylinders instead of five.
 *   
 *    In both cases cylinders are orientedd with their axes parallel
 *   to the z axis, as this means that the euclidean coordinates of 
 *   walker position do not have to be rotated into the cylindrical
 *   frame and thus greatly simplifies the geometrical calculations.
 *   
 * @author matt m.hall@cs.ucl.ac.uk
 *
 */
public class ParallelCylinderSubstrate extends CylinderSubstrate{
	
	/** logging object */
	Logger logger = Logger.getLogger(this.getClass().getName());
	
	/** dimensionality of space */
	private final int D=DiffusionSimulation.D;
	
	/** spacing of cylinders */
	private final double R;
	
	/** size of repeating cell */
	private final double[] l;
	
	/** cylinder orientation */
	private final double[] V;
	
	/** space to store substrate coordinates */
	private final double[] subsCoords= new double[D];
	
	/** array of three zeros */
	private final double[] noOffset= new double[]{0.0, 0.0, 0.0}; 
	
	// labels for packing configs
	public static final int SQUARE = 1;    // square packing
	public static final int HEX = 2;       // hexagonal packing
	public static final int DISTRIBUTED=3; // distributed radii (handled in DistributedRadiusCylinderSubstrate.java)
	
	/** packing in individual instance */
	private final int packing;

	/** type of cylinder on the substrate */
	private final CylType cylType;
	
	/**
	 * construct a lattice of parallel simple cylinders of
	 * given radius, orientation and spacing
	 * 
	 * @param R cylinder spacing
	 * @param r radius of cylinders
	 * @param packing cylinder packing (HEX or SQUARE)
	 */
	public ParallelCylinderSubstrate(SimulationParams simParams){
		
		super(new double[]{2*SimulationParams.sim_R, 
		                   2*SimulationParams.sim_R, 
		                   2*SimulationParams.sim_R}, simParams, (SimulationParams.cylinderType==CylType.NESTED));
		
		double p= simParams.getP();
		
		Cylinder[] cylinder;
		double[] P= new double[D];
		
		if(SimulationParams.sim_R<2.0*SimulationParams.sim_r){
			logger.warning("Cylinder spacing on parallel cylinder lattice (R="
					+SimulationParams.sim_R+") is less than cylinder radius (r="
					+SimulationParams.sim_r+")");
		}
		

		this.R=SimulationParams.sim_R;
	        
		double r=SimulationParams.sim_r;
	        
		this.V= new double[D];
		for(int i=0; i<D; i++){
			P[i]=R/2.0;
			this.V[i]=0.0;
		}
		V[D-1]=1.0;
		
		this.packing=SimulationParams.sim_cyl_pack;
		
		//this.cylType= CylType.BASIC;
		this.cylType= SimulationParams.cylinderType;
		
		if(packing==SQUARE){
		
			double cylArea=Math.PI*r*r;
			double sqArea= R*R;
			
			double V_I= cylArea/sqArea;
			
			logger.info("Constructing square-packed parallel cylinder substrate. radius= "+r+" separation= "+R);			
			logger.info("intracellular vol frac = "+V_I);

			
			
			/* square packing needs only one cylinder
			 * but boundaries require that we put an 
			 * extra 8 to form the neighbouring cells
			 * 
			 *  o o o
			 *  o o o
			 *  o o o
			 *  
			 */
			cylinder= new Cylinder[9];

			for(int i=0; i<3; i++){
				P[0]= (i)*R;
				
				for(int j=0; j<3; j++){
				    P[1]= (j)*R;
				    
					//cylinder[3*i+j]=new BasicCylinder(P, r, p);
				    cylinder[3*i+j]= CylinderFactory.getCylinder(P,  r,  p,  cylType);
				}
			}
			
			
			// and cell size is the same as cylinder spacing
			this.l= new double[]{2*R, 2*R, 2*R};
		}
		else if(packing==HEX){
			
			int counter=0;
			
			/**
			 * hexagonal packing needs four cylinders:
			 *   o
			 *  o o
			 *   o 
			 * but periodic boundaries mean that we need
			 * an additional four:
			 * o o o
			 *  o o
			 * o o o
			 */
			cylinder= new Cylinder[8];
			
			l= new double[D];
			
			double k=Math.sqrt(3.0/4.0);
			
			// cell size is sqrt(3/4)R
			l[0]=2*R;
			l[1]=2*k*R;
			l[2]=2*R;
			
			super.setSubstrateDims(l);
			
			double V_I= (Math.PI*r*r)/(k*R*R);
			
			logger.info("constructing hexagonally packed cylinder substrate, radius= "+r+" separation= "+R);
			logger.info("intracellular vol frac = "+V_I);
	
			// (0,0)
			P[0]=0.0;
			P[1]=0.0;
			cylinder[counter++]= CylinderFactory.getCylinder(P,  r,  p,  cylType);//new BasicCylinder(P, r, p);
            
			// (R,0)
			P[0]=R;
			cylinder[counter++]= CylinderFactory.getCylinder(P,  r,  p,  cylType);//new BasicCylinder(P, r, p);
            
			// (2R,0)
			P[0]=2*R;
			cylinder[counter++]= CylinderFactory.getCylinder(P,  r,  p,  cylType);//new BasicCylinder(P, r, p);
            
			// (R/2,kR)
			P[0]=R/2;
			P[1]=k*R;
            cylinder[counter++]= CylinderFactory.getCylinder(P,  r,  p,  cylType);//new BasicCylinder(P, r, p);
			
            // (3R/2,kR)
			P[0]=3*R/2;
            cylinder[counter++]= CylinderFactory.getCylinder(P,  r,  p,  cylType);//new BasicCylinder(P, r, p);
			
            // (0,2kR)
            P[0]=0.0;
            P[1]=2*k*R;
            cylinder[counter++]= CylinderFactory.getCylinder(P,  r,  p,  cylType);//new BasicCylinder(P, r, p);

            // (R,2kR)
            P[0]=R;
            cylinder[counter++]= CylinderFactory.getCylinder(P,  r,  p,  cylType);//new BasicCylinder(P, r, p);

            // (2R,2kR)
            P[0]=2*R;
            cylinder[counter++]= CylinderFactory.getCylinder(P,  r,  p,  cylType);//new BasicCylinder(P, r, p);
	
		}
		else{
			cylinder=null;
			throw new LoggedException("packing label = "+packing+" not understood, please use square ("+SQUARE+") or hexagonal ("+HEX+")");
		}
		
		
		if(cylType==CylType.NESTED){
			ArrayList<Cylinder> allCyls=new ArrayList<Cylinder>();
			
			for(int i=0; i< cylinder.length; i++){
				allCyls.addAll(((NestedCylinder)cylinder[i]).allCylinders());
			}
			
			cylinder= new Cylinder[allCyls.size()];
			
			for(int i=0; i<cylinder.length; i++){
				cylinder[i]= allCyls.get(i);
			}
		}
		
		
		setCylinders(cylinder);
		
		int[] n= new int[]{300, 300,1};
				
		initialiseSpatialOptimisation(n);
		
	}
	
	/**
	 * constructor for thick-walled cylinders. exactly the same as the
	 * above but with more parameters for the myelinated cylinders.
	 * 
	 * @param cylType which type of cylinder to use
	 * @param simParams simulation parameters object
	 */
	public ParallelCylinderSubstrate(CylType type, SimulationParams simParams){
		
		super(new double[]{4*SimulationParams.sim_L, 
		                   4*SimulationParams.sim_L, 
		                   4*SimulationParams.sim_L}, simParams, false);
		
		double R= SimulationParams.sim_R;
		double r1= SimulationParams.sim_cyl_r1;
		double r2= SimulationParams.sim_r;
		double D1= SimulationParams.sim_cyl_D1;
		double D2= SimulationParams.sim_cyl_D2;
		double p= SimulationParams.sim_p;
		int packing= SimulationParams.sim_cyl_pack;
		
		
		
		Cylinder[] cylinder;
		double[] P= new double[D];
		
		if(R<2.0*r2){
			logger.warning("Cylinder spacing on parallel cylinder lattice (R="
					+R+") is less than cylinder radius (r="+r2+")");
		}
		
		
		this.V= new double[D];
		for(int i=0; i<D; i++){
			P[i]=R/2.0;
			this.V[i]=0.0;
		}
		V[D-1]=1.0;
		
		this.R=R;
		
		this.packing=packing;
		
		this.cylType=type;
		
		if(packing==SQUARE){
		
			double cylArea=Math.PI*r2*r2;
			double sqArea= R*R;
			
			double V_I= cylArea/sqArea;
			
			logger.info("Constructing square-packed parallel myelinated cylinder substrate. radius= "+r2+" separation= "+R);			
			logger.info("intracellular vol frac = "+V_I);

			
			
			/* square packing needs only one cylinder
			 * but boundaries require that we put an 
			 * extra 8 to form the neighbouring cells
			 * 
			 *  o o o
			 *  o o o
			 *  o o o
			 *  
			 */
			cylinder= new Cylinder[9];

			for(int i=0; i<3; i++){
				//P[0]= (i-1)*R;
			    P[0]= (i)*R;
			    
				for(int j=0; j<3; j++){
					//P[1]= (j-1)*R;
				    P[1]= (j)*R;
				    
					//cylinder[3*i+j]=new MyelinatedCylinder(P, r1, r2, D1, D2, p);
				    cylinder[3*i+j]= CylinderFactory.getCylinder(P, simParams);
				}
			}
			
			
			// and cell size is the same as cylinder spacing
			this.l=new double[]{R, R, R};
		}
		else if(packing==HEX){
			
			int counter=0;
			
			/**
			 * hexagonal packing needs four cylinders:
			 *   o
			 *  o o
			 *   o 
			 * but periodic boundaries mean that we need
			 * an additional four:
			 * o o o
			 *  o o
			 * o o o
			 */
			cylinder= new Cylinder[8];
			
			double k=Math.sqrt(3.0/4.0)*R;
			// cell size is 2*sqrt(3/4)R
			l= new double[D];
			l[0]=2*R;
			l[1]=2*k;
			l[2]=2*R;
			
			super.setSubstrateDims(l);
			
			double V_I= (Math.PI*r2*r2)/(l[1]*R);
			
			logger.info("constructing hexagonally packed myelinated cylinder substrate, radius= "+r2+" separation= "+R);
			logger.info("intracellular vol frac = "+V_I);

			P[0]=0.0;
            P[1]=0.0;
            //cylinder[counter++]= new MyelinatedCylinder(P, r1, r2, D1, D2, p);
            cylinder[counter++]= CylinderFactory.getCylinder(P, simParams);
            
            P[0]=R;
            //cylinder[counter++]= new MyelinatedCylinder(P, r1, r2, D1, D2, p);
            cylinder[counter++]= CylinderFactory.getCylinder(P, simParams);
            
            P[0]=2*R;
            //cylinder[counter++]= new MyelinatedCylinder(P, r1, r2, D1, D2, p);
            cylinder[counter++]= CylinderFactory.getCylinder(P, simParams);
            
            P[0]=R/2;
            P[1]=k;
            //cylinder[counter++]= new MyelinatedCylinder(P, r1, r2, D1, D2, p);
            cylinder[counter++]= CylinderFactory.getCylinder(P, simParams);
            
            P[0]=3*R/2;
            //cylinder[counter++]= new MyelinatedCylinder(P, r1, r2, D1, D2, p);
            cylinder[counter++]= CylinderFactory.getCylinder(P, simParams);
            
            P[0]=0.0;
            P[1]=2*k;
            //cylinder[counter++]= new MyelinatedCylinder(P, r1, r2, D1, D2, p);
            cylinder[counter++]= CylinderFactory.getCylinder(P, simParams);
            
            P[0]=R;
            //cylinder[counter++]= new MyelinatedCylinder(P, r1, r2, D1, D2, p);
            cylinder[counter++]= CylinderFactory.getCylinder(P, simParams);
            
            P[0]=2*R;
            //cylinder[counter++]= new MyelinatedCylinder(P, r1, r2, D1, D2, p);
            cylinder[counter++]= CylinderFactory.getCylinder(P, simParams);
			
		}
		else{
			cylinder=null;
			throw new LoggedException("packing label = "+packing+" not understood, please use square ("+SQUARE+") or hexagonal ("+HEX+")");
		}
		
		setCylinders(cylinder);
		
	}
	
	
	/**
	 * @return size of a unit cell. for hex packing, returns the long axis length
	 */
	public double[] getSubstrateSize() {
		return l;
	}

	/**
	 * @return  coord of the centre of a unit cell
	 */
	public double getPeakCoord() {
		return l[0]/2.0;
	}

	/** 
	 * overrides getDiffusivityAt method in Substrate
	 * so that we can do spatial variation in diffusivity.
	 * 
	 * @param walkerPos position to check in world coords
	 * 
	 * @return diffusivity at that location
	 */
	public final double getDiffusivityAt(double[] walkerPos){
		
		double[] subsCoords= new double[D]
		                                ;
		getSubstrateCoords(walkerPos, noOffset, subsCoords);
			
		for(int i=0; i<cylinder.length; i++){
			if(cylinder[i].inside(subsCoords)){
				return cylinder[i].getDiffusivityAt(subsCoords);
			}
		}
		
		return CL_Initializer.DIFF_CONST;
	}
	
	// extras
	/** 
	 * gets the cylinder spacing
	 * same as short axis length in hex packing cell. 
	 * 
	 * @return spacing between cylinders
	 */
	public double getCylinderSpacing(){
		return R;
	}

	/**
	 * gets the packing type
	 * 
	 * @return square or hexagonal
	 */
	public int getPacking(){
		return packing;
	}
	
	public String getPackingString(){
		if(packing==SQUARE){
			return new String("square");
		}
		else if (packing==HEX){
			return new String("hex");
		}
		else{
			return new String("unknown! code ="+packing);
		}
	}

	/** 
	 * transforms the walker coords into the cell
	 * 
	 * @param r position in non-periodic space
	 * @param offset the offset vector from walker position
	 * 
	 * @return position in cell (periodic) space
	 */
	public void getSubstrateCoords(double[] r, double[] offset, double[] subsCoords){
			
		for(int i=0; i<D; i++){
			double pos= r[i]+offset[i];
			double windingNum= Math.floor(pos/l[i]);
			subsCoords[i]= pos-windingNum*l[i];
		}
		
	}
	
	/** test code for barrier crossing with cylinders
	 * 
	 * instantiates a square lattice with R=2.0, r=0.5 and a 
	 * walker at r=(1.4,1) with a step of s=(0.2,0) which
	 * should take it across the barrier at (1.5, 1.0) 
	 * giving t=0.1
	 * 
	 * 
	 * @param args
	 */
	public static void main(String[] args){
		
		int D= DiffusionSimulation.D;
		
		int tmax=10;
		
		double[] r0 = new double[]{1.4, 1.0, 0.0};
		Walker walker = new Walker(r0);
		double[] s= new double[]{0.2, 0.0, 0.0};
		
		double[] normal = new double[D];
		double[] d = new double[1];
		
	      SimulableScheme scheme;
	      try{
	          URI uri= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
	          
	          String path= uri.getPath();
	          
	          scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
	      }  
	      catch(URISyntaxException urise){
	          throw new LoggedException(urise);
	      }

		SimulationParams simParams= new SimulationParams(1, 1000, 1.0, 1, 
				SubstrateFactory.SubstrateType.CYL_1_FIXED, StepType.FIXEDLENGTH,
				1.0, scheme);
		
		
		System.err.println("testing geometry directly");
		// first off, test the parallel cylinder substrate directly
		ParallelCylinderSubstrate pcs= new ParallelCylinderSubstrate(simParams);
		
		System.err.println("substrate type is "+pcs.getPackingString());
		System.err.println("cylinder spacing is "+pcs.getCylinderSpacing());
		
		boolean[] in= new boolean[1];
		double[] offset= new double[]{0.0, 0.0, 0.0};
		double[] p= new double[1];
		
		boolean crosses=false;
		try {
			crosses = pcs.crossesMembrane(walker, offset, s, normal, d, false, 0.2, in, p, false, null);
		} catch (StepRejectedException e) {
			e.printStackTrace();
		}
		System.err.println("crosses is "+crosses+" normal=("+normal[0]+","+normal[1]+","+normal[2]+") d="+d[0]);
		
		System.err.println("\ntesting geometry in substrate framework");

		// now instantiate a substrate and test step amendment
		Substrate substrate= SubstrateFactory.getSubstrate(SubstrateFactory.SubstrateType.CYL_1_INFLAM, simParams);

		
		try {
			substrate.crossesMembrane(walker, offset, s, normal, d, false, 0.2, in, p, false, null);
		} catch (StepRejectedException e) {
			e.printStackTrace();
		}
		System.err.println("crosses is "+crosses+" normal=("+normal[0]+","+normal[1]+","+normal[2]+") d="+d[0]);

		double[] toBarrier= new double[D];
		double[] reflected= new double[D];
		double[] unreflected= new double[D];
		
		substrate.testAmendment(walker, s, normal, d, 0.2, toBarrier, reflected, unreflected);
		
		System.err.println("toBarrier= ("+toBarrier[0]+","+toBarrier[1]+","+toBarrier[2]+")");
		System.err.println("reflected= ("+reflected[0]+","+reflected[1]+","+reflected[2]+")");
		System.err.println("unreflected= ("+unreflected[0]+","+unreflected[1]+","+unreflected[2]+")");
		
	}
	
}

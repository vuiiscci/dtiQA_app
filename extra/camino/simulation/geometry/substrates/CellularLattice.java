/* CelllarLattice.java created on 28-Nov-2005
 * (simulation)
 * 
 * author: Matt Hall (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation.geometry.substrates;


import imaging.RectGradSteTanScheme;
import imaging.SimulableScheme;

import java.io.FileWriter;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Logger;

import misc.LoggedException;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.Walker;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.geometry.elements.Triangle;
import simulation.geometry.substrates.SubstrateFactory.SubstrateType;
import tools.CL_Initializer;

/**
 *
 * Abstract class that implements a general cubic lattice geometry
 * initialisation of the geometry is left abstract so that inheriting 
 * classes may implement any pattern they choose.
 * 
 * Cells on the lattice have a well-defined spatial extend (they are 
 * cubes) and if defined as occupied (or active. these terms are used 
 * synonymously) then the cube has membranes on each side. A membrane
 * allows a walker to pass through it with probability p or elastically
 * reflects it with prob 1-p.
 * 
 * This lattice has periodic boundary conditions.
 * 
 * @author Matt Hall (m.hall@cs.ucl.ac.uk)
 *
 */
public abstract class CellularLattice extends Substrate{

    /** logging object */
    private Logger logger= Logger.getLogger(this.getClass().getName());
 
    /** the cell size */
    protected final double l;
    
    /** the lattice size */
    protected final int L;
    
    /** dimension of the substrate */
    protected final int D=DiffusionSimulation.D;
    
    /** the occupation flag array */
    protected boolean[] occupied=null;
    
    /** the length of the occupation array */
    protected int occupiedLength;
    
    /** tolerence for floating point equality */
    public static final double TINYNUM= 5E-14;
    
    /** membrane permeability */
    private final double p;
    
    public CellularLattice(SimulationParams simParams){
        
    	super(simParams, new double[]{SimulationParams.sim_l*SimulationParams.sim_L,
    	                              SimulationParams.sim_l*SimulationParams.sim_L,
    	                              SimulationParams.sim_l*SimulationParams.sim_L});
    	
        this.l=SimulationParams.sim_l;
        this.L=(int)SimulationParams.sim_L;
        
        this.p=simParams.getP();
        
        occupiedLength=1;
        for(int i=0; i<D; i++){
            occupiedLength*=L;
        }
        
        occupied=new boolean[occupiedLength];
    }
    
    /** abstract method for initialising cell activation. active
     * cells have membranes, whereas inactive cells do not. Thus
     * many of the substrate properties are contained in the pattern
     * of lattice activation.
     * 
     * should be called from the inheriting class's constructor.
     *
     */
    public abstract void initLattice();
    
    
    /** returns the linear sequence index of a set of D cell
     *  This method appkies boundaries internally but DOES NOT
     *  CHANGE THE VALUES IN THE CELL ARRAY. This is to prevent 
     *  all the boundary-distance nonsense later on.
     *  
     * 
     *  @param cell vector of cell indeices in D dimensions
     *  @return index in 1D sequence
     */    
    private int getOccupiedArrayIndex(int[] cell){
        
        int index=0;
        
        // 1D case is easy
        if(D==1){
            int cellNum=cell[0];
            
            if(cellNum<0){
                cellNum+=L;
            }
            else if(cellNum>=L){
                cellNum-=L;
            }
            
            return cell[0];
        }
        
        // nD case -- apply boundaries and get occupation
        for(int i=0; i<D; i++){
            int cellNum=cell[i]%L;

            if(cellNum<0){
                cellNum+=L;
            }
            else if(cellNum>=L){
                cellNum-=L;
            }
            
            index+=cellNum*((int)Math.pow(L, i));
        }
        
        return index;
    }

    /** gets the cell coords from the spatia location
     * 
     * @param r spatial position
     * 
     * @return array of cell coords
     */
    public int[] getCell(double[] r){
        
        int[] cell=new int[D];
        //int[] windingNumber=getWindingNumbers(r);
        
        //double latticeSize=getSubstrateSize();
        
        for(int i=0; i<D; i++){
        	//double offset = windingNumber[i]*latticeSize;
            //cell[i]=(int)Math.floor((r[i]-offset)/l);
            cell[i]=(int)Math.floor(r[i]/l);
        }
        
        return cell;
    }
    
    /** check if given cell is occupied on the lattice
     * 
     * @param cell cell coords
     * @return true or false;
     */
    public boolean getCellOccupation(int[] cell){
        int index= getOccupiedArrayIndex(cell);
        
        return occupied[index];
    }
    
    /** checks if a step will take a walker across a membrane or not
     *  by comparing the cell coords before and after the step.
     * 
     *  if they're the same, then no crossing otherwise we've crossed.
     * 
     * @param walker the walker making the step
     * @param offset displacement from walker position to check against 
     * @param stepVector the step vector
     * @param normal the normal of the surface
     * @param d space to store disance of barrier from origin 
     * @param skipCurrent flag to skip barrier we're sitting on top of
     * 
     * @see simulation.geometry.Geometry#crossesMembrane(simulation.dynamics.Walker, double[], double[])
     */
    public boolean crossesMembrane(Walker walker, double[] offset, double[] stepVector, 
    						double[] normal, double[] d, boolean skipCurrent, 
    						double origLength, boolean[] in, double[] p, boolean report, FileWriter debugWriter) {
        
        boolean crosses=false;
                
        double[] stepPos= new double[D];
        double modStep=0.0;
        
        double[] skippingNormal=null;
        double skippingDistance=0.0;

        double[] walkerPos= new double[D];
        
        for(int i=0; i<D; i++){
        	walkerPos[i]=walker.r[i]+offset[i];
            stepPos[i]=walkerPos[i]+stepVector[i];
            modStep+=stepVector[i]*stepVector[i];
        }
        
        
        modStep=Math.sqrt(modStep);
        
        if(modStep/origLength <= TINYNUM){
        	return false;
        }

        int[] currentCell=getCell(walker.r);
        int[] nextCell=getCell(stepPos);
        
        int numDiffs=0;
        
        for(int i=0; i<D; i++){
            if(currentCell[i]!=nextCell[i]){               
                crosses=true;
                numDiffs++;
            }
        }
        
        
        
        

        if(crosses){
            boolean currentCellOcc=getCellOccupation(currentCell);
            boolean newCellOcc=getCellOccupation(nextCell);
            
            
            if(currentCellOcc||newCellOcc){
                int index=0;

                if(skipCurrent){
                	skippingDistance= d[0];

                	skippingNormal= new double[D];
                	for(int i=0; i<D; i++){
                		skippingNormal[i]=normal[i];
                	}
                }
                
                if(numDiffs==1){
                	// if number of differences is 1, get the normal
                	index=calcBarrierNormal(currentCell, nextCell, normal);
                	
                }
                else{
                	// otherwise more than one barrier is crossed
                	// so we need to work out which one is closest
                	int[] trialNextCell=new int[D];
                	double[] trialNormal= new double[D];
                	double trialDist=0.0;
                	boolean trialCellOcc;
                	
                	int closest=-1;
                	double minDist=Double.MAX_VALUE;
                	
                	// try all the possible barriers and find the closest
                	for(int i=0; i<D; i++){
                		if(currentCell[i]!=nextCell[i]){
                			double normalDotPos=0.0;
                			double normalDotStep=0.0;

                			for(int j=0; j<D; j++){
                				
                				if(j==i){
                					trialNextCell[j]=nextCell[j];
                					trialNormal[j]=1.0;
                					normalDotPos=walker.r[j];
                					normalDotStep=stepVector[j];
                				}
                				else{
                					trialNextCell[j]=currentCell[j];
                					trialNormal[j]=0.0;
                				}
                			}
                			trialCellOcc=getCellOccupation(trialNextCell);
                			
                			// if the current cell and trial cell are both
                			// unoccupied then there's no barrier!
                			if(!(currentCellOcc||trialCellOcc)){
                				continue;
                			}
                			
                			double planeFromOrig;
                			if(trialNextCell[i]>currentCell[i]){
                				planeFromOrig= trialNextCell[i]*l;
                			}
                			else{
                				planeFromOrig= currentCell[i]*l;
                			}
                			
                			double distFromWalker=normalDotPos - planeFromOrig;
                			double cosTheta = normalDotStep/modStep;
                			
                			double distToBarrier=distFromWalker/cosTheta;

                        	if(skipCurrent){  // this means we're sitting on a barrier
                       		
                        		double relStepLengthToBarrier = distToBarrier/origLength;
                        		
                        		// check if we are very close to current barrier
                        		// and that this barrier is the right distance from origin
                        		if((Math.abs(relStepLengthToBarrier)<=TINYNUM)&&(Math.abs(d[0]-skippingDistance)<=TINYNUM)){
                        			boolean skipping=true;
                        			
                        			// if so, check the barrier normal against the skipping normal
                        			for(int j=0; j<D; j++){
                        				if(Math.abs(normal[j]-skippingNormal[j])>TINYNUM){
                        					skipping=false;
                        					break;
                        				}
                        			}
                        			
                        			// if skipping this barrier, skip on...
                        			if(skipping){
                        				continue;
                        			}
                        		}                       		
                        	}
                        	
                        	
                        	if(distFromWalker<minDist){
                				minDist=distToBarrier;
                				closest=i;
                			}
                		}
                	}
                	
                	// this catches the case where none of the possible trials leads to a crossing
                	// e.g. crossing two cells, one is ignored, other is between unoccupied cells
                	if(closest==-1){
                		return false;
                	}
                	
                	// construct normal and distance from orig for closest plane
                	for(int i=0; i<D; i++){
                		if(i==closest){
                			normal[i]=1.0;
                			index=i;
                		}
                		else{
                			normal[i]=0.0;
                			nextCell[i]=currentCell[i];
                		}
                	}                	
                }
                
                
                
                // get plane's distance from origin
                if(currentCell[index]>nextCell[index]){
                    d[0]=currentCell[index]*l;
                }
                else{
                    d[0]=nextCell[index]*l;
                }
                
                
                double pointPlaneDist = Math.abs(normal[index]*walker.r[index] - d[0]);
        		double stepLength=0.0;
        		for(int i=0; i<D; i++){
        			stepLength+=stepVector[i]*stepVector[i];
        		}
        		stepLength=Math.sqrt(stepLength);

                // check if we should be ignoring this barrier
            	if(skipCurrent){  // this means we're sitting on a barrier

            		double normalDotStep=0.0;
            		for(int i=0; i<D; i++){
            			normalDotStep+=normal[i]*stepVector[i];
            		}
                
            		double cosTheta = normalDotStep/stepLength;
                               		
            		double relStepLengthToBarrier = pointPlaneDist/(cosTheta*origLength);
            		
            		// check if we are very close to current barrier
            		// and that this barrier is the right distance from origin
            		if((Math.abs(relStepLengthToBarrier)<=TINYNUM)&&(Math.abs(d[0]-skippingDistance)<=TINYNUM)){
            			boolean skipping=true;
            			
            			// if so, check the barrier normal against the skipping normal
            			for(int j=0; j<D; j++){
            				if(Math.abs(normal[j]-skippingNormal[j])>TINYNUM){
            					skipping=false;
            					break;
            				}
            			}
            			
            			// if skipping this barrier, skip on...
            			if(skipping){
            				return false;
            			}
            		}                       		
            	}
                
                
                
                // if the following is true, then something is very wrong...
                if(pointPlaneDist>stepLength){
                	logger.severe("step=("+stepVector[0]/l+", "+stepVector[1]/l+", "+stepVector[2]/l+")    cell "+currentCell[0]+", "+currentCell[1]+", "+currentCell[2]);
                	logger.severe("new pos = ("+stepPos[0]/l+", "+stepPos[1]/l+", "+stepPos[2]/l+")    cell "+nextCell[0]+", "+nextCell[1]+", "+nextCell[2]);
                	throw new RuntimeException("distance of walker at ("
                			+walker.r[0]/l+", "+walker.r[1]/l+", "+walker.r[2]/l+") from plane at "
                			+d[0]/l+" index "+index+" is "+pointPlaneDist/l+" > "+stepLength/l);
                }

                p[0]=this.p;

                return true;
            }
        }
        
        return false;
    }
    
    /** returns the size of the substrate
     * @see simulation.geometry.Geometry#getSubstrateSize()
     */
    public final double[] getSubstrateSize() {
        return new double[]{L*l, L*l, L*l};
    }

    public Collection<Triangle> getTriangles(){
    	
    	int n= (int)Math.ceil(L/l);

    	Collection<Triangle> tris= new ArrayList<Triangle>(n*n*n*12);    	
    	
    	int[] cellInd= new int[D];
    	
    	for(cellInd[0]=0; cellInd[0]<n; cellInd[0]++){
    		for(cellInd[1]=0; cellInd[1]<n; cellInd[1]++){
    			for(cellInd[2]=0; cellInd[2]<n; cellInd[2]++){
    				if(getCellOccupation(cellInd)){
    					tris.addAll(trianglesForCubeAt(cellInd));
    				}
    			}
    		}
    	}
    	
    	return tris;
    }
    
    private final Collection<Triangle> trianglesForCubeAt(int[] cell){
    	
    	double[] b = new double[]{l*cell[0], l*cell[1], l*cell[2]};
    	double[] t= new double[]{l*(cell[0]+1), l*(cell[1]+1), l*(cell[2]+1)};
    	
    	Collection<Triangle> cube= new ArrayList<Triangle>(12);
    	
    	
    	// left
    	double[] v1= new double[]{b[0], b[1], b[2]};
    	double[] v2= new double[]{b[0], t[1], t[2]};
    	double[] v3= new double[]{b[0], t[1], b[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));
    	
    	v1= new double[]{b[0], t[1], t[2]};
    	v2= new double[]{b[0], b[1], b[2]};
    	v3= new double[]{b[0], b[1], t[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));


    	// right
    	v1= new double[]{t[0], b[1], b[2]};
    	v2= new double[]{t[0], t[1], t[2]};
    	v3= new double[]{t[0], t[1], b[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));
    	
    	v1= new double[]{t[0], t[1], t[2]};
    	v2= new double[]{t[0], b[1], b[2]};
    	v3= new double[]{t[0], b[1], t[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));
    	
    	
    	// top
    	v1= new double[]{t[0], t[1], b[2]};
    	v2= new double[]{b[0], t[1], b[2]};
    	v3= new double[]{b[0], t[1], t[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));

    	v1= new double[]{t[0], t[1], b[2]};
    	v2= new double[]{b[0], t[1], t[2]};
    	v3= new double[]{t[0], t[1], t[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));

    	
    	// bottom
    	v1= new double[]{t[0], b[1], b[2]};
    	v2= new double[]{b[0], b[1], b[2]};
    	v3= new double[]{b[0], b[1], t[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));

    	v1= new double[]{t[0], b[1], b[2]};
    	v2= new double[]{b[0], b[1], t[2]};
    	v3= new double[]{t[0], b[1], t[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));
    	
    	
    	// front
    	v1= new double[]{t[0], t[1], t[2]};
    	v2= new double[]{b[0], t[1], t[2]};
    	v3= new double[]{b[0], b[1], t[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));

    	v1= new double[]{t[0], t[1], t[2]};
    	v2= new double[]{b[0], b[1], t[2]};
    	v3= new double[]{t[0], b[1], t[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));

    	
    	// back
    	v1= new double[]{t[0], t[1], b[2]};
    	v2= new double[]{b[0], t[1], b[2]};
    	v3= new double[]{b[0], b[1], b[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));

    	v1= new double[]{t[0], t[1], b[2]};
    	v2= new double[]{b[0], b[1], b[2]};
    	v3= new double[]{t[0], b[1], b[2]};
    	
    	cube.add(new Triangle(v1, v2, v3, 0.0));

    	
    	return cube;
    }

    /** applies periodic boundary conditions to walker positions
     * 
     * @param walker the array of walkers walker
     */
    public void applyBoundaries(Walker[] walker){

/*        double substrateSize=getSubstrateSize();
        
        for(int i=0; i<walker.length; i++){
            for(int j=0; j<D; j++){
                if(walker[i].r[j]<0.0){
                    walker[i].r[j]+=substrateSize;
		    walker[i].r0[j]+=substrateSize;
                }
                if(walker[i].r[j]>=substrateSize){
                    walker[i].r[j]-=substrateSize;
		    walker[i].r0[j]-=substrateSize;
                }
            }
        }*/
    }
    
    /** initialiser. does nothing.
     * 
     */
    public void init(){
    	
    }
    
    
    /**
     * constructs the array of "winding numbers" for walker positions.
     * This is loosely the number of times their paths have wrapped around
     * the lattice in each direction. It is used to map positions off the 
     * lattice back onto it in a periodic fashion.
     * 
     * @param r euclidean position in space
     * 
     * @return array of D winding numbers
     */
    public int[] getWindingNumbers(double[] r){
    	
    	int[] windingNumbers = new int[D];
    	
    	double[] latticeSize= getSubstrateSize();
    	    	
    	for(int i=0; i<D; i++){
    		windingNumbers[i]= (int)Math.floor(r[i]/latticeSize[i]);
    	}

    	return windingNumbers;
    }
    
    /** wrapper for winding nimber routine using walker instead of position
     * 
     * @param walker walker whose winding number we want to find
     * @return winding number array
     */
    public int[] getWindingNumbers(Walker walker){
    	return getWindingNumbers(walker.r);
    }
    
    
    /**
     * the reason for this is to avoid initialising everyone on top of a 
     * membrane on lattices with even L. That really screws up the step
     * amending code.
     * 
     * @return a coord for the peak vector
     */
    public double getPeakCoord(){
        
        if(L%2==1){
            return l*(double)L/2.0;
        }
        else{
            return l*(double)(L+1)/2.0;
        }
    }
    
    /** calculates the normal between two adjacent cells
     * 
     * @param currentCell coords of starting cell
     * @param nextCell coords of finishing cell
     * @param normal space to store the normal
     *
     */
    private int calcBarrierNormal(int[] currentCell, int[] nextCell, double[] normal){
    	
    	int index=-1;
    	int numDiffs=0;
    	
    	for(int i=0; i<D; i++){
    		if(currentCell[i]!=nextCell[i]){
    			normal[i]=1.0;
    			index=i;
    			numDiffs++;
    		}
    		else{
    			normal[i]=0.0;
    		}
		}
    	
    	
    	if(numDiffs==1){
    		return index;
    	}
    	else{
        	for(int i=0; i<D; i++){
        		normal[i]/=Math.sqrt((double)numDiffs);
        	}

    		return -numDiffs;
    	}
    }
    
    public static void testGetCell(){
    	
    	SimulableScheme scheme;
    	try{
            URI uri= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= uri.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }
        
        
        
    	SimulationParams.sim_l=1.0;
    	SimulationParams.sim_L=20;
    	
    	SimulationParams simParams = new SimulationParams(1, 1000, 
                0.0, SimulationParams.UNIFORM, SubstrateType.CELL_STRIPED, 
                StepType.FIXEDLENGTH, 1.0, scheme);

    	
    	
    	CellularLattice lattice = new StripedCellularLattice(simParams);
    	
    	double[] r;
    	int[] cell;
    	
    	// test cell coords
    	r= new double[]{10.0, 10.0, 10.0};
    	cell = lattice.getCell(r);
    	
    	System.err.println("("+r[0]+", "+r[1]+", "+r[2]+") -- > ("
    			+cell[0]+","+cell[1]+","+cell[2]+")");

    
    	r= new double[]{-10.0, 10.0, -10.0};
    	cell = lattice.getCell(r);
    	
    	System.err.println("("+r[0]+", "+r[1]+", "+r[2]+") -- > ("
    			+cell[0]+","+cell[1]+","+cell[2]+")");
    	
    	r= new double[]{10.0, 30.0, 10.0};
    	cell = lattice.getCell(r);
    	
    	System.err.println("("+r[0]+", "+r[1]+", "+r[2]+") -- > ("
    			+cell[0]+","+cell[1]+","+cell[2]+")");

    	r= new double[]{11.0, -22.5, 33.1};
    	cell = lattice.getCell(r);
    	
    	System.err.println("("+r[0]+", "+r[1]+", "+r[2]+") -- > ("
    			+cell[0]+","+cell[1]+","+cell[2]+")");

    	r= new double[]{-1.1, 0.7, -103.12};
    	cell = lattice.getCell(r);
    	
    	System.err.println("("+r[0]+", "+r[1]+", "+r[2]+") -- > ("
    			+cell[0]+","+cell[1]+","+cell[2]+")");

}
    
    
    public static void testWindingNumbers(){
    	
    	int tmax=10;

        SimulableScheme scheme;
        try{
            URI uri= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= uri.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }

        SimulationParams.sim_l=1.0;
        SimulationParams.sim_L=20;
    	
    	SimulationParams simParams = new SimulationParams(1, tmax, 
                0.0, SimulationParams.UNIFORM, SubstrateType.CELL_STRIPED,
                StepType.FIXEDLENGTH, 1.0, scheme);


    	
    	// construct lattice object
    	CellularLattice lattice = new StripedCellularLattice(simParams);
    	
    	// and a walker or two
    	Walker walker= new Walker(new double[] {0.0, 10.0, 0.0});
    	
    	// also need winding numbers array
    	int[] windingNumbers;
    	
    	System.err.println("lattice size 20, cell size 1.0");
    	
    	// test the winding number routine -- single direction
    	windingNumbers= lattice.getWindingNumbers(walker);
    	
    	System.err.println("(0.0, 10.0, 0.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");
 
    	walker = new Walker(new double[]{0.0, 30.0, 0.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(0.0, 30.0, 0.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	walker = new Walker(new double[]{0.0, -10.0, 0.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(0.0, -10.0, 0.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	walker = new Walker(new double[]{0.0, -30.0, 0.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(0.0, -30.0, 0.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	// two directions
    	walker = new Walker(new double[]{10.0, 10.0, 0.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(10.0, 10.0, 0.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	walker = new Walker(new double[]{10.0, -10.0, 0.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(10.0, -10.0, 0.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	walker = new Walker(new double[]{30.0, -30.0, 0.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(30.0, -30.0, 0.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	walker = new Walker(new double[]{-30.0, -30.0, 0.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(-30.0, -30.0, 0.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	// three directions
    	
    	walker = new Walker(new double[]{10.0, 10.0, 10.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(10.0, 10.0, 10.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	
    	walker = new Walker(new double[]{10.0, 10.0, -10.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(10.0, 10.0, -10.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	
    	
    	walker = new Walker(new double[]{-10.0, 10.0, -10.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(-10.0, 10.0, -10.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	
    	walker = new Walker(new double[]{10.0, 30.0, 10.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(10.0, 30.0, 10.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");

    	
    	walker = new Walker(new double[]{100.0, -30.0, -10.0});
    	
    	windingNumbers = lattice.getWindingNumbers(walker);

    	System.err.println("(100.0, -30.0, -10.0) --> ("+windingNumbers[0]+", "
    			+windingNumbers[1]+", "+windingNumbers[2]+")");


    	
    }

    public static void testBarrierReflection(){
    	
    	System.err.println("testing barrier reflection code");
    	
    	int tmax =100000;
    	
    	int D=DiffusionSimulation.D;
    
    	SimulationParams.sim_l=2.0;
    	SimulationParams.sim_L=2;
    	
        SimulableScheme scheme;
        try{
            URI uri= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= uri.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }
    	
    	SimulationParams simParams = new SimulationParams(1000, tmax, 0.0, SimulationParams.UNIFORM, 
    	        SubstrateFactory.SubstrateType.CELL_STRIPED, StepType.FIXEDLENGTH, 
    	        10.0, scheme);
    	
    	Substrate substrate = SubstrateFactory.getSubstrate(SubstrateType.CELL_STRIPED, simParams);
    	
    	System.err.println("testing single crossing");
    	Walker walker= new Walker(new double[]{1.5, 1.2, 1.5});
    	
    	double[] step= new double[]{1.0, -1.0, 0.0};
    	
    	double[] normal= new double[D];
    	
    	double[] d = new double[1];
    	
    	double[] toBarrier = new double[D];
    	double[] amended = new double[D];
    	double[] unamended = new double[D];
    	
    	boolean[] in = new boolean[1];
    	
    	double[] offset= new double[]{0.0, 0.0, 0.0};
    	
    	double[] p= new double[1];
    	
    	// test step amendment
    	System.err.println("\t... membrane crossing detection");
    	boolean singleCrossesMembrane=false;
		try {
			singleCrossesMembrane = substrate.crossesMembrane(walker, offset, step, normal, d, false, Math.sqrt(2.0), in, p, false, null);
		} catch (StepRejectedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    	System.err.println("* single crossing flag is "+singleCrossesMembrane);
    	System.err.println("* normal between ("+walker.r[0]+","+walker.r[1]+","+walker.r[2]+")"
    						+"and ("+(walker.r[0]+step[0])+","+(walker.r[1]+step[1])+","+(walker.r[2]+step[2])
    						+"is ("+normal[0]+","+normal[1]+","+normal[2]+") distance "+d[0]);
    	
    	System.err.println("\t... step amendment");
    	substrate.testAmendment(walker, step, normal, d, Math.sqrt(2.0), toBarrier, amended, unamended);
    	
    	System.err.println("* step to barrier is ("+toBarrier[0]+","+toBarrier[1]+","+toBarrier[2]+")");
    	System.err.println("* amended step is ("+amended[0]+","+amended[1]+","+amended[2]+")");
    	System.err.println("* unamended step is "+unamended[0]+","+unamended[1]+","+unamended[2]+")");
    	
    	// test barrier skipping 
    	System.err.println("\t... barrier skipping");
    	walker.makeStep(toBarrier);
    	
    	boolean unreflectedCrossing=false;
		try {
			unreflectedCrossing = substrate.crossesMembrane(walker, offset, unamended, normal, d, true, Math.sqrt(2.0), in, p, false, null);
		} catch (StepRejectedException e) {
			e.printStackTrace();
		} 
		
    	boolean reflectedCrossing=false;
		try {
			reflectedCrossing = substrate.crossesMembrane(walker, offset, amended, normal, d, true, Math.sqrt(2.0), in, p, false, null);
		} catch (StepRejectedException e) {
			e.printStackTrace();
		}
    	
    	
    	System.err.println("* unflected crossing flag is "+unreflectedCrossing+" ("
    			+walker.r[0]+","+walker.r[1]+","+walker.r[2]+") -> ("
    			+(walker.r[0]+unamended[0])+","+(walker.r[1]+unamended[1])+","+(walker.r[2]+unamended[2])+")");
    	System.err.println("* reflected crossing flag is "+reflectedCrossing+ "("
    			+walker.r[0]+","+walker.r[1]+","+walker.r[2]+") -> ("
    			+(walker.r[0]+amended[0])+","+(walker.r[1]+amended[1])+","+(walker.r[2]+amended[2])+")");

    	System.err.println("\t... testing multiple barrier reflection");
    	walker= new Walker(new double[]{1.5,2.6,1.5});
    	
    	boolean multipleCrossing=false;
		try {
			multipleCrossing = substrate.crossesMembrane(walker, offset, step, normal, d, false, Math.sqrt(2.0), in, p, false, null);
		} catch (StepRejectedException e) {
			e.printStackTrace();
		}
    	
    	System.err.println("* multiple crossing flag is "+multipleCrossing);
    	System.err.println("* first normal encountered is ("+normal[0]+","+normal[1]+","+normal[2]+") distance "+d[0]+" units from origin");
    	
    	substrate.testAmendment(walker, step, normal, d, Math.sqrt(2.0), toBarrier, amended, unamended);
    	
    	System.err.println("* unflected crossing flag is "+unreflectedCrossing+" ("
    			+walker.r[0]+","+walker.r[1]+","+walker.r[2]+") -> ("
    			+(walker.r[0]+unamended[0])+","+(walker.r[1]+unamended[1])+","+(walker.r[2]+unamended[2])+")");
    	System.err.println("* reflected crossing flag is "+reflectedCrossing+ "("
    			+walker.r[0]+","+walker.r[1]+","+walker.r[2]+") -> ("
    			+(walker.r[0]+amended[0])+","+(walker.r[1]+amended[1])+","+(walker.r[2]+amended[2])+")");
    	
    	walker.makeStep(toBarrier);
    	
    	try {
			unreflectedCrossing= substrate.crossesMembrane(walker, offset, unamended, normal, d, true, Math.sqrt(2.0), in, p, false, null);
		} catch (StepRejectedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	try {
			reflectedCrossing= substrate.crossesMembrane(walker, offset, amended, normal, d, true, Math.sqrt(2.0), in, p, false, null);
		} catch (StepRejectedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    	System.err.println("* at barrier...");
    	System.err.println("* unreflected crossing is "+unreflectedCrossing);
    	System.err.println("* reflected crossing is "+reflectedCrossing);
    	
    }
    
    
    
    public static void testBoundaryCrossing(){

        // initialise a simualtion on a cellular lattice with
        // impenetrable bounaries
        
        SimulableScheme scheme;
        try{
            URI uri= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= uri.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }

        SimulationParams.sim_l=1.0;
        SimulationParams.sim_L=20;
        SimulationParams.sim_stripethickness=1;
        
        
        
        SimulationParams simParams= new SimulationParams(1000, 100000, 0.0, SimulationParams.UNIFORM, 
                SubstrateType.CELL_STRIPED, StepType.FIXEDLENGTH, 
                10.0, scheme);

        // the step length (normally set from diffusion constant)
        double[] stepParams={0.1};
        
        simParams.setStepParams(stepParams);
        
        DiffusionSimulation diffSim= new DiffusionSimulation(simParams, scheme);
  
        diffSim.nextVoxel();
        
        double meanSquareDisp = diffSim.getMeanSquareDisplacement();
        
        //assertTrue(meanSquareDisp<=2.0/Math.sqrt(2.0));
        
    }
    
    public static void testPeriodicBoundaries(){
        double[] r0= new double[] {0.0, 0.0, 0.0};
        double[] step= new double[] {-0.1, 0.0, 0.0};
        
        Walker[] walker = new Walker[1];
        
        walker[0]=new Walker(r0);
        
        SimulableScheme scheme;
        try{
            URI uri= DiffusionSimulation.class.getResource("/test/bmx7_ED.scheme1").toURI();
            
            String path= uri.getPath();
            
            scheme= (SimulableScheme)RectGradSteTanScheme.readScheme(path);
        }
        catch(URISyntaxException urise){
            throw new LoggedException(urise);
        }

    	SimulationParams simParams = new SimulationParams(1, 1000, 
                0.0, SimulationParams.UNIFORM, SubstrateType.CELL_STRIPED,
                StepType.FIXEDLENGTH, 1.0, scheme);


    	SimulationParams.sim_l=1.0;
    	SimulationParams.sim_L=1;
        
        CellularLattice lattice = new StripedCellularLattice(simParams);
        
        lattice.applyBoundaries(walker);
        
        /*
        assertTrue(step[0]>0.0);
        assertTrue(step[0]<1.0);
        */
        
    }

	public boolean intracellular(Walker walker){

		if(getCellOccupation(getCell(walker.r))){
			return true;
		}
		else{
			return false;	
		}
		
	}

    
    public static void main(String[] args) {

    	System.err.println("Camino diffusion simulation engine test code\n");
    	testBarrierReflection();
    	
    }
}

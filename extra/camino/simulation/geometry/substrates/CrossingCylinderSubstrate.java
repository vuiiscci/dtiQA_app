package simulation.geometry.substrates;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.logging.Logger;

import numerics.RealMatrix;
import numerics.Rotations;
import numerics.Vector3D;
import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.Walker;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.dynamics.exceptions.TooDamnCloseException;
import simulation.geometry.elements.BasicCylinder;
import simulation.geometry.elements.Cylinder;
import simulation.geometry.elements.Triangle;


/**
 * substrate will two principle directions.
 * cylinders are square-packed and packed into
 * alternating sheets:
 * 
 * o o o
 * =====
 * o o o
 * 
 * default angle between populations 90 degrees.
 * 
 * The implementation here differs slightly from 
 * the parallel cylinders in that the unit cell 
 * needs to be larger. Here the cell is 2L*2L*2L
 * (cylinder sep L) so that we can cover two rows
 * of cylinders.
 * 
 * This work by having two cylinders aligned to the z axis.
 * the substrate coordinate system for the cell is split
 * into two: the top half is rotated about he y axis 
 * with respect tot he lower, which is unrotated.
 * 
 * Walker positions are rotated with respect to the cylinder 
 * in the top coord system and edges are dealt with 
 * accrordingly, eliminating the edge match issue.
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class CrossingCylinderSubstrate extends CylinderSubstrate {

	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** dimensionality of space */
	private final int D= DiffusionSimulation.D;
	
	/** crossing angle */
	private final double cAngle;
	
	/** cylinder separation */
	private final double R;
	
	/** rotation matrix into crossing coord system */
	private final RealMatrix invRotMat;
	
	/** rotation matrix out of crossing coord system */
	private final RealMatrix rotMat;
	
	/**  has the walker being checked for intersection been rotated? */
	private boolean rotated= false;
	
	/** height of central plane if going from lower to upper */
	private final double xMidUp;
	
	/** height of central plane if going from upper to lower */
	private final double xMidDown;
	
	/** 
	 * constructor for default crossing angle.
	 * 
	 * @param R cylinder separation
	 * @param r cylinder radius
	 */
	public CrossingCylinderSubstrate(SimulationParams simParams){
		super(new double[]{SimulationParams.sim_R, 
		                   2.0*SimulationParams.sim_R, 
		                   SimulationParams.sim_R}, simParams, false);
		
		double R=SimulationParams.sim_R;
		double r=SimulationParams.sim_r;
		
		SimulationParams.sim_l=SimulationParams.sim_R*2.0;
		
		if(R<=2.0*r){
			logger.warning("Cylinder separation on crossing cylinder substrate is less than twice the radius.");
			logger.warning("overlapping or abutting cylinders on this substrate will behave unpredictably");
		}
		
		Cylinder[] cylinder= new BasicCylinder[2];
		
		double p= simParams.getP();
		
		// set crossing angle
		this.cAngle= SimulationParams.sim_cAngle;

		double[] yAxis= new double[]{0.0, 1.0, 0.0};
		
		// construct crossing population cylinder axis
		rotMat= Rotations.getRotMat(yAxis, -cAngle);
		invRotMat= Rotations.getRotMat(yAxis, cAngle);
		
		// set cylinder separation
		this.R= R;
		
		// set central plane heights
		this.xMidUp= R*(1.0+1.0/a);
		this.xMidDown= R*(1.0-1.0/a);
		
		// bottom row (z-axis)
		cylinder[0]= new BasicCylinder(new double[]{R/2, R/2, 0.0}, r, p);
		cylinder[1]= new BasicCylinder(new double[]{R/2, 3*R/2, 0.0}, r, p);

		super.setCylinders(cylinder);
		
	}
	
	/**
	 * override the get substrate coords method from substrate
	 * 
	 */
	public void getSubstrateCoords(double[] walkerPos, double[] offset, double[] newPos){
	    
	    // get substrate coords
	    super.getSubstrateCoords(walkerPos, offset, newPos);
	    rotated=false;
	    
	    if(newPos[1]>=R){
	        // if in the top half of the substrate, rotate the point
	        // and recheck substrate coords
	        double[] tempPos= new double[D];
        
	        // must include offset in rotation
	        for(int i=0; i<D; i++){
	            tempPos[i]= walkerPos[i]+offset[i];
	        }
        
	        // rotate into upper coord system
	        tempPos= Rotations.transformPoint(invRotMat, tempPos);
	        
	        // recheck substrate coords
	        super.getSubstrateCoords(tempPos, new double[]{0.0, 0.0, 0.0}, newPos);
	        
	        // 
	        rotated=true;
	    }
	    
	    
	}
	
	/**
	 * overrides the default step mapper to take account of the rotated patch 
	 * in the coordinate system. if the spin position needs to be rotated, this 
	 * will rotate the step as well, otherwise it'll just return the original step
	 * 
	 * @param walker the walker
	 * @param offset displacement from current position from which step is taken
	 * @param step the step we want to make
	 * 
	 * @return appropriately rotated or original step
	 */
	public double[] mapStepIntoSubstrate(Walker walker, double[] offset, double[] step){
	    
	    double[] newPos= new double[D];
	    
	    // check position for rotation
	    getSubstrateCoords(walker.r, offset, newPos);
	    
	    if(rotated){
	        // if we're in the rotated part of the substrate rotate the step
	        return Rotations.transformPoint(invRotMat, step);
	    }
	    
	    // otherwise just return the original step
	    return step;
	}
	
	
	/**
	 * if the step has been rotated, this undoes the rotation
	 * 
	 * @param subsCoords current spin coords in substrate frame
	 * @param offset offset from current position
	 * @param step vector to un-rotate
	 * 
	 * @return unrotated step or original is no rotation has occurred
	 */
	public double[] unmapStepFromSubstrate(double[] subsCoords, double[] offset, double[] step){
	    
	    if(rotated){
	        
	        return Rotations.transformPoint(rotMat, step);
	        
	    }
	    
	    return step;
	}
	
	/** 
     * checks if a walker's step will take it across a membrane or not.
     * this just involves checking every cylinder in turn. The cylinder
     * crosses() method should fill in the blanks where necessary.
     * 
     * @param walker the walker to check
     * @param offset offset from current walker position for beginning of step
     * @param step the step being made
     * @param normal space to store the normal if a barrier is crossed
     * @param d space to store barrier distance from origin dotted with normal
     * @param skipCurrent flag indicating we're sitting on a barrier and should 
     *        ignore the closest one.
     * @param originLength the original length of the step vector;
     * 
     * @return true or false
	 * @throws StepRejectedException 
     */
    public boolean crossesMembrane(Walker walker, double[] rawOffset, double[] rawStep,
            double[] normal, double[] d, boolean skipCurrent, double origLength, 
            boolean[] in, double[] p, boolean report, FileWriter debugWriter) throws StepRejectedException {
        
        double len= 0.0;
        double[] walkerPos= new double[D];
        double[] intDist= new double[1];
        
        boolean[] intIn= new boolean[1];
        double[] intP= new double[1];
        
        for(int i=0; i<rawStep.length; i++){
            len += rawStep[i]*rawStep[i];
        }
        len=Math.sqrt(len);
        
        if(len/origLength<=1E-14){
            return false;
        }

        double[] offset;
        double[] step;
        getSubstrateCoords(walker.r, rawOffset, walkerPos);

        if(rotated){
            // need to rotate the step
            step=Rotations.transformPoint(invRotMat, rawStep);
            offset=new double[D];//Rotations.transformPoint(invRotMat, rawOffset);
        }
        else{
            step=rawStep;
            offset=rawOffset;
        }
        
        
        
        int numIntersections=0;
        boolean skip=skipCurrent;
        
        for(int i=0; i<cylinder.length; i++){
            
            if(skipCurrent){
                if(i==lastCrossed){
                    intIn[0]=in[0];
                    skip=true;
                }
                else{
                    intIn[0]=cylinder[i].inside(walkerPos);
                    skip= false;
                }               
            }
            
            boolean crosses=false;
            try{
            	crosses=cylinder[i].crosses(walkerPos, step, normal, d, skip, origLength, intDist, intIn, intP, walker.R);
            }
            catch(TooDamnCloseException tdce){
            	throw new StepRejectedException(tdce.getMessage());
            }
            	
            // check the cylinders
            if(crosses){

                shortestDist[numIntersections]=intDist[0];
                intOriginDist[numIntersections]=d[0];
                for(int j=0; j<D; j++){
                    intNormal[numIntersections][j]=normal[j];
                }
                intInside[numIntersections]=intIn[0];
                cylCrossed[numIntersections]=i;

                cyl_p[numIntersections]=intP[0];
                
                
                isBoundary[numIntersections]=false;
                
                numIntersections++;
                
            }
            
        }
        
        // check if we've intersected the boundaries of the substrate
        if(checkBoundaryIntersection(walkerPos, offset, step, normal, d, skipCurrent, origLength, intDist, intIn, intP)){
            shortestDist[numIntersections]=intDist[0];
            intOriginDist[numIntersections]=d[0];
            for(int j=0; j<D; j++){
                intNormal[numIntersections][j]=normal[j];
            }
            intInside[numIntersections]=intIn[0];
            
            // in this case we haven't crossed a cylinder, so must check them all next time
            cylCrossed[numIntersections]=-1; 

            cyl_p[numIntersections]=intP[0];
            
            
            isBoundary[numIntersections]=true;
        
            numIntersections++;
        }
        
        
        // check if we've intersected the boundary of the coordinate systems
        if(checkMidPlaneIntersection(walkerPos, offset, step, normal, d, skipCurrent, origLength, intDist, intIn, intP)){
            shortestDist[numIntersections]=intDist[0];
            intOriginDist[numIntersections]=d[0];
            for(int j=0; j<D; j++){
                intNormal[numIntersections][j]=normal[j];
            }
            intInside[numIntersections]=intIn[0];
            
            // in this case we haven't crossed a cylinder, so must check them all next time
            cylCrossed[numIntersections]=-1; 

            cyl_p[numIntersections]=intP[0];
            
            
            isBoundary[numIntersections]=true;
            
            numIntersections++;
        }
        
        
  
        if(numIntersections==1){
            // if there's one interesection everything is  
            // as it should be so just copy the distance
            // and normal over and return true
            d[0]=intOriginDist[0];
            for(int i=0; i<D; i++){
                normal[i]=intNormal[0][i];
            }
            in[0]=intInside[0];
            
            lastCrossed=cylCrossed[0];
            
            p[0]=cyl_p[0];
            
            
            intersectsBoundary= isBoundary[0];
            
            return true;
        }
        
        if(numIntersections>1){
            
            // if there are more than intersections
            // we must pick the one that would happen
            // first, and set up the distance and 
            // normal accordingly
            double closestDist=Double.MAX_VALUE;
            int closest=-1;
            
            // find closest interesection to the walker
            for(int i=0; i<numIntersections; i++){
                if(shortestDist[i]<closestDist){
                    closestDist=shortestDist[i];
                    closest=i;
                    lastCrossed=cylCrossed[i];
                }
            }
            
            // find the second closest
            double secondClosestDist= Double.MAX_VALUE;
            for(int i=0; i<numIntersections; i++){
            	if(i==closest){
            		continue;
            	}

            	if(shortestDist[i]<secondClosestDist){
            		secondClosestDist= shortestDist[i];
            	}
            }

            // check that there aren't two barrier within walker radius
            if((secondClosestDist-closestDist)*len<=walker.R){
            	throw new StepRejectedException("Two barriers within walker radius. rejecting step");
            }
            
            if(closest==-1){
                return false;
            }
            
            // set the distance and normal
            d[0]=intOriginDist[closest];
            
            for(int i=0; i<D; i++){
                normal[i]=intNormal[closest][i];
            }
            in[0]=intInside[closest];
            
            p[0]=cyl_p[closest];
            
            
            intersectsBoundary= isBoundary[closest];
            
            return true;
        }
        
        // in this case there are no intersections
        return false;
    }
	
	
	
	/**
	 * similar to the boundary intersection checking in substrate,
	 * this check is a step will cross the xy-plane half way up the
	 * substrate -- this is the border between the two coordinate
	 * systems an needs to be checked in addition to everything 
	 * else. 
	 * 
	 * this check itself is almost identical to the substrate edge
	 * check. there is a small offset to the position of the plane 
	 * that eliminates numerical nonsense using the same constant a
	 * that's used in the substrate edge check.
	 * 
	 * @param walkerPos the position of the walker in substrate coords
     * @param offset current offset from starting point
     * @param stepVector proposed step vector
     * @param normal space to store boundary normal
     * @param d distance to boundary
     * @param skipCurrent should we skip the current boundary
     * @param origLength original step length
     * @param in initially intracellular?
     * @param p permeability of boundary (always one here)
     * 
     * @return true is boundary crossed
	 */
	protected boolean checkMidPlaneIntersection(double[] walkerPos, double[] offset, double[] stepVector,
            double[] normal, double[] d, boolean skipCurrent, double origLength, double[] intDist,  
            boolean[] in, double[] p){
	
	        
	    if(walkerPos[1]>=R){
	        if(walkerPos[1]+stepVector[1]<R){
	            // step crosses from upper half to lower
	            double dist=(xMidDown-walkerPos[1])/stepVector[1];
	            d[0]= xMidDown;
	            p[0]= 1.0;
	            normal[0]=0;
	            normal[1]=1;
	            normal[2]=0;
	            intDist[0]=dist;
	            
	            return true;
	        }
	    }
	    else{
	        if(walkerPos[1]+stepVector[1]>=R){
                // step crosses from lower half to upper
                double dist=(xMidUp-walkerPos[1])/stepVector[1];
                d[0]= xMidUp;
                p[0]= 1.0;
                normal[0]=0;
                normal[1]=1;
                normal[2]=0;
                intDist[0]=dist;
                
                return true;
            }
	    }
	    
	    // if we're here, the step isn't crossing the central line
	    return false;
	}
	
	
	/**
	 * returns the peak coordinate (the place where the
	 * spike sits on each axis for delta-peaked initial
	 * conditions)
	 * 
	 * @return cylinder sep (half the size of the unit cell)
	 * 
	 */
	public double getPeakCoord() {
		
		return R;
	}

	/**
	 * return the size of the unit cell
	 * 
	 * @return 2 x cylinder sep
	 */
	public double[] getSubstrateSize() {
		
		return new double[] {R, 2.0*R, R};
	}

	
    /**
     * get the rotation matrix
     * 
     * @return elements of rotMat in OpenGL order
     */
    public double[] getRotMat(){
        
        double[]  rotMatArr = new double[D*D];
        
        for(int i=0; i<D; i++){
            for(int j=0; j<D; j++){
                int index= D*j +i;
                
                rotMatArr[index]= rotMat.entry(i, j);
            }
        }
            
        return rotMatArr;
    }

    /**
     * get the inverse rotation matrix
     * 
     * @return elements of rotMat in OpenGL order
     */
    public double[] getInvRotMat(){
        
        double[]  invRotMatArr = new double[D*D];
        
        for(int i=0; i<D; i++){
            for(int j=0; j<D; j++){
                int index= D*j +i;
                
                invRotMatArr[index]= invRotMat.entry(i, j);
            }
        }
            
        return invRotMatArr;
    }
    
    /**
     * override get triangles method to account for rotation of upper cylinder
     * 
     */
    public Collection<Triangle> getTriangles(){
		
		Collection<Triangle> tris= new ArrayList<Triangle>(subsObj.length);
		
		// unrotated cylinders
		tris.addAll(cylinder[0].getTriangles());
		double[] pos= cylinder[0].getPosition();
		pos[0]+=R;
		Cylinder cyl1= new BasicCylinder(pos, cylinder[0].getRadius(), 0.0);
		tris.addAll(cyl1.getTriangles());
		pos[0]-=2*R;
		Cylinder cyl2= new BasicCylinder(pos, cylinder[0].getRadius(), 0.0);
		tris.addAll(cyl2.getTriangles());
		
		pos[1]+=2*R;
		cyl1= new BasicCylinder(pos, cylinder[0].getRadius(), 0.0);
		tris.addAll(cyl1.getTriangles());
		pos[0]+=R;
		cyl1= new BasicCylinder(pos, cylinder[0].getRadius(), 0.0);
		tris.addAll(cyl1.getTriangles());
		
		pos[0]+=R;
		cyl1= new BasicCylinder(pos, cylinder[0].getRadius(), 0.0);
		tris.addAll(cyl1.getTriangles());
		
        Collection<Triangle> upperCylTris= cylinder[1].getTriangles();
        
        for(Iterator<Triangle> triIt= upperCylTris.iterator(); triIt.hasNext(); ){
        	
        	Triangle tri= triIt.next();
        	
        	double[][] vert= new double[3][];
        	
        	for(int i=0; i<3; i++){
        		vert[i]= tri.getVertex(i);
        		
        		vert[i]= Rotations.transformPoint(invRotMat, vert[i]);
        	}
        	
        	tris.add(new Triangle(vert[0], vert[1], vert[2], 0.0));
        	
        }
        
        pos= cylinder[1].getPosition();
        pos[0]+=R;
		cyl1= new BasicCylinder(pos, cylinder[1].getRadius(), 0.0);
		upperCylTris= cyl1.getTriangles();
        
        for(Iterator<Triangle> triIt= upperCylTris.iterator(); triIt.hasNext(); ){
        	
        	Triangle tri= triIt.next();
        	
        	double[][] vert= new double[3][];
        	
        	for(int i=0; i<3; i++){
        		vert[i]= tri.getVertex(i);
        		
        		vert[i]= Rotations.transformPoint(invRotMat, vert[i]);
        	}
        	
        	tris.add(new Triangle(vert[0], vert[1], vert[2], 0.0));
        	
        }
		
		
		
		pos[0]-=2*R;
		cyl2= new BasicCylinder(pos, cylinder[1].getRadius(), 0.0);
		upperCylTris= cyl2.getTriangles();
        
        for(Iterator<Triangle> triIt= upperCylTris.iterator(); triIt.hasNext(); ){
        	
        	Triangle tri= triIt.next();
        	
        	double[][] vert= new double[3][];
        	
        	for(int i=0; i<3; i++){
        		vert[i]= tri.getVertex(i);
        		
        		vert[i]= Rotations.transformPoint(invRotMat, vert[i]);
        	}
        	
        	tris.add(new Triangle(vert[0], vert[1], vert[2], 0.0));
        	
        }


        
           	
    	return tris; 
        
	}

}

package simulation.geometry.substrates;


import imaging.DW_Scheme;
import imaging.SimulableScheme;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.logging.Logger;

import misc.DT;
import misc.LoggedException;

import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.Walker;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.dynamics.exceptions.StepRejectedException;
import simulation.dynamics.exceptions.TooDamnCloseException;
import simulation.geometry.PLYreader;
import simulation.geometry.elements.SubstrateObject;
import simulation.geometry.elements.Triangle;
import simulation.geometry.substrates.SubstrateFactory.SubstrateType;
import tools.CL_Initializer;

public class MeshSubstrate extends Substrate {

	/** dimensionality of space */
	private static final int D= DiffusionSimulation.D;
	
	/** triangles forming the object mesh */
	private final Collection<Triangle> triangles;
	
	/** fraction of substrate size that is contained in central voxel */
	private final double voxelSizeFrac= SimulationParams.sim_voxelSizeFrac;
	
	/** bottomLeft corner of central voxel */
	private final double[] voxelMin;
	
	/** topRight corner of central voxel */
	private final double[] voxelMax;
	
	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** index of the triangle to skip, if we're skipping */
	private Triangle skipping= null;
	
	/** size of (cubic) unit cell around object */
	//private final double[] sep;
	
	/** length scale that largest vertex displacement vector in object is scaled to */
	private final double size;
	
	/** the ray vector that is used as a step in the intracellular routine */
	private final double[] ray;

	/** map of triangles to their exception coplanar lists */
	private final HashMap<Triangle, ArrayList<Triangle>> coplanarMap;
	
	/** flag indicating if this is a convenient convex hull surface or not */
	private final boolean convexHull;
	
	/** 
	 * read a substrate from the named PLY file and scale it to the 
	 * specified length. 
	 * 
	 * currently assumes that the object in the PLY file is centred
	 * on the origin. this is simple enough to change. eventually
	 * i'll add an option to specify better configurations such as 
	 * multiple objects, positions and rotations.
	 * 
	 * the object will be scaled so that the longest displacement 
	 * vector to a vertex from the origin takes on the length given
	 * in the scale parameter. 
	 * 
	 * @param fname name of PLY file to read
	 * @param scale length to scale object to
	 */
	public MeshSubstrate(SimulationParams simParams){
		
		super(simParams, new double[] {SimulationParams.sim_R, 
		                               SimulationParams.sim_R, 
		                               SimulationParams.sim_R});
		
		

		
		String fname= SimulationParams.sim_plyfile;
		
		// store the object size
		this.size=Double.MAX_VALUE;

		// assemble the ray vector -- parallel to an axis and slightly larger than object
		this.ray=new double[]{0.0, Double.MAX_VALUE, 0.0};
		
		logger.info("reading object mesh from '"+fname+"'");
		triangles =PLYreader.readPLYfile(fname, 1.0, simParams.getP());
		
		this.convexHull= PLYreader.closedSurface;
		
		logger.info("getting coplanar map");
		this.coplanarMap= PLYreader.getCoplanarListMap();
		
		// find centre of mass of all vertices
		logger.info("finding centre of mass");
        double[] m= new double[D];
        int count=0;
		
        Iterator<Triangle> triIt;
		for(triIt= triangles.iterator(); triIt.hasNext();){
		    
		    Triangle triangle = triIt.next();
		    for(int i=0; i<3; i++){
		        double[] vertex= triangle.getVertex(i);
		        for(int j=0; j<D; j++){
		            m[j]+=vertex[j];
		        }
		        
		        count++;
		    }
		}
		
		// normalise by number of vertices
		for(int i=0; i<D; i++){
		    m[i]/=count;
		}
		logger.info("centre of mass = ("+m[0]+","+m[1]+","+m[2]+")");
		
		// top and bottom corners of box
		logger.info("finding top and bottom corners of mesh");
		double[] bottomLeft= new double[]{Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE};
		double[] topRight= new double[]{-Double.MAX_VALUE, -Double.MAX_VALUE, -Double.MAX_VALUE};

		for(triIt= triangles.iterator(); triIt.hasNext(); ){
            
            // for each triangle
            Triangle triangle= triIt.next();

            for(int k=0; k<3; k++){
                // for each vertex
                double[] vert= triangle.getVertex(k);
            
                for(int i=0; i<D; i++){
                    // is this a new record in either direction?
                    if(vert[i]<bottomLeft[i]){
                        bottomLeft[i]=vert[i];
                    }
                    
                    if(vert[i]>topRight[i]){
                        topRight[i]=vert[i];
                    }
                }
            }
	    }
		
		logger.info("Mesh dimension: "+(topRight[0]-bottomLeft[0])+","+(topRight[1]-bottomLeft[1])+","+(topRight[2]-bottomLeft[2]));
		
		/*
		 * check if user has specified a new mesh-separation vector and adjust
		 * the bottom-left and top-right coords accordingly
		 */
		if(SimulationParams.sim_mesh_sep!= null){
			
			for(int i=0; i<D; i++){
				if(SimulationParams.sim_mesh_sep[i]>0.0){
					
					double L= topRight[i]-bottomLeft[i];
					
					logger.info("mesh separation "+SimulationParams.sim_mesh_sep[i]+" in direction "+i);
					
					if(SimulationParams.sim_mesh_sep[i]<=L){
						logger.warning("Specified mesh separation in direction "+i+" is smaller than mesh dimension " +L +
								". Outer portions of the mesh will be ignored.");
					}
					
					double halfDiff= (SimulationParams.sim_mesh_sep[i]-L)/2;
					
					bottomLeft[i]-=halfDiff;
					topRight[i]+=halfDiff;
					
				}
				
			}

			logger.info("voxel bundaries adjusted to specified mesh separation.");
			logger.info("outer and inner voxel dimensions will be based on specified separations");
			
		}
		
		logger.info("moving origin to bottom left of box, reinitialising triangles");
		/* 
		 * subtract bottom left coord off all vertex positions
		 * and reinitialise the intersection checking vectors
		 * in all triangles
		 */
		for(triIt=triangles.iterator(); triIt.hasNext(); ){
		
		    Triangle triangle= triIt.next();

		    for(int k=0; k<3; k++){
		        double[] vert= triangle.getVertex(k);
		        
		        for(int i=0; i<D; i++){
		            vert[i]-=bottomLeft[i];
		        }
	            triangle.setVertex(k, vert);
		    }
		    
		    triangle.initVectors();
		}
		
		logger.info("constructing array for spatial optimisation");
		SubstrateObject[] newTriangles= new Triangle[triangles.size()];
		int triInd=0;
		for(triIt= triangles.iterator(); triIt.hasNext(); triInd++){
		    
		    newTriangles[triInd]= triIt.next();
		    
		}
		
		/* substrate object array in superclass needs to contain 
		 * our processed triangles before we can initialise 
		 * spatial optimisation
		 */
		super.subsObj= newTriangles;
		
		
		// set dimensions of outer and central voxel
		voxelMin= new double[D];
		voxelMax= new double[D];
		
		for(int i=0; i<D; i++){
		    // outer voxel (substrate size)
		    this.L[i]=topRight[i]-bottomLeft[i];
		    super.L[i]= this.L[i];
		    
		    // central region from which data is generated
		    double centre= L[i]/2;
		    double voxSize= voxelSizeFrac*L[i];
		    
		    voxelMin[i]= centre-voxSize/2;
		    voxelMax[i]= centre+voxSize/2;
		    
		}

		// reinitialise values for boundary intersection
		initBoundaryIntersectionArrays();
		
		// TODO: investigate grid sizes
		// initialise spatial optimisation
		int[] n= new int[]{10, 10, 5};
		logger.info("initialising spatial optimisation. n=("+n[0]+","+n[1]+","+n[2]+")");
		initialiseSpatialOptimisation(n);
		
		logger.info("mesh processing complete.");
		
	}
	
	
	
	public boolean crossesMembrane(Walker walker, double[] offset, double[] stepVector,
			double[] normal, double[] d, boolean skipCurrent, double origLength, 
			boolean[] in, double[] p, boolean report, FileWriter debugWriter) throws StepRejectedException {
		
		
		double nearest=2.0;
		
		boolean crossing=false;
		
		double[] tempNormal= new double[D];
		double[] tempD= new double[1];
		double[] intDist= new double[1];
		double[] tempP= new double[1];
		
		
		
		// map walker position into unit cell
		getSubstrateCoords(walker.r, offset, subsCoords);

        initCandidates(walker, offset, stepVector);
        
        Triangle toSkip=null;
        
        int count =0;
		// check interactions with all triangles
		//while(triIt.hasNext()){
		while(moreCandidates()){
			// get the next triangle
			Triangle tri=(Triangle)nextCandidate();

			//ignore null objects
			if(tri==null){
				continue;
			}
			
			
			// check skipping
			if(skipCurrent){
				if(tri==skipping){
					continue;
				}				
			}
			
			
			count++;
			
			boolean crosses=false;
			try{
				crosses=tri.crosses(subsCoords, stepVector, tempNormal, tempD, false, origLength, intDist, null, tempP, walker.R);
			}
			catch(TooDamnCloseException tdce){
				throw new StepRejectedException(tdce.getMessage());
			}
			
			// check intersection
			if(crosses){
				crossing=true;
				
				// if crossing is nearer than any previous one
				// log the triangle's details
				if(intDist[0]<nearest){
					nearest=intDist[0];             // set the arclength for future comparisons
					for(int j=0; j<D; j++){			// copy the normal over
						normal[j]=tempNormal[j];
					}
					d[0]=tempD[0];					// set the distance to interaction point
					toSkip= tri;                  // set the skipping triangle
	                p[0]=tri.getPermeability(0);    // set the permeability
				}
			}
			
			// increment triangle counter
		}

		// check intersection with substrate boundaries
		if(checkBoundaryIntersection(subsCoords, offset, stepVector, tempNormal, tempD, false, origLength, intDist, null, tempP)){
		    crossing=true;
            // just like before, if crossing is nearer than 
            // any previous one log the triangle's details
            if(intDist[0]<nearest){
                nearest=intDist[0];             // set the arclength for future comparisons
                for(int j=0; j<D; j++){         // copy the normal over
                    normal[j]=tempNormal[j];
                }
                d[0]=tempD[0];                  // set the distance to interaction point
                toSkip= null;                 // skipping triangle is null
                p[0]=1.0;                       // completely permeable
            }
		}
		
		
		if(!crossing){
			skipping=null;
		}
		
		if(toSkip!=null){
		    skipping=toSkip;
		}
		
		return crossing;
	}

	/**
	 * return peak coord, used as the centre of the object
	 */
	public double getPeakCoord() {
		
	    return L[0]/2.0;
	    
	}

	/**
	 * returns the size of the unit cell
	 */
	public double[] getSubstrateSize() {

	    return L;
	}

	/** 
	 * not used, does nothing.
	 */
	public void init() {

	}
	
	/**
	 * 6/11/2012: maps from the two versions don't agree. This is due to
	 * rounding differences in the way in which triangles are associated
	 * with bounding boxes, but this APPEARS not to be important. 
	 * 
	 * Voxel/Object maps from this code for hexagonal cylinder meshes
	 * do not leak, and the displacement distribution of spins on an 
	 * asparagus mesh with 500000 triangles APPEARS ok. Cautiously optimistic!
	 * 
	 * overrides the spatial optimisation initialiser in Substrate
	 * with a method optimised for triangles. 
	 * 
	 * This generates a warning due to the instantiation of an array 
	 * of generified ArrayLists without an explicit unchecked cast.
	 * This is because of some subtleties in Java's implementation
	 * of generics. Any solution to this would involve messing around
	 * with reflection and frankly is not worth it. It could also be
	 * dealt with via @suppressWarnings but this is only legal since
	 * Java 1.7 and consequently would cause some issues for anyone 
	 * using an earlier compliance level. In light of this, I've left 
	 * the code as it is and the warning in place.
	 * 
	 * Testing on large-scale meshes show a 200-fold speed increase over
	 * generic version. Not too shabby!
	 * 
	 * @param n number of subvoxels in x, y & z directions.
	 */
	public void initialiseSpatialOptimisation(int[] n){
		
        logger.info("Using triangle-optimised spatial optimisation algorithm");
        
        int numSubVoxels= initSpatialOptArrays(n);
        
        double[] xmin= new double[D];
        double[] xmax= new double[D];
        
        
		ArrayList<Triangle>[] intersectingObjects= new ArrayList[numSubVoxels]; 
        
        // loop over all objects
        for(int o=0; o<subsObj.length; o++){
        	
        	// in this case, we know that all objects are triangles
        	Triangle tri=(Triangle)subsObj[o];
        	
        	// initialise min and max arrays with zeroth vertex
        	double[] v= tri.getVertex(0);
        	for(int i=0; i<D; i++){
        		xmin[i]=v[i];
        		xmax[i]=v[i];
        	}
        	
        	// find min and max extent of triangle (minimal axis-aligned b-box)
        	for(int i=1; i<3; i++){
        		v=tri.getVertex(i);
        		for(int j=0; j<D; j++){
        			if(v[j]<xmin[j]){
        				xmin[j]=v[j];
        			}
        			if(v[j]>xmax[j]){
        				xmax[j]=v[j];
        			}
        		}
        	}
        	
        	// get indices of min and max subvoxels
        	int[] nmin= new int[D];
        	int[] nmax= new int[D];
        	
        	boolean same=true;
        	
        	for(int i=0; i<D; i++){
        		nmin[i]=(int)Math.floor(xmin[i]/s[i]);
        		nmax[i]=(int)Math.floor(xmax[i]/s[i]);
        		
        		// mild hack! clamp cell ranges to avoid rounding error
        		// in extreme-valued vertices
        		if(nmin[i]==-1){
        			nmin[i]=0;
        		}
        		if(nmax[i]==n[i]){
        			nmax[i]=n[i]-1;
        		}
        		
        		if(nmin[i]!=nmax[i]){
        			same=false;
        		}
        	}
        	
        	// if min and max subvoxels are the same, this is the only voxel intersected
        	// and therefore contains the triangle completely. In this case there's no
        	// need to perform the intersection check, we know that this is the one and
        	// only, hence just add it to the map and continue.
        	if(same){
        		
        		int index= getSubVoxelIndex(nmin[0], nmin[1], nmin[2]);
    			
        		//if(index==40){
    			//	System.err.println("catch case");
    			//}
        		
        		if(intersectingObjects[index]==null){
        			intersectingObjects[index]= new ArrayList<Triangle>();
        		}
        		
        		intersectingObjects[index].add(tri);
        		
        		continue;
        	}
        	
        	// if the triangles isn't completely contained in one subvoxel, we need to
        	// do an explicit intersection check across all subvoxels in the region
        	// defined by nmin and nmax (NB- less-equal not less-than on loops)
        	for(int i=nmin[0]; i<=nmax[0]; i++){
        		for(int j=nmin[1]; j<=nmax[1]; j++){
        			for(int k= nmin[2]; k<=nmax[2]; k++){
        				
        				// get the linear index for the subvoxel indices
        				int index= getSubVoxelIndex(i, j, k);
        				
        				//if(index==40){
        				//	System.err.println("catch case");
        				//}
        				
        				// get the coords for the limiting corners of the subvoxel
                        double[] bottomLeft= getBottomLeft(i, j, k);
                        double[] topRight= getUpperRight(i, j, k);
                        
                        // perform intersection check
                        if(tri.intersectsCubicRegion(bottomLeft, topRight)){
                        	
                        	// is this the first intersection for this subvoxel?
                        	if(intersectingObjects[index]==null){
                    			intersectingObjects[index]= new ArrayList<Triangle>();
                    		}
                        	
                            intersectingObjects[index].add(tri);                            
                        }
        			}
        		}
        	}
        }
        
        // we've now checked all objects for subvoxel intersection and built lists
        // of intersecting triangles for each subvoxel. Now need to prse these into 
        // arrays of SubstrateObjects and make the actual map.
        for(int i=0; i<intersectingObjects.length; i++){
        	
        	// no intersecting objects? set map to null and continue
        	if(intersectingObjects[i]==null){
        		voxToObjects[i]=null;
        		continue;
        	}
        	
        	// otherwise copy the contents of the intersecting object ArrayList to the fixed-length array
        	voxToObjects[i]= new SubstrateObject[intersectingObjects[i].size()];
            for(int o=0; o<intersectingObjects[i].size(); o++){
                voxToObjects[i][o]= intersectingObjects[i].get(o);
            }
        }
        
        // finally, set the spatially optimised flag to true, and we're done.
        spatialOptInitialised=true;

	}
	
	
	
	
	/**
	 * performs the ray test to check if a point is inside the
	 * object or not. This test can cope with non-convex objects
	 * (like sharks) and is almost as simple as the dot product 
	 * test used for convex objects.
	 * 
	 * We project a ray in an arbitrary direction and count the
	 * number of intersections with the surface of the object.
	 * If this number if odd, the point is inside the object.
	 * 
	 * This is only well-defined for closed surfaces.
	 * any holes will cause it to get it wrong some of the time.
	 * 
	 * @param walker the walker whose position to test.
	 * 
	 * @return true if inside object, false if outside
	 */
	public boolean intracellular(Walker walker) {
		
	    if(convexHull){
    		// ArrayList of coplanar triangles to
    		// skip after intersections are found
    		ArrayList<Triangle> skips=new ArrayList<Triangle>();
    		
    		final double[] offset= new double[]{0.0, 0.0, 0.0};
    		final double[] subsCoords= new double[D];
    		
    		
    		// map walker's position into unit cell
    		//double[] subsCoords=getSubstrateCoords(walker.r, offset);
    		getSubstrateCoords(walker.r, offset, subsCoords);
    		
    		// space to store the geometric stuff
    		// even though we ignore it here, it still
    		// needs to go somewhere
    		double[] tempNorm= new double[D];
    		double[] tempD= new double[1];
    		double[] intDist= new double[1];
    		double[] p= new double[1];
    		
    		// number of intersections
    		int n=0;
    				
    		// triangle iterator
    		Iterator<Triangle> triIt=triangles.iterator();
    		
    		// loop over all triangles
    		while(triIt.hasNext()){
    			// get triangle
    			Triangle triangle=(Triangle)triIt.next();
    
    			if(skips.contains(triangle)){
    				continue;
    			}
    			
    			boolean crosses=false;
    			try{
    				crosses=triangle.crosses(subsCoords, ray, tempNorm, tempD, false, 2.0*size, intDist, null, p, walker.R);
    			}
    			catch(TooDamnCloseException tdce){
    				throw new LoggedException(tdce);
    			}
    			
    			// check for intersection with ray
    			if(crosses){
    				// if yes, add coplanar triangles to the skipping array
    				ArrayList<Triangle> newSkips=coplanarMap.get(triangle);
    				
    				if(newSkips!=null){
    					skips.addAll(newSkips);
    				}
    				
    				// increment the intersection counter
    				n++;
    			}
    			
    		}
    		
    		// if n is odd, we're inside
    		if(n%2==1){
    		    
    		    //System.err.println("closed surface mesh intracellular returning true. n= "+n+" ray= ("+ray[0]+","+ray[1]+","+ray[2]+")");
    			return true;
    		}
    		
    		// otherwise outside
            //System.err.println("closed surface mesh intracellular returning false. n= "+n+" ray= ("+ray[0]+","+ray[1]+","+ray[2]+")");
    		return false;
	    }
	    else{
	        // if not a convex hull there's no way to perform this test
	        // reliably in general, so we won't even try.
	        return false;
	    }
	}
	
	
	/**
	 * check is a given walker is in the voxel. if a walker
	 * is in the central region after being mapped into 
	 * substrate coords then this returns true, otherwise false
	 * and the walker won't be included in data synthesis
	 * 
	 * overrides the default voxel check from root substrate
	 * class.
	 * 
	 * @param pos walker pos in global coords
	 * 
	 * @return true if in central voxel, otherwise false
	 */
	public boolean voxelContains(double[] pos){
	    
	    final double[] nullOffset= new double[]{0.0, 0.0, 0.0};
	    
	    final double[] subsCoords= new double[D];
	    
	    //double[] subsCoords= getSubstrateCoords(pos, nullOffset);
	    getSubstrateCoords(pos, nullOffset, subsCoords);

	    // check if each of the coords is outside the central voxel
	    for(int i=0; i<D; i++){
	        if(subsCoords[i]<voxelMin[i]){
	            return false;
	        }
	        
	        if(subsCoords[i]>voxelMax[i]){
	            return false;
	        }
	    }
	    
	    // if we've got here, we're ok.
	    return true;
	}
	
	private static void testIntracellular(){
        int tmax=100;
        
        SimulableScheme scheme= null;
        
             
        SimulationParams simParams= new SimulationParams(1, 1000, 0.0, 
                SimulationParams.UNIFORM, SubstrateType.TRI_PLY_MESH, StepType.FIXEDLENGTH, 
                1.5, scheme);
        
        SimulationParams.sim_plyfile=new String("cube.ply");
        
        // construct substrate from file
        MeshSubstrate meshSubs= new MeshSubstrate(simParams);
        
        // fish out the triangles
        Collection<Triangle> mesh= meshSubs.triangles;
        
        // and let's have a look at them...
        System.err.println("read "+mesh.size()+" triangles");
        
        int i=0;
        Iterator<Triangle> meshIt= mesh.iterator();
        while(meshIt.hasNext()){
            Triangle triangle= (Triangle)meshIt.next();
            
            System.err.println(++i+" "+triangle);
        }
        
        // now let's have a look at the intracellular routine
        double midPt= meshSubs.getPeakCoord();
        
        double[] pos= new double[]{midPt, midPt, midPt};
        
        Walker walker= new Walker(pos);
        
        System.err.println("midpoint is "+midPt);
        
        System.err.println("intracellular(midpoint) is "+meshSubs.intracellular(walker));
        
        double size= meshSubs.size;
        
        System.err.println("mesh size is "+size);
        
        pos[0]+=1.1*size;
        
        walker = new Walker(pos);
        
        System.err.println("extra cellular point is ("+pos[0]+","+pos[1]+","+pos[2]+")");
        System.err.println("intracellular is "+meshSubs.intracellular(walker));
	    
	}
	
	
	/**
	 * returns the list of triangles. This method is used by the visualiser.
	 */
	public Collection<Triangle> getTriangles(){
		return triangles;
	}
	
	
	/** 
	 * test ply substrate reading and intra cellular routine
	 * 
	 * @param args
	 */
	public static void main(String[] args){
		
	}
	
}

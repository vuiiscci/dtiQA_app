package simulation.geometry;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.logging.Logger;

import simulation.DiffusionSimulation;
import simulation.geometry.elements.Triangle;

import misc.LoggedException;

public class PLYreader {

    /** logging object */
	private static Logger logger = Logger.getLogger("simulation.geometry.PLYreader");
	
	/** is the object a closed surface (specified in file header) */
	public static boolean closedSurface= false;
	
	/**
	 * during the read process we construct a list of triangles
	 * that are subdivisions of a single polygon (eg. the two
	 * triangles that make up a square face on a cube). These
	 * are used to prevent the same face being multiply counted
	 * when intersections are with the common border of the 
	 * triangles.
	 * 
	 * This is only useful when doing the interior/exterior check
	 * and is not referenced during the normal intersection check.
	 */
	private static HashMap<Triangle, ArrayList<Triangle>> coplanarMap=null;
	
	/**
	 * static method to read a named PLY format file, parse the 
	 * contents and return a <code>Collection</code> of 
	 * <code>Triangle</code> objects.
	 * 
	 * Currently this only works with ascii PLY files and will
	 * ignore anything other than vertex and connectivity 
	 * information.
	 * 
	 * 
	 * @param fname
	 * @return collection of triangles read from the file
	 */
	public static final Collection<Triangle> readPLYfile(String PLYname, double scale, double p){

		/** dimensionality of space */
		final int D= DiffusionSimulation.D;
		
		/** somewhere to put the triangles */
		Collection<Triangle> triangle;
		
		/** file reader capable of getting one line at a time */
		BufferedReader plyFile;
		
		/** string to store each line read for file */
		String line;
		
		/** string to store tokens in each line */
		String nextToken=null;
		
		/** tokeniser to break line into elements */
		StringTokenizer tokeniser;
		
		/** number of vertices in model */
		int numVertices=0;
		
		/** number of faces (in general not same as number of triangles) */
		int numFaces=0;
		
		/** array of vertex coordinates */
		double[][] vertex;
		
		
		// open file
		try{
			plyFile= new BufferedReader(new FileReader(PLYname));
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}
		
		// reader the first line
		try{
			// check the first line of the PLY file
			line= plyFile.readLine();
			
			if(!line.equalsIgnoreCase("ply")){
				throw new LoggedException(PLYname+" is not a PLY file");
			}
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}
		
		// read the header
		try{
			logger.info("reading PLY header");
			line=plyFile.readLine();
			while(line!=null){
				tokeniser=new StringTokenizer(line, " ");
				while(tokeniser.hasMoreTokens()){
					nextToken= tokeniser.nextToken();
					if(nextToken.equalsIgnoreCase("element")){
						nextToken= tokeniser.nextToken();
						if(nextToken.equalsIgnoreCase("vertex")){
							nextToken= tokeniser.nextToken();
							numVertices= Integer.parseInt(nextToken);
						}
						if(nextToken.equalsIgnoreCase("face")){
							nextToken= tokeniser.nextToken();
							numFaces= Integer.parseInt(nextToken);
						}
					}
					if(nextToken.equalsIgnoreCase("comment")){
					    nextToken= tokeniser.nextToken();
					    if(nextToken.equalsIgnoreCase("closed")){
					        nextToken= tokeniser.nextToken();
					        if(nextToken.equalsIgnoreCase("surface")){
					            logger.info("object forms a closed surface");
					            PLYreader.closedSurface= true;
					        }
					    }
					    
					}
					if(nextToken.equalsIgnoreCase("end_header")){
						logger.info("end of header");
						break;
					}
				}
				if(nextToken.equalsIgnoreCase("end_header")){
					break;
				}
				line=plyFile.readLine();
			}
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}
		
		logger.info("reading mesh data");
		
		// allocate vertex array
		vertex=new double[numVertices][D];
		
		// read vertices
		try{
			for(int i=0; i<numVertices; i++){
				line=plyFile.readLine();	
			
				tokeniser= new StringTokenizer(line, " ");
				for(int j=0; j<D; j++){
					nextToken=tokeniser.nextToken();
					vertex[i][j]=Double.parseDouble(nextToken)/scale;
				}
			}
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}
		
		triangle= new ArrayList<Triangle>(numFaces);
		
		// construct triangles
		int numTriangles=0;
		
		coplanarMap= new HashMap<Triangle, ArrayList<Triangle>>();
		
		
		try{
			for(int face=0; face<numFaces; face++){
				numTriangles++;
			
				line=plyFile.readLine();
				
				tokeniser=new StringTokenizer(line, " ");
				nextToken= tokeniser.nextToken();
				
				int numVerts= Integer.parseInt(nextToken);
				
				if(numVerts<3){
					logger.warning("found a polygon with "+numVerts+" vertices. skipping.");
					continue;
				}
				
				if(numVerts==3){
					// trianglular faces are easy...
					double[][] verts= new double[numVerts][D];

					for(int i=0; i<numVerts; i++){
						nextToken=tokeniser.nextToken();
						int ind= Integer.parseInt(nextToken);
						
						for(int j=0; j<D; j++){
							verts[i][j]= vertex[ind][j];
						}
					}
					triangle.add(new Triangle(verts[0], verts[1], verts[2], p));
				}
				else if(numVerts==4){
					// quad faces can be split into two triangles
					int[] ind= new int[numVerts];
					
					// parse all the indices
					for(int i=0; i<numVerts; i++){
						nextToken=tokeniser.nextToken();
						ind[i]= Integer.parseInt(nextToken);
					}
					
					// increment the numTriangles counter again
					numTriangles++;
					
					double[][] verts= new double[numVerts][D];
					
					for(int i=0; i<numVerts; i++){
						for(int j=0; j<D; j++){
							verts[i][j]=vertex[ind[i]][j];
						}
					}
					
					// add triangles to collection
					Triangle tri1= new Triangle(verts[0], verts[1], verts[2], p);
					Triangle tri2= new Triangle(verts[0], verts[2], verts[3], p);
					triangle.add(tri1);
					triangle.add(tri2);
					
					// assemble coplanar lists
					ArrayList<Triangle> copTri1= new ArrayList<Triangle>();
					ArrayList<Triangle> copTri2= new ArrayList<Triangle>();
					
					// add other triangle to lists
					copTri1.add(tri2);
					copTri2.add(tri1);
					
					// add triangle-list pairs to coplanar map
					coplanarMap.put(tri1, copTri1);
					coplanarMap.put(tri2, copTri2);
				}
				else{
					// general convex polygon. star-triangulate.
					int[] ind= new int[numVerts];
					
					// parse indices
					for(int i=0; i<numVerts; i++){
						nextToken=tokeniser.nextToken();
						ind[i]=Integer.parseInt(nextToken);
					}
					
					// find the centre of the polygon (algebraic mean of vertices)
					double[] centre= new double[]{0.0, 0.0, 0.0};
					for(int i=0; i<numVerts; i++){
						for(int j=0; j<D; j++){
							centre[j]+=vertex[ind[i]][j];
						}
					}
					for(int j=0; j<D; j++){
						centre[j]/=numVerts;
					}
					
					numTriangles+=numVerts-1;
					
					// compile vertex list
					double[][] verts= new double[numVerts][D];
					for(int i=numVerts; i>0; i--){
						for(int j=0; j<D; j++){
							verts[i-1][j]=vertex[ind[i-1]][j];
						}
					}
					
					ArrayList<Triangle> subTris= new ArrayList<Triangle>();
					
					// each triangle goes (vert[i], vert[(i+1)%numVerts], centre)
					for(int i=numVerts; i>0; i--){
						// construct triangle
						Triangle newTri=new Triangle(verts[i-1], verts[i%numVerts], centre, p);
						
						// add to set of all triangles
						triangle.add(newTri);
						
						// add to set of subdivisions in this polygon
						subTris.add(newTri);
					}
					
					// make coplanar lists
					Iterator triIt= subTris.iterator();
					while(triIt.hasNext()){							// loop over all triangles in polygon
						Triangle thisTri=(Triangle)triIt.next();    // get current triangle
						ArrayList<Triangle> copList= new ArrayList<Triangle>();		// list of all other triangles

						Iterator otherTriIt=subTris.iterator();			// loop over triangles again
						while(otherTriIt.hasNext()){                        // include all triangles except current one
							Triangle thatTri=(Triangle)otherTriIt.next();
							if(thatTri==thisTri){
								continue;
							}
							copList.add(thatTri);
						}
						
						// add triangle-coplanar list pair to map
						coplanarMap.put(thisTri, copList);
					}
				}		
			}
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}

		logger.info("read and constructed "+numFaces+" faces using "+numTriangles+" triangles");

		try{
			plyFile.close();
		}
		catch(IOException ioe){
			throw new LoggedException(ioe);
		}
		
		return triangle;
	}
	
	/**
	 * returns the map of triangle it triangle indices of triangle
	 * subdivisions that should be ignored in an interior check if the
	 * original triangle is intersected by the ray.
	 * 
	 * note that it is assumed that subdivided polygons are coplanar.
	 * 
	 * @return map of Triangle to ArrayList of triangle indices of other 
	 * 			triangles to be ignored.
	 * 
	 * @throws LoggedException if the map has not been assembled prior
	 *                         to calling this function, a LoggedException
	 *                         is thrown.
	 */
	public static final HashMap<Triangle, ArrayList<Triangle>> getCoplanarListMap() throws LoggedException{
		
		if(coplanarMap==null){
			throw new LoggedException("coplanar triangle map has not been initialised before being fetched");
		}

		return coplanarMap;
	
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// test ply reader with a simple file
		String fname= new String("shark.ply");
		
		// scale factor for large objects 
		double scale= 1.0;
		
		// test file reading
		System.err.print("reading "+fname+"... ");
		Collection<Triangle> triangles= PLYreader.readPLYfile(fname, scale, 0.0);
		System.err.println("done\n");
		System.err.flush();
		
		// leaf through the list
		Iterator triIt= triangles.iterator();
		int i=0;
		
		System.err.println("read "+triangles.size()+" triangles:");
		
		while(triIt.hasNext()){
			
			Triangle triangle= (Triangle)triIt.next();
			System.err.println((i++)+" "+triangle);
			
		}
		
		System.err.println("all done");
	}

}

/**
 * The elements sub-package contains everything that is
 * used to make up a substrate. This includes geometric 
 * primitives like triangles and cylinders but also
 * potentially compound objects like a triangle mesh
 * or sub-mesh.
 * 
 * Objects designed to be queried for intersection or other
 * dynamical properties implements the SubstrateObject
 * interface, others that are used for different other reasons
 * do not implement the interface. An example of this second 
 * class of object is Quad, which is used as part of bounding
 * boxes rather than an actual substrate element.
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 * 
 */
package simulation.geometry.elements;

import java.util.ArrayList;

import simulation.dynamics.exceptions.TooDamnCloseException;

/**
 * substrate objects are sub-units of a substrate geometry
 * such as a cylinder, a triangle or a complete mesh. A substrate object
 * can potentially contain quite a lot of internal structure,
 * including other substrate objects invisible to the outside
 * substrate and any spatial partitioning system it needs.
 * 
 * The key thing is that is knows how to check for barrier
 * collisions and is able to return material properties at
 * various locations inside it's volume.
 * 
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public interface SubstrateObject {

	/**
	 * checks if a given step will cross a membrane in the object
	 * 
	 * @param walkerPos position of walker in space
	 * @param step step it wants to make
	 * @param normal space to return surface normal
	 * @param d space to return distance to interaction point in normal dir
	 * @param skipCurrent should we skip the barrier the walker is currently sitting on?
	 * @param origLength original length of step
	 * @param intDist proportion of step before interaction point
	 * @param in are we initially inside or outside the object?
	 * 
	 * @return true or false. if true, geometric quatities are provided
	 * 
	 * TODO: the permiability for the membrane should also be returned
	 *       when crossing
	 */
	public boolean crosses(double[] walkerPos, double[] step,
			double[] normal, double[] d, boolean skipCurrent, double origLength, 
			double[] intDist, boolean[] in, double[] p, double walkerRad) throws TooDamnCloseException;
	
	
	/**
	 * returns the diffusivity at the given location. this
	 * allows for arbitrary spatial variation.
	 * 
	 * @param subsCoords coords mapped onto substrate
	 * @return diffusivity
	 */
	public double getDiffusivityAt(double[] subsCoords);
	
	/**
	 * returns the oriented bounding box for the object.
	 * This is mostly ised during initialisation of the 
	 * partitioning tree.
	 * 
	 * @return bounding box optimally oriented bounding box
	 */
	public BoundingBox getBoundingBox();
	
	/**
	 * checks if a step will intersect the bounding box for the 
	 * object.
	 * 
	 * @param subsCoords
	 * @param step
	 * @return true or false. nothing else needed until the main check
	 */
	public boolean boundingBoxIntersects(double[] subsCoords, double[] step);
	
	/** checks if a point is inside or outside the object 
	 * 
	 * @param location in substrate coords
	 * 
	 * @return true or false
	 */
	public boolean inside(double[] subsCoords);
	
	/** 
	 * checks if the object intersects a cubic region. Note that this is different
	 * from the boundingBoxIntersects() method (@see boundingBoxIntersects(double[], double[]))
	 * in that instead of checking if a give step intersects the bounding box for this object,
	 * it check if the object itself intersets an axis-aligned cubic region.
	 * 
	 * @param bottomLeft
	 * @param topRight
	 * 
	 * @return true or false
	 */
	public boolean intersectsCubicRegion(double[] bottomLeft, double[] topRight);
	
	/**
	 * returns the shortest distance from the object to a given point
	 * 
	 * @param r0 the position to check
	 */
	public double getDistanceFrom(double[] r0);

	
	/**
	 * @return position in space
	 */
	public double[] getPosition();
	
	
	/**
	 * returns a set of Triangle objects so that the object can be turned into a mesh
	 * or visualised easily.
	 * 
	 * @return ArrayList<Triangle> of triangles for object. 
	 */
	public ArrayList<Triangle> getTriangles();
}

package tractography;

import java.io.*;

import numerics.*;

/**
 * <dl>
 * <dt>Purpose: To hold a collection of <code>Tract</code> objects.
 * <BR><BR>
 * 
 * <dt>Description:
 * <dd> This class is a container for a variable number of <code>Tract</code> objects. 
 * Bounds checking is performed as <code>Tract</code>'s are added, and the capacity of the 
 * <code>TractCollection</code> will increase as necessary.
 *
 * </dl>
 *
 * @version $Id$
 * @author  Philip Cook
 * @see tractography.Tract
 * 
 */
public final class TractCollection {

    /** Number of tracts currently in this collection. */
    private int numberOfTracts; 
    
    /** Array of Tracts. */
    private Tract[] tracts; 
    
    /** Number of Tracts this object can hold. */
    private int capacity; 
    
    /** Percentage by which to increase capacity when array is full. */
    private double growBy; 

    /** Total number of points in this collection. */
    private int totalPoints;
    

    /** Increase capacity of this TractCollection. */
    private void increaseCapacity() {
	
	capacity += (int)((double)capacity + (double)capacity * growBy / 100.0);
	
	Tract[] newTracts = new Tract[capacity];
	
	for (int i = 0; i < numberOfTracts; i++) {
	    newTracts[i] = tracts[i];
	}
	
	tracts = newTracts;
	
    }

   
    /** Default constructor. Capacity is 100 <code>Tract</code>'s, grow factor is 100%. */
    public TractCollection() {
	growBy = 100.0;
	capacity = 100;
	tracts = new Tract[capacity];
	numberOfTracts = 0;
	totalPoints = 0;
    }

    /** Construct a TractCollection with a specified capacity and growth factor.
     * @param initialCapacity the initial number of <code>Tract</code>'s the collection can hold.
     * @param growFactor the percentage to increase capacity by when this 
     * <code>TractCollection</code> is full.
     */
    public TractCollection(int initialCapacity, double growFactor) {
	growBy = growFactor;
	capacity = initialCapacity;
	tracts = new Tract[capacity];
	numberOfTracts = 0;
	totalPoints = 0;
    }

    
    /** @return the number of tracts in this collection. */
    public int numberOfTracts() {
	return numberOfTracts;
    }


    /** @return the number of points in this collection. */
    public int totalPoints() {
	return totalPoints;
    }
    

    /** 
     * Adds a tract to this collection.
     */
    public void addTract(Tract t) {
	tracts[numberOfTracts] = t;
	numberOfTracts++;
	totalPoints += t.numberOfPoints();
	if (numberOfTracts == capacity) {
	    increaseCapacity();
	}
    }

    
    /** 
     * Adds all tracts from another <code>TractCollection</code> to this <code>TractCollection</code>.
     */
    public void addTractCollection(TractCollection tc) {

	if (tc.numberOfTracts + numberOfTracts >= capacity) {
	    capacity = tc.numberOfTracts + numberOfTracts + 1;
	    
	
	    Tract[] newTracts = new Tract[capacity];
	    
	    for (int i = 0; i < numberOfTracts; i++) {
		newTracts[i] = tracts[i];
	    }
	
	    int oldNumberOfTracts = numberOfTracts;
	    
	    for (int i = 0; i < tc.numberOfTracts; i++) {
		newTracts[oldNumberOfTracts + i] = tc.tracts[i];
	    }
	    
	    numberOfTracts = oldNumberOfTracts + tc.numberOfTracts;
	
	    tracts = newTracts;
	
	}
	else {

	    for (int i = 0; i < tc.numberOfTracts; i++) {
		tracts[numberOfTracts + i] = tc.tracts[i];
	    }
	    
	    numberOfTracts += tc.numberOfTracts;

	}
	
	
    }
    

    /**
     * @return A string representation of this collection, of the form:
     * <p>
     * <code>
     * Class name
     * Total Number Of Points
     * Number of Tracts
     * </code>
     * @see #print()
     */
    public String toString() {
	return "tractography.TractCollection\nTotal points: " + 
	    totalPoints + "\nTotal Tracts: " + numberOfTracts + "\n";
    }


    /** 
     * Get a Tract. For speed, this will not return a defensive copy. 
     * The client is trusted to not modify the tracts by adding points
     * @param n the Tract index.
     * @return the nth tract.
     */
    public Tract getTract(int n) {
	if (n < 0 || n > numberOfTracts) {
	    throw new IndexOutOfBoundsException("Illegal Tract index: " + n);
	}
	else {
	    return tracts[n];
	}
    }
    
   
    /** 
     * Replace a Tract in this collection. 
     * @param t the Tract to be inserted into this collection.
     * @param n the index of the Tract to be replaced.
     */
    public void replaceTract(Tract t, int n) {
	
	if (n < 0 || n > numberOfTracts) {
	    return;
	}
	else {
	    totalPoints -= tracts[n].numberOfPoints();
	    tracts[n] = t;
	    totalPoints += t.numberOfPoints();
	}
   
    }

 
    /** 
     * Remove small Tracts from a collection. Seeds in isotropic regions will return a 
     * <code>Tract</code> containing one point
     * to be returned. This method removes Tracts below a specified length from the collection.
     * @param original the TractCollection to be optimised.
     * @param threshold the minimum number of points that a Tract may contain.
     * @return a TractCollection containing all Tracts which contain as many or more points than 
     * <code>threshold</code>.
     */
    public static TractCollection removeShortTracts(TractCollection original, int threshold) {

	TractCollection optimised = new TractCollection(original.numberOfTracts + 1, 10.0);

	for (int i = 0; i < original.numberOfTracts; i++) {
	    if (original.getTract(i).numberOfPoints() >= threshold) {
		optimised.addTract( original.getTract(i) );
	    }
	}

	return optimised;

    }

    
    /**
     * Writes tracts so that they can be read later with <code>readTracts</code>.
     *
     */
    public void writeRawTracts(DataOutputStream dout) throws IOException {
	
	for (int fibre = 0; fibre < numberOfTracts(); fibre++) {
	    
	    tracts[fibre].writeRaw(dout);
	   
	}
	
	
    }

}

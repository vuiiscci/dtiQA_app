package tractography;


/**
 * Interface for a collection of distinct regions of interest. 
 *
 * @author Philip Cook
 * @version $Id$
 */
public interface ROI_Collection {


    /**
     * @return the number of regions in this collection.
     * 
     */
    public int numberOfRegions();


    /**
     * Get a specific ROI.
     *
     * @param index the index of the required region, 
     * where <code>-1 < index < numberOfRegions() </code>.
     */
    public RegionOfInterest getRegion(int index);


    /**
     * Get all regions in an array.
     *
     */
    public RegionOfInterest[] getAllRegions();
    
}

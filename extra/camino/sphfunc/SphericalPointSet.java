package sphfunc;

import misc.SphericalPoints;
import tools.*;

/**
 * <dl>
 * <dt>Purpose: Implements a storage class for sets of points on the sphere
 * that can read arrays out of the files downloaded from Hardin Sloane and
 * Smith.
 * <dd>
 * 
 * <dt>Description: All html is removed from the files and the number of points
 * is added as the first line. The points are then stored as x y z coordinates
 * on consecutive lines.
 * <dd>
 * 
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id: SphericalPointSet.java,v 1.3 2005/08/18 11:16:06
 *         ucacmgh Exp $
 *  
 */
public class SphericalPointSet {

    
    /**
     * The data array containing the points is publicly available, it has
     * dimensions [number of points][3] and contains the x y and z coordinates
     * of each point.
     */
    public double[][] data;

    /**
     * Constructor reads the data from the file and stores it in data.
     */
    public SphericalPointSet(String filename) {

        FileInput f = new FileInput(filename);
        int noPoints = f.readInteger();

        data = new double[noPoints][3];

        for (int i = 0; i < noPoints; i++) {
            for (int j = 0; j < 3; j++) {
                data[i][j] = f.readDouble();
            }
        }

        f.close();
    }

}

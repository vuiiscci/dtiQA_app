package sphfunc;

import misc.*;

/**
 * <dl>
 * <dt>Purpose: Implements a simple database of icosahedral spherical codes of
 * increasing degree.
 * <dd>
 * 
 * <dt>Description: The arrays of points are stored in text files and only read
 * in when required. Then they are stored in memory.
 * 
 * <dd>
 * 
 * 
 * </dl>
 * 
 * @author Danny Alexander
 *  
 * @version $Id$
 */
public class ISCodes {

    /**
     * Directory containing the point set files. 
     */
    private static String directory = "/ISCodes/";


    /**
     * Contains the names of the files containing the point sets.
     */
    private static String[] filenames = { "isc01082.txt", "isc01922.txt", "isc03002.txt",
            "isc04322.txt", "isc05882.txt", "isc08672.txt", "isc12002.txt",
            "isc15872.txt", "isc20282.txt", "isc27002.txt", "isc32672.txt",
            "isc38882.txt", "isc48002.txt", "isc55472.txt", "isc63482.txt",
            "isc72032.txt", "isc78032.txt" };

    /**
     * Contains the names of the files containing the point sets.
     */
    private static String[] filenamesForMaxEnt = { "isc01082.txt", "isc01922.txt",
            "isc04322.txt", "isc08672.txt", "isc15872.txt", "isc32672.txt" };

    /**
     * Array of references to point set objects.
     */
    private static SphericalPointSet[] ptSets = new SphericalPointSet[filenames.length];

    /**
     * Array of references to point set objects used in MEFitterTest.
     */
    private static SphericalPointSet[] ptSetsForME = new SphericalPointSet[filenamesForMaxEnt.length];

    /**
     * Returns the i-th point set.
     */
    public static SphericalPointSet getPointSet(int i) throws ISCodesException {

        if (i < 0 || i > filenames.length) {
            throw new ISCodesException("Point set out of range.");
        }

        if (ptSets[i] == null) {
            java.net.URL url = ISCodes.class.getResource(directory + filenames[i]);
            if(url == null) {
                throw new LoggedException("Could not find file: " + directory + filenames[i]);
            }
                
            ptSets[i] = new SphericalPointSet(url.getFile());
        }

        return ptSets[i];
    }

    /**
     * Returns the i-th point set in a smaller set of pointsets chosen for doing
     * integration in <code>MEFitterTest</code>.
     */
    public static SphericalPointSet getPointSetForMaxEnt(int i) throws ISCodesException {
        if (i < 0 || i > filenamesForMaxEnt.length) {
            throw new ISCodesException("Max. Ent. point set out of range.");
        }

        if (ptSetsForME[i] == null) {
            java.net.URL url = ISCodes.class.getResource(directory
                    + filenamesForMaxEnt[i]);
            if(url == null) {
                throw new LoggedException("Could not find file: " + directory + filenamesForMaxEnt[i]);
            }
            ptSetsForME[i] = new SphericalPointSet(url.getFile());
        }

        return ptSetsForME[i];
    }

    /**
     * Returns the number of available point sets.
     */
    public static int getNoPointSets() {
        return filenames.length;
    }

    /**
     * Returns the number of available point sets for <code>MEFitterTest</code>.
     */
    public static int getNoPointSetsForMaxEnt() {
        return filenamesForMaxEnt.length;
    }

}


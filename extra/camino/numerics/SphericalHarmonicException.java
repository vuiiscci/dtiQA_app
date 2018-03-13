package numerics;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Exception class used by the SphericalHarmonics class.
 * 
 * <dt>Description:
 * 
 * <dd>Standard extension of the Exception class.
 * 
 * </dl>
 * 
 * @version $Id: SphericalHarmonicException.java,v 1.3 2005/08/18 11:12:21
 *          ucacmgh Exp $
 * @author Danny Alexander
 *  
 */
public class SphericalHarmonicException extends Exception {

    SphericalHarmonicException() {
        super();
    }

    SphericalHarmonicException(String s) {
        super(s);
    }
}


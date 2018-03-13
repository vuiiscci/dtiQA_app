package numerics;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Exception class used by the IncompleteBeta class.
 * 
 * <dt>Description:
 * 
 * <dd>Standard extension of the Exception class.
 * 
 * </dl>
 * 
 * @version $Id: SphericalHarmonicException.java,v 2.1 2004/03/16 12:12:44 ucacd
 *          xa Exp $
 * @author Danny Alexander
 *  
 */
public class IncompleteBetaException extends Exception {
    IncompleteBetaException() {
        super();
    }

    IncompleteBetaException(String s) {
        super(s);
    }
}


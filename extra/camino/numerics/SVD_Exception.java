package numerics;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Exception class used in the SVD method of RealMatrix.
 * 
 * <dt>Description:
 * 
 * <dd>Standard extension of the Exception class.
 * 
 * </dl>
 * 
 * @version $Id$
 * @author Danny Alexander
 *  
 */
public class SVD_Exception extends Exception {

    SVD_Exception() {
        super();
    }

    public SVD_Exception(String s) {
        super(s);
    }
}


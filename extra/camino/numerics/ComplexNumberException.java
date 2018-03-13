package numerics;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Exception for Complex class.
 * 
 * <dt>Description:
 * 
 * <dd>Standard extension of Exception class.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: ComplexNumberException.java,v 1.3 2005/08/18 11:12:22 ucacmgh
 *          Exp $
 *  
 */
public class ComplexNumberException extends Exception {
    ComplexNumberException() {
        super();
    }

    ComplexNumberException(String s) {
        super(s);
    }
}


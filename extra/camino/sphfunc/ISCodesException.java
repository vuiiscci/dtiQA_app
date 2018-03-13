package sphfunc;

/**
 * <dl>
 * <dt>Purpose: Exception class used in ISCodes database.
 * <dd>
 * 
 * <dt>Description:
 * <dd>Exception class used in ISCodes database.
 * 
 * </dl>
 * 
 * $Id$
 * 
 * @author Danny Alexander
 *  
 */
public class ISCodesException extends Exception {
    
    ISCodesException() {
        super();
    }

    ISCodesException(String s) {
        super(s);
    }
}


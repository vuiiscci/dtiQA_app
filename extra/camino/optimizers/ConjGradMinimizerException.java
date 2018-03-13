package optimizers;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Exception class used in ConjGradMinimizer.
 * 
 * <dt>Description:
 * 
 * <dd>Standard extension of the Exception class.
 * 
 * </dl>
 * 
 * @version $Id: ConjGradMinimizerException.java,v 1.3 2005/08/18 11:13:57
 *          ucacmgh Exp $
 * @author Danny Alexander
 *  
 */
public class ConjGradMinimizerException extends Exception {

    ConjGradMinimizerException() {
        super();
    }

    ConjGradMinimizerException(String s) {
        super(s);
    }
}


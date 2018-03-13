package optimizers;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Exception class used in Marquardt Minimiser.
 * 
 * <dt>Description:
 * 
 * <dd>Specifically for going over the maximum number of iterations.
 * 
 * </dl>
 * 
 * @version $Id: MarquardtMinimiserException.java,v 1.3 2005/08/18 11:13:57
 *          ucacmgh Exp $
 * @author Danny Alexander
 *  
 */
public class MarquardtMinimiserNonConvergenceException extends MarquardtMinimiserException {

    public MarquardtMinimiserNonConvergenceException() {
        super();
    }

    public MarquardtMinimiserNonConvergenceException(String s) {
        super(s);
    }
}


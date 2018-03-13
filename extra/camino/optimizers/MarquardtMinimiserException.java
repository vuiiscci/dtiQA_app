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
 * <dd>Standard extension of the Exception class.
 * 
 * </dl>
 * 
 * @version $Id: MarquardtMinimiserException.java,v 1.3 2005/08/18 11:13:57
 *          ucacmgh Exp $
 * @author Danny Alexander
 *  
 */
public class MarquardtMinimiserException extends MinimizerException {

    public MarquardtMinimiserException() {
        super();
    }

    public MarquardtMinimiserException(String s) {
        super(s);
    }
}


package optimizers;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General Exception class for Minimizer objects.
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
public class MinimizerException extends Exception {

    public MinimizerException() {
        super();
    }

    public MinimizerException(String s) {
        super(s);
    }
}


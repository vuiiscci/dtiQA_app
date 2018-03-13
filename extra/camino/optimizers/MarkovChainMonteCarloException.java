package optimizers;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Exception class used in MCMC Minimiser.
 * 
 * <dt>Description:
 * 
 * <dd>Standard extension of the Exception class.
 * 
 * </dl>
 * 
 * @version $Id:$
 * @author Aidos Abzhanov
 *  
 */
public class MarkovChainMonteCarloException extends MinimizerException {
	
	/**
	 * Default constructor
	 */
    public MarkovChainMonteCarloException() {
        super();
    }
	
	/**
     * Constructor
     * 
     * @param params
     *            The exception value.
     */	
    public MarkovChainMonteCarloException(String s) {
        super(s);
    }
}

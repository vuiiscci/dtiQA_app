package numerics;

/**
 * A mathematical function that returns a single value from a single variable.
 *
 * @version $Id$
 * @author  Philip Cook
 *
 * 
 */
public abstract class SVAnalyticalFunction {  

    public final boolean differentiable;
    
    /**
     * @param differentiable true if a derivative is defined.
     */
    public SVAnalyticalFunction(boolean differentiable) {
	this.differentiable = differentiable;
    }

    /**
     * @return the value of the function for the given variable value.
     * 
     */ 
    public abstract double evaluate(double var);
    

    /**
     * @return the derivative of this function. If the derivative is not defined, 
     * this method should throw a java.lang.UnsupportedOperationException.
     *
     */
    public SVAnalyticalFunction differentiated() {
	throw new java.lang.UnsupportedOperationException("No derivative defined for equation " + this);
    }

    /**
     *
     * @return a LaTeX string describing this function
     */
    public abstract String toString();
    
    
}

package numerics;

/**
 * 
 * The exception thrown due to errors in calculating the error function
 * 
 * @author Matt Hall m.hall@cs.ucl.ac.uk
 * 
 * 
 *  
 */
public class ErrorFunctionException extends Exception {

    /**
     *  
     */
    public ErrorFunctionException() {
        super();
    }

    /**
     * @param message
     */
    public ErrorFunctionException(String message) {
        super(message);
    }

    /**
     * @param cause
     */
    public ErrorFunctionException(Throwable cause) {
        super(cause);
    }

    /**
     * @param message
     * @param cause
     */
    public ErrorFunctionException(String message, Throwable cause) {
        super(message, cause);
    }

}
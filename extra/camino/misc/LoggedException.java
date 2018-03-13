package misc;

import java.util.logging.*;

/**
 * Superclass for a RuntimeException that is automatically logged at the SEVERE level 
 * upon construction. Also contains methods to log other exceptions. All log messages 
 * contain a full stack trace.
 *
 * @author Philip Cook 
 * @version $Id$
 *
 */
public class LoggedException extends RuntimeException {

    private static Logger logger = Logger.getLogger("camino.misc.LoggedException");

    public static boolean exceptionLogging = true;

    
    /**
     * Turns exception logging on or off. By default, exception logging is on.
     *
     * @param log if true, exception logging is on. If false, exceptions will not be logged, though
     * they will still kill the JVM if they are not caught. 
     *
     */
    public static void setExceptionLogging(boolean log) {
	exceptionLogging = log;
    }

    /**
     * Create an exception with an empty message.
     */
    public LoggedException() {
        super();
	if (exceptionLogging) {
	    logExceptionSevere(this, Thread.currentThread().getName());
	}
    }

    /**
     * @param cause the cause of this exception. The message is <code>cause.toString()</code>.
     */
    public LoggedException(Throwable cause) {
      super(cause);
        
      if (exceptionLogging) {
	  logExceptionSevere(this, Thread.currentThread().getName());
      }
        
    }

    /**
     * @param cause the cause of this exception.
     * @param message the cause is not added to the message by this constructor. 
     */
    public LoggedException(String message, Throwable cause) {
      super(message, cause);
        
      	if (exceptionLogging) {
	    logExceptionSevere(this, Thread.currentThread().getName());
	}

    }
    
    /**
     * @param message a message for this exception.
     */
    public LoggedException(String message) {
        super(message);
        
      	if (exceptionLogging) {
	    logExceptionSevere(this, Thread.currentThread().getName());
	}
    }


    /**
     * Logs an exception at the SEVERE level, also logs an additional message.
     *
     */
    public static void logExceptionSevere(Exception e, String message, String threadName) {
        logException(e, threadName, Level.SEVERE);
        logger.log(Level.SEVERE, message);
    }



    /**
     * Logs an exception at the WARNING level, also logs an additional message.
     *
     */
    public static void logExceptionWarning(Exception e, String message, String threadName) {
        logException(e, threadName, Level.WARNING);     
        logger.log(Level.WARNING, message);

    }


    /**
     * Logs an exception at the SEVERE level.
     *
     */
    public static void logExceptionSevere(Exception e, String threadName) {
        logException(e, threadName, Level.SEVERE);
    }

    /**
     * Logs an exception at the WARNING level.
     *
     */
    public static void logExceptionWarning(Exception e, String threadName) {
        logException(e, threadName, Level.WARNING);
    }


    /**
     * Logs an exception.
     *
     * @param level the <code>Level</code> to log the exception at.
     * @see java.util.logging.Level.
     */
    public static void logException(Exception e, String threadName, Level level) {

	// thread and exception message 
        String stackString = "Exception in thread \"" + 
            threadName + "\" " + e.getClass() + ": " + e.getMessage() + "\n";
        
        StackTraceElement[] stackTrace = e.getStackTrace();
	// stack elements
        for (int i = 0; i < stackTrace.length; i++) {
            stackString += "\tat " + stackTrace[i] + "\n";
        }

        logger.log(level, stackString);
        
    }


}

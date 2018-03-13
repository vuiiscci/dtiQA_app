package simulation.dynamics.exceptions;

import java.util.logging.Level;

import misc.LoggedException;

/**
 * Exception thrown when pathological steps are encountered by
 * simulation update. These are logged at level FINEST to avoid them getting
 * into the default logs but can be captured if need be by lowering the log
 * level in the logger.properties file.
 * 
 * This doesn't inherit from LoggedException because it is checked and 
 * should be caught when thrown. Logged exceptions are unchecked, plus
 * always log at SEVERE. This logs at FINEST.
 * 
 * 
 * @author gmorgan
 *
 */
public class StepRejectedException extends Exception {

	public StepRejectedException(){
		if(LoggedException.exceptionLogging){
			LoggedException.logException(this, Thread.currentThread().getName(), Level.FINEST);
		}
	}
	
	public StepRejectedException(String message){	
		if(LoggedException.exceptionLogging){
			LoggedException.logException(this, Thread.currentThread().getName(), Level.FINEST);
		}
	}
	
	
}

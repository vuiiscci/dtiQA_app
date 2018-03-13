package simulation.dynamics.exceptions;

import java.util.logging.Level;

import misc.LoggedException;

/**
 * exception thrown when a given step does not cross a barrier but means that
 * the walker will be less than a walker radius from the surface. Approaching 
 * this close is rare, but leads to leaks due to inaacuracies in floating point
 * arithemtic.
 * 
 * @author matt hall (matt.hall@ucl.ac.uk)
 *
 */
public class TooDamnCloseException extends Exception {

	public TooDamnCloseException(){
		if(LoggedException.exceptionLogging){
			LoggedException.logException(this, Thread.currentThread().getName(), Level.FINEST);
		}
	}
	
	public TooDamnCloseException(String message){	
		if(LoggedException.exceptionLogging){
			LoggedException.logException(this, Thread.currentThread().getName(), Level.FINEST);
		}
	}
	
}

package apps;

import data.OutputManager;
import tools.CL_Initializer;

/** 
 * first stab at a base-level executable class. this is 
 * inherited by any camino command entrypoint
 *
 * @author matt, shahrum.
 *
 */
 public abstract class Executable{

//	public final OutputManager om;

 
	/**
	 * constructor taking commandline args
	 */
	 public Executable(String[] args){

		initDefaultVals();
		initOptions(args);

		//om = new OutputManager();
		
		initVariables();
	}
	

	/**
	 * Initialise the default values of class-level fields, as they would be at declaration.
	 * <p> 
         * Class variables must be initialized in this method, because we call initOptions(String[]) and initVariables() 
         * from this class constructor. Variable declarations in the subclass aren't initialized at this time.
         * 
         * Therefore, avoid class variable declarations like "private int x = 1;", do "private int x;" then set "x=1" 
         * in this method.
         *
	 */
	 public abstract void initDefaultVals();
	
	/**
	 * default commandline parsing and initialisation
	 */
	 public void initOptions(String[] args){
	 
		// Parse the command line arguments
	    CL_Initializer.CL_init(args);
	    CL_Initializer.checkParsing(args);
	}
		
	/**
	 * abstract initialisation method for commands
	 */
	 public abstract void initVariables();
	 
	 /**
	  * abstract method execution for commands
	  */
	 public abstract void execute(OutputManager om);
	 
}
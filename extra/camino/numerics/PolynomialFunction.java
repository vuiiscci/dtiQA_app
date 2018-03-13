package numerics;

import data.*;
import misc.*;

import java.io.*;
import java.util.ArrayList;

/**
 * Represents a polynomial function of any degree.
 * 
 * @author Philip Cook
 * @version $Id$
 */
public abstract class PolynomialFunction {


    public PolynomialFunction() {

    }

    
    /**
     * Evaluate the function. See subclasses for the order of the variables in the array
     * <code>vars</code>.
     *
     * @return the value of the function with the given parameters.
     *
     */
    public abstract double evaluate(double[] vars);



    /**
     * Reads a function from a file. 
     *
     * @throws LoggedException if the function cannot be read.
     */
    public static PolynomialFunction[] readFunctions(String file) {

        VoxelOrderDataSource data = new VoxelOrderDataSource(file, 1, "double");

        int numVariables = (int)data.nextVoxel()[0];

	if (numVariables == 1) {
	    return PolynomialFunction1D.readFunctions(file);
	}
	else if (numVariables == 2) {
	    return PolynomialFunction2D.readFunctions(file);
	}
	else {
	    throw new LoggedException("Cannot create function with " + numVariables + 
				      " variables");
	}


    }


    /**
     * @return the number of variables in this function.
     */
    public abstract int variables();


}

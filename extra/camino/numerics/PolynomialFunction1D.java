package numerics;

import data.*;
import misc.*;

import java.io.*;
import java.util.ArrayList;

/**
 * Represents a polynomial function of any degree with one variable {x}.
 * 
 * @author Philip Cook
 * @version $Id$
 */
public class PolynomialFunction1D extends PolynomialFunction {


    private final int order;

    private final double[] coeffs;
    


    /**
     * Construct a function with a specified order and associated coefficients.
     * @param coefficients in the order 1, x, x^2...x^maxOrd
     */
    public PolynomialFunction1D(int maxOrd, double[] coefficients) {

        order = maxOrd;

        coeffs = new double[coefficients.length];

        System.arraycopy(coefficients, 0, coeffs, 0, coefficients.length);


    }


    public double evaluate(double[] vars) {

        int counter = 0;

        double answer = 0.0;

        // inefficient use of Math.pow
        for (int i = 0; i <= order; i++) {
	    answer += Math.pow(vars[0],i) * coeffs[i];
        }

        return answer;
        
    }
    
       

    /**
     * Reads a function from a file. 
     *
     * @throws LoggedException if the function cannot be read.
     */
    public static PolynomialFunction1D[] readFunctions(String file) {

        VoxelOrderDataSource data = new VoxelOrderDataSource(file, 1, "double");

        int numVariables = (int)data.nextVoxel()[0];

	if (numVariables != 1) {
	    throw new LoggedException("Tried to create function with 1 variable, but " + 
				      "number of variables in function is " + numVariables);
	}

        int order = (int)data.nextVoxel()[0];

        int numCoeffs = order + 1;

        ArrayList<PolynomialFunction1D> list = new ArrayList<PolynomialFunction1D>();

        while (data.more()) {

            double[] coeffs = new double[numCoeffs];
            
            for (int i = 0; i < numCoeffs; i++) {
                coeffs[i] = data.nextVoxel()[0];
            }

            list.add(new PolynomialFunction1D(order, coeffs));
            
        }
               
        PolynomialFunction1D[] array = new PolynomialFunction1D[list.size()];

        list.toArray(array);
        
        return array;

    }



    public String toString() {

        String equation = "";
        int counter = 0;

        for (int i = 0; i <= order; i++) {
	    equation += coeffs[counter] + " * x^" + i + " + ";
        }

        // trim trailing plus
        equation = equation.substring(0, equation.length() - 3);

        return equation;
        
    }


    public int variables() {
	return 1;
    }


}

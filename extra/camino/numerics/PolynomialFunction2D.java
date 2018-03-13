package numerics;

import data.*;
import misc.*;

import java.io.*;
import java.util.ArrayList;

/**
 * Represents a polynomial function of any degree with two variables {x, y}.
 * 
 * @author Philip Cook
 * @version $Id$
 */
public class PolynomialFunction2D extends PolynomialFunction {


    private final int order;

    private final double[] coeffs;
    


    /**
     * Construct a function with a specified order and associated coefficients.
     * @param coefficients in the order 1, y, y^2...y^maxOrd, x, xy, xy^(maxOrd - 1)...x^maxOrd
     */
    public PolynomialFunction2D(int maxOrd, double[] coefficients) {

        order = maxOrd;

        coeffs = new double[coefficients.length];

        System.arraycopy(coefficients, 0, coeffs, 0, coefficients.length);


    }


    public double evaluate(double[] vars) {

        int counter = 0;

        double answer = 0.0;

        // inefficient use of Math.pow
        for (int i = 0; i <= order; i++) {
            for (int j = 0; j <= order; j++) {
                if (i + j <= order) {
                    answer += Math.pow(vars[0],i) * Math.pow(vars[1], j) * coeffs[counter];
                    counter++;
                }
            }
        }

        return answer;
        
    }
    
       

    /**
     * Reads a function from a file. 
     *
     * @throws LoggedException if the function cannot be read.
     */
    public static PolynomialFunction2D[] readFunctions(String file) {

        VoxelOrderDataSource data = new VoxelOrderDataSource(file, 1, "double");

        int numVariables = (int)data.nextVoxel()[0];

	if (numVariables != 2) {
	    throw new LoggedException("Tried to create function with 2 variables, but " + 
				      "number of variables in function is " + numVariables);
	}

        int order = (int)data.nextVoxel()[0];

        int numCoeffs = (order+1)*(order+2) / 2;

        ArrayList<PolynomialFunction2D> list = new ArrayList<PolynomialFunction2D>();

        while (data.more()) {

            double[] coeffs = new double[numCoeffs];
            
            for (int i = 0; i < numCoeffs; i++) {
                coeffs[i] = data.nextVoxel()[0];
            }

            list.add(new PolynomialFunction2D(order, coeffs));
            
        }
               
        PolynomialFunction2D[] array = new PolynomialFunction2D[list.size()];

        list.toArray(array);
        
        return array;

    }



    public String toString() {

        String equation = "";
        int counter = 0;

        for (int i = 0; i <= order; i++) {
            for (int j = 0; j <= order; j++) {
                if (i + j <= order) {
                    equation += coeffs[counter] + " * x^" + i + " * y^" + j + " + ";
                    counter++;
                }
            }
        }

        // trim trailing plus
        equation = equation.substring(0, equation.length() - 3);

        return equation;
        
    }


    public int variables() {
	return 2;
    }




}

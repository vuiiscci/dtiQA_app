package numerics;

/**
 * Solves nonlinear equations by the Newton-Raphson method.
 * <p> 
 * The method tries to evaluate a function of the form f(x) - constant = 0. 
 * The user supplies the constant and an initial guess of x, x_0. 
 * The point where the tangent to f(x_0), f\prime(x_0), crosses the x-axis is taken to be
 * closer to the solution. Therefore, x_1 = x_0 - \frac{f(x_0)}{f\prime(x_0)}. 
 * The process is continued until \mid\frac{ x_{i+1} - x_i }{ x_{i+1} }\mid \times 100 \leq maxEpsilon, 
 * where maxEpsilon is a user supplied parameter.
 * 
 * 
 * @version $Id$
 * @author  Philip Cook
 * @see sphereDistFit.SVAnalyticalFunction
 *
 * 
 */
public class NewtonRaphsonSolver {

    
    /**
     * Solve the function f(x) - c = 0.
     * @param f must be a differentiable function.
     * @param c the known value of f(x) for unknown x.
     * @param x0 the initial guess for x.
     * @param maxEpsilon the maximum acceptable error (see class comment for definition).
     *
     * @return x when the solution has converged to within epsilon
     *
     *@throws ConvergenceException if after 1000 iterations, convergence has not been reached.
     */
    public static double solve(SVAnalyticalFunction f, double c, double x0, double maxEpsilon) 
	throws ConvergenceException {
	
	
	SVAnalyticalFunction fPrime = f.differentiated();

	double epsilon = Double.MAX_VALUE;

	int counter = 0;

	double xi = x0;

	double nextX = 0.0;

	while (epsilon > maxEpsilon) {
	    
	    double fpxi = fPrime.evaluate(xi);

	    // x_{i+1} = x_{i} - \frac{ f(x_{i}) }{ f\prime(x_{i}) }. 	    
	    nextX = xi - (f.evaluate(xi) - c) / fpxi;

	    if (Double.isNaN(nextX)) {
		throw new ConvergenceException("xi is NaN, probable divergence");
	    }
            else if (Double.isInfinite(nextX)) {
		throw new ConvergenceException("xi is Infinite, probable divergence");
	    }

	    // \mid\frac{ x_{i+1} - x_i }{ x_{i+1} }\mid \times 100
	    epsilon = Math.abs( (nextX - xi) / (nextX) ) * 100.0;

	    if (counter++ > 9999) {
		throw new ConvergenceException("Couldn't converge -- try a different x0");
	    }

	    xi = nextX;

	}

	return xi;
	    
    }





}

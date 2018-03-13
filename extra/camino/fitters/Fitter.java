package fitters;

import optimizers.*;
import imaging.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>General base class for objects that fit a model to data.
 * Derived classes specify which model and how to fit it.
 * 
 * <dt>Description:
 * 
 * <dd> Specifically, derived classes implement a Codec that specifies
 * the mapping from model parameters to optimized parameters.  They
 * create a specific Minimizer object to do the fitting and implement
 * a routine to find a starting point for the optimization.
 * 
 * </dl>
 *
 * @author  Danny
 * @version $Id$
 *
 */
public abstract class Fitter {

    protected DW_Scheme scheme; 
    protected Codec cod; 
    protected Minimizer minimizer;


    /**
     * Fits the model to one voxel's worth of data.
     *
     * @param data The measurements to fit to, which must come from
     * the imaging protocol in the Scheme object passed to the
     * constructor.
     *
     * @return An array of candidate solutions returned by the
     * minimizer.  The array has size Lx(N+1) where N is the number of
     * model parameters and L is the number of candidate solutions the
     * minimizer computes.  Each row contains the N model parameter
     * values followed by the value of the objective function for that
     * combination of parameters.
     */
    public double[][] fit(double[] data) throws MinimizerException {
        minimizer.setMeasurements(data);
        double[] startPoint = getStartPoint(data);
        minimizer.setInitParams(cod.modelToOpt(startPoint));
        minimizer.minimise();
        return minimizer.getSolutions();
    }
    /**
     * Fits the model to one voxel's worth of data.
     *
     * @param data The measurements to fit to, which must come from
     * the imaging protocol in the Scheme object passed to the
     * constructor.
     * @param gradAdj an adjustment to thee gradient used in the scheme.
     *
     * @return An array of candidate solutions returned by the
     * minimizer.  The array has size Lx(N+1) where N is the number of
     * model parameters and L is the number of candidate solutions the
     * minimizer computes.  Each row contains the N model parameter
     * values followed by the value of the objective function for that
     * combination of parameters.
     */
    public double[][] fit(double[] data, double[] gradAdj) throws MinimizerException {
        minimizer.adjustScheme(gradAdj);
        return fit(data);
    }

    /**
     * Returns the number of parameters the fitter returns for
     * each voxel.  By default, this is the number of solutions
     * returned times the number of model parameters plus two
     * (for the exit code and value of the objective function).
     *
     * @return The number of values returned per voxel.
     */
    public int getNumValuesPerRun() {
        return (cod.getNumModelParams()+2)*minimizer.getNumSolutions();
    }


    /**
     * Returns the number of parameters the fitter returns for each
     * individual solution (there may be lots of these per run/voxel).
     * By default, this is the number of model parameters plus two
     * (for the exit code and value of the objective function).
     *
     * @return The number of values returned per solution.
     */
    public int getNumValuesPerSolution() {
        return (cod.getNumModelParams()+2);
    }


    /**
     * Returns the first index of the theta-phi pair of parameters
     * encoding a single fiber orientation.  Theta and phi are assumed
     * consecutive and in that order.  If the model contains no
     * orientation parameter, return 0, which is the default.  The
     * value accounts for the exit code placed in position 0.
     * Thus if vals is an array corresponding to one row of the matrix
     * returned by Fitter.fit, then theta = vals[getDirectionIndex()]
     * and phi = vals[getDirectionIndex()+1].
     *
     * @return The first index of the direction parameter.
     */
    public int getDirectionIndex() {
        return cod.getDirectionIndex() + 1;
    }

    
    /**
     * Determines a starting set of parameters from the data.
     *
     * @param data The set of measurements to fit to.
     *
     * @return The starting model parameter values.
     */
    protected abstract double[] getStartPoint(double[] data);


}


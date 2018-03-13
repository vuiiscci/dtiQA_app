package tractography;

import misc.*;
import numerics.*;
import tools.*;


/**
 * Computes a statistic from streamlines, either from the streamline itself or by the values of a 
 * scalar image along the streamline. The whole streamline may or a subset of it may be used in the
 * computation of the statistic. 
 * 
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class TractStatisticFilter {


    // dimensions of seed space
    private final int[] dataDims;
	
    // mm
    private final double[] voxelDims;

    private ScalarImage scalars = null;

    private String tractStat = "mean";

    private boolean interpolated = false;
    
    // true if we have scalars
    private boolean tractScalarStatistic;


    /**
     * Construct an filter that uses the scalar image to define the tract statistic.
     * 
     * @param scalars a ScalarImage. The interpolation scheme of the image may be changed
     * if interpolation is enabled for the filter. The image should be in the
     * same space as the tracts.
     *
     */
    public TractStatisticFilter(ScalarImage scalars) {
	this.scalars = scalars;
	dataDims = scalars.getDataDims();
	voxelDims = scalars.getVoxelDims();
	tractScalarStatistic = true;
    }


    /**
     * For Matlab.
     * @param data the scalar data.
     * @param voxelDims voxel dimensions in mm.
     */
    public TractStatisticFilter(double[][][] scalarData, double[] voxelDims) {
	this.scalars = new ScalarImage(scalarData, voxelDims);
	dataDims = scalars.getDataDims();
	this.voxelDims = scalars.getVoxelDims();
	tractScalarStatistic = true;
    }


    /**
     * Constructor for filters that do not use scalars.
     */
    public TractStatisticFilter(int[] dataDims, double[] voxelDims) {
	this.dataDims = dataDims;
	this.voxelDims = voxelDims;
	tractScalarStatistic = false;
    }
    

    /**
     * Sets the statistic to be computed for each tract.
     *
     * @param stat either: none, mean, min, max, median, var, meanvar (mean and variance) or length.
     */
    public void setTractStatistic(String stat) {
	tractStat = stat;
    }


    /**
     * If interpolation is enabled, the scalar image is
     * linearly interpolated.
     *
     */
    public void setInterpolate(boolean interp) {
	interpolated = interp;

	if (interpolated && scalars != null) {
	    scalars.setInterpolation("linear");
	}
    }


    /**
     * Process a tract.
     *
     * @return the statistic of interest for the tract.
     *
     */
    public double[] processTract(Tract t) {

	double[] stat = null;
    
        if (tractScalarStatistic) {

	    // values for each point
	    double[] values = null;

	    // weights for use in calculation of mean only
	    double[] weights = null;

	    Point3D[] points = t.getPoints();
	    
	    values = scalars.valuesAt(points);
	    
	    
	    if (tractStat.equals("none")) {
		stat = values;
	    }
	    else if (tractStat.equals("mean")) {
		if (weights != null) {
		    stat = new double[1];
		    stat[0] = ArrayOps.weightedMean(values, weights);
		}
		else {
		    stat = new double[1];
		    stat[0] = ArrayOps.mean(values);
		}
	    }
	    else if (tractStat.equals("max")) {
		stat = new double[1];
		stat[0] = ArrayOps.max(values);
	    }
	    else if (tractStat.equals("min")) {
		stat = new double[1];
		stat[0] = ArrayOps.min(values);
	    }
	    else if (tractStat.equals("median")) {
		stat = new double[1];
		stat[0] = ArrayOps.median(values);
	    }
	    else if (tractStat.equals("sum")) {
		stat = new double[1];
		stat[0] = ArrayOps.sum(values);
	    }
	    else if (tractStat.equals("var")) {

		stat = new double[1];

		if (weights != null) {
		    stat[0] = 
			ArrayOps.weightedVar(values, weights, ArrayOps.weightedMean(values, weights));
		}
		else {
		    stat[0] = ArrayOps.var(values, ArrayOps.mean(values));
		}
	    }
	    else if (tractStat.equals("meanvar")) {

		stat = new double[2];
		
		if (weights != null) {
		    
		    stat[0] = ArrayOps.weightedMean(values, weights);
		    
		    stat[1] = 
			ArrayOps.weightedVar(values, weights, stat[0]);
		}
		else {
		    stat[0] = ArrayOps.mean(values);
		    stat[1] = ArrayOps.var(values, stat[0]);
		}
	    }
	    else {
		throw new LoggedException("Unsupported tract statistic: " + tractStat);
	    }

	}
	else {
	    
	    // statistic does not depend on scalars

	    stat = new double[1];

	    if (tractStat.equals("length")) {
		stat[0] = t.length();
	    }
	    else if (tractStat.equals("endpointsep")) {
		stat[0] = t.endpointSeparation();
	    }
	    else {
		throw new LoggedException("Unsupported tract statistic " + 
					  "(may be supported with a scalar image): " + tractStat);
	    }


	}


	return stat;

    }



}

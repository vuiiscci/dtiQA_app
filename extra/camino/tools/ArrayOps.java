package tools;

import misc.*;

import java.util.Arrays;
import java.text.DecimalFormat;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Simple operations on arrays.
 * 
 * <dt>Description:
 * 
 * <dd>Implements various common operations on arrays.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class ArrayOps {

    public static final DecimalFormat decimalDF = new DecimalFormat("0.00000E00"); 
    public static final DecimalFormat integerDF = new DecimalFormat("0000"); 

    /**
     * Returns the largest value in the array.
     * 
     * @param arr
     *            The array.
     * 
     * @return The maximum
     */
    public static double max(double[] arr) {
        double max = arr[0];
        for (int i = 1; i < arr.length; i++) {
            max = (arr[i] > max) ? arr[i] : max;
        }

        return max;
    }

    
    /**
     * Returns the smallest value in the array.
     * 
     * @param arr
     *            The array.
     * 
     * @return The minimum
     */
    public static double min(double[] arr) {
        double min = arr[0];
        for (int i = 1; i < arr.length; i++) {
            min = (arr[i] < min) ? arr[i] : min;
        }

        return min;
    }


    /**
     * Returns the mean value in the array.
     * 
     * @param arr
     *            The array.
     * 
     * @return The mean, or 0.0 if the array has zero length.
     */   
    public static double mean(double[] arr) {

	if (arr.length == 0) {
	    return 0.0;
	}

	double sum = 0.0;

	for (int i = 0; i < arr.length; i++) {
	    sum += arr[i];
	}

	return sum / (double)(arr.length);

    }



    /**
     * Returns the sum of the values in the array.
     * 
     * @param arr
     *            The array.
     * 
     * @return The sum, or 0.0 if the array has zero length.
     */   
    public static double sum(double[] arr) {

	if (arr.length == 0) {
	    return 0.0;
	}

	double sum = 0.0;

	for (int i = 0; i < arr.length; i++) {
	    sum += arr[i];
	}

	return sum;

    }




    /**
     * Returns the mean value in the array.
     * 
     * @param arr
     *            The array.
     * @param weights
     *            An array of non-negative weights, the same length as <code>arr</code>.
     * 
     * @return The mean, or 0.0 if the array has zero length.
     */   
    public static double weightedMean(double[] arr, double[] weights) {

	if (arr.length == 0) {
	    return 0.0;
	}
	if (arr.length != weights.length) {

	    throw new LoggedException("Number of elements must match number of weights");
	}

	double weightedSum = 0.0;
	double norm = 0.0;

	for (int i = 0; i < arr.length; i++) {
	    weightedSum += arr[i] * weights[i];
	    norm += weights[i]; 
	    
	    if (weights[i] < 0.0) {

		throw new LoggedException("Weights must be non-negative");
	    }
	}

	if (norm == 0.0) {
	    return 0.0;
	}

	return weightedSum / norm;

    }



    /**
     * Returns the geometric mean value in the array.
     * 
     * @param arr
     *            The array.
     * 
     * @return The geometric mean, or 0.0 if the array has zero length. 
     */   
    public static double geoMean(double[] arr) {

	if (arr.length == 0) {
	    return 0.0;
	}

        double gm = 1.0;
        for (int i = 0; i < arr.length; i++) {
            gm *= arr[i];
        }
        gm = Math.pow(gm, 1.0 / (double) arr.length);
	
        return gm;
    }


    
    /**
     * Returns the median value in the array. <code>NaN</code> and <code>Infinity</code> are
     * handled in the same way as in the Java <code>Arrays.sort</code> method. If there are an 
     * even number of values in the array, the method returns the mean of the middle two values.
     * 
     * @param arr
     *            The array.
     * 
     * @return The median, or 0.0 if the array has zero length.
     */   
    public static double median(double[] arr) {

	if (arr.length == 0) {
	    return 0.0;
	}
	
	int length = arr.length;

	double[] copy = new double[length];

	System.arraycopy(arr, 0, copy, 0, length);

	java.util.Arrays.sort(copy);

	if (length % 2 == 0) {
	    return (copy[length / 2 - 1] + copy[length / 2]) / 2.0;
	}
	else {
	    return copy[length / 2];
	}
    }


    
    
    /**
     * Returns the variance of the values in the array.
     * 
     * @param arr
     *            The array.
     * 
     * @param mean the mean value.
     * 
     * @return The variance \sum_i (arr[i] - mean)^2 / (arr.length -1), or 0.0 if arr.length < 2.
     */  
    public static double var(double[] arr, double mean) {

	if (arr.length < 2) {
	    return 0.0;
	}
	
	double var = 0.0;

	for (int i = 0; i < arr.length; i++) {
	    var += (arr[i] - mean) * (arr[i] - mean);
	}
	
	var = var / (arr.length - 1.0);

	return var;
    }


    /**
     * Returns the weighted variance of the values in the array.
     * 
     * @param arr
     *            The array.
     * 
     * @param wMean the weighted mean value.
     * 
     * @return The variance \sum_i w_i (arr[i] - wMean)^2 / [(arr.length - 1) mean(w)] , or 0.0 if arr.length < 2.
     */  
    public static double weightedVar(double[] arr, double[] weights, double wMean) {

	if (arr.length < 2) {
	    return 0.0;
	}
	
	double var = 0.0;

	double meanWeight = 0.0;

	for (int i = 0; i < arr.length; i++) {
	    var += weights[i] * (arr[i] - wMean) * (arr[i] - wMean);

	    meanWeight += weights[i];

	    if (weights[i] < 0.0) {
		throw new LoggedException("Weights must be non-negative");
	    }

	}
	
	meanWeight = meanWeight / (double)arr.length;

	if (meanWeight == 0.0) {
	    return 0.0;
	}
	
	var = var / (meanWeight  * (arr.length - 1.0));

	return var;
    }


    /**
     * Returns the median absolute deviation of the array.
     *
     * @return MAD = median(abs(arr - median))
     */
    public static double mad(double[] arr, double median) {
	if (arr.length == 0) {
	    return 0.0;
	}
	
        double[] absDev = new double[arr.length];

        for (int i = 0; i < arr.length; i++) {
            absDev[i] = Math.abs(arr[i] - median);
        }

        return median(absDev);

    }
    

    /**
     * Converts an array to a string representation with a specified number of decimal places
     */
    public static String toString(double[] data) {
        
        StringBuffer buf = new StringBuffer();
        
        for (int i = 0; i < data.length; i++) {
            buf.append("\t[");
            buf.append(integerDF.format(i));
            buf.append("]\t");
            buf.append(decimalDF.format(data[i]));
            buf.append("\n");
        }

        return buf.toString();
    }

}


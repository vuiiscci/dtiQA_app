package numerics;

import java.lang.Math;
import java.lang.Cloneable;
import java.io.Serializable;

/**
 * Provides a symmetric matrix in lower triangular format. This is for memory efficient storage and it doesn't do any special math. 
 * RealMatrix should be used to do matrix manipulations.
 * 
 * @author  Philip Cook
 *  
 */
public class SymmetricMatrix  {


    /**
     * The size of the matrix
     */
    private final int size;
    
    
    /**
     * The array of matrix entries, where entries[i] is an array of length i+1. That is, you may 
     * access any element [i][j] where j &lt;= i. 
     */
    private final double[][] entries;


    /**
     * Make a symmetric matrix of size s.
     */
    public SymmetricMatrix(int s) {
        size = s;
        entries = new double[size][];

        for (int i = 0; i < size; i++) {
            entries[i] = new double[i+1];
        }
    }



    /**
     * Gets element [i][j] == [j][i].
     *
     */
    public double get(int i, int j) {
        
        if (  !( (i > -1) && (i < size) ) || !( (j > -1) && (j < size) )  ) {
            throw new IllegalArgumentException("Invalid matrix index " + i + ", " + j );
        }


        if (j > i) {
            int tmp = i;
            i = j;
            j = tmp;
        }

        return entries[i][j];
    }


    /**
     * Sets element [i][j] and [j][i].
     *
     */
    public void set(int i, int j, double value) {
        
        if (  !( (i > -1) && (i < size) ) || !( (j > -1) && (j < size) )  ) {
            throw new IllegalArgumentException("Invalid matrix index " + i + ", " + j );
        }


        if (j > i) {
            int tmp = i;
            i = j;
            j = tmp;
        }

        entries[i][j] = value;
    }


    /**
     * Add a value to a particular entry.
     *
     */
    public void add(int i, int j, double value) {
   
        if (  !( (i > -1) && (i < size) ) || !( (j > -1) && (j < size) )  ) {
            throw new IllegalArgumentException("Invalid matrix index " + i + ", " + j );
        }


        if (j > i) {
            int tmp = i;
            i = j;
            j = tmp;
        }

        entries[i][j] += value;

    }


    /**
     * Multiply a particular entry by a scalar.
     *
     */
    public void scale(int i, int j, double value) {
   
        if (  !( (i > -1) && (i < size) ) || !( (j > -1) && (j < size) )  ) {
            throw new IllegalArgumentException("Invalid matrix index " + i + ", " + j );
        }


        if (j > i) {
            int tmp = i;
            i = j;
            j = tmp;
        }
        
        entries[i][j] *= value;
    }


    /**
     * Multiply all entries by a scalar.
     *
     */
    public void scale(double value) {
   
        for (int i = 0; i < size; i++) {
            for (int j = 0; j <= i; j++) {
                entries[i][j] *= value;
            }
        }
        
    }
    

    /**
     * @return a RealMatrix version of this object.
     */
    public RealMatrix toRealMatrix() {

        RealMatrix r = new RealMatrix(size, size);

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < i; j++) {
                r.entries[i][j] = entries[i][j];
                r.entries[j][i] = entries[i][j];
            }
            r.entries[i][i] = entries[i][i];
        }

        return r;

    }


    /**
     * The size of this matrix, which is the number of rows (and columns).
     */
    public int size() {
        return size;
    }

}
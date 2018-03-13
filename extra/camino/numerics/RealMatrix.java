package numerics;

import java.lang.Math;
import java.lang.Cloneable;
import java.io.Serializable;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Matrix of real numbers.
 * 
 * <dt>Description:
 * 
 * <dd>Provides methods to generate, access and manipulate matrices.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class RealMatrix implements Serializable, Cloneable {

    private static final long serialVersionUID = -495648904583884654L;

    /**
     * The number of rows in the matrix.
     */
    protected int r;

    /**
     * The number of columns in the matrix.
     */
    protected int c;

    /**
     * The array of matrix entries.
     */
    public double[][] entries;

    /**
     * Default constructor creates an empty matrix.
     */
    public RealMatrix() {
        r = 0;
        c = 0;
        entries = null;
    }

    /**
     * Creates a matrix full of zeros with the specified size.
     * 
     * @param rows
     *            Number of rows in the matrix.
     * 
     * @param columns
     *            Number of columns in the matrix.
     */
    public RealMatrix(int rows, int columns) {
        r = rows;
        c = columns;
        entries = null;
        entries = new double[r][c];
//         for (int i = 0; i < r; i++)
//             for (int j = 0; j < c; j++)
//                 entries[i][j] = 0.0;
    }


    /**
     * Creates a matrix using the supplied array as the entries. Note that no bounds checking will be
     * done, so ensure that the array represents a square matrix. 
     * 
     * @param e the matrix in array form, where there are <code>e.length</code> rows and 
     * <code>e[0].length</code> columns.
     *         
     */
    public RealMatrix(double[][] e) {
        r = e.length;
        c = e[0].length;
        entries = e;
    }


    public Object clone() {
        RealMatrix cl = new RealMatrix(this.r, this.c);
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                cl.entries[i][j] = this.entries[i][j];
        return cl;
    }

    /**
     * Returns a square identity matrix.
     * 
     * @param rows
     *            Size of the matrix.
     * 
     * @return A rows-by-rows identity matrix.
     */
    public static RealMatrix identity(int rows) {
        RealMatrix id = new RealMatrix(rows, rows);
        for (int i = 0; i < rows; i++) {
            id.setEntry(i, i, 1.0);
        }

        return id;
    }

    /**
     * Returns a single entry of the matrix.
     * 
     * @param row
     *            The row of the entry.
     * 
     * @param column
     *            The column of the entry
     * 
     * @return entry (row, column)
     */
    public double entry(int row, int column) {
        if (row < 0 || row >= r || column < 0 || column >= c)
                throw new RuntimeException("RealMatrix entry out of range! row= " + row
                        + " num rows= " + r + " column = " + column + " num columns = "
                        + c);
        return entries[row][column];
    }

    /**
     * Sets a matrix entry.
     * 
     * @param row
     *            The row of the entry.
     * 
     * @param column
     *            The column of the entry.
     * 
     * @param val
     *            The new value of the entry.
     */
    public void setEntry(int row, int column, double val) {
        if (row < 0 || row >= r || column < 0 || column >= c)
                throw new RuntimeException("RealMatrix entry out of range! row= " + row
                        + " num rows= " + r + " column = " + column + " num columns = "
                        + c);
        entries[row][column] = val;
    }

    /**
     * Returns the number of rows.
     * 
     * @return The number of rows.
     */
    public int rows() {
        return r;
    }

    /**
     * Returns the number of columns.
     * 
     * @return The number of columns.
     */
    public int columns() {
        return c;
    }

    /**
     * Switches two columns of the matrix.
     * 
     * @param a
     *            Index of first column to swap.
     * 
     * @param b
     *            Index of second column to swap.
     */
    public void swapColumns(int a, int b) {
        double temp = 0;
        for (int i = 0; i < r; i++) {
            temp = entries[i][a];
            entries[i][a] = entries[i][b];
            entries[i][b] = temp;
        }
    }

    /**
     * Switches two rows of the matrix.
     * 
     * @param a
     *            Index of first row to swap.
     * 
     * @param b
     *            Index of second row to swap.
     */
    public void swapRows(int a, int b) {
        double temp = 0;
        for (int i = 0; i < c; i++) {
            temp = entries[a][i];
            entries[a][i] = entries[b][i];
            entries[b][i] = temp;
        }
    }

    /**
     * Multiplies matrices.
     * 
     * @param m
     *            The matrix to postmultiply by.
     * 
     * @return The result of the multiplication.
     */
    public RealMatrix product(RealMatrix m) {
        RealMatrix prod = new RealMatrix(r, m.columns());
        if (c != m.rows()) {
            System.err.println("WARNING: Cannot multiply matrices (columns in this =" + c
                    + ", rows in given = " + m.rows() + ") - zero matrix returned");
            return prod;
        }
        for (int i = 0; i < r; i++)
            for (int j = 0; j < m.columns(); j++) {
                double value = 0;
                for (int k = 0; k < c; k++)
                    value += entries[i][k] * m.entries[k][j];
                prod.entries[i][j] = value;
            }
        return prod;
    }



    /**
     * Adds matrices.
     * 
     * @param m
     *            Matrix to add.
     * 
     * @return The result of the addition.
     */
    public RealMatrix add(RealMatrix m) {
        RealMatrix sum = new RealMatrix(r, c);
        if (r != m.rows() || c != m.columns()) {
            System.err.println("Cannot add matrices of different sizes [this is a (" + r
                    + "x" + c + "), given a (" + m.rows() + "x" + m.columns()
                    + ")]	- zero matrix returned.");
            return sum;
        }
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                sum.entries[i][j] = entries[i][j] + m.entries[i][j];
        return sum;
    }


    /**
     * Subtracts matrices.
     * 
     * @param m
     *            Matrix to subtract from this one.
     * 
     * @return The result of the subtraction.
     */
    public RealMatrix sub(RealMatrix m) {
        RealMatrix sum = new RealMatrix(r, c);
        if (r != m.rows() || c != m.columns()) {
            System.err.println("Cannot subtract matrices of different sizes [this is a (" + r
                    + "x" + c + "), given a (" + m.rows() + "x" + m.columns()
                    + ")]	- zero matrix returned.");
            return sum;
        }
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                sum.entries[i][j] = entries[i][j] - m.entries[i][j];
        return sum;
    }



    /**
     * Negates the matrix.
     * 
     * @return The result of the negation.
     */
    public RealMatrix negate() {
        RealMatrix neg = new RealMatrix(r, c);
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                neg.entries[i][j] = -entries[i][j];
        return neg;
    }

    /**
     * Multiplies the matrix by a scalar.
     * 
     * @param scalar
     *            The scalar multiple.
     */
    public void scale(double scalar) {
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                entries[i][j] *= scalar;
    }

    /**
     * Multiplies the matrix by a scalar and returns the result leaving
     * the original unchanged.
     * 
     * @param scalar
     *            The scalar multiple.
     *
     * @return scaled matrix
     */
    public RealMatrix scalarMult(double scalar) {
        RealMatrix sc = new RealMatrix(r, c);
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                sc.entries[i][j] = entries[i][j]*scalar;

        return sc;
    }

    /**
     * Inverts a matrix inefficiently by computing the determinant and adjoint.
     * 
     * @return The inverse matrix.
     */
    public RealMatrix inverse() {
        RealMatrix inv = new RealMatrix(r, c);
        if (r != c) {
            System.err
                    .println("Cannot invert non-square RealMatrix - NULL RealMatrix returned");
            return inv;
        }
        inv = adjoint();
        inv = inv.transpose();
        inv.scale(1.0 / det());
        return inv;
    }

    /**
     * Computes the adjoint matrix.
     * 
     * @return The adjoint.
     */
    public RealMatrix adjoint() {
        RealMatrix adj = new RealMatrix(r, c);
        if (r != c) {
            System.err
                    .println("Cannot adjoint non-square RealMatrix - NULL RealMatrix returned");
            return adj;
        }
        if (r == 1) {
            adj = (RealMatrix) this.clone();
            return adj;
        }
        int ctr1 = -1;
        for (int i = 0; i < r; i++) {
            int ctr2 = -ctr1;
            for (int j = 0; j < c; j++) {
                //Construct sub RealMatrix.

                RealMatrix subm = new RealMatrix(r - 1, c - 1);
                int rctr = 0;
                for (int i2 = 0; i2 < r; i2++) {

                    int cctr = 0;
                    for (int j2 = 0; j2 < c; j2++) {
                        if (i2 != i && j2 != j) {
                            subm.entries[rctr][cctr] = entries[i2][j2];
                            cctr += 1;
                        }
                    }
                    if (i2 != i) rctr += 1;
                }
                adj.entries[i][j] = ((double) ctr2) * subm.det();
                ctr2 = -ctr2;
            }
            ctr1 = -ctr1;
        }
        return adj;
    }

    /**
     * Computes the determinant. Beware that this is O(N!) for non-triangular matrices.
     * 
     * @return The determinant.
     */
    public double det() {
        if (r != c) {
            System.err
                    .println("Cannot find determinant of non-square RealMatrix - zero returned");
            return 0;
        }
        if (r == 1) return entries[0][0];

        //If the RealMatrix is triangular can find the determinant efficiently.

        if (triangular()) {
            double d = entries[0][0];
            for (int i = 1; i < r; i++)
                d *= entries[i][i];
            return d;
        }
        double d = 0;
        int ctr = 1;
        for (int i = 0; i < c; i++) {
            RealMatrix subm = new RealMatrix(r - 1, c - 1);
            for (int i2 = 0; i2 < r - 1; i2++) {
                int cctr = 0;
                for (int j = 0; j < c; j++) {
                    if (j != i) {
                        subm.entries[i2][cctr] = entries[i2 + 1][j];
                        cctr += 1;
                    }
                }
            }
            d += ((double) ctr) * entries[0][i] * subm.det();
            ctr = -ctr;
        }
        return d;
    }

    /**
     * Computes the trace of the matrix.
     * 
     * @return The trace.
     */
    public double trace() {
        double tr = 0;
        for (int i = 0; i < r && i < c; i++)
            tr += entries[i][i];
        return tr;
    }


    
    /**
     * Computes the transpose.
     * 
     * @return The transpose of the matrix.
     */
    public RealMatrix transpose() {
        RealMatrix tr = new RealMatrix(c, r);
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                tr.entries[j][i] = entries[i][j];
        return tr;
    }

    /**
     * Does the jacobi transformation. Adapted from numerical recipes in C, pg
     * 467.
     * 
     * @return Array of matrices. The first element is the diagonal matrix of
     *         eigenvalues. The second is the array of eigenvectors. The columns
     *         of the second matrix are the eigenvectors.
     */
    public RealMatrix[] jacobi() {

        if (!(symmetric())) {
            throw new RuntimeException(
                    "RealMatrix is not symmetric, cannot perform Jacobi transformation");
        }

        //Initialise evecs to the identity and diag to this.

        RealMatrix evecs = new RealMatrix(r, r);
        for (int i = 0; i < r; i++)
            evecs.entries[i][i] = 1.0;

        RealMatrix diag = new RealMatrix(r, r);
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++)
                diag.entries[i][j] = entries[i][j];

        //Typical matrices require 6-10 sweeps to achieve convergence
        //to machine precision. A maximum of 50 is allowed.

        for (int i = 0; i < 50; i++) {

            //Find sum of the off-diagonal elements.

            double sm = 0.0;
            for (int j = 0; j < r; j++)
                for (int k = 0; k < j; k++)
                    sm += Math.abs(diag.entries[j][k]);

            //Default return that relies on quadratic convergence to
            //machine underflow.

            if (sm == 0.0) {
                RealMatrix[] results = new RealMatrix[2];
                results[0] = diag;
                results[1] = evecs;

                return results;
            }

            //In the first three sweeps, rotation is carried out only
            //if the selected element is greater than some
            //threshold. Thereafter it is carried out anyway.

            double thresh = 0.0;
            if (i < 4) thresh = (0.2 * sm) / ((double) (r * r));

            for (int j = 0; j < r - 1; j++)
                for (int k = j + 1; k < r; k++) {
                    //After four sweeps, skip the rotation if the
                    //off-diagonal element is small.

                    double g = 100.0 * Math.abs(diag.entries[j][k]);
                    if (i > 4
                            && (double) (Math.abs(diag.entries[j][j]) + g) == (double) Math
                                    .abs(diag.entries[j][j])
                            && (double) (Math.abs(diag.entries[k][k]) + g) == (double) Math
                                    .abs(diag.entries[k][k])) {
                        diag.entries[j][k] = 0.0;
                        diag.entries[k][j] = 0.0;
                    }
                    else if (Math.abs(diag.entries[j][k]) > thresh) {
                        double h = diag.entries[k][k] - diag.entries[j][j];
                        double t = 0;
                        if ((double) (Math.abs(h) + g) == (double) Math.abs(h))
                            t = diag.entries[j][k] / h;
                        else {
                            double theta = 0.5 * h / (diag.entries[j][k]);
                            t = 1.0 / (Math.abs(theta) + Math.sqrt(1.0 + theta * theta));
                            if (theta < 0.0) t = -t;
                        }
                        double cs = 1.0 / Math.sqrt(1.0 + t * t);
                        double sn = t * cs;
                        RealMatrix rot = new RealMatrix(r, r);
                        for (int p = 0; p < r; p++)
                            rot.entries[p][p] = 1.0;
                        rot.entries[j][j] = cs;
                        rot.entries[k][k] = cs;
                        rot.entries[j][k] = sn;
                        rot.entries[k][j] = -sn;
                        diag = ((rot.transpose()).product(diag)).product(rot);
                        evecs = evecs.product(rot);
                    }
                }
        }
        //System.err.println("Too many iterations in Jacobi routine");
        //System.exit(0);
        throw new RuntimeException("Too many iterations in Jacobi routine");

        //RealMatrix[] rP = null;
        //return rP;
    }

    /**
     * Tests whether the matrix is symmetric.
     * 
     * @return The result of the test.
     */
    public boolean symmetric() {
        if (r != c) return false;
        for (int i = 0; i < r; i++)
            for (int j = 0; j < i; j++)
                if (entries[i][j] != entries[j][i]) return false;
        return true;
    }

    /**
     * Tests whether a matrix is triangular.
     * 
     * @return The result of the test.
     */
    public boolean triangular() {
        if (r != c) return false;
        if (r == 1) return true;
        boolean isu = true;
        boolean isl = true;
        for (int i = 1; i < r; i++)
            for (int j = 0; j < i; j++) {
                isu = isu && (entries[i][j] == 0);
                isl = isl && (entries[j][i] == 0);
            }
        return (isu || isl);
    }

    /**
     * Does LU inversion.
     * 
     * @return The inverse matrix.
     */
    public RealMatrix luinvert() {
        if (r != c) {
            throw new RuntimeException("Cannot LU-decompose non-square RealMatrix");
        }
        RealMatrix inv = new RealMatrix(r, r);

        RealMatrix tempmat = new RealMatrix(r, c);
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                tempmat.entries[i][j] = entries[i][j];

        int[] indx = new int[r];
        double[] col = new double[r];
        tempmat.ludecomp(indx);
        for (int j = 0; j < r; j++) {
            for (int i = 0; i < r; i++)
                col[i] = 0.0;
            col[j] = 1.0;
            tempmat.lubksub(indx, col);
            for (int i = 0; i < r; i++)
                inv.entries[i][j] = col[i];
        }
        return inv;
    }

    public String toString() {
        String str = "";
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++)
                str += entries[i][j] + "\t";
            str += "\n";
        }

        return str;
    }

    /**
     * Does LU decomposition. Adapted from Numerical Recipes.
     * 
     * @param indx
     *            Array of indexes.
     */
    private void ludecomp(int[] indx) {
        RealMatrix vv = new RealMatrix(r, 1);
        double big;
        double dum;
        int imax = 0;
        double temp;
        double sum;
        for (int i = 0; i < r; i++) {
            big = 0.0;
            for (int j = 0; j < r; j++)
                if ((temp = Math.abs(entries[i][j])) > big) big = temp;
            if (big == 0.0) {
                throw new RuntimeException("Cannot LU-decompose singular RealMatrix");
            }
            vv.entries[i][0] = 1.0 / big;
        }
        for (int j = 0; j < r; j++) {
            for (int i = 0; i <= j; i++) {
                sum = entries[i][j];
                for (int k = 0; k < i; k++)
                    sum -= entries[i][k] * entries[k][j];
                entries[i][j] = sum;
            }
            big = 0.0;
            for (int i = j + 1; i < r; i++) {
                sum = entries[i][j];
                for (int k = 0; k < j; k++)
                    sum -= entries[i][k] * entries[k][j];
                entries[i][j] = sum;
                if ((dum = vv.entries[i][0] * Math.abs(sum)) >= big) {
                    big = dum;
                    imax = i;
                }
            }
            if (j != imax) {
                for (int k = 0; k < r; k++) {
                    dum = entries[imax][k];
                    entries[imax][k] = entries[j][k];
                    entries[j][k] = dum;
                }
                vv.entries[imax][0] = vv.entries[j][0];
            }
            indx[j] = imax;
            if (j != r - 1) {
                dum = 1.0 / entries[j][j];
                for (int i = j + 1; i < r; i++)
                    entries[i][j] *= dum;
            }
        }
    }

    /**
     * Does the back substitution for LU decomposition.
     */
    private void lubksub(int[] indx, double[] b) {
        double sum;
        int ii = 0;
        for (int i = 0; i < r; i++) {
            int ip = indx[i];
            sum = b[ip];
            b[ip] = b[i];
            if (ii != 0)
                for (int j = ii - 1; j < i; j++)
                    sum -= entries[i][j] * b[j];
            else if (sum != 0) ii = i + 1;
            b[i] = sum;
        }
        for (int i = r - 1; i >= 0; i--) {
            sum = b[i];
            for (int j = i + 1; j < r; j++)
                sum -= entries[i][j] * b[j];
            b[i] = sum / entries[i][i];
        }
    }

    /**
     * Computes the singular value decomposition U, S and V of the matrix A, so
     * that A = U S V^T.
     * 
     * @return {U, S, V}.
     */
    public RealMatrix[] svd() throws SVD_Exception {

        // Copy the array of this matrix. Note that the NRC routines
        // start counting from 1 in their arrays.
        double[][] uMat = new double[r + 1][c + 1];
        for (int i = 1; i <= r; i++) {
            for (int j = 1; j <= c; j++) {
                uMat[i][j] = entries[i - 1][j - 1];
            }
        }

        // Initialize the v matrix.
        double[][] vMat = new double[c + 1][c + 1];

        // And the vector of singular values.
        double[] sMat = new double[c + 1];

        // Call the numerical recipes routine.
        dsvdcmp(uMat, r, c, sMat, vMat);

        // Now create the actual matrix objects.
        RealMatrix[] results = new RealMatrix[3];

        // U
        results[0] = new RealMatrix(uMat.length - 1, uMat[0].length - 1);
        for (int i = 0; i < uMat.length - 1; i++) {
            for (int j = 0; j < uMat[0].length - 1; j++) {
                results[0].setEntry(i, j, uMat[i + 1][j + 1]);
            }
        }

        // S
        results[1] = new RealMatrix(sMat.length - 1, sMat.length - 1);
        for (int i = 0; i < sMat.length - 1; i++) {
            results[1].setEntry(i, i, sMat[i + 1]);
        }

        // V
        results[2] = new RealMatrix(vMat.length - 1, vMat[0].length - 1);
        for (int i = 0; i < vMat.length - 1; i++) {
            for (int j = 0; j < vMat[0].length - 1; j++) {
                results[2].setEntry(i, j, vMat[i + 1][j + 1]);
            }
        }

        return results;

    }

    /**
     * Numerical recipes routine for singular value decomposition.
     * 
     * @param a
     *            The matrix to decompose, which gets replaced by matrix U.
     * 
     * @param m
     *            The number of rows.
     * 
     * @param n
     *            The number of columns.
     * 
     * @param w
     *            Array to fill with singular values.
     * 
     * @param v
     *            Array to fill with matrix V.
     */
    protected void dsvdcmp(double[][] a, int m, int n, double[] w, double[][] v)
            throws SVD_Exception {
        int flag, i, its, j, jj, k, l, nm;
        double anorm, c, f, g, h, s, scale, x, y, z;

        // Java wants these initialized.
        l = 0;
        nm = 0;

        double[] rv1 = new double[n + 1];
        g = scale = anorm = 0.0;
        for (i = 1; i <= n; i++) {
            l = i + 1;
            rv1[i] = scale * g;
            g = s = scale = 0.0;
            if (i <= m) {
                for (k = i; k <= m; k++)
                    scale += Math.abs(a[k][i]);
                // if (scale) {
                if (scale != 0) {
                    for (k = i; k <= m; k++) {
                        a[k][i] /= scale;
                        s += a[k][i] * a[k][i];
                    }
                    f = a[i][i];
                    //g = -SIGN(sqrt(s),f);
                    g = -(f >= 0.0 ? Math.sqrt(s) : -Math.sqrt(s));
                    h = f * g - s;
                    a[i][i] = f - g;
                    for (j = l; j <= n; j++) {
                        for (s = 0.0, k = i; k <= m; k++)
                            s += a[k][i] * a[k][j];
                        f = s / h;
                        for (k = i; k <= m; k++)
                            a[k][j] += f * a[k][i];
                    }
                    for (k = i; k <= m; k++)
                        a[k][i] *= scale;
                }
            }
            w[i] = scale * g;
            g = s = scale = 0.0;
            if (i <= m && i != n) {
                for (k = l; k <= n; k++)
                    scale += Math.abs(a[i][k]);
                // if (scale) {
                if (scale != 0.0) {
                    for (k = l; k <= n; k++) {
                        a[i][k] /= scale;
                        s += a[i][k] * a[i][k];
                    }
                    f = a[i][l];
                    // g = -SIGN(sqrt(s),f);
                    g = -(f >= 0.0 ? Math.sqrt(s) : -Math.sqrt(s));
                    h = f * g - s;
                    a[i][l] = f - g;
                    for (k = l; k <= n; k++)
                        rv1[k] = a[i][k] / h;
                    for (j = l; j <= m; j++) {
                        for (s = 0.0, k = l; k <= n; k++)
                            s += a[j][k] * a[i][k];
                        for (k = l; k <= n; k++)
                            a[j][k] += s * rv1[k];
                    }
                    for (k = l; k <= n; k++)
                        a[i][k] *= scale;
                }
            }
            // anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
            anorm = (anorm > (Math.abs(w[i]) + Math.abs(rv1[i]))) ? anorm : (Math
                    .abs(w[i]) + Math.abs(rv1[i]));
        }
        for (i = n; i >= 1; i--) {
            if (i < n) {
                // if (g) {
                if (g != 0.0) {
                    for (j = l; j <= n; j++)
                        v[j][i] = (a[i][j] / a[i][l]) / g;
                    for (j = l; j <= n; j++) {
                        for (s = 0.0, k = l; k <= n; k++)
                            s += a[i][k] * v[k][j];
                        for (k = l; k <= n; k++)
                            v[k][j] += s * v[k][i];
                    }
                }
                for (j = l; j <= n; j++)
                    v[i][j] = v[j][i] = 0.0;
            }
            v[i][i] = 1.0;
            g = rv1[i];
            l = i;
        }
        // for (i=IMIN(m,n);i>=1;i--) {
        for (i = (m < n ? m : n); i >= 1; i--) {
            l = i + 1;
            g = w[i];
            for (j = l; j <= n; j++)
                a[i][j] = 0.0;
            // if (g) {
            if (g != 0.0) {
                g = 1.0 / g;
                for (j = l; j <= n; j++) {
                    for (s = 0.0, k = l; k <= m; k++)
                        s += a[k][i] * a[k][j];
                    f = (s / a[i][i]) * g;
                    for (k = i; k <= m; k++)
                        a[k][j] += f * a[k][i];
                }
                for (j = i; j <= m; j++)
                    a[j][i] *= g;
            }
            else
                for (j = i; j <= m; j++)
                    a[j][i] = 0.0;
            ++a[i][i];
        }
        for (k = n; k >= 1; k--) {
            for (its = 1; its <= 30; its++) {
                flag = 1;
                for (l = k; l >= 1; l--) {
                    nm = l - 1;
                    if ((double) (Math.abs(rv1[l]) + anorm) == anorm) {
                        flag = 0;
                        break;
                    }
                    if ((double) (Math.abs(w[nm]) + anorm) == anorm) break;
                }
                // if (flag) {
                if (flag != 0) {
                    c = 0.0;
                    s = 1.0;
                    for (i = l; i <= k; i++) {
                        f = s * rv1[i];
                        rv1[i] = c * rv1[i];
                        if ((double) (Math.abs(f) + anorm) == anorm) break;
                        g = w[i];
                        h = dpythag(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = -f * h;
                        for (j = 1; j <= m; j++) {
                            y = a[j][nm];
                            z = a[j][i];
                            a[j][nm] = y * c + z * s;
                            a[j][i] = z * c - y * s;
                        }
                    }
                }
                z = w[k];
                if (l == k) {
                    if (z < 0.0) {
                        w[k] = -z;
                        for (j = 1; j <= n; j++)
                            v[j][k] = -v[j][k];
                    }
                    break;
                }
                if (its == 30)
                        throw new SVD_Exception("no convergence in 30 dsvdcmp iterations");
                x = w[l];
                nm = k - 1;
                y = w[nm];
                g = rv1[nm];
                h = rv1[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                g = dpythag(f, 1.0);
                // f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                f = ((x - z) * (x + z) + h * ((y / (f + (f >= 0 ? g : -g))) - h)) / x;
                c = s = 1.0;
                for (j = l; j <= nm; j++) {
                    i = j + 1;
                    g = rv1[i];
                    y = w[i];
                    h = s * g;
                    g = c * g;
                    z = dpythag(f, h);
                    rv1[j] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y *= c;
                    for (jj = 1; jj <= n; jj++) {
                        x = v[jj][j];
                        z = v[jj][i];
                        v[jj][j] = x * c + z * s;
                        v[jj][i] = z * c - x * s;
                    }
                    z = dpythag(f, h);
                    w[j] = z;
                    // if (z) {
                    if (z != 0.0) {
                        z = 1.0 / z;
                        c = f * z;
                        s = h * z;
                    }
                    f = c * g + s * y;
                    x = c * y - s * g;
                    for (jj = 1; jj <= m; jj++) {
                        y = a[jj][j];
                        z = a[jj][i];
                        a[jj][j] = y * c + z * s;
                        a[jj][i] = z * c - y * s;
                    }
                }
                rv1[l] = 0.0;
                rv1[k] = f;
                w[k] = x;
            }
        }
    }

    /**
     * Computes (a^2 + b^2)^(1/2) without destructive overflow or underflow.
     * 
     * @param a
     *            a
     * 
     * @param b
     *            b
     * 
     * @return (a^2 + b^2)^(1/2)
     */
    protected double dpythag(double a, double b) {
        double absa, absb;
        absa = Math.abs(a);
        absb = Math.abs(b);
        if (absa > absb)
            return absa * Math.sqrt(1.0 + (absb / absa) * (absb / absa));
        else
            return (absb == 0.0 ? 0.0 : absb
                    * Math.sqrt(1.0 + (absa / absb) * (absa / absb)));
    }

    /* 
     * Calculates the pseudo inverse of the matrix using SVD
     * svThresh set the threshold
     */
    public RealMatrix pseudoInv(double svThresh) 
    {
    	// Get the singular value decomposition.
	RealMatrix[] svd = null;
	try 
	{
	    svd = this.svd();
	}
	catch (Exception e) 
	{
            //System.err.println(this);
	    throw new RuntimeException(e);
	}

	// Find the maximum singular value.
	double maxSV = svd[1].entry(0, 0);
	for (int i = 0; (i < svd[1].rows() && i < svd[1].columns()); i++) 
	{
	    if (svd[1].entry(i, i) > maxSV) 
	    {
		maxSV = svd[1].entry(i, i);
	    }
	}

	// Invert the singular values.
	for (int i = 0; (i < svd[1].rows() && i < svd[1].columns()); i++) 
	{
	    double curSV = svd[1].entry(i, i);
	    svd[1].setEntry(i, i, (curSV > maxSV / svThresh) ? (1.0 / curSV) : 0.0);
	}

	// Get the pseudo inverse and return.
	return svd[2].product(svd[1].transpose()).product(svd[0].transpose());
    }

    /* 
     * Calculates the pseudo inverse of the matrix using SVD
     * Calls RealMatrix pseudoInv(RealMatrix M, double svThresh) using a threshold of 1.0E12
     */
    public RealMatrix pseudoInv() 
    {
        return pseudoInv(1.0E12);
    }
    
    public boolean isSingular()
    {    	
    	for (int c=0; c < this.c; c++)
    	{
    		boolean res=true;
        	for (int r=0; r < this.r; r++)
        	{
        		res=res && this.entries[r][c]==0;        		
        	}
        	
        	if (res)
        	{
        		return true;
        	}
    	}
    	return false;
    }
    
    
}


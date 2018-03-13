package numerics;

import java.text.*;

/**
 * <dl>
 * <dt>Purpose: Container for eigenvectors and eigenvalues.
 * <dd>
 * <BR><BR>
 * <dt>Description: 
 * <dd> A simple container class to hold eigenvectors and eigenvalues.
 * Whoever constructs this object should make sure that eigenvector[i] corresponds to eigenvalue[i].
 *
 * @author Philip Cook
 * @version $Id$
 * 
 */
public class EigenSystem3D {

    private static DecimalFormat df = new DecimalFormat("0.00000E00");


    public final Vector3D[] eigenvectors;
    public final double[] eigenvalues;

    public EigenSystem3D(Vector3D[] vectors, double[] values) {

	if ( (vectors.length != 3) || (vectors.length != values.length) ) {
	    throw new java.lang.IllegalArgumentException("Number of eigenvalues must match number of eigenvectors");
	}
	
	
	eigenvectors = vectors;
	eigenvalues = values;

    }

    /**
     * Convert a sorted eigen system where the eigenvalues are (in decreasing order)
     * seig([0][0], [0][1], [0][2]). 
     * The eigenvectors should be (seig[1][n], seig[2][n], seig[3][n]).
     *
     */
    public EigenSystem3D(double[][] seig) {

	// bounds check
	eigenvectors = new Vector3D[3];
	eigenvalues = new double[3];

	eigenvectors[0] = new Vector3D(seig[1][0], seig[2][0], seig[3][0]);
	eigenvectors[1] = new Vector3D(seig[1][1], seig[2][1], seig[3][1]);
	eigenvectors[2] = new Vector3D(seig[1][2], seig[2][2], seig[3][2]);
	
	eigenvalues[0] = seig[0][0];
	eigenvalues[1] = seig[0][1];
	eigenvalues[2] = seig[0][2];
    }


    /**
     * Sorts eigen system by eigenvalue, greatest first. 
     */
   public static EigenSystem3D sort(Jama.EigenvalueDecomposition eig) {
	
	double[][] v = eig.getV().getArray();
	
	double[] lams = eig.getRealEigenvalues();

	if (lams.length != 3) {
	    throw new java.lang.IllegalArgumentException("eig is not from a 3x3 matrix");
	}

	int minInd = 0;
	int maxInd = 0;
	for(int j = 1; j < 3; j++) {
	    if(lams[j] >= lams[maxInd]) {
		maxInd = j;
	    }
	    if(lams[j] < lams[minInd]) {
		minInd = j;
	    }
	}

	int midInd = 3 - minInd - maxInd;


	Vector3D[] evecs = new Vector3D[3];
	double[] evals = new double[3];

	evals[0] = lams[maxInd];
	evals[1] = lams[midInd];
	evals[2] = lams[minInd];

	evecs[0] = new Vector3D(v[0][maxInd], v[1][maxInd], v[2][maxInd]);
	evecs[1] = new Vector3D(v[0][midInd], v[1][midInd], v[2][midInd]);
	evecs[2] = new Vector3D(v[0][minInd], v[1][minInd], v[2][minInd]);

	return new EigenSystem3D(evecs, evals);

    }




    /**
     * Sorts eigen system by eigenvalue, greatest first.
     */
    public static EigenSystem3D sort(RealMatrix A) {

	if (A.rows() != 3 || A.columns() != 3) {
	    throw new IllegalArgumentException("Method requires 3x3 matrix");
	}

	RealMatrix[] eig = A.jacobi();
       
	int minInd = 0;
	int maxInd = 0;
       
	for(int j = 1; j < 3; j++) {
	    if(eig[0].entries[j][j] >= eig[0].entries[maxInd][maxInd]) {
		maxInd = j;
	    }
	    if(eig[0].entries[j][j] < eig[0].entries[minInd][minInd]) {
		minInd = j;
	    }
	}
       
	int midInd = 3 - minInd - maxInd;
	
	Vector3D[] evecs = new Vector3D[3];
	double[] evals = new double[3];

	evecs[0] = new Vector3D(eig[1].entries[0][maxInd], eig[1].entries[1][maxInd], 
				eig[1].entries[2][maxInd]);

	evecs[1] = new Vector3D(eig[1].entries[0][midInd], eig[1].entries[1][midInd], 
				eig[1].entries[2][midInd]);

	evecs[2] = new Vector3D(eig[1].entries[0][minInd], eig[1].entries[1][minInd], 
				eig[1].entries[2][minInd]);

	evals[0] = eig[0].entries[maxInd][maxInd];
	evals[1] = eig[0].entries[midInd][midInd];
	evals[2] = eig[0].entries[minInd][minInd];
    
	return new EigenSystem3D(evecs, evals);
         
    }


    /**
     * String representation is 3x4 matrix of eigenvalues and eigenvectors.
     */
    public String toString() {
	
	String s = "";

	for (int i = 0; i < 3; i++) {
	    s = s + df.format(eigenvalues[i]) + " " + df.format(eigenvectors[i].x) + " " + 
		df.format(eigenvectors[i].y) + " " + df.format(eigenvectors[i].z);

	    s = s + "\n";
	}

	return s;
    }

}

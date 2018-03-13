package numerics;

import java.util.Random;



/**
 * Angular Central Gaussian Distribution on a 3D sphere. Gives p(x) = |A|^{-1/2}(x^T A x)^{-3/2}.
 *
 * @see section 9.4.4 of Mardia and Jupp ("Directional Statistics", Wiley, 2000)
 * @version $Id$
 * @author Philip Cook
 */
public final class ACG_Distribution implements AxialDistribution {


    private final double sigma1;
    private final double sigma2;
    private final double sigma3;

 
    // matrix that diagonalizes A
    private final RealMatrix T;

    private final RealMatrix invA;

    // normalization constant
    private final double normC;
    
    private final Random ran;


    // don't want anyone to use this
    private ACG_Distribution() {

	T = null;
	invA = null;

	sigma1 = 0.0;
	sigma2 = 0.0;
	sigma3 = 0.0;

	normC = 1.0;
	ran = null;

    }


    /**
     * @param A the covariance matrix of a Gaussian distribution N(0,A).
     * @param r a Random number generator used for the generation of samples.
     */
    public ACG_Distribution(RealMatrix A, Random r) {

	invA = A.inverse();

	normC = 1.0 / Math.sqrt(A.det());

	// System.out.println(A);

	RealMatrix[] eig = A.jacobi();

	int minInd = 0;
	int maxInd = 0;

	for(int j=1; j<3; j++) {
	    if(eig[0].entries[j][j] >= eig[0].entries[maxInd][maxInd]) {
		maxInd = j;
	    }
	    if(eig[0].entries[j][j] < eig[0].entries[minInd][minInd]) {
		minInd = j;
	    }
	}

	int midInd = 3 - minInd - maxInd;
	
	sigma1 = Math.sqrt(eig[0].entries[maxInd][maxInd]);
	sigma2 = Math.sqrt(eig[0].entries[midInd][midInd]);
	sigma3 = Math.sqrt(eig[0].entries[minInd][minInd]);

	RealMatrix tmp = new RealMatrix(3,3);

	tmp.entries[0][0] = eig[1].entries[0][maxInd];
	tmp.entries[1][0] = eig[1].entries[1][maxInd];
	tmp.entries[2][0] = eig[1].entries[2][maxInd];

	tmp.entries[0][1] = eig[1].entries[0][midInd];
	tmp.entries[1][1] = eig[1].entries[1][midInd];
	tmp.entries[2][1] = eig[1].entries[2][midInd];

	tmp.entries[0][2] = eig[1].entries[0][minInd];
	tmp.entries[1][2] = eig[1].entries[1][minInd];
	tmp.entries[2][2] = eig[1].entries[2][minInd];

	T = tmp;

	ran = r;

    }



    /**
     * @param evecs the eigenvectors of A, ordered by the associated eigenvalue, with
     * the greatest first.
     * @param sigmaSqs eigenvalues of A, in descending order, such that <code>sigmaSqs[i]</code> 
     * corresponds to <code>evecs[i]</code>.
     * @param r a Random number generator used for the generation of samples.
     *
     */
    public ACG_Distribution(Vector3D[] evecs, double[] sigmaSqs, Random r) {

	RealMatrix[] eigEigT = new RealMatrix[3];
	
	for (int i = 0; i < 3; i++) {
	    
	    eigEigT[i] = new RealMatrix(3,3);

	    eigEigT[i].entries[0][0] = evecs[i].x * evecs[i].x;
	    eigEigT[i].entries[0][1] = evecs[i].x * evecs[i].y;
	    eigEigT[i].entries[0][2] = evecs[i].x * evecs[i].z;
	    eigEigT[i].entries[1][0] = eigEigT[i].entries[0][1];
	    eigEigT[i].entries[1][1] = evecs[i].y * evecs[i].y;
	    eigEigT[i].entries[1][2] = evecs[i].y * evecs[i].z;
	    eigEigT[i].entries[2][0] = eigEigT[i].entries[0][2];
	    eigEigT[i].entries[2][1] = eigEigT[i].entries[1][2];
	    eigEigT[i].entries[2][2] = evecs[i].z * evecs[i].z;

	    eigEigT[i].scale(sigmaSqs[i]);
	}

	RealMatrix A = eigEigT[0].add(eigEigT[1].add(eigEigT[2]));

	//	System.out.println(A);

	invA = A.inverse();

	normC = 1.0 / Math.sqrt(A.det());

	sigma1 = Math.sqrt(sigmaSqs[0]);
	sigma2 = Math.sqrt(sigmaSqs[1]);
	sigma3 = Math.sqrt(sigmaSqs[2]);

	RealMatrix tmp = new RealMatrix(3,3);

	tmp.entries[0][0] = evecs[0].x;
	tmp.entries[1][0] = evecs[0].y;
	tmp.entries[2][0] = evecs[0].z;

	tmp.entries[0][1] = evecs[1].x;
	tmp.entries[1][1] = evecs[1].y;
	tmp.entries[2][1] = evecs[1].z;

	tmp.entries[0][2] = evecs[2].x;
	tmp.entries[1][2] = evecs[2].y;
	tmp.entries[2][2] = evecs[2].z;

	T = tmp;

	ran = r;


    }

    /**
     * Gets the pdf of the distribution.
     *
     * @param x a unit vector on the unit sphere.
     */
    public double pdf(Vector3D x) {
 	double xTinvAx = x.x * 
 	    (x.x * invA.entries[0][0] + x.y * invA.entries[1][0] + x.z * invA.entries[2][0]) + 
 	    x.y * (x.x * invA.entries[0][1] + x.y * invA.entries[1][1] + x.z * invA.entries[2][1]) +
 	    x.z * (x.x * invA.entries[0][2] + x.y * invA.entries[1][2] + x.z * invA.entries[2][2]);


	return normC * Math.pow(xTinvAx, -1.5);
    }


    /**
     * Gets the log pdf of the distribution.
     *
     * @param x a unit vector on the unit sphere.
     */
    public double logPDF(Vector3D x) {
	
	double xTinvAx = x.x * 
 	    (x.x * invA.entries[0][0] + x.y * invA.entries[1][0] + x.z * invA.entries[2][0]) + 
 	    x.y * (x.x * invA.entries[0][1] + x.y * invA.entries[1][1] + x.z * invA.entries[2][1]) +
 	    x.z * (x.x * invA.entries[0][2] + x.y * invA.entries[1][2] + x.z * invA.entries[2][2]);

 	return Math.log(normC) + -1.5 * Math.log(xTinvAx);

    }

    


    /**
     *
     * @return the next vector from the distribution.
     * 
     */
    public Vector3D nextVector() {
	// generate 3D vector x distributed as N(0, A). Then x.normalized() is distributed as ACD(A)
	
	double xx = ran.nextGaussian() * sigma1;
	double yy = ran.nextGaussian() * sigma2;
	double zz = ran.nextGaussian() * sigma3;

	double sx = T.entries[0][0] * xx + T.entries[0][1] * yy + T.entries[0][2] * zz;
	double sy = T.entries[1][0] * xx + T.entries[1][1] * yy + T.entries[1][2] * zz;
	double sz = T.entries[2][0] * xx + T.entries[2][1] * yy + T.entries[2][2] * zz;

	return new Vector3D(sx, sy, sz).normalized();
	
    }

    /**
     * @return String representation of this object.
     */
    public String toString() {
	return new String("ACG_Distribution, sigma1 == " + sigma1 + 
	    ", sigma2 == " + sigma2 + ", sigma3 == " + sigma3);
    }


    /**
     * @return the eigenvectors of A
     */
    public Vector3D[] eigenvectors() {
	return new Vector3D[] {
		new Vector3D(T.entries[0][0], T.entries[1][0], T.entries[2][0]),
		new Vector3D(T.entries[0][1], T.entries[1][1], T.entries[2][1]),
		new Vector3D(T.entries[0][2], T.entries[1][2], T.entries[2][2])
	    };

    } 


    /**
     * @return the square root of the eigenvalues of A.
     */ 
    public double[] sigmas() {
	return new double[] {sigma1, sigma2, sigma3};
    }
}

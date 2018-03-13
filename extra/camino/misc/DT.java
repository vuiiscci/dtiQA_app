package misc;

import numerics.*;

import java.util.logging.*;


/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Encapsulates a diffusion tensor.
 * 
 * <dt>Description:
 * 
 * <dd>Stores the six components and allows common diffusion tensor operations.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class DT implements Cloneable {

    private static Logger logger = Logger.getLogger("camino.misc.DT");

    /**
     * The six components of the diffusion tensor.
     */
    protected final double dxx, dxy, dxz, dyy, dyz, dzz;


    /**
     * Constructor requires the six components.
     * 
     * @param Dxx
     *            The xx component of the diffusion tensor.
     * 
     * @param Dxy
     *            The xy component of the diffusion tensor.
     * 
     * @param Dxz
     *            The xz component of the diffusion tensor.
     * 
     * @param Dyy
     *            The yy component of the diffusion tensor.
     * 
     * @param Dyz
     *            The yz component of the diffusion tensor.
     * 
     * @param Dzz
     *            The zz component of the diffusion tensor.
     */
    public DT(double Dxx, double Dxy, double Dxz, double Dyy, double Dyz, double Dzz) {
        dxx = Dxx;
        dxy = Dxy;
        dxz = Dxz;
        dyy = Dyy;
        dyz = Dyz;
        dzz = Dzz;
    }


    /**
     * Constructor from an array containing the six components.
     * 
     * @param dtComps
     *            {Dxx, Dxy, Dxz, Dyy, Dyz, Dzz}
     */
    public DT(double[] dtComps) {
        dxx = dtComps[0];
        dxy = dtComps[1];
        dxz = dtComps[2];
        dyy = dtComps[3];
        dyz = dtComps[4];
        dzz = dtComps[5];
    }



    


    /**
     * Contracts the diffusion tensor by a vector.
     * 
     * @param q
     *            The vector to contract by.
     * 
     * @return q^T D q.
     */
    public double contractBy(double[] q) {
        return q[0] * (dxx * q[0] + dxy * q[1] + dxz * q[2]) + q[1]
                * (dxy * q[0] + dyy * q[1] + dyz * q[2]) + q[2]
                * (dxz * q[0] + dyz * q[1] + dzz * q[2]);
    }

    public RealMatrix toMatrix() {
        RealMatrix r = new RealMatrix(3,3);

        r.entries[0][0] = dxx;
        r.entries[0][1] = dxy;
        r.entries[0][2] = dxz;
        r.entries[1][1] = dyy;
        r.entries[1][2] = dyz;
        r.entries[2][2] = dzz;
        r.entries[1][0] = dxy;
	r.entries[2][0] = dxz;
	r.entries[2][1] = dyz;

        return r;
    }

    /**
     * Multiplies the vector by the diffusion tensor.
     * 
     * @param q
     *            The vector to multiply.
     * 
     * @return D q.
     */
    public double[] multiply(double[] q) {
        double[] Dq = new double[3];
        Dq[0] = dxx * q[0] + dxy * q[1] + dxz * q[2];
        Dq[1] = dxy * q[0] + dyy * q[1] + dyz * q[2];
        Dq[2] = dxz * q[0] + dyz * q[1] + dzz * q[2];
        return Dq;
    }

    public Object clone() {
        return new DT(dxx, dxy, dxz, dyy, dyz, dzz);
    }


    /**
     * Computes the determinant.
     * 
     * @return The determinant of the diffusion tensor.
     */
    public double determinant() {
        return (-dxz * dxz * dyy + dxy * dxz * dyz + dxy * dxz * dyz - dxx * dyz * dyz
                - dxy * dxy * dzz + dxx * dyy * dzz);
    }


    /**
     * Computes the trace of the diffusion tensor.
     * 
     * @return trace(D).
     */
    public double trace() {
        return dxx + dyy + dzz;
    }


    /**
     * Compute the fractional anisotropy of the diffusion tensor.
     * 
     * @return FA.
     */
    public double fa() {
        double modSq = trace() / 3.0;
        modSq *= modSq;

        double stsp = dxx * dxx + 2.0 * dxy * dxy + 2.0 * dxz * dxz + dyy * dyy + 2.0
                * dyz * dyz + dzz * dzz;

        double fracAnis = 0.0;
        if (stsp != 0.0) {
            
            // may be negative due to precision error if the tensor is close to isotropic
            double faSq = 3.0 * (1.0 - 3.0 * modSq / stsp) / 2.0;

            if (faSq >= 0.0) {
                fracAnis = Math.sqrt(faSq);
            }
            else if (faSq > -1E-9) {
                fracAnis = 0.0;
            }
            else {
                logger.warning("Spurious FA value : " + fracAnis + " for tensor\n" + toString());
            }
        }

        return fracAnis;
    }


    /**
     * Applies a similarity transform to the diffusion tensor using the linear
     * transformation matrix supplied. This applied the inverse transformation
     * to method transform.
     * 
     * @param trans
     *            A 3x3 linear transformation matrix.
     * 
     * @return The transformed tensor: trans^T D trans.
     */
    public DT iTransform(RealMatrix trans) {

        double a = trans.entry(0, 0);
        double b = trans.entry(0, 1);
        double c = trans.entry(0, 2);
        double d = trans.entry(1, 0);
        double e = trans.entry(1, 1);
        double f = trans.entry(1, 2);
        double g = trans.entry(2, 0);
        double h = trans.entry(2, 1);
        double i = trans.entry(2, 2);

        return doTrans(a, b, c, d, e, f, g, h, i);
    }


    /**
     * Called by <code>iTransform<code> and <code>transform</code> to avoid code duplication.
     */
    protected DT doTrans(double a, double b, double c, double d, double e, double f,
            double g, double h, double i) {

        // Do the transformation directly using result computed in
        // mathematica.
        double newDxx = a * a * dxx + 2 * a * d * dxy + d * d * dyy + a * dxz * g + a
                * dxz * g + 2 * d * dyz * g + dzz * g * g;
        double newDxy = a * b * dxx + b * d * dxy + a * dxy * e + d * dyy * e + b * dxz
                * g + dyz * e * g + a * dxz * h + d * dyz * h + dzz * g * h;
        double newDxz = a * c * dxx + c * d * dxy + a * dxy * f + d * dyy * f + c * dxz
                * g + dyz * f * g + a * dxz * i + d * dyz * i + dzz * g * i;
        double newDyy = b * b * dxx + 2 * b * dxy * e + dyy * e * e + b * dxz * h + b
                * dxz * h + 2 * dyz * e * h + dzz * h * h;
        double newDyz = b * c * dxx + c * dxy * e + b * dxy * f + dyy * e * f + c * dxz
                * h + dyz * f * h + b * dxz * i + dyz * e * i + dzz * h * i;
        double newDzz = c * c * dxx + 2 * c * dxy * f + dyy * f * f + c * dxz * i + c
                * dxz * i + 2 * dyz * f * i + dzz * i * i;

        return new DT(newDxx, newDxy, newDxz, newDyy, newDyz, newDzz);
    }


    /**
     * Applies a similarity transform to the diffusion tensor using the linear
     * transformation matrix supplied.
     * 
     * @param trans
     *            A 3x3 linear transformation matrix.
     * 
     * @return The transformed tensor: trans D trans^T.
     */
    public DT transform(RealMatrix trans) {

        double a = trans.entry(0, 0);
        double d = trans.entry(0, 1);
        double g = trans.entry(0, 2);
        double b = trans.entry(1, 0);
        double e = trans.entry(1, 1);
        double h = trans.entry(1, 2);
        double c = trans.entry(2, 0);
        double f = trans.entry(2, 1);
        double i = trans.entry(2, 2);

        return doTrans(a, b, c, d, e, f, g, h, i);
    }


    /**
     * Computes the inverse diffusion tensor.
     * 
     * @return The diffusion tensor with components of the inverse diffusion
     *         tensor matrix.
     */
    public DT inverse() {
        double det = determinant();

        double invDxx = (-dyz * dyz + dyy * dzz) / det;
        double invDxy = (dxz * dyz - dxy * dzz) / det;
        double invDxz = (-dxz * dyy + dxy * dyz) / det;
        double invDyy = (-dxz * dxz + dxx * dzz) / det;
        double invDyz = (dxy * dxz - dxx * dyz) / det;
        double invDzz = (-dxy * dxy + dxx * dyy) / det;

        return new DT(invDxx, invDxy, invDxz, invDyy, invDyz, invDzz);
    }


    /**
     * Scales a diffusion tensor.
     * 
     * @param scalar
     *            The scaling factor.
     * 
     * @return The scaled tensor.
     */
    public DT scale(double scalar) {
        return new DT(scalar * dxx, scalar * dxy, scalar * dxz, scalar * dyy, scalar
                * dyz, scalar * dzz);
    }

    /**
     * Returns an array containing the components of the diffusion tensor.
     * 
     * @return {Dxx, Dxy, Dxz, Dyy, Dyz, Dzz}
     */
    public double[] getComponents() {
        double[] comps = new double[6];
        comps[0] = dxx;
        comps[1] = dxy;
        comps[2] = dxz;
        comps[3] = dyy;
        comps[4] = dyz;
        comps[5] = dzz;

        return comps;
    }


    /**
     * Eigensystem computed using the <code>Jama.EigenvalueDecomposition</code> class.
     * 
     * @return A matrix containing the eigenvalues and eigenvectors. ([0][0],
     *         [0][1], [0][2]) are the eigenvalues. ([1][0], [2][0], [3][0]) is
     *         the eigenvector corresponding to eigenvalue [0][0].
     */
    public double[][] eigenSystem() {

        Jama.Matrix matrix = new Jama.Matrix(3, 3);

        matrix.set(0, 0, dxx);
        matrix.set(0, 1, dxy);
        matrix.set(1, 0, dxy);
        matrix.set(1, 1, dyy);
        matrix.set(0, 2, dxz);
        matrix.set(2, 0, dxz);
        matrix.set(2, 2, dzz);
        matrix.set(1, 2, dyz);
        matrix.set(2, 1, dyz);

        try {
            Jama.EigenvalueDecomposition e = new Jama.EigenvalueDecomposition(matrix);
            
            double[] eigenValues = e.getRealEigenvalues();
            
            Jama.Matrix V = e.getV();
            
            double[][] eigenSystem = new double[4][3];
            double[][] aV = V.getArray();

            for (int i = 0; i < 3; i++) {
                eigenSystem[0][i] = eigenValues[i];
                for (int j = 0; j < 3; j++) {
                    eigenSystem[j + 1][i] = aV[j][i];
                }
            }
            
            return eigenSystem;
        }
        catch (Exception e) {
            // JAMA can freak out if there is an overflow
            LoggedException.logExceptionWarning(e, "Bad tensor data, returning default values, program will continue.\n Tensor was :\n" + toString(), Thread.currentThread().getName());

	               

            double[][] eigenSystem = new double[4][3];

            // unit axes for vectors
            // zero eigenvalues
            eigenSystem[1][0] = 1.0;
            eigenSystem[1][1] = 1.0;
            eigenSystem[1][2] = 1.0;

            return eigenSystem;
            
        }
    }

 
    /**
     * Returns a sorted eigensystem, with the largest eigenvalue and
     * corresponding eigenvector in column 0.
     * 
     * @return Array containing the eigenvalues and eigenvectors ordered by
     *         eigenvalue. Element [0][0] is the largest eigenvalue, ([1][0],
     *         [2][0], [3][0]) the corresponding eigenvector. Element [0][2] is
     *         the smallest eigenvalue.
     */
    public double[][] sortedEigenSystem() {
        double[][] eigenSys = eigenSystem();

        int minInd = 0;
        int maxInd = 0;
        for (int j = 1; j < 3; j++) {
            if (eigenSys[0][j] >= eigenSys[0][maxInd]) {
                maxInd = j;
            }
            if (eigenSys[0][j] < eigenSys[0][minInd]) {
                minInd = j;
            }
        }

        int midInd = 3 - minInd - maxInd;

        double[][] sEigenSys = new double[4][3];
        for (int j = 0; j < 4; j++) {
            sEigenSys[j][0] = eigenSys[j][maxInd];
            sEigenSys[j][1] = eigenSys[j][midInd];
            sEigenSys[j][2] = eigenSys[j][minInd];
        }

        return sEigenSys;
    }


    /**
     * Reorients the tensor using the PPD algorithm.
     *
     * @param jac
     *            The jacobian of the transformation.
     *
     * @return
     *            The reoriented DT.
     */
    public DT ppd(RealMatrix jac) {
        
        // Get the eigensystem of the DT.
        double[][] eig = sortedEigenSystem();
        
        // Get first two eigenvectors and eigenvalues in matrices
        RealMatrix evals = new RealMatrix(3,3);
        RealMatrix e1 = new RealMatrix(3,1);
        RealMatrix e2 = new RealMatrix(3,1);
        for(int i=0; i<3; i++) {
            evals.entries[i][i] = eig[0][i];
            e1.entries[i][0] = eig[i+1][0];
            e2.entries[i][0] = eig[i+1][1];
        }

        //System.err.println("Evals\n" + evals);
        //System.err.println("E1\n" + e1);
        //System.err.println("E2\n" + e2);

        // Compute the transformed eigenvectors.
        RealMatrix n1 = jac.product(e1);
        RealMatrix n2 = jac.product(e2);
        //System.err.println("N1\n" + n1);
        //System.err.println("N2\n" + n2);

        // Normalize them.
        double modN1 = Math.sqrt(n1.entries[0][0]*n1.entries[0][0]
                                 + n1.entries[1][0]*n1.entries[1][0]
                                 + n1.entries[2][0]*n1.entries[2][0]);
        double modN2 = Math.sqrt(n2.entries[0][0]*n2.entries[0][0]
                                 + n2.entries[1][0]*n2.entries[1][0]
                                 + n2.entries[2][0]*n2.entries[2][0]);

        for(int i=0; i<3; i++) {
            n1.entries[i][0] = n1.entries[i][0]/modN1;
            n2.entries[i][0] = n2.entries[i][0]/modN2;
        }
        //System.err.println("N1norm\n" + n1);
        //System.err.println("N2norm\n" + n2);

        // Make new orthogonal eigensystem.
        RealMatrix evecTrans = new RealMatrix(3,3);
        
        // The transformed major eigenvector is n1.
        for(int i=0; i<3; i++) {
            evecTrans.entries[i][0] = n1.entries[i][0];
        }

        // Compute the projection of n2 onto the plane perpendicular
        // to n1.
        double n1n2 = n1.entries[0][0]*n2.entries[0][0]
            + n1.entries[1][0]*n2.entries[1][0]
            + n1.entries[2][0]*n2.entries[2][0];

        RealMatrix newE2 = n2.add(n1.scalarMult(-n1n2));

        // Renormalize and add to matrix
        double modNewE2 = Math.sqrt(newE2.entries[0][0]*newE2.entries[0][0]
                                 + newE2.entries[1][0]*newE2.entries[1][0]
                                 + newE2.entries[2][0]*newE2.entries[2][0]);
        for(int i=0; i<3; i++) {
            evecTrans.entries[i][1] = newE2.entries[i][0]/modNewE2;
        }

        // Finally the new E3 is orthogonal to the new E1 and E2.
        evecTrans.entries[0][2] = evecTrans.entries[1][0]*evecTrans.entries[2][1] - evecTrans.entries[1][1]*evecTrans.entries[2][0];
        evecTrans.entries[1][2] = evecTrans.entries[2][0]*evecTrans.entries[0][1] - evecTrans.entries[2][1]*evecTrans.entries[0][0];
        evecTrans.entries[2][2] = evecTrans.entries[0][0]*evecTrans.entries[1][1] - evecTrans.entries[0][1]*evecTrans.entries[1][0];

        //System.err.println("evecTrans\n" + evecTrans);

        // Reconstruct the transformed diffusion tensor.
        RealMatrix dTrans = evecTrans.product(evals.product(evecTrans.transpose()));

        //System.err.println("dTrans\n" + dTrans);

        return new DT(dTrans.entries[0][0], dTrans.entries[0][1], dTrans.entries[0][2], dTrans.entries[1][1], dTrans.entries[1][2], dTrans.entries[2][2]);
    }


    /**
     * Returns the principal direction.
     * 
     * @return {x, y, z}.
     */
    public double[] getPD() {

        double[][] s = sortedEigenSystem();
        double[] pd = new double[3];
        pd[0] = s[1][0];
        pd[1] = s[2][0];
        pd[2] = s[3][0];

        return pd;
    }


    /**
     * Returns an array containing the three eigenvalues of the diffusion tensor.
     * 
     * @return Unsorted array of eigenvalues.
     */
    public double[] eigenValues() {
        double[][] eSys = eigenSystem();
        double[] eVals = new double[3];
        for (int i = 0; i < 3; i++) {
            eVals[i] = eSys[0][i];
        }

        return eVals;
    }


    /**
     *
     * Construct a DT from an eigen system of the format produced by eigenSystem().
     *
     */
    public static DT dtFromEig(double[][] eig) {

        // D == l1 * (e1 * e1T) + l2 * (e2 * e2T) + l3 * (e3 * e3T)

        // long-winded way avoids computing inverse

        RealMatrix e1 = new RealMatrix(new double[][] { {eig[1][0]}, {eig[2][0]}, {eig[3][0]} });
        RealMatrix e2 = new RealMatrix(new double[][] { {eig[1][1]}, {eig[2][1]}, {eig[3][1]} });
        RealMatrix e3 = new RealMatrix(new double[][] { {eig[1][2]}, {eig[2][2]}, {eig[3][2]} });

        RealMatrix e1e1T = e1.product(e1.transpose()).scalarMult(eig[0][0]);
        RealMatrix e2e2T = e2.product(e2.transpose()).scalarMult(eig[0][1]);
        RealMatrix e3e3T = e3.product(e3.transpose()).scalarMult(eig[0][2]);

        RealMatrix sum = e1e1T.add(e2e2T).add(e3e3T);

        return new DT(sum.entries[0][0], sum.entries[0][1], sum.entries[0][2], sum.entries[1][1], sum.entries[1][2], sum.entries[2][2]);

    }



    /**
     * Returns a string containing the diffusion tensor matrix.
     */
    public String toString() {
        String s = dxx + " " + dxy + " " + dxz + "\n";
        s += dxy + " " + dyy + " " + dyz + "\n";
        s += dxz + " " + dyz + " " + dzz + "\n";

        return s;
    }

   
}


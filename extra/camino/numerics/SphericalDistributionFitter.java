package numerics;


import Jama.*;

/**
 * Superclass for classes that fit parametric distributions to vectors or axes on the sphere.
 *
 * @author Philip Cook
 * @version $Id$
 */
public abstract class SphericalDistributionFitter {


    /**
     * @return the normalized scatter matrix of the samples.
     * 
     */
    public static RealMatrix tBar(Vector3D[] sampleVecs) {
	Matrix tBar = new Matrix(3,3);

	for (int i = 0; i < sampleVecs.length; i++) {
	    
	    Matrix m = sampleVecs[i].toJamaMatrix();
	    
	    tBar.plusEquals( m.times(m.transpose()).times(1.0 / sampleVecs.length) );
	}

	return new RealMatrix(tBar.getArray());

    }


    /**
     *
     * @return the sorted eigen system of the normalized scatter matrix of the sampleVecs. 
     *
     */
    public static EigenSystem3D tBarEigenSystem(RealMatrix tBar) {
	return EigenSystem3D.sort(new Matrix(tBar.entries).eig());
    }
	  


    /**
     *
     * @return the sorted eigen system of the normalized scatter matrix of the sampleVecs. 
     *
     */
    public static EigenSystem3D tBarEigenSystem(Vector3D[] sampleVecs) {
	
	Matrix tBar = new Matrix(3,3);

	for (int i = 0; i < sampleVecs.length; i++) {
	    
	    Matrix m = sampleVecs[i].toJamaMatrix();
	    
	    tBar.plusEquals( m.times(m.transpose()).times(1.0 / sampleVecs.length) );
	}

	return EigenSystem3D.sort(tBar.eig());

    }


    /**
     * Gets the eigen system of the scatter matrix from an array.
     *
     * @param tBarComps {txx, txy, txz, tyy, tyz, tzz} (upper triangular) or 
     * {txx, txy, tyy, txz, tyz, tzz} (lower triangular).
     * @param upperTriangular true if the comps are in upper triangular format.
     * @return the sorted eigen system of the normalized scatter matrix.
     *
     */
    public static EigenSystem3D tBarEigenSystem(double[] tBarComps, boolean upperTriangular) {

	if (upperTriangular) {
	    return tBarEigenSystem(tBarComps);
	}
	else {
	    double[] tmp = new double[6];

	    tmp[0] = tBarComps[0];
	    tmp[1] = tBarComps[1];
	    tmp[2] = tBarComps[3];
	    tmp[3] = tBarComps[2];
	    tmp[4] = tBarComps[4];
	    tmp[5] = tBarComps[5];

	    return tBarEigenSystem(tmp);
	}
    }

    /**
     * Gets the eigen system of the scatter matrix from an array.
     *
     * @param tBarComps {txx, txy, txz, tyy, tyz, tzz}.
     *
     * @return the sorted eigen system of the normalized scatter matrix.
     *
     */
    public static EigenSystem3D tBarEigenSystem(double[] tBarComps) {
	    Matrix tBar = new Matrix(3,3);

	    double mxx = tBarComps[0];
	    double mxy = tBarComps[1];
	    double mxz = tBarComps[2];
	    double myy = tBarComps[3];
	    double myz = tBarComps[4];
	    double mzz = tBarComps[5];

	    tBar.set(0,0, mxx);
	    tBar.set(1,1, myy);
	    tBar.set(2,2, mzz);

	    tBar.set(0,1, mxy);
	    tBar.set(1,0, mxy);

	    tBar.set(0,2, mxz);
	    tBar.set(2,0, mxz);

	    tBar.set(1,2, myz);
	    tBar.set(2,1, myz);
	    
	    return EigenSystem3D.sort(tBar.eig());

    }

}

package numerics;



import java.util.Arrays;

/**
 * This class implements several routines for fitting an ACG distribution. 
 * <p> The routines come from Tyler, "Statistical analysis for the angular central gaussian distribution", Biometrika, 74:579-590, (1984).
 *
 * @author Philip Cook
 * @version $Id$
 */
public class ACG_Fitter extends SphericalDistributionFitter {


    // don't want this to be instantiable
    private ACG_Fitter() {

    }

   
    /**
     *
     * Finds the covariance matrix A of the ACG
     * distribution.
     *
     * @param sampleVecs samples from the distribution. 
     *
     * @return the solution that maximises the likelihood of the <code>sampleVecs</code>.
     *
     */
    public static RealMatrix findA(Vector3D[] sampleVecs) {
	
	RealMatrix lastAHat = RealMatrix.identity(3);

	double diffLastIteration = Double.MAX_VALUE;

        int iterations = 0;

	while (diffLastIteration > 1E-10 && iterations < 10000) {

	    RealMatrix aHat = new RealMatrix(3,3);
	    RealMatrix lastAHatInv = lastAHat.inverse();

	    double norm = 0.0;

	    for (int i = 0; i < sampleVecs.length; i++) {
		double xTaHatInvX = 
		    sampleVecs[i].x * (sampleVecs[i].x * lastAHatInv.entries[0][0] + 
				       sampleVecs[i].y * lastAHatInv.entries[1][0] + 
				       sampleVecs[i].z * lastAHatInv.entries[2][0]) + 
		    sampleVecs[i].y * (sampleVecs[i].x * lastAHatInv.entries[0][1] + 
				       sampleVecs[i].y * lastAHatInv.entries[1][1] + 
				       sampleVecs[i].z * lastAHatInv.entries[2][1]) +
		    sampleVecs[i].z * (sampleVecs[i].x * lastAHatInv.entries[0][2] + 
				       sampleVecs[i].y * lastAHatInv.entries[1][2] + 
				       sampleVecs[i].z * lastAHatInv.entries[2][2]);

// 		RealMatrix x = sampleVecs[i].toRealMatrix();
// 		RealMatrix xT = x.transpose();

// 		double xTaHatInvX = xT.product(lastAHatInv).product(x).entries[0][0];

  		RealMatrix xxT = new RealMatrix(3,3);
// 		RealMatrix xxT = x.product(xT);

		xxT.entries[0][0] = sampleVecs[i].x * sampleVecs[i].x;
		xxT.entries[0][1] = sampleVecs[i].x * sampleVecs[i].y;
		xxT.entries[0][2] = sampleVecs[i].x * sampleVecs[i].z;
		xxT.entries[1][0] = xxT.entries[0][1];
		xxT.entries[1][1] = sampleVecs[i].y * sampleVecs[i].y;
		xxT.entries[1][2] = sampleVecs[i].y * sampleVecs[i].z;
		xxT.entries[2][0] = xxT.entries[0][2];
		xxT.entries[2][1] = xxT.entries[1][2];
		xxT.entries[2][2] = sampleVecs[i].z * sampleVecs[i].z;

 		xxT.scale(1.0 / xTaHatInvX);

 		aHat = aHat.add(xxT);
		
 		norm += 1.0 / xTaHatInvX;

 	    }

 	    aHat.scale(3.0 / norm);

	    diffLastIteration = Math.abs(aHat.entries[0][0] - lastAHat.entries[0][0]) + 
		Math.abs(aHat.entries[0][1] - lastAHat.entries[0][1]) +
		Math.abs(aHat.entries[0][2] - lastAHat.entries[0][2]) +
		Math.abs(aHat.entries[1][1] - lastAHat.entries[1][1]) + 
		Math.abs(aHat.entries[1][2] - lastAHat.entries[1][2]) + 
		Math.abs(aHat.entries[2][2] - lastAHat.entries[2][2]);

	    lastAHat = aHat;

            iterations++;
	}

	return lastAHat;

    }




}

package tractography;

import numerics.*;
import java.util.Random;

/**
 *
 * Doesn't randomize; used for testing.
 * 
 */
public class PICoNonRandomizer implements PICoRandomizer {


    protected final TractographyImage image;


    public PICoNonRandomizer(TractographyImage im) {
	image = im;
    }


    /** 
     * Get the next randomized PDs for this voxel. 
     * Subsequent calls with the same arguments will return the same vectors, 
     * until #resetRandomization() is called.
     * 
     *
     */
    public Vector3D[] getRandomizedPDs(int i, int j, int k) {
	return image.getPDs(i,j,k);
    }

  

    /**
     * Always returns 0.5
     */
    public double pdf(int i, int j, int k, int pdIndex, Vector3D v) {
	return 0.5;
    }


}



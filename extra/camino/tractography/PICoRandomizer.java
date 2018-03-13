package tractography;

import numerics.*;
import java.util.Random;

/**
 *
 * Provides samples from the fibre orientation distribution function (PDF). 
 * There may be multiple fibres in each voxel, and hence multiple PDFs. 
 *
 * @version $Id $
 * @author  Philip Cook
 * 
 * 
 */
public interface PICoRandomizer {

    /**
     * Randomizes the PDs in the voxel (i,j,k) and returns them.
     */
    public Vector3D[] getRandomizedPDs(int i, int j, int k);

  

    
    /**
     * Evaluates the PDF for a vector at voxel (i,j,k).
     *
     * @param v a unit vector.
     * @param pdIndex the index of the PDF in the voxel.
     * 
     * @return the value of the PDF evaluated at vector <code>v</code>.
     */
    public double pdf(int i, int j, int k, int pdIndex, Vector3D v);


}

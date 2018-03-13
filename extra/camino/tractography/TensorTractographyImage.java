package tractography;

import data.*;
import imaging.*;
import misc.DT;
import numerics.*;


/**
 * Interface for tractography images
 *
 * @version $Id$
 * @author  Philip Cook
 * 
 */
public interface TensorTractographyImage extends TractographyImage {

   
    /**
     * Gets the diffusion tensors in some voxel. The order matches the order of vectors from 
     * #TractographyImage.getPDs(int, int, int).
     *
     */
    public DT[] getDTs(int i, int j, int k);


    /**
     * @return the mixing parameters in this voxel. The ordering matches the order of tensors from 
     * #getDTs(int, int, int)
     */
    public double[] getMix(int i, int j, int k);

}

package tractography;

import misc.*;
import numerics.*;


/**
 * 
 * Interface for interpolators that handle diffusion tensor data.
 *
 * @author Philip Cook
 * @version $Id$
 */
public interface TensorInterpolator extends ImageInterpolator {


    /** 
     * Gets the DT at some point. 
     * 
     * @param point the point in mm to interpolate at. 
     * @param previousDirection the direction of the previous tracking step. This
     * will be used with multi-tensor images to decide which tensors to use.
     * 
     * @return the interpolated tensor.
     */
    public DT getDT(Point3D point, Vector3D previousDirection);



}

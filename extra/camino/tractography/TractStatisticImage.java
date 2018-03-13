package tractography;

import misc.*;
import numerics.*;
import tools.*;


/**
 * Gathers statistics from a TractStatisticFilter and creates an image using the tracts.
 * A single value is associated with each streamline, the type of calue is set by the
 * <code>setTractStatistic</code> method. The values from all streamlines are combined into
 * the output image.
 * <p>
 * <b>Tract-based scalar images</b>
 * <BR>
 * Tract-based scalar images combine streamlines with a scalar image (for example, FA). The streamlines
 * should be defined in the same space as the scalar image. The vector value of the scalar data is evaluated at
 * each point along the streamline; a scalar function (min, max, mean, median, var, sum) of this data vector is 
 * computed and then asssigned to the image. 
 * </p>
 * <p>
 * <b>Tract-based images without scalars</b>
 * <BR>
 * Without a scalar image, the tract statistic is derived directly from the streamlines. Currently, the only
 * available statistic is the length. 
 * </p>
 * <p>
 * <b>Creating the image</b>
 * <BR>
 * The statistic from each tract is assigned to the seed voxel of the tract, or optionally to every 
 * voxel intersected by the tract. The <code>setImageStatistic</code> and <code>setCountIntersect</code>
 * methods determines how the tract values are combined into the output image.
 * values
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class TractStatisticImage {


    // dimensions of seed space
    private final int[] dataDims;
	
    // mm
    private final double[] voxelDims;

    private ScalarImage scalars = null;

    private String imageStat = "mean";

    private TractStatisticFilter tractStat = null;

    private VoxelwiseStatisticalImage image = null;

    private boolean seedVoxelOnly = true;


    /**
     * Construct an image that uses the scalar image to define the tract statistic.
     * 
     * @param scalarImg a ScalarImage. The interpolation scheme of the image may be changed
     * if interpolation is enabled for the TractStatisticImage. The image should be in the
     * same space as the tracts.
     *
     */
    public TractStatisticImage(ScalarImage scalarImg) {
        scalars = scalarImg;
	dataDims = scalars.getDataDims();
	voxelDims = scalars.getVoxelDims();
	tractStat = new TractStatisticFilter(scalars);
    }


    /**
     * For Matlab.
     * @param data the scalar data.
     * @param voxelDims voxel dimensions in mm.
     */
    public TractStatisticImage(double[][][] scalarData, double[] voxelDims) {
	this.scalars = new ScalarImage(scalarData, voxelDims);
	dataDims = scalars.getDataDims();
	this.voxelDims = scalars.getVoxelDims();
	tractStat = new TractStatisticFilter(scalarData, voxelDims);
    }


    /**
     * Constructs an image without scalars. 
     * @param dataDims image dimensions.
     * @param voxelDims voxel dimensions, in mm.
     */
    public TractStatisticImage(int[] dataDims, double[] voxelDims) {
	this.dataDims = dataDims;
	this.voxelDims = voxelDims;
	tractStat = new TractStatisticFilter(dataDims, voxelDims);
	tractStat.setTractStatistic("length");
    }
    

    /**
     * Sets the tract statistic.
     * @param stat one of the scalar statistics supported by <code>TractStatisticFilter</code>. Vector
     * statistics, such as "meanvar" are not supported by this image (the first element will be used).
     */
    public void setTractStatistic(String stat) {
	tractStat.setTractStatistic(stat);
    }

    
    /**
     * Sets the image statistic, which determines how multiple values in each voxel are combined. The
     * computation of mean and variance differs depending on whether the tract statistic is added to the 
     * seed voxel of the tract, or to all voxels intersected by the tract. 
     * 
     * 
     * @param stat one of: mean, min, max, median, var.
     */
    public void setImageStatistic(String stat) {

	imageStat = stat;

	if ( imageStat.equals("median") || imageStat.equals("var") ) {
	    image = new SparseVectorImage(dataDims, voxelDims);
	}
	else {
	    image = new DynamicScalarImage(dataDims, voxelDims);
	}

    }



    /**
     * By default, the image is not interpolated. 
     *
     * @param interp if true, use linear interpolation of the scalar image and treat the tracts
     * as interpolated. 
     *
     * @see TractStatisticFilter
     */
    public void setInterpolate(boolean interp) {
	tractStat.setInterpolate(interp);
    }



    /**
     * By default, the tract statistic is only added to the voxel where the tract was seeded. If this
     * is called with a <code>true</code> parameter, each tract statistic is added to all voxels 
     * intersected by the tract.
     *
     *
     */
    public void setCountIntersect(boolean countIntersect) {
	seedVoxelOnly = !countIntersect;
    }


    /**
     * Add all of the tracts in the collection to the image.
     */
    public void processTracts(TractCollection tc) {
	for (int i = 0; i < tc.numberOfTracts(); i++) {
	    processTract(tc.getTract(i));
	}
    }



    /**
     * Processes a single tract. 
     */
    public void processTract(Tract t) {
        
    
	double tractValue = tractStat.processTract(t)[0];

	if (seedVoxelOnly) {
	    Point3D seed = t.getPoint(t.seedPointIndex());

            int x = (int)(seed.x / voxelDims[0]);
            int y = (int)(seed.y / voxelDims[1]);
            int z = (int)(seed.z / voxelDims[2]);
	    
	    image.addValue(x,y,z, tractValue); 
	}
	else {
	    
	    Voxel[] voxels = t.toVoxelList(voxelDims[0], voxelDims[1], voxelDims[2]).getVoxels();
	    
	    for (int v = 0; v < voxels.length; v++) {
		image.addValue(voxels[v].x, voxels[v].y, voxels[v].z, tractValue);
	    }
	    
	}
	
    }

 

    /**
     * Gets the tract-based statistical image.
     *
     */
    public double[][][] getImageStatistic() {
        
        return image.getVoxelStatistic(imageStat);
    }



}

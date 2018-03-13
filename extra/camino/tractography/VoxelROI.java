package tractography;

import imaging.ImageHeader;
import misc.LoggedException;
import numerics.*;

import java.io.*;
import java.util.*;


/**
 * A free form ROI for tractography. The ROI consists of labelled regions.
 * The regions are labelled as integers and may take any positive value.
 * 
 *
 * @author Philip Cook
 * @version $Id$
 */
public final class VoxelROI implements ROI_Collection {

    private final PointListROI[] regions;

    private final int[] regionLabels;

    private final int numberOfRegions;

 

    /**
     * Create an ROI from integer-valued regions in an image. Seed points are placed at the centre
     * of each voxel, in physical space. They are translated to points in Camino space according to the
     * diffusion image header.
     *
     * @param roiFile the image containing the seed regions.
     *
     * @param diffusionSpace the header describing the space of the diffusion image. The physical space
     * of the seed and diffusion images must be aligned, but they can be in different voxel space.
     *  
     */
    public VoxelROI(String roiFile, ImageHeader diffusionSpace) {

        ImageHeader header = null;

        try {
            header = ImageHeader.readHeader(roiFile);
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }

        RealMatrix voxToPhys = header.getVoxelToPhysicalTransform();

        RealMatrix physToDiffusionVox = diffusionSpace.getPhysicalToVoxelTransform();

        boolean seedSpaceIsDiffusion = header.sameSpace(diffusionSpace);

        double[][][] roiVol = header.readSingleVolumeData();


        int xDataDim = header.xDataDim();
        int yDataDim = header.yDataDim();
        int zDataDim = header.zDataDim();

        HashMap<Integer, ArrayList<Point3D>> roiMap = new HashMap<Integer, ArrayList<Point3D>>(100, 0.75f);

	for (int k = 0; k < zDataDim; k++) {
	    for (int j = 0; j < yDataDim; j++) {
		for (int i = 0; i < xDataDim; i++) {
                    
                    if (roiVol[i][j][k] > 0.0) {
                        Integer label = new Integer((int)roiVol[i][j][k]);

                        ArrayList<Point3D> subROI = null;

                        if (roiMap.containsKey(label)) {
                            subROI = roiMap.get(label);
                        }
                        else {
                            subROI = new ArrayList<Point3D>(100);
                            
                            roiMap.put(label, subROI);
                        }


                        Point3D vox = new Point3D(i,j,k);

                        if (!seedSpaceIsDiffusion) {
                            // voxel coordinates map to voxel center via physical space transform
                            vox = new Point3D(i, j, k).transform(voxToPhys).transform(physToDiffusionVox);
                        }

                        // Camino space is indexed from voxel corner, hence the offset
                        Point3D seed = new Point3D( (vox.x + 0.5) * diffusionSpace.xVoxelDim(), 
                                                    (vox.y + 0.5) * diffusionSpace.yVoxelDim(), 
                                                    (vox.z + 0.5) * diffusionSpace.zVoxelDim());

                        subROI.add(seed);

                    }
                }
	    }
	}


        Integer[] keys = roiMap.keySet().toArray(new Integer[roiMap.size()]);

        Arrays.sort(keys);

        numberOfRegions = keys.length;

        regionLabels = new int[numberOfRegions];

        regions = new PointListROI[numberOfRegions];

        for (int i = 0; i < keys.length; i++) {
            Point3D[] points = roiMap.get(keys[i]).toArray(new Point3D[roiMap.get(keys[i]).size()]);
            
            regions[i] = new PointListROI(points, keys[i].intValue());

            regionLabels[i] = keys[i].intValue();
        }
      
    }





    /**
     * @return the number of labeled regions within this ROI.
     */
    public int numberOfRegions() {
	return numberOfRegions;
    }


 
    /**
     * Get a specific ROI, by its region label
     *
     * @return the ROI, or <code>null</code> no region with this label exists.
     */
    public RegionOfInterest getRegion(int label) {
	int index = Arrays.binarySearch(regionLabels, label);

        if (index > -1) {
            return regions[index];
        }
        else {
            throw new LoggedException("No such ROI label " + label);
        }
    }



    public RegionOfInterest[] getAllRegions() {

        RegionOfInterest[] defCopy = new RegionOfInterest[numberOfRegions];

        System.arraycopy(regions, 0, defCopy, 0, numberOfRegions);

        return defCopy;

    }

}

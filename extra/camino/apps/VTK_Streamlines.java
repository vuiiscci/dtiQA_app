package apps;


import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;
import tractography.*;

import java.util.logging.*;
import java.io.*;

/**
 * 
 * Converts streamlines from raw format to VTK polylines (in binary form). 
 *
 * @author Philip Cook
 * @version $Id$
 * 
 */
public class VTK_Streamlines {



    /**
     * Output manager
     */
    private static OutputManager om;


    public static void main(String[] args) {


	CL_Initializer.inputModel = "raw";

	CL_Initializer.CL_init(args);

        om = new OutputManager();
        DataOutputStream out = om.getOutputStream();

	int xDataDim = CL_Initializer.dataDims[0];
	int yDataDim = CL_Initializer.dataDims[1];
	int zDataDim = CL_Initializer.dataDims[2];
	    
	double xVoxelDim = Math.abs(CL_Initializer.voxelDims[0]);
	double yVoxelDim = Math.abs(CL_Initializer.voxelDims[1]);
	double zVoxelDim = Math.abs(CL_Initializer.voxelDims[2]);

	// raw data used for scalar image
	double[][][] scalarVol = null;

	// provides interpolated values of scalarVol
	ScalarImage scalarImage = null;

	String seedFile = null;
	String scalarFile = null;
	String targetFile = null;

	// colour (single scalar) is identical for each point in the streamline
	// streamline colour depends on the ROI index of the seed
	boolean colourBySeedROI = false;

	// colour (single scalar) is identical for each point in the streamline
	// streamline colour depends on the ROI index of the first target a streamline hits
	boolean colourByTargetROI = false;

	// colour (RGB triplet) depends on local orientation of streamline
	boolean colourByOrientation = false;

	// colour (single scalar) depends on local scalar value 
	boolean colourByScalar = false;

	// single scalar "colours" are just scalar values, the LUT maps them to colours

	
	// for colour by scalar only; don't want to interpolate targets
	boolean interpolateScalars = false;
	
	
	for (int i = 0; i < args.length; i++) {

	    // image args
	    if (args[i].equals("-scalarfile")) {
		scalarFile = args[i + 1];
		colourByScalar = true;
		CL_Initializer.markAsParsed(i, 2);
	    } 
            else if (args[i].equals("-targetfile")) {
                targetFile = args[i + 1];
		colourByTargetROI = true;
                CL_Initializer.markAsParsed(i,2);
	    }
	    else if (args[i].equals("-seedfile")) {
		seedFile = args[i + 1];
		colourBySeedROI = true;
		CL_Initializer.markAsParsed(i, 2);
	    } 
	    else if (args[i].equals("-colourorient")) {
		colourByOrientation = true;
		CL_Initializer.markAsParsed(i);
	    } 
	    if (args[i].equals("-interpolatescalars") || args[i].equals("-interpolate")) {
		interpolateScalars = true;
		CL_Initializer.markAsParsed(i);
	    }

	}

	CL_Initializer.checkParsing(args);

        // means do we output scalars, which may be defined from the tract orientation or from an image
	boolean haveScalars = colourByScalar || colourByTargetROI || colourBySeedROI || colourByOrientation;


	String colourSourceFile = null;

	if (colourByScalar) {
	    colourSourceFile = scalarFile;
	}
	if (colourByTargetROI) {
	    colourSourceFile = targetFile;
	}
	if (colourBySeedROI) {
	    colourSourceFile = seedFile;
	}

        RealMatrix physToVox = RealMatrix.identity(4);

	if (colourSourceFile != null) {
            CL_Initializer.headerTemplateFile = colourSourceFile;
            CL_Initializer.initInputSpaceAndHeaderOptions();
            physToVox = CL_Initializer.headerTemplate.getPhysicalToVoxelTransform();


            xDataDim = CL_Initializer.dataDims[0];
            yDataDim = CL_Initializer.dataDims[1];
            zDataDim = CL_Initializer.dataDims[2];
	    
            xVoxelDim = Math.abs(CL_Initializer.voxelDims[0]);
            yVoxelDim = Math.abs(CL_Initializer.voxelDims[1]);
            zVoxelDim = Math.abs(CL_Initializer.voxelDims[2]);
            
	    ImageHeader ih = null;
	    
	    try {
		ih = ImageHeader.readHeader(colourSourceFile);
	    }
	    catch (IOException e) {
		throw new LoggedException("Can't read file " + colourSourceFile, e);
	    }
    
	    scalarVol = ih.readSingleVolumeData();

	    if (colourByScalar) {
		scalarImage = new ScalarImage(scalarVol, new double[] {xVoxelDim, yVoxelDim, zVoxelDim});

		if (interpolateScalars) {
		    scalarImage.setInterpolation("linear");
		}
	    }
	}
	
        // Read tracts in physical space, transform them if we need to later
	TractSource tractSource = new TractSource(CL_Initializer.inputFile, null);
	
	// collect all tracts, then convert to polydata
	TractCollection tc = new TractCollection(1000, 100.0);

	while (tractSource.more()) {
	    tc.addTract(tractSource.nextTract());
	}

	try {

	    // VTK header
	    String hdr = "# vtk DataFile Version 3.0\nCamino tracts\nBINARY\nDATASET POLYDATA\n";

	    hdr = hdr + "POINTS " + tc.totalPoints() + " float\n";
	
	    out.write(hdr.getBytes("US-ASCII"));

	    // now output all the points
            
            for (int i = 0; i < tc.numberOfTracts(); i++) {

                // Output points directly in physical space
		Tract t = new Tract(tc.getTract(i));

		for (int p = 0; p < t.numberOfPoints(); p++) {
		    Point3D point = t.getPoint(p);

		    out.writeFloat((float)point.x);
		    out.writeFloat((float)point.y);
		    out.writeFloat((float)point.z);
		}
	    }

	    // now write lines, ie which points belong to which tract
	    String lines = "\nLINES " + tc.numberOfTracts() + " " + ( tc.numberOfTracts() + tc.totalPoints() ) + "\n";

	    out.write(lines.getBytes("US-ASCII"));	

	    // lines reference point ID, where each point has a unique ID
	    int pointCounter = 0;

	    for (int i = 0; i < tc.numberOfTracts(); i++) {
		Tract t = tc.getTract(i);
	    
		out.writeInt(t.numberOfPoints());

		for (int p = 0; p < t.numberOfPoints(); p++) {
		    out.writeInt(pointCounter);
		    pointCounter++;
		}
	    }	

	    String cells = "CELL_DATA " + tc.numberOfTracts() + "\nPOINT_DATA " + tc.totalPoints() + "\n";

	    out.write(cells.getBytes("US-ASCII"));	

	    if (haveScalars) {

		// now need to write scalar data
		// SCALARS dataName dataType numComp 
		// LOOKUP_TABLE tableName 
	
		String scalars = "";

		if (colourByOrientation) {
		    // unsigned char type assumed
		    scalars = "COLOR_SCALARS RGB_scalars 3\n";
		}
		else if (colourByTargetROI) {
		    scalars = "SCALARS target_scalars float 1\nLOOKUP_TABLE default\n";
		}
		else if (colourBySeedROI) {
		    scalars = "SCALARS seed_scalars float 1\nLOOKUP_TABLE default\n";
		}
		else if (colourByScalar) {
		    scalars = "SCALARS custom_scalars float 1\nLOOKUP_TABLE default\n";		
		}

		out.write(scalars.getBytes("US-ASCII"));	



		// write scalar values as binary
		for (int i = 0; i < tc.numberOfTracts(); i++) {

                    // copy in case we need to put this into Camino space
		    Tract t = new Tract(tc.getTract(i));
		
		    if (colourByOrientation) {

			Point3D point = null;
			Point3D next = null;

			// if Tract has 1 point, it is coloured red
			Vector3D orientation = new Vector3D(1.0, 0.0, 0.0);

			// color point 0 by direction from 0 -> 1
			for (int p = 0; p < t.numberOfPoints() - 1; p++) {
			
			    point = t.getPoint(p);
			    next = t.getPoint(p + 1);
			    
			    orientation = new Vector3D(point, next).normalized();

			    out.writeByte( (byte)(255 * Math.abs(orientation.x)) );
			    out.writeByte( (byte)(255 * Math.abs(orientation.y)) );
			    out.writeByte( (byte)(255 * Math.abs(orientation.z)) );
			}
			
			// colour last point by orientation from (numPoints - 2) -> (numPoints - 1)
			out.writeByte( (byte)(255 * Math.abs(orientation.x)) );
			out.writeByte( (byte)(255 * Math.abs(orientation.y)) );
			out.writeByte( (byte)(255 * Math.abs(orientation.z)) );
			
			
		    }
		    else {
                        t.transformToCaminoSpace(physToVox, xVoxelDim, yVoxelDim, zVoxelDim);
                    }


                    if (colourByTargetROI) {
                        
     			VoxelList voxelList = t.toVoxelList(xVoxelDim, yVoxelDim, zVoxelDim);
		    
			Voxel[] voxels = voxelList.getVoxels();
		    
			int voxelSeedIndex = voxelList.seedPointIndex();

		    
			// first target hit in each direction
			int upwardTarget = 0;
			int downwardTarget = 0;

			upward: 
			for (int v = voxelSeedIndex; v < voxels.length; v++) {
			    int voxelTarget = (int)scalarVol[voxels[v].x][voxels[v].y][voxels[v].z];
			
			    if (voxelTarget > 0) {
				upwardTarget = voxelTarget;
				break upward;
			    }
			}

			downward: 
			for (int v = voxelSeedIndex; v >= 0; v--) {
			    int voxelTarget = (int)scalarVol[voxels[v].x][voxels[v].y][voxels[v].z];
			
			    if (voxelTarget > 0) {
				downwardTarget = voxelTarget;
				break downward;
			    }
			}
		    
			// if only hit the target in one direction, colour whole streamline with target ID
			if (upwardTarget > 0 && downwardTarget == 0) {
			    downwardTarget = upwardTarget;
			}
			else if (upwardTarget == 0 && downwardTarget > 0) {
			    upwardTarget = downwardTarget;
			}

			// write values
			for (int p = 0; p < t.seedPointIndex(); p++) {
			    out.writeFloat((float)downwardTarget);
			}
			for (int p = t.seedPointIndex(); p < t.numberOfPoints(); p++) {
			    out.writeFloat((float)upwardTarget);
			}
		    
		    } // if colourByTarget
		    else if (colourBySeedROI) {

			Point3D seedPoint = t.getPoint(t.seedPointIndex());

			int xSeedVoxel = (int)(seedPoint.x / xVoxelDim);
			int ySeedVoxel = (int)(seedPoint.y / yVoxelDim);
			int zSeedVoxel = (int)(seedPoint.z / zVoxelDim);

			int seedROI = (int)scalarVol[xSeedVoxel][ySeedVoxel][zSeedVoxel];

			for (int p = 0; p < t.numberOfPoints(); p++) {
			    out.writeFloat((float)seedROI);
			}
		    
		    }
		    else if (colourByScalar) {
			for (int p = 0; p < t.numberOfPoints(); p++) {

			    Point3D point = t.getPoint(p);
			
			    out.writeFloat((float)scalarImage.valueAt(point));
			}
		    }
		
	    
		}	
	    }
	
	    om.close();
	}
	catch (IOException e) {
	    throw new LoggedException(e);
	}
	
    }




}

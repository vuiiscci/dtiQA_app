package tractography;

import data.*;
import imaging.ImageHeader;
import misc.*;
import numerics.*;

import java.io.*;
import java.util.logging.*;
import java.util.zip.*;

/**
 * Provides an interface for input of Tract objects.
 *
 * @author Philip Cook
 * @version $Id$
 */
public final class TractSource {


    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.tractography.TractSource");

    private final DataInputStream din;

    private Tract nextTract = null;

    private boolean reachedEndOfFile;
    private boolean noMoreTracts;
  
    private final RealMatrix physToVoxTransform;

    private final boolean transformToCaminoSpace;

    private final double xVoxelDim;
    private final double yVoxelDim;
    private final double zVoxelDim;

    /**
     * Tract source reads tracts in raw format.
     *
     *
     */
    public TractSource(String filename) {
        this(filename, null);
    }


    /**
     * Tract source reads tracts in raw format, and translates them from physical to Camino space 
     * for the image described by the header.
     *
     *
     */
    public TractSource(String filename, ImageHeader header) {

	int bufferSize = OutputManager.FILEBUFFERSIZE;

        if (header != null) {
            transformToCaminoSpace = true;
            physToVoxTransform = header.getPhysicalToVoxelTransform();

            xVoxelDim = header.xVoxelDim();
            yVoxelDim = header.yVoxelDim();
            zVoxelDim = header.zVoxelDim();

        }
        else {
            transformToCaminoSpace = false;
            physToVoxTransform = null;

            xVoxelDim = 0.0;
            yVoxelDim = 0.0;
            zVoxelDim = 0.0;
        }

	try {
	    if (filename == null) {
		din = new DataInputStream(new BufferedInputStream(System.in, bufferSize));
        	logger.info("reading data from standard input");
	    }
	    else {
		
		if (filename.endsWith(".gz")) {
		    FileInputStream fin = new FileInputStream(filename);
		    
		    din = 
			new DataInputStream(new GZIPInputStream(new BufferedInputStream(fin, bufferSize)));
		}
		else {
		    FileInputStream fin = new FileInputStream(filename);
		    
		    din = new DataInputStream(new BufferedInputStream(fin, bufferSize));
		    
		}
	    }
	}
	catch (IOException e) {
	    throw new LoggedException(e);
	}

	init();

	
    }


    private void init() {
        reachedEndOfFile = false;
        noMoreTracts = false;

        // Read in the first voxel ready for the first data request.
        try {
            readNextTract();
	    
	    if (reachedEndOfFile) {
		noMoreTracts = true;
	    }

        }
        catch (Exception e) {
            throw new LoggedException(e);
        }

    }

    private void readNextTract() {
        nextTract = readTractFromRaw();
    }
    

    private Tract readTractFromRaw() {
	
	boolean gotWholeTract = true;

	try {
            
            // (PAC) Would be more efficient to just read all points then construct tract directly
            // Need to profile to assess priority 


	    Tract zeroToSeed = new Tract(100, 100.0); // points 0 to seedPointIndex
	    Tract seedToEnd = new Tract(100, 100.0); // points seedPointIndex to end
	
	    int numPoints = (int)din.readFloat();  // should get EOF here after last tract

	    gotWholeTract = false;
	
	    int seedPointIndex = (int)din.readFloat();
	
	    Point3D[] points = new Point3D[numPoints];
	
	    for (int p = 0; p < numPoints; p++) {
		double x = din.readFloat();
		double y = din.readFloat();
		double z = din.readFloat();
	    
		points[p] = new Point3D(x,y,z);
	    }
	    gotWholeTract = true;
	
	    for (int p = seedPointIndex; p >= 0; p--) {
		zeroToSeed.addPoint(points[p]);
	    }
	
	    for (int p = seedPointIndex; p < numPoints; p++) {
		seedToEnd.addPoint(points[p]);
	    }
	
	    seedToEnd.joinTract(zeroToSeed);

	    return seedToEnd;

	}
	catch (Exception e) {
	    if (e instanceof EOFException) {
		if (!gotWholeTract) {
		    throw new LoggedException("EOF before whole tract was read");
		}
		reachedEndOfFile = true;
		return  null;
	    }
	    else {
		throw new LoggedException("Unexpected " + e + " . Check that input is in RAW format");
	    }
	}
	
    }


    
    /**
     * 
     *
     *  @throws LoggedException if more tracts are read than the file contains.
     */
    public Tract nextTract() {

	if (noMoreTracts) {
	    throw new LoggedException("No more tracts in input");
	}

	Tract current = nextTract;
	
	readNextTract();

	if (reachedEndOfFile) {
            noMoreTracts = true;
        }

	if (transformToCaminoSpace) {
            current.transformToCaminoSpace(physToVoxTransform, xVoxelDim, yVoxelDim, zVoxelDim);
        }

        return current;
    }


    /**
     * @return true if there are more tracts. When this returns false, further calls
     *
     */
    public boolean more() {
	return !noMoreTracts;
    }

  
}

package apps;

import java.io.*;
import java.util.logging.Logger;

import tools.*;
import data.*;
import sphfunc.*;
import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Generates a bitmap showing the spherical function in each voxel
 * of a region of the input data.
 * 
 * <dt>Description:
 * 
 * <dd>A bitmap image is divided into equal boxes for each voxel in
 * the region to be displayed. In each voxel, we compute the value of
 * the function for every element of a spherical point set. We scale
 * each point by the function value and scale all the points to fit in
 * a box with vertices (1, 1, 1) and (-1, -1, -1).
 * 
 * We then project the points to 2D and add a translation
 * corresponding to the voxel coordinates. We scale the points to fit
 * in a fraction of the image box allocated to that voxel. For each
 * point, we increment the corresponding image array element.
 * 
 * Finally we output the image as a byte array (.gray format).
 * 
 * </dl>
 * 
 * @author Danny Alexander
 *
 * $Id$
 *  
 */
public class SphFuncBitMap extends Executable{

    public SphFuncBitMap(String[] args){
    	super(args);
    }
	
    private static Logger logger = Logger.getLogger("camino.apps.SphFuncBitMap");

    /**
     * Buffer size in input and output stream.
     */
    public static int BUFFERSIZE;// = 1000000;

    // Parameters of the bitmap image.

    /**
     * Number of pixels in x-direction of each minifig.
     */
    private static int miniFigSizeX;// = 15;

    /**
     * Number of pixels in y-direction of each minifig.
     */
    private static int miniFigSizeY;// = 15;

    /**
     * Number of pixels separating minifigs in x-direction.
     */
    private static int miniFigSeparationX;// = 1;

    /**
     * Number of pixels separating minifigs in y-direction.
     */
    private static int miniFigSeparationY;// = 1;

    /**
     * Total pixels per input voxel in x-direction.
     */
    private static int pixelsPerVoxelX;// = miniFigSizeX + miniFigSeparationX;

    /**
     * Total pixels per input voxel in y-direction.
     */
    private static int pixelsPerVoxelY;// = miniFigSizeY + miniFigSeparationY;

    /**
     * Orientation of the slice. 0 indicates an axial slice; 1 a coronal
     * slice; 2 a sagittal slice.
     *
     * This is obsolete now.  Use PROJ array instead.
     */
    private static int ORIENTATION;// = 0;

    /**
     * The ORIENTATION variable is old and the projection is better
     * specified using the PROJ array.  The first element specifies
     * which projection of the spherical functions go along the image
     * x-axis, the second the image y-axis.  1 means project along
     * the x direction, 2 along y, 3 along z.  Negate the element to
     * reverse the direction of any axis.
     */
    private static int[] PROJ = {1, 2};

    /**
     * Defines the colour coding.  1 corresponds to x, 2 to y and 3 to
     * z.  The positions in the array correspond to red, green and
     * blue, respectively.  So the default is to code y with red, x
     * with green and z with blue.
     */
    private static int[] COLCODE = {2, 1, 3};

    //Parameters of plot

    /**
     * Index of the point set to use.
     */
    private static int pointSetInd;// = 0;

    /**
     * Normalizing diffusivity for plotting diffusion tensors.
     */
    private static double maxD;// = 1.5E-9;

    /**
     * The program raises the spherical function value by this power
     * to emphasize or deemphasize peaks.
     */
    private static double powerScale;// = 1.0;

    //Parameters of the input data.
 
    /**
     * x-size of input data array.
     */
    private static int xSize;// = 128;

    /**
     * y-size of input data array.
     */
    private static int ySize;// = 128;

    /**
     * The array that holds the bitmap.
     */
    private static char[][][] imArray;

    /**
     * Array holding the colour for the sphfunc icons.
     */
    private static char[] col = {255, 255, 255};

    /**
     * Flag indicating whether to colour code the direction of
     * each plotted point.
     */
    private static boolean colCodeDirections;

    //Defines the bounding box of the region of interest.
    private int xMin;
    private int xMax;
    private int yMin;
    private int yMax;

    /**
     * Output manager
     */
    //    private static OutputManager om;

    //Specifies the interval between voxels to be displayed.
    //When set to one, every voxel is displayed. If n, then
    //only every n-th voxel (in each axis).
    private int interval;

    //If set, a small point set (elec120) is
    //used to get a coarse picture.
    private boolean smallPointSet;

    // Specify the type of normalization of the spherical functions
    // to do.  If this is set to true, the display uses min-max
    // normalization as in Tuch's displays.  Otherwise it just
    // scales each function to the range [0, 1].
    private boolean minmaxnorm;

    //Specify the type of the output image.
    private boolean grayscale;

    // Filename for the image to load in and use as background
    private String backDropFile;

    // String specifying the kind of interpolation to use in the
    // backdrop image.
    private String backInterp;

    // A spherical function may be negative in certain directions.
    // We ignore points where the function is negative and do
    // not plot them.
    //tmp commented out       logger.info("Ignoring points with radius less than zero");

		
    public void initDefaultVals(){
        BUFFERSIZE = 1000000;
        miniFigSizeX = 15;
        miniFigSizeY = 15;
        miniFigSeparationX = 1;
        miniFigSeparationY = 1;
        ORIENTATION = 0;
        //int[] PROJ={1,2};
        int[] COLCODE = {2, 1, 3};
        pointSetInd = 0;
        maxD = 1.5E-9;
        powerScale = 1.0;
        xSize = 128;
        ySize = 128;
        char[] col = {255, 255, 255};
        pixelsPerVoxelX = miniFigSizeX + miniFigSeparationX;
        pixelsPerVoxelY = miniFigSizeY + miniFigSeparationY;
        xMin = -1;
        xMax = -1;
        yMin = -1;
        yMax = -1;
        interval = 1;
        smallPointSet = false;
        minmaxnorm = false;
        grayscale = true;
        backDropFile = null;
        backInterp = "bilinear";
        //logger.info("Ignoring points with radius less than zero");	
    }
	
	
    public void initOptions(String[] args) {

        // Get general command line arguments.
        CL_Initializer.inputDataType = "double";
        CL_Initializer.CL_init(args);

        if(CL_Initializer.schemeFile != null)
            CL_Initializer.initImagingScheme();
        
        // Check the max ent control points are known
        if (CL_Initializer.inputModel.equals("maxent"))
            if(CL_Initializer.mePointSet < 0 && CL_Initializer.schemeFile == null) {
                throw new LoggedException("Maximum entropy control points must be specified by a schemefile or -mepointset");
            }
            else {
                CL_Initializer.initMaxEnt();
            }
        

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-pointset")) {
                pointSetInd = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
            if (args[i].equals("-orientation")) {
                ORIENTATION = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
            /*			System.err.println(PROJ[0]);
                                System.err.println(PROJ[1]);*/
            if (args[i].equals("-projection")) {
                PROJ[0] = Integer.parseInt(args[i + 1]);
                PROJ[1] = Integer.parseInt(args[i + 2]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
                CL_Initializer.markAsParsed(i + 2);
            }
            if (args[i].equals("-colcode")) {
                COLCODE[0] = Integer.parseInt(args[i + 1]);
                COLCODE[1] = Integer.parseInt(args[i + 2]);
                COLCODE[2] = Integer.parseInt(args[i + 3]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
                CL_Initializer.markAsParsed(i + 2);
                CL_Initializer.markAsParsed(i + 3);
            }
            if (args[i].equals("-xsize")) {
                xSize = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
            if (args[i].equals("-ysize")) {
                ySize = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
            if (args[i].equals("-minifigsize")) {
                miniFigSizeX = Integer.parseInt(args[i + 1]);
                miniFigSizeY = Integer.parseInt(args[i + 2]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
                CL_Initializer.markAsParsed(i + 2);
            }
            if (args[i].equals("-minifigseparation")) {
                miniFigSeparationX = Integer.parseInt(args[i + 1]);
                miniFigSeparationY = Integer.parseInt(args[i + 2]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
                CL_Initializer.markAsParsed(i + 2);
            }
            if (args[i].equals("-box")) {
                xMin = Integer.parseInt(args[i + 1]);
                xMax = Integer.parseInt(args[i + 2]);
                yMin = Integer.parseInt(args[i + 3]);
                yMax = Integer.parseInt(args[i + 4]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
                CL_Initializer.markAsParsed(i + 2);
                CL_Initializer.markAsParsed(i + 3);
                CL_Initializer.markAsParsed(i + 4);
            }
            if (args[i].equals("-coarse")) {
                smallPointSet = true;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-minmaxnorm")) {
                minmaxnorm = true;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-iconcol")) {
                col[0] = (char)Integer.parseInt(args[i + 1]);
                col[1] = (char)Integer.parseInt(args[i + 2]);
                col[2] = (char)Integer.parseInt(args[i + 3]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
                CL_Initializer.markAsParsed(i + 2);
                CL_Initializer.markAsParsed(i + 3);
                grayscale = false;
            }
            if (args[i].equals("-dircolcode")) {
                colCodeDirections = true;
                grayscale = false;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-backdrop")) {
                backDropFile = args[i+1];
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i+1);
            }
            if (args[i].equals("-backdropinterp")) {
                backInterp = args[i+1];
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i+1);
            }
            if (args[i].equals("-maxd")) {
                maxD = Double.parseDouble(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
            if (args[i].equals("-powerscale")) {
                powerScale = Double.parseDouble(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
            if (args[i].equals("-interval")) {
                interval = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i+1);
            }
        }

        CL_Initializer.checkParsing(args);
        /*        System.err.println(PROJ[0]);
                  System.err.println(PROJ[1]);*/
        // om = new OutputManager();
    }
	
    public void initVariables(){
        if(CL_Initializer.inputModel.equals("dt") ||
           CL_Initializer.inputModel.equals("twotensor") ||
           CL_Initializer.inputModel.equals("threetensor") ||
           CL_Initializer.inputModel.equals("multitensor"))

            CL_Initializer.initTensorDataSource();
        else
            CL_Initializer.initSphFuncDataSource();
		
        // Figure out the range of voxels to include.
        if (xMax <= 0 || xMax > xSize - 1) {
            xMax = xSize - 1;
        }
        if (yMax <= 0 || yMax > ySize - 1) {
            yMax = ySize - 1;
        }
        if (xMin <= 0 || xMin > xSize) {
            xMin = 0;
        }
        if (yMin <= 0 || yMin > ySize) {
            yMin = 0;
        }

        //Construct the bitmap array.
        pixelsPerVoxelX = miniFigSizeX + miniFigSeparationX;
        pixelsPerVoxelY = miniFigSizeY + miniFigSeparationY;
        int xVoxRes = (xMax + 1 - xMin) / interval;
        int yVoxRes = (yMax + 1 - yMin) / interval;
        int imSizeX = xVoxRes * pixelsPerVoxelX;
        int imSizeY = yVoxRes * pixelsPerVoxelY;

        imArray = new char[imSizeX][imSizeY][3];
    }
	
    public void execute(OutputManager om){	

        // Read in the backdrop image.
        if(backDropFile != null) {
            char[][] backDropImage = readBackDrop(backDropFile, xSize, ySize);

            // Initialize the output array with the backdrop.
            initOutputArray(backDropImage, interval, xMax, yMin, pixelsPerVoxelX, pixelsPerVoxelY, backInterp);
        }

        // Now create icons in each pixel from the input data.
        try {

            // Load in the point sets.
            SphericalPointSet[] ps = new SphericalPointSet[ISCodes
                                                           .getNoPointSetsForMaxEnt()];
            for (int i = 0; i < ps.length; i++) {
                ps[i] = ISCodes.getPointSetForMaxEnt(i);
            }

            for (int x = 0; x < xSize; x++) {
                for (int y = 0; y < ySize; y++) {

                    // Read in the coefficients.
                    double[] coeffs = CL_Initializer.data.nextVoxel();

                    // Check for background voxels and plot nothing
                    // there.
                    if(coeffs[0] >= 0) {

                        // Construct a spherical function using the
                        // coefficients.
                        SphericalFunction s = null;
                        DT[] ds = null;
                        if (CL_Initializer.inputModel.equals("sh")) {
                            s = new EvenSHS(coeffs, CL_Initializer.maxOrder);
                        }
                        else if (CL_Initializer.inputModel.equals("dt")) {
                            DT d = new DT(coeffs[2], coeffs[3], coeffs[4], coeffs[5], coeffs[6], coeffs[7]);
                            ds = new DT[1];
                            ds[0] = d;
                        }
                        else if(CL_Initializer.inputModel.equals("twotensor")) {
                            ds = new DT[2];

                            // Need to weight by volume fractions and
                            // then normalize by max volume fraction.
                            // Accomplish this just by scaling each DT
                            // by appropriate factors.
                            double maxVF = coeffs[3]>coeffs[10]?coeffs[3]:coeffs[10];
                            double sc1 = coeffs[3]/maxVF;
                            double sc2 = coeffs[10]/maxVF;

                            ds[0] = new DT(sc1*coeffs[4], sc1*coeffs[5], sc1*coeffs[6], sc1*coeffs[7], sc1*coeffs[8], sc1*coeffs[9]);
                            ds[1] = new DT(sc2*coeffs[11], sc2*coeffs[12], sc2*coeffs[13], sc2*coeffs[14], sc2*coeffs[15], sc2*coeffs[16]);
                        }
                        else if(CL_Initializer.inputModel.equals("threetensor")) {
                            ds = new DT[3];

                            double maxVF = coeffs[3];
                            if(coeffs[10]>maxVF)
                                maxVF = coeffs[10];
                            if(coeffs[17]>maxVF)
                                maxVF = coeffs[17];

                            double sc1 = coeffs[3]/maxVF;
                            double sc2 = coeffs[10]/maxVF;
                            double sc3 = coeffs[17]/maxVF;

                            ds[0] = new DT(sc1*coeffs[4], sc1*coeffs[5], sc1*coeffs[6], sc1*coeffs[7], sc1*coeffs[8], sc1*coeffs[9]);
                            ds[1] = new DT(sc2*coeffs[11], sc2*coeffs[12], sc2*coeffs[13], sc2*coeffs[14], sc2*coeffs[15], sc2*coeffs[16]);
                            ds[2] = new DT(sc3*coeffs[18], sc3*coeffs[19], sc3*coeffs[20], sc3*coeffs[21], sc3*coeffs[22], sc3*coeffs[23]);
                        }
                        else if(CL_Initializer.inputModel.equals("multitensor")) {
                            // Element 2 of the voxel array is the number of
                            // fitted components.
                            ds = new DT[(int)coeffs[2]];

                            // Find all the scaling factors and their max.
                            double[] sc = new double[(int)coeffs[2]];
                            double maxVF = coeffs[3];
                            for(int i=0; i<sc.length; i++) {
                                sc[i] = coeffs[3+7*i];
                                if(sc[i]>maxVF)
                                    maxVF = sc[i];
                            }
                            for(int i=0; i<sc.length; i++)
                                sc[i] = sc[i]/maxVF;                           
                            
                            for(int i=0; i<ds.length; i++) {
                                ds[i] = new DT(sc[i]*coeffs[4+7*i], sc[i]*coeffs[5+7*i], sc[i]*coeffs[6+7*i], sc[i]*coeffs[7+7*i], sc[i]*coeffs[8+7*i], sc[i]*coeffs[9+7*i]);
                            }
                        }
                        else if (CL_Initializer.inputModel.equals("maxent")) {
                            s = new MaxEntProfile(coeffs, CL_Initializer.kernelParams);
                        }
                        else {
                            s = new TuchRBF_Sum(coeffs);
                        }
                        
                        if (x >= xMin && y >= yMin && x <= xMax && y <= yMax
                            && (x % interval) == 0 && (y % interval) == 0) {
                            
                            //Get the point set for the plot.
                            double[][] pointSet = getPointSet(smallPointSet, ps);
                            
                            //Plot the function on the bitmap array.
                            char[][][] miniFig;
                            if(CL_Initializer.inputModel.equals("dt") ||
                               CL_Initializer.inputModel.equals("twotensor") ||
                               CL_Initializer.inputModel.equals("threetensor") ||
                               CL_Initializer.inputModel.equals("multitensor")) {
                                miniFig = getMiniFigXY_DTs(ds, pointSet);
                            }
                            else {
                                miniFig = getMiniFigXY(s, pointSet, minmaxnorm);
                            }
                                
                            addMiniFig((x - xMin) / interval, (y - yMin) / interval, miniFig);
                            
                        } // if in image range.
                    } // if background.
                } //y
            } //x
        }
        catch (Exception e) {
            
            LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
        }

        //Finally output the image array in binary form.
        System.err.println("Image size: " + imArray[0].length + "x" + imArray.length);
        outputBitMap(om.getOutputStream(), grayscale);

        // Tidy up.
        om.close();
    }


    /**
     * Reads in the backdrop image from the specified file.
     *
     * @param backDropFile Name of file containing the backdrop image.
     *
     * @param xSize Number of pixels in x direction.
     *
     * @param ySize Number of pixels in y direction.
     *
     * @return Array containing the image scaled to the range [0,255].
     */
    private static char[][] readBackDrop(String backDropFile, int xSize, int ySize) {
        logger.info("Loading " + backDropFile + ".  Expecting type " + CL_Initializer.inputDataType + ".");

        double[][] dBackdrop = new double[xSize][ySize];
        DataSource vin = ExternalDataSource.getDataSource(backDropFile, 1, CL_Initializer.inputDataType);

        double maxBD = 0.0;
        double minBD = 0.0;
        for(int i=0; i<xSize; i++) {
            for(int j=0; j<ySize; j++) try {
                    dBackdrop[i][j] = (double)(vin.nextVoxel()[0]);
                    maxBD = (i==0)?dBackdrop[i][j]:((dBackdrop[i][j] > maxBD) ? dBackdrop[i][j] : maxBD);
                    minBD = (i==0)?dBackdrop[i][j]:((dBackdrop[i][j] < minBD) ? dBackdrop[i][j] : minBD);
                
                }
                catch(Exception e) {
                    logger.warning("Backdrop image smaller than expected.  Continuing...");
                    break;
                }
        }

        // Normalize and convert to char array.  We also flip in the
        // x direction for consistency with the minifig array.
        char[][] cBackdrop = new char[xSize][ySize];
        if(maxBD>minBD) {
            for(int i=0; i<xSize; i++) {
                for(int j=0; j<ySize; j++) {
                    cBackdrop[i][j] = (char)(255.0*(dBackdrop[i][j] - minBD)/(maxBD - minBD));
                }
            }
        }

        return cBackdrop;
    }
                

    /**
     * Initializes the output array by putting in the backdrop image.
     *
     * @param backDropImage The backdrop image array.
     *
     * @param interval The interval between voxels to display in the output.
     *
     * @param xMax The maximum x coordinate to display.
     *
     * @param yMin The minimum y coordinate to display.
     *
     * @param pixelsPerVoxelX Minifig x size.
     *
     * @param pixelsPerVoxelY Minifig y size.
     *
     * @param backInterp String specifying the kind of interpolation.
     * Either "nn" (for nearest neighbour) or "bilinear".
     */
    private static void initOutputArray(char[][] backDropImage, int interval, int xMax, int yMin, int pixelsPerVoxelX, int pixelsPerVoxelY, String backInterp) {
        for(int i=0; i<imArray.length; i++) {
            for(int j=0; j<imArray[0].length; j++) {

                // Compute the corresponding indices to the backdrop image.
                // We need to account for the fact that the array is
                // flipped in the x direction
                double bdx = xMax - ((double)((i+pixelsPerVoxelX/2)*interval)/(double)pixelsPerVoxelX) + 1;
                double bdy = ((double)((j-pixelsPerVoxelY/2)*interval)/(double)pixelsPerVoxelY + yMin);

                //System.err.println(i + " " + j + " " + bdx + " " + bdy + " " + (int)(bdx + 0.5) + " " + (int)(bdy + 0.5));

                if(backInterp.equals("nn")) {
                    // Use nearest neighbour interpolation.
                    for(int k=0; k<3; k++) {
                        imArray[i][j][k] = backDropImage[(int)bdx][(int)bdy];
                    }
                }
                else {
                    // Use bilinear.
                    int ix = (int)bdx;
                    int iy = (int)bdy;
                    int ipx = ((ix + 1)<backDropImage.length)?(ix + 1):ix;
                    int ipy = ((iy + 1)<backDropImage[0].length)?(iy + 1):iy;

                    double dx = bdx - (double)ix;
                    double dy = bdy - (double)iy;

                    double gl = (1-dx)*(1-dy)*backDropImage[ix][iy]
                        + dx*(1-dy)*backDropImage[ipx][iy]
                        + (1-dx)*dy*backDropImage[ix][ipy]
                        + dx*dy*backDropImage[ipx][ipy];

                    for(int k=0; k<3; k++) {
                        imArray[i][j][k] = (char)gl;
                    }
                }
            }
        }
    }                
        

    /**
     * Returns the point set to be used for the plot.
     *
     * @param smallPointSet A Boolean indicating whether the -coarse
     * option was specified.
     *
     * @param ps List of point sets available.
     *
     * @return The point set to use.
     */
    private static double[][] getPointSet(boolean smallPointSet, SphericalPointSet[] ps) {

        //Just use the smallest point set.
        double[][] pointSet = ps[pointSetInd].data;
        if (smallPointSet) {
            pointSet = SphericalPoints.getElecPointSet(120);
        }

        return pointSet;
    }

    /**
     * Creates a minifig bitmap for a single spherical function.
     *
     * @param s The spherical function to plot.
     *
     * @param pointSet The point set to use to generate the plot.
     *
     * @param minmaxnorm Boolean indicating whether or not to do
     * minmax scaling.
     *
     * @return An RGB minifig containing the plot.
     */
    private static char[][][] getMiniFigXY(SphericalFunction s, double[][] pointSet, boolean minmaxnorm) {

        char[][][] miniFig = new char[miniFigSizeX][miniFigSizeY][3];

        double[] rs = new double[pointSet.length];
        double maxR = 0.0;
        double minR = 0.0;
        for (int i = 0; i < pointSet.length; i++) {

            rs[i] = s.getRadius(pointSet[i][0], pointSet[i][1], pointSet[i][2]);
            // Do the power scaling if required.
            if (powerScale != 1) {
                rs[i] = Math.pow(rs[i], powerScale);
            }

            maxR = (rs[i] > maxR) ? rs[i] : maxR;
            minR = (i==0)?rs[i]:((rs[i] < minR) ? rs[i] : minR);
        }

	// Since we ignore points with negative radius, truncate
	// minR at zero.
	minR = (minR<0.0)?0.0:minR;

        for (int i = 0; i < pointSet.length; i++) {
            double r = rs[i] / maxR;
	    if(minmaxnorm) {
		r = (maxR>minR)?(rs[i] - minR)/(maxR - minR):0.0;
	    }

            //Map the coordinates to the range
            //[0, 1].
            //Also here, we swap the x and y
            //coordinates to make the PAS orientation
            //consistent with the image.
            if (r > 0) {
                double yPt = pointSet[i][0] * r;
                double xPt = pointSet[i][1] * r;
                double zPt = pointSet[i][2] * r;

                addPointToMiniFig(xPt, yPt, zPt, miniFig);
            }
            else {
                if(r < 0) logger.info("Radius less than zero.  Ignoring.");
            }
        } //i

        return miniFig;
    }

    /**
     * Creates a minifig bitmap for a single diffusion tensor.
     *
     * @param ds List of diffusion tensors
     *
     * @param pointSet The point set used to generate the plot.
     *
     * @return An RGB minifig containing the plot.
     */
    private static char[][][] getMiniFigXY_DTs(DT[] ds, double[][] pointSet) {

        char[][][] miniFig = new char[miniFigSizeX][miniFigSizeY][3];

        for(int j=0; j<ds.length; j++) {

            DT d = ds[j];

            double[][] newPts = new double[pointSet.length][3];
            for (int i = 0; i < pointSet.length; i++) {
                newPts[i] = d.multiply(pointSet[i]);
            }

            for (int i = 0; i < pointSet.length; i++) {

                // Map the coordinates to a prespecified range.
                // Also here, we swap the x and y
                // coordinates to make the function orientation
                // consistent with the image.
                double yPt = newPts[i][0] / maxD;
                double xPt = newPts[i][1] / maxD;
                double zPt = newPts[i][2] / maxD;

                addPointToMiniFig(xPt, yPt, zPt, miniFig);
            }
        }

        return miniFig;
    }

    /**
     * Adds a point to the image array.
     *
     * @param xPt The x coordinate of the point on the spherical
     * function.
     *
     * @param yPt The y coordinate of the point on the spherical
     * function.
     *
     * @param zPt The z coordinate of the point on the spherical
     * function.
     * 
     * @param miniFit The minifig array to plot the point in.
     */
    private static void addPointToMiniFig(double xPt, double yPt, double zPt,
                                          char[][][] miniFig) {

        //Project the three dimensional point to the two dimensions of the
        //image.
        double xPt2D = xPt;
        double yPt2D = yPt;
        //System.err.println(PROJ[0]);
        // The projection is specified in the PROJ array.
        if(Math.abs(PROJ[0]) == 1) {
            xPt2D = (PROJ[0]>0)?xPt:-xPt;
        }
        else if(Math.abs(PROJ[0]) == 2) {
            xPt2D = (PROJ[0]>0)?yPt:-yPt;
        }
        else if(Math.abs(PROJ[0]) == 3) {
            xPt2D = (PROJ[0]>0)?zPt:-zPt;
        }
        if(Math.abs(PROJ[1]) == 1) {
            yPt2D = (PROJ[1]>0)?xPt:-xPt;
        }
        else if(Math.abs(PROJ[1]) == 2) {
            yPt2D = (PROJ[1]>0)?yPt:-yPt;
        }
        else if(Math.abs(PROJ[1]) == 3) {
            yPt2D = (PROJ[1]>0)?zPt:-zPt;
        }

        // This remains for compatibility with old versions but
        // the -orientation option has been replaced by the more
        // general -projection, which sets the PROJ array instead.
        if (ORIENTATION == 1) {
            xPt2D = zPt;
            yPt2D = yPt;
        }
        if (ORIENTATION == 2) {
            xPt2D = zPt;
            yPt2D = -xPt;
        }
        if (ORIENTATION == 3) {
            xPt2D = -zPt;
            yPt2D = yPt;
        }
        if (ORIENTATION == 4) {
            xPt2D = -zPt;
            yPt2D = -xPt;
        }
        if (ORIENTATION == 5) {
            xPt2D = yPt;
            yPt2D = xPt;
        }
        if (ORIENTATION == 6) {
            xPt2D = zPt;
            yPt2D = yPt;
        }

        //Now find the pixel coordinate corresponding to this
        //point.
        int xMF = (int) (((double) miniFigSizeX) * (xPt2D + 1.0) / 2.0);
        int yMF = (int) (((double) miniFigSizeY) * (yPt2D + 1.0) / 2.0);

        //Check we have not gone out of range of the mini-fig.
        if(xMF >= 0 && yMF >= 0 && xMF < miniFigSizeX && yMF < miniFigSizeY) {

            // Get the RGB colour for this pixel
            char[] rgb = getPixelColour(xPt, yPt, zPt);
            
            //Set the minifig pixel.
            for(int k=0; k<3; k++) {
                // We choose the brightest possible colour in each pixel.
                // This is a short cut to speed up processing where we
                // choose the max in each channel separately.
                miniFig[xMF][yMF][k] = (rgb[k]>miniFig[xMF][yMF][k])?rgb[k]:miniFig[xMF][yMF][k];
            }
        }
    }

    /**
     * Adds a minifig to the main image.
     *
     * @param x The x-position of the voxel in the array of minifigs.
     * 
     * @param y The y-position of the voxel in the array of minifigs.
     * 
     * @param miniFig The minifig to add at position x y.
     */
    private static void addMiniFig(int x, int y, char[][][] miniFig) {
        int baseX = imArray.length - pixelsPerVoxelX * (x + 1);
        int baseY = pixelsPerVoxelY * y;

        for (int i = 0; i < miniFig.length; i++) {
            for (int j = 0; j < miniFig[0].length; j++) {

                // Only add points in the sphfunc icon.
		if(miniFig[i][j][0] > 0 || miniFig[i][j][1] > 0 || miniFig[i][j][2] > 0) 
                    for(int k=0; k<3; k++) {
                        imArray[baseX + i][baseY + j][k] = miniFig[i][j][k];
                    }
            }
        }
    }

    /**
     * Returns the colour for a pixel depending on the 3D coordinates
     * that project to it.
     * 
     * @param x The x coordinate of the spherical function point.
     * 
     * @param y The y coordinate of the spherical function point.
     * 
     * @param z The z coordinate of the spherical function point.
     * 
     * @return The RGB vector.
     */
    private static char[] getPixelColour(double x, double y, double z) {
        char[] rgb = new char[3];

        // For consistency with image orientation, x and y need to
        // be transposed.
        if(colCodeDirections) {
            if(COLCODE[0] == 1) {
                rgb[0] = (char)(Math.abs(x)*255.0);
            }
            else if(COLCODE[0] == 2) {
                rgb[0] = (char)(Math.abs(y)*255.0);
            }
            else {
                rgb[0] = (char)(Math.abs(z)*255.0);
            }
            
            if(COLCODE[1] == 1) {
                rgb[1] = (char)(Math.abs(x)*255.0);
            }
            else if(COLCODE[1] == 2) {
                rgb[1] = (char)(Math.abs(y)*255.0);
            }
            else {
                rgb[1] = (char)(Math.abs(z)*255.0);
            }
            
            if(COLCODE[2] == 1) {
                rgb[2] = (char)(Math.abs(x)*255.0);
            }
            else if(COLCODE[2] == 2) {
                rgb[2] = (char)(Math.abs(y)*255.0);
            }
            else {
                rgb[2] = (char)(Math.abs(z)*255.0);
            }
            
        }
        else {
            rgb[0] = col[0];
            rgb[1] = col[1];
            rgb[2] = col[2];
        }

        return rgb;
    }


    /**
     * Outputs the bitmap image to a file.
     * 
     * @param out The output stream for the data.
     * 
     * @param grayscale A boolean indicating whether the output image
     * is grayscale of colour.
     */
    public static void outputBitMap(DataOutputStream out, boolean grayscale) {
        try {

            for (int i = 0; i < imArray.length; i++) {
                for (int j = 0; j < imArray[0].length; j++) {
                    char outputValue = (char) imArray[i][j][0];

                    //If the grayscale flag is not specified, we make the
                    //image binary.
                    int maxChar = 255;
		    if(grayscale) {
                        // Plot back on white to save ink!
			outputValue = (char) (maxChar - (int) outputValue);
			out.writeByte(outputValue);
		    }
		    else {
			for(int k=0; k<3; k++) {
			    outputValue = (char)imArray[i][j][k];
			    out.writeByte(outputValue);
			}
		    }
                }
            }
        }
        catch (Exception e) {
            
            LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
        }
    }


}


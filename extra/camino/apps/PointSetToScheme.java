package apps;

import data.*;
import imaging.*;
import misc.*;
import tools.*;

import java.util.*;
import java.util.logging.Logger;
import java.util.regex.*;
import java.io.*;
import java.text.*;

/**
 * Converts point sets to scheme files, given additional imaging parameters.
 *
 * @author Philip Cook
 * @version $Id$
 */
public class PointSetToScheme extends Executable{
    
	private static DecimalFormat df;
    private static DecimalFormat bValDF;
    
    private static final Logger logger = Logger.getLogger("apps.PointSetToScheme");
    
    public PointSetToScheme(String[] args){
        super(args);
    }
    
	// flags to flip directions
    private boolean flipX;// = false;
    private boolean flipY;// = false;
    private boolean flipZ;// = false;
	
    // number of zero measurements to prepend onto scheme
    private int M;// = 0;
        
    // number of times all measurements are acquired
    private int numScans;// = 1;
        
    // if true, repeat measurements are stored together
    private boolean interleave;// = false;
	
    // params expected for each scheme type
	private int bVectorParams;// = 1;
    private int rectGradST_Params;// = 4;
	private int rectGradTRSE_Params;// = 8;

	// Some scanners define b-value as |g|^2 * \beta
    // where \beta is what they CLAIM the b-value is.
    // If you use normalized gradient directions, you need to increase b accordingly
    // to make the same ADC calculation.
    private boolean useGradMod;

	private double[] dwPars;

    private String schemeVersion;// = B_VectorScheme.VERSION;		
	
	public void initDefaultVals(){
		flipX = false;
		flipY = false;
		flipZ = false;       
		
		M = 0;
        numScans = 1;
        interleave = false;
		
		// params expected for each scheme type
        bVectorParams = 1;
        rectGradST_Params = 4;
        rectGradTRSE_Params = 8;
        useGradMod = false;

        // insist on UK style numbers, no comma decimal points
        Locale.setDefault(Locale.UK);
        
        // now that the locale is set, the decimalformats can be constructed
        df = new DecimalFormat("   0.000000;  -0.000000");
        bValDF = new DecimalFormat("   0.000E00");
                
        
		double[] dwPars = new double[0];
		schemeVersion = B_VectorScheme.VERSION;		
    }
    
    public void initOptions(String[] args){
        CL_Initializer.CL_init(args);
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-bvalue")) {
                dwPars = new double[1];
                
                dwPars[0] = Double.parseDouble(args[i+1]);
                
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-usegradmod")) {
                useGradMod = true;
                logger.info("Gradient direction magnitude will be incorporated into b-values");
                CL_Initializer.markAsParsed(i);
            }
            else if (args[i].equals("-addzeromeas")) {
                M = Integer.parseInt(args[i+1]);
                
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-flipx")) {
                flipX = true;
                CL_Initializer.markAsParsed(i);
            }
            else if (args[i].equals("-flipy")) {
                flipY = true;
                CL_Initializer.markAsParsed(i);
                
            }
            else if (args[i].equals("-flipz")) {
                flipZ = true;
                CL_Initializer.markAsParsed(i);
            }
            else if (args[i].equals("-numscans")) {
                numScans = Integer.parseInt(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            else if (args[i].equals("-interleave")) {
                interleave = true;
                CL_Initializer.markAsParsed(i);
            }
            else if (args[i].equals("-version")) {
                schemeVersion = args[i+1];
                CL_Initializer.markAsParsed(i,2);
            }
            else if (args[i].equals("-dwparams")) {
                int counter = 0;
                
                double[] tmp = new double[100];
                
                while (i + counter + 1 < args.length && !args[i+counter+1].startsWith("-")) {
                    tmp[counter] = Double.parseDouble(args[i+counter+1]);
                    counter++;
                }
                
                dwPars = new double[counter];
                
                System.arraycopy(tmp, 0, dwPars, 0, counter);
                
                CL_Initializer.markAsParsed(i, counter+1);
                
            }
            
        }
        
        CL_Initializer.checkParsing(args);
    }
    
	public void initVariables(){
 
	}
	
	public void execute(OutputManager om){	    
//	public static void main(String[] args) {
                
        // can scale b-values with grad dir magnitude, but it is unclear how to
        // incorporate it into other scheme versions
        if ( useGradMod && !schemeVersion.equals(B_VectorScheme.VERSION) ) {
            throw new LoggedException("Can only use gradient magnitude with BVECTOR scheme version");
        }        
        try {
            
            if (CL_Initializer.inputFile == null) {
                logger.info("reading data from standard input");
            }
            
            // points not normalized
            double[][] points = readPoints(CL_Initializer.inputFile, flipX, flipY, flipZ);
            
            int N = points.length;
            
            // prepend zero measurements
            double[][] withZeros = new double[N + M][3];
            
            for (int i = 0; i < N; i++) {
                System.arraycopy(points[i], 0, withZeros[i+M], 0, 3);
            }
            
            points = withZeros;
            
            int numMeas = points.length;
            
            double[][] withRepeats = new double[numScans * numMeas][3];
            
            if (numScans > 1) {
                
                if (interleave) {
                    for (int m = 0; m < numMeas; m++) {
                        for (int s = 0; s < numScans; s++) {
                            System.arraycopy(points[m], 0, withRepeats[m * numScans + s], 0, 3);
                        }
                    }
                    
                }
                else {
                    for (int s = 0; s < numScans; s++) {
                        for (int m = 0; m < numMeas; m++) {
                            System.arraycopy(points[m], 0, withRepeats[s * numMeas + m], 0, 3);
                        }
                    }
                }
                
                points = withRepeats;
            }
            
            String scheme = null;
            
            if (schemeVersion.equals(B_VectorScheme.VERSION)) {
                if (dwPars.length != bVectorParams) {
                    throw new LoggedException("Scheme version is " + schemeVersion + ", expected " +
                    bVectorParams + " param, got " + dwPars.length);
                }
                scheme = pointSetToB_VectorScheme(points, useGradMod, dwPars);
            }
            else if (schemeVersion.equals(RectGradSteTanScheme.VERSION)) {
                
                if (dwPars.length != rectGradST_Params) {
                    throw new LoggedException("Scheme version is " + schemeVersion + ", expected " +
                    rectGradST_Params + " params, got " + dwPars.length);
                }
                scheme = pointSetToRectGradSteTanScheme(points, dwPars);
            }
            else if (schemeVersion.equals(RectGradTRSE_Scheme.VERSION)) {
                
                if (dwPars.length != rectGradTRSE_Params) {
                    throw new LoggedException("Scheme version is " + schemeVersion + ", expected " +
                    rectGradTRSE_Params + " params, got " + dwPars.length);
                }
                scheme = pointSetToRectGradTRSE_Scheme(points, dwPars);
            }
            
            
            if (OutputManager.outputFile == null) {
                System.out.print(scheme);
            }
            else {
                
                FileOutput out = new FileOutput(OutputManager.outputFile);
                
                out.writeString(scheme);
                
                out.close();
            }
            
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }     
    }
    
    
    /**
     * Converts a point set and imaging parameters to a String representation of a
     * <code>RectGradStejskalTannerScheme</code>.
     *
     * @param points the gradient directions, including zero gradients.
     * @param dwPars diffusion weighting parameters: {|G|, DELTA, delta, TE}. All of these
     * parameters except TE will be set to zero for any measurement where the gradient direction is zero.
     *
     * @see imaging.RectGradSteTanScheme
     */
    public static String pointSetToRectGradSteTanScheme(double[][] points, double[] dwPars) {
        
        double modG = dwPars[0];
        double bigDel = dwPars[1];
        double smallDel = dwPars[2];
        double TE = dwPars[3];
        
        int numMeas = points.length;
        
        double[] mod = getPointModulus(points);
        
        boolean pointsNormalized = true;
        
        for (int i = 0; i < numMeas; i++) {
            if (mod[i] > 0.0 && !(Math.abs(1.0 - mod[i]) < 1E-5)) {
                pointsNormalized = false;
            }
        }
        
        if (!pointsNormalized) {
            logger.warning("Some gradient directions do not have unit length. Directions normalized in output.");
        }
        
        points = normalizePoints(points);
        
        StringBuffer sb = new StringBuffer();
        
        sb.append("#Stejskal-Tanner scheme with rectangular gradient pulses\n");
        sb.append("#g_x\tg_y\tg_z\t|g|\tDELTA\tdelta\tTE\n");
        sb.append("VERSION: " + RectGradSteTanScheme.VERSION);
        sb.append("\n");
        
        for (int i = 0; i < numMeas; i++) {
            
            sb.append(df.format(points[i][0]));
            sb.append(df.format(points[i][1]));
            sb.append(df.format(points[i][2]));
            
            if (points[i][0] == 0.0 && points[i][1] == 0.0 && points[i][2] == 0.0) {
                sb.append(df.format(0.0));
                sb.append(df.format(0.0));
                sb.append(df.format(0.0));
                sb.append(df.format(TE));
            }
            else {
                sb.append(df.format(modG));
                sb.append(df.format(bigDel));
                sb.append(df.format(smallDel));
                sb.append(df.format(TE));
            }
            
            sb.append("\n");
            
        }
        
        return sb.toString();
        
    }
    
    
    /**
     * Converts a point set and imaging parameters to a String representation of a
     * <code>B_VectorScheme</code>.
     *
     * @param points - the gradient directions, including zero gradients.
     * @param useGradMod - if true, use the gradient direction modulus to scale the b-value.
     * @param dwPars diffusion weighting parameters: {b}. The b-value will be set to zero
     * for any measurement with a zero gradient direction.
     *
     * @see imaging.B_VectorScheme
     */
    public static String pointSetToB_VectorScheme(double[][] points, boolean useGradMod, double[] dwPars) {
        
        double bValue = dwPars[0];
        
        int numMeas = points.length;
        
        StringBuffer sb = new StringBuffer();
        
        sb.append("#B-vector scheme. Contains gradient directions and b-values\n");
        sb.append("#g_x\tg_y\tg_z\tb\n");
        
        sb.append("VERSION: " + B_VectorScheme.VERSION + "\n");
        
        double[] mod = getPointModulus(points);
        
        points = normalizePoints(points);
        
        for (int i = 0; i < numMeas; i++) {
            
            sb.append(df.format(points[i][0]));
            sb.append(df.format(points[i][1]));
            sb.append(df.format(points[i][2]));
            
            if (points[i][0] == 0.0 && points[i][1] == 0.0 && points[i][2] == 0.0) {
                sb.append(bValDF.format(0.0));
            }
            else if (useGradMod) {
                sb.append(bValDF.format(bValue * mod[i] * mod[i]));
            }
            else {
                sb.append(bValDF.format(bValue));
            }
            
            sb.append("\n");
        }
        
        
        return sb.toString();
    }
    
    
    /**
     * Converts a point set and imaging parameters to a String representation of a
     * <code>RectGradTRSE_Scheme</code>.
     *
     * @param points the gradient directions, including zero gradients.
     * @param dwPars diffusion weighting parameters: {|G|, delta1, t_delta1, delta2, t_delta2,
     * delta3, t_tdelta4, TE}. All of these parameters except TE will be set to zero for any measurement
     * where the gradient direction is zero.
     *
     * @see imaging.RectGradTRSE_Scheme
     */
    public static String pointSetToRectGradTRSE_Scheme(double[][] points, double[] dwPars) {
        
        double modG = dwPars[0];
        double t_del1 = dwPars[1];
        double del1 = dwPars[2];
        double t_del2 = dwPars[3];
        double del2 = dwPars[4];
        double del3 = dwPars[5];
        double t_del4 = dwPars[6];
        double te = dwPars[7];
        
        int numMeas = points.length;
        
        double[] mod = getPointModulus(points);
        
        boolean pointsNormalized = true;
        
        for (int i = 0; i < numMeas; i++) {
            if (mod[i] > 0.0 && !(Math.abs(1.0 - mod[i]) < 1E-5)) {
                pointsNormalized = false;
            }
        }
        
        if (!pointsNormalized) {
            logger.warning("Some gradient directions do not have unit length. Directions normalized in output.");
        }
        
        points = normalizePoints(points);
        
        StringBuffer sb = new StringBuffer();
        
        sb.append("#Twice-refocused scheme with rectangular gradient pulses\n");
        sb.append("#g_x\tg_y\tg_z\t|g|\tt_del1\tdel1\tt_del2\tdel2\tdel3\tt_del4\tTE\n");
        sb.append("VERSION: " + RectGradTRSE_Scheme.VERSION);
        sb.append("\n");
        
        for (int i = 0; i < numMeas; i++) {
            
            sb.append(df.format(points[i][0]));
            sb.append(df.format(points[i][1]));
            sb.append(df.format(points[i][2]));
            
            if (points[i][0] == 0.0 && points[i][1] == 0.0 && points[i][2] == 0.0) {
                sb.append(df.format(0.0));
                sb.append(df.format(0.0));
                sb.append(df.format(0.0));
                sb.append(df.format(0.0));
                sb.append(df.format(0.0));
                sb.append(df.format(0.0));
                sb.append(df.format(0.0));
                sb.append(df.format(te));
            }
            else {
                sb.append(df.format(modG));
                sb.append(df.format(t_del1));
                sb.append(df.format(del1));
                sb.append(df.format(t_del2));
                sb.append(df.format(del2));
                sb.append(df.format(del3));
                sb.append(df.format(t_del4));
                sb.append(df.format(te));
            }
            
            sb.append("\n");
        }
        
        return sb.toString();
        
    }
    
    
    /**
     * Reads points in the format <BR> <code>
     * numPoints <BR>
     * x y z <BR>
     * x y z <BR>
     * ... <BR> </code>
     * or in the format <BR> <code>
     * numPoints <BR>
     * x <BR>
     * y <BR>
     * z <BR>
     * x <BR>
     * ... <BR> </code>
     * from a text file. If the points are not unit vectors, they are normalized.
     *
     * @param filename the name of the file to read from. If <code>null</code>, the method reads
     * from the standard input.
     */
    public static double[][] readPoints(String filename) throws IOException {
        return readPoints(filename, false, false, false);
    }
    
    
    /**
     * Reads points in the format <BR> <code>
     * numPoints <BR>
     * x y z <BR>
     * x y z <BR>
     * ... <BR> </code>
     * or in the format <BR> <code>
     * numPoints <BR>
     * x <BR>
     * y <BR>
     * z <BR>
     * x <BR>
     * ... <BR> </code>
     * or
     * <code>
     * x y z <BR>
     * x y z <BR>
     * ... <BR> </code>
     * from a text file. Does not normalize points to unity, call <code>normalizePoints</code> to do that.
     *
     * @param filename the name of the file to read from. If <code>null</code>, the method reads
     * from the standard input.
     */
    public static double[][] readPoints(String filename, boolean flipX,
    boolean flipY, boolean flipZ) throws IOException {
        
        
        // handle input a little differently because it's text
        Vector<String> lines = new Vector<String>();
        
        Scanner fileScanner = null;
        
        if (filename == null) {
            fileScanner = new Scanner(System.in);
        }
        else {
            fileScanner = new Scanner(new File(filename));
        }
        
        // Read in the file line by line.
        fileScanner.useDelimiter("\r\n|\n");
        
        // Store all the lines in a vector.
        while(fileScanner.hasNext()) {
            
            // ignore empty lines - avoids problems with trailing newlines
            
            String next = fileScanner.next();
            
            if (next.length() > 0) {
                lines.add(next);
            }
            
        }
        
        fileScanner.close();
        
        
        // now figure out format of input file
        // if there is a header line, the first entry must be the number of points
        boolean firstLineIsHdr = false;
        
        boolean oneNumberPerLine = false;
        
        int numberOfPoints = 0;
        
        try {
            String[] firstLineElements = lines.elementAt(0).trim().split("\\s+");
            
            numberOfPoints = Integer.parseInt(firstLineElements[0]);
            
            // does this number correspond to the number of points?
            if ( numberOfPoints == ((lines.size() - 1) / 3) || numberOfPoints == (lines.size() - 1) ) {
                firstLineIsHdr = true;
            }
            else {
                // if the first line does not start with the number of points, it must be a point
                
                if (firstLineElements.length != 3) {
                    throw new LoggedException("Unknown point set format");
                }
                
                numberOfPoints = lines.size();
            }
            
        }
        catch (NumberFormatException e) {
            // if number of points is not first, then there must be three numbers per line
            // and number of lines == number of points
            numberOfPoints = lines.size();
        }
        
        
        if (firstLineIsHdr) {
            try {
                Double.parseDouble(lines.elementAt(1).trim());
                oneNumberPerLine = true;
            }
            catch (NumberFormatException e) {
                oneNumberPerLine = false;
            }
        }
        
        double[][] points = new double[numberOfPoints][3];
        
        if (oneNumberPerLine) {
            for (int i = 0; i < numberOfPoints; i++) {
                for (int j = 0; j < 3; j++) {
                    double component = Double.parseDouble(lines.elementAt(i * 3 + j + 1));
                    points[i][j] = component;
                }
            }
        }
        else {
            
            int offset = firstLineIsHdr ? 1 : 0;
            
            for (int i = 0; i < numberOfPoints; i++) {
                
                String[] comps = lines.elementAt(i+offset).trim().split("\\s+");
                
                if (comps.length > 3) {
                    throw new LoggedException("Unknown point set format");
                }
                
                for (int c = 0; c < 3; c++) {
                    double comp = Double.parseDouble(comps[c]);
                    points[i][c] = comp;
                }
            }
        }
        
        
        for (int i = 0; i < numberOfPoints; i++) {
            if (flipX && points[i][0] != 0.0) {
                points[i][0] = -points[i][0];
            }
            if (flipY && points[i][1] != 0.0) {
                points[i][1] = -points[i][1];
            }
            if (flipZ && points[i][2] != 0.0) {
                points[i][2] = -points[i][2];
            }
        }
        
        return points;
        
    }
    
    
    
    /**
     * Gets the modulus of the vector from (0, 0, 0) to the point.
     *
     */
    public static double[] getPointModulus(double[][] points) {
        
        int numberOfPoints = points.length;
        
        double[] mod = new double[numberOfPoints];
        
        // normalize points
        for (int i = 0; i < numberOfPoints; i++) {
            mod[i] = Math.sqrt(points[i][0] * points[i][0] + points[i][1] * points[i][1] +
            points[i][2] * points[i][2]);
        }
        
        return mod;
    }
    
    
    /**
     * Normalizes points to unit length.
     *
     */
    public static double[][] normalizePoints(double[][] points) {
        
        int numberOfPoints = points.length;
        
        double[][] normPoints = new double[numberOfPoints][3];
        
        // normalize points
        for (int i = 0; i < numberOfPoints; i++) {
            double mod = Math.sqrt(points[i][0] * points[i][0] + points[i][1] * points[i][1] +
            points[i][2] * points[i][2]);
            
            // don't worry about small rounding errors
            if (Math.abs(1.0 - mod) > 1E-4 && mod > 0.0) {
                points[i][0] = points[i][0] / mod;
                points[i][1] = points[i][1] / mod;
                points[i][2] = points[i][2] / mod;
            }
            
            // get rid of -0.0
            for (int j = 0; j < 3; j++) {
                if (points[i][j] == 0.0) {
                    points[i][j] = Math.abs(points[i][j]);
                }
            }
            
        }
        
        return points;
    }
    
    
}
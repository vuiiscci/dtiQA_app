package apps;

import java.io.*;
import java.text.*;
import java.util.logging.Logger;

import misc.LoggedException;

import tools.CL_Initializer;

import data.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Converts binary data to text.
 * 
 * <dt>Description:
 * 
 * <dd>Reads in binary data with format specified on the command line, although
 * assumed big-endian, and writes out the values as numbers in text line by
 * line.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 *  
 */
public class BinaryToText extends Executable{

	public BinaryToText(String[] args){
		super(args);
	}
	
    /**
     * the logger object for message handling
     */
    private static Logger logger = Logger.getLogger("camino.apps.BinaryToText");

    /**
     * The type of the input data.
     */
    //private static String inputDataType = "float";
	private String inputDataType;

    /**
     * Input data file. If left null, the program reads data from the standard
     * input.
     */
    //private static String inputFile = null;
	private String inputFile;

    /**
     * The stream from which to read the input.
     */
    //private static DataInputStream in = new DataInputStream(System.in);
	private static DataInputStream in;

    // This is the number of decimal places to print in scientific
    // notation.
    private int decimalPlaces;

    private String zeros;

    /** If true, annotate output with line numbers / voxel indexing */
    private boolean annotate;

    /** Components per voxel - if annotate is true and this is > 0, label output by voxel */
    private int components;
	
	public void initDefaultVals(){
		inputDataType = "float";
		inputFile = null;
		in = new DataInputStream(System.in);
		decimalPlaces = 6;
		zeros = "";
                annotate = false;
		components = 0;
	}	
	
	public void initOptions(String[] args){
//    public static void main(String[] args) {
        // parse arguments
        CL_Initializer.CL_init(args);


        // Parse the command line arguments
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-inputdatatype")) {
                inputDataType = args[i + 1];
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
            if (args[i].equals("-dp")) {
                decimalPlaces = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i);
                CL_Initializer.markAsParsed(i + 1);
            }
	    if (args[i].equals("-linenum")) {
		annotate = true;
	       CL_Initializer.markAsParsed(i);
	    }
	    if (args[i].equals("-components")) {
		annotate = true;
		components = Integer.parseInt(args[i + 1]);
                CL_Initializer.markAsParsed(i,2);
	    }
        }

        // check all arguments parsed
        CL_Initializer.checkParsing(args);
	}
	
	public void initVariables() {
	
	}
	
	public void execute(OutputManager om){
        // Set up the scientific notation format.
        for (int i = 0; i < decimalPlaces; i++) {
            zeros += "0";
        }

        DecimalFormat sciForm = new DecimalFormat("0." + zeros + "E00;-0." + zeros
                + "E00");

	DecimalFormat lineNumForm = new DecimalFormat("0000");

	int lineNum = 0;

        Object[] o = new Object[1];
		while (true) {
	    
            readNext(inputDataType, o);

	    if (annotate) {
		// if no components given, just output line numbers
		if (components == 0) {
		    System.out.print(lineNumForm.format(lineNum + 1) + " : ");
		}
		else {
		    if (lineNum % components == 0) {
			System.out.println("--- Voxel : " + lineNumForm.format(lineNum / components + 1) + " ---");
		    }
		    System.out.print(lineNumForm.format(lineNum % components + 1) + " : ");
		}
	    }

            if (inputDataType.equals("char") || inputDataType.equals("byte")
                    || inputDataType.equals("short") || inputDataType.equals("int")
                    || inputDataType.equals("long")) {

                System.out.println(MessageFormat.format("{0,number, 0;-0}", o));
            }
            else if (inputDataType.equals("float") || inputDataType.equals("double")) {

                System.out.println(sciForm.format(((Number) o[0]).doubleValue()));
            }

	    lineNum++;

        }
    }

    /**
     * Reads the next value and places it in the first element of array <code>o</code>.
     *
     * @param inputDataType as defined in {@link data.ExternalDataSource}.
     * 
     */
    public static void readNext(String inputDataType, Object[] o) {

        try {

            if (inputDataType.equals("byte")) {
                o[0] = new Byte(in.readByte());
            }
            else if (inputDataType.equals("char")) {
                int b = in.readByte();
                o[0] = new Integer(b >= 0 ? b : 256 + b);
            }
            else if (inputDataType.equals("short")) {
                o[0] = new Short(in.readShort());
            }
            else if (inputDataType.equals("int")) {
                o[0] = new Integer(in.readInt());
            }
            else if (inputDataType.equals("long")) {
                o[0] = new Long(in.readLong());
            }
            else if (inputDataType.equals("float")) {
                o[0] = new Float(in.readFloat());
            }
            else if (inputDataType.equals("double")) {
                o[0] = new Double(in.readDouble());
            }

        }
        catch (Exception e) {

            // End of stream reached. Tidy up and exit.
            try {
                in.close();
            }
            catch (Exception e2) {
                LoggedException.logExceptionWarning(e2, Thread.currentThread().getName());
                
            }

            System.exit(0);
        }
    }

}

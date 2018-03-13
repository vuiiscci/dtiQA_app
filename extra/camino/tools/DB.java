package tools;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Contains printIfVerbose used for debugging.
 * 
 * <dt>Description:
 * 
 * <dd>Used to print output from a program at specified levels of verbosity.
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id$
 *  
 */
public class DB {

    /**
     * Specifies the level of verbosity.
     */
    private static int VERBOSITY = 0;

    /**
     * Output file stream.
     */
    private static FileOutput fOut;

    /**
     * Specifies the level of verbosity for an instance of DB.
     */
    private int VERBOSITY_NS = 0;

    /**
     * Output file stream for an instance of DB.
     */
    private FileOutput fOut_ns;

    /**
     * Default constructor.
     */
    public DB() {
    }

    /**
     * Constructs a DB object given a verbosity level.
     */
    public DB(int v) {
        VERBOSITY_NS = v;
    }

    /**
     * Constructs a DB object given a verbosity level and filename for output.
     */
    public DB(int v, String fName) {
        VERBOSITY_NS = v;
        fOut_ns = new FileOutput(fName);
    }

    /**
     * Prints the object if vLevel is less than or equal to the specified
     * verbosity. If an output file has been specified the output goes there,
     * otherwise it goes to the standard output.
     */
    public static void printIfVerbose(Object o, int vLevel) {
        if (VERBOSITY >= vLevel) {
            if (fOut != null)
                try {
                    fOut.writeString(o.toString());
                    fOut.writeString("\n");
                }
                catch (Exception e) {
                    System.err.println("WARNING: " + e);
                }
            else {
                System.err.println(o);
            }
        }
    }

    /**
     * Prints the object if vLevel is less than or equal to the specified
     * verbosity. If an output file has been specified the output goes there,
     * otherwise it goes to the standard output. This method uses the output
     * stream and verbosity level associated with a particular instance of the
     * class.
     */
    public void printIfVerboseInst(Object o, int vLevel) {
        if (VERBOSITY_NS >= vLevel) {
            if (fOut_ns != null)
                try {
                    fOut_ns.writeString(o.toString());
                    fOut_ns.writeString("\n");
                }
                catch (Exception e) {
                    System.err.println("WARNING: " + e);
                }
            else {
                System.err.println(o);
            }
        }
    }

    /**
     * Sets the verbosity level.
     */
    public static void setVerbosity(int v) {
        VERBOSITY = v;
    }

    /**
     * Sets the verbosity level for an instance of the class.
     */
    public void setVerbosityInst(int v) {
        VERBOSITY_NS = v;
    }

    /**
     * Sets the output filename.
     */
    public static void setOutputFile(String fName) {
        close();
        fOut = new FileOutput(fName);
    }

    /**
     * Sets the output filename for an instance of the class.
     */
    public void setOutputFileInst(String fName) {
        closeInst();
        fOut_ns = new FileOutput(fName);
    }

    /**
     * Returns the current verbosity level.
     */
    public static int getVerbosity() {
        return VERBOSITY;
    }

    /**
     * Returns the current verbosity level for an instance of the class.
     */
    public int getVerbosityInst() {
        return VERBOSITY_NS;
    }

    /**
     * Closes the output file.
     */
    public static void close() {
        if (fOut != null) {
            fOut.close();
        }
    }

    /**
     * Closes the output file for an instance of the class.
     */
    public void closeInst() {
        if (fOut_ns != null) {
            fOut_ns.close();
        }
    }

}
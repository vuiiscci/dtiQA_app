package tools;

import java.io.*;

import misc.LoggedException;

/**
 * A simple input class to read values from a file of characters. If
 * any file errors occur, methods in this class will display an error
 * message and terminate the program.  <p> Be aware that this class
 * will return a value when it encounters EOF. For example, if you
 * call readDouble() and there is no more data, the method returns
 * 0.0. Further calls will result in a LoggedException being
 * thrown. This is a potential source of bugs because the file could
 * be one line shorter than expected and you'd never know.
 * 
 * @version $Id$
 * @author Graham Roberts
 */
public class FileInput {

    /**
     * Construct FileInput object given a file name.
     */
    public FileInput(String fname) {
        filename = fname;
        try {
            reader = new BufferedReader(new FileReader(filename));
        }
        catch (FileNotFoundException e) {
            throw new LoggedException("WARNING: Can't open file: " + filename);
        }
    }

    /**
     * Close the file when finished
     */
    public void close() {
        try {
            reader.close();
        }
        catch (IOException e) {
            throw new LoggedException("Can't close file: " + filename);
        }
    }

    /**
     * Return true if the end of file has been reached.
     */
    public boolean eof() {
        return eof;
    }

    /**
     * Read an int value from file.
     */
    public final synchronized int readInteger() {
        //
        //  The place to hold the data received from the file so
        //  that we can parse it to find the int value.
        //

        if (eof) {
            throw new LoggedException("Tried to read past end of file");
        }


        String input = "";
        int val = 0;

        //
        //  Read a String of characters from the file.
        //
        try {
            input = reader.readLine();
            if (input == null) {
                eof = true;
            }
        }
        catch (IOException e) {
            throw new LoggedException("readInteger failed for file: " + filename);
        }
        //
        //  Parse the String to construct an int value.
        //
        try {
            val = Integer.parseInt(input);
        }
        catch (NumberFormatException e) {
            throw new LoggedException(e);
        }
        return val;
    }

    /**
     * Read a long value from file.
     */
    public final synchronized long readLong() {
        //
        //  The place to hold the data received from the file so
        //  that we can parse it to find the long value.
        //

        if (eof) {
            throw new LoggedException("Tried to read past end of file");
        }

        String input = "";
        long val = 0L;
        //
        //  Read a String of characters from the file.
        //
        try {
            input = reader.readLine();
            if (input == null) {
                eof = true;
            }
        }
        catch (IOException e) {
            throw new LoggedException("readLong failed for file: " + filename);
        }
        //
        //  Parse the String to construct a long value.
        //
        try {
            val = Long.parseLong(input);
        }
        catch (NumberFormatException e) {
            throw new LoggedException(e);
        }
        return val;
    }

    /**
     * Read a double value from file.
     */
    public final synchronized double readDouble() {
        //
        //  The place to hold the data received from the file so
        //  that we can parse it to find the double value.
        //

        if (eof) {
            throw new LoggedException("Tried to read past end of file");
        }

        String input = "";
        double val = 0.0D;
        //
        //  Read a String of characters from the file.
        //
        try {
            input = reader.readLine();
            if (input == null) {
                eof = true;
            }
        }
        catch (IOException e) {
            throw new LoggedException("readDouble failed for file: " + filename);
        }
        //
        //  Parse the String to construct a double value.
        //
        try {
	    val = Double.parseDouble(input);
	}
        catch (NumberFormatException e) {
            throw new LoggedException(e);
        }
        return val;
    }

    /**
     * Read a float value from file.
     */
    public final synchronized float readFloat() {
        //
        //  The place to hold the data received from the file so
        //  that we can parse it to find the float value.
        //

        if (eof) {
            throw new LoggedException("Tried to read past end of file");
        }

        String input = "";
        float val = 0.0F;
        //
        //  Read a String of characters from the file.
        //
        try {
            input = reader.readLine();
            if (input == null) {
                eof = true;
            }
        }
        catch (IOException e) {
            throw new LoggedException("readFloat failed for file: " + filename);
        }
        //
        //  Parse the String to construct a float value.
        //
        try {
            val = Float.parseFloat(input); 
	}
        catch (NumberFormatException e) {
            throw new LoggedException(e);
        }
        return val;
    }

    /**
     * Read a char value from file.
     */
    public final synchronized char readCharacter() {
        //
        //  No need to parse anything, just get a character and return
        //  it..
        //
        if (eof) {
            throw new LoggedException("Tried to read past end of file");
        }

        char c = ' ';
        try {
            int n = reader.read();
            if (n == -1) {
                eof = true;
            }
            else {
                c = (char) n;
            }
        }
        catch (IOException e) {
            throw new LoggedException("readCharacter failed for file: " + filename);
        }
        return c;
    }

    /**
     * Read an String value from file.
     */
    public final synchronized String readString() {
        //
        //  No need to parse anything, just get a string and return
        //  it..
        //

        if (eof) {
            throw new LoggedException("Tried to read past end of file");
        }

        String s = "";
        try {
            s = reader.readLine();
            if (s == null) {
                eof = true;

            }
        }
        catch (IOException e) {
            throw new LoggedException("readString failed for file: " + filename);
        }
        return s;
    }


    /**
     * Instance variables to store file details.
     */
    private String filename = "";

    private BufferedReader reader = null;

    private boolean eof = false;
}

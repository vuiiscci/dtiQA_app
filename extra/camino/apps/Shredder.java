package apps;

import java.io.*;
import misc.LoggedException;
import java.util.logging.Logger;

import data.*;

/**
 *
 * Extracts periodic chunks from binary data. Call with command line:
 * <p> java apps.Shredder &ltoffset &rt &lt chunkSize &rt &lt space
 * &rt <p> Shredder makes an initial offset of <code>offset</code>
 * bytes. It then reads and outputs <code>chunkSize</code> bytes,
 * skips <code>space</code> bytes, and repeats until there is no more
 * input.  <p> If the <code>chunkSize</code> is negative, chunks of
 * size |<code>chunkSize</code>| are read and the byte ordering of
 * each chunk is reversed. If byte swapping is required, the
 * <code>chunkSize</code> should correspond to the size of the data
 * type.
 *
 * @author Danny Alexander
 * @author Philip Cook
 * @version $Id$
 */
public class Shredder extends Executable{
  
	public Shredder(String[] args){
    	super(args);
    }
	
    private static Logger logger = Logger.getLogger("camino.apps.Shredder");
    
	public static int FILEBUFFERSIZE;
	public static DataInputStream in;
	public static DataOutputStream out;
    private static int skipChunkSize;
    private static byte[] skipChunk;
	private int offset;
    private int chunkSize;
	private int space;
	private byte[] chunk;
	private boolean reverse;


    public void initDefaultVals(){
	
		FILEBUFFERSIZE = data.ExternalDataSource.FILEBUFFERSIZE;
    
		in = new DataInputStream(new BufferedInputStream(System.in, FILEBUFFERSIZE));
    
		out = new DataOutputStream(new BufferedOutputStream(System.out, FILEBUFFERSIZE));
    
    //  chunk size for skip method
		skipChunkSize = FILEBUFFERSIZE;
    
		skipChunk = null;
    
	}	
	
//    public static void main(String[] args) {
	public void initOptions(String[] args){
        
        if (args.length != 3) {
            System.err
            .println("Usage: java Shredder <offset> <chunk size> <space>.\n\n");
            System.exit(0);
        }
        
        offset = Integer.parseInt(args[0]);
        chunkSize = Integer.parseInt(args[1]);
        space = Integer.parseInt(args[2]);
    }
     

	public void initVariables(){
		skipChunk = new byte[skipChunkSize];
                       
        reverse = false;
        if (chunkSize < 0) {
            reverse = true;
            chunkSize = -chunkSize;
        }
        
        chunk = new byte[chunkSize]; // stuff we want to output        
	}
	
	public void execute(OutputManager om){	
    
        if (chunkSize == 0) {
            throw new LoggedException("Can't read zero byte chunks");
        }
 		if ( !skip(offset) ) {
            throw new LoggedException("Could not skip initial offset of " + offset + " bytes");
        }
	
        while (getChunk(chunk)) { // get chunk bytes
            
            if (reverse) {
                arrFlip(chunk);
            }
            
            output(chunk);
            
            if (!skip(space)) { // try to skip
                break;
            }
        }
        
        // Close output stream.
        try {
            out.close();
        }
        catch (Exception e) {
            throw new LoggedException(e);
        }
        
        
    }
    
    
    /**
     * Reads in the next chunk from the standard input into the array
     * chunk and returns a boolean indicating the success of the
     * operation.
     *
     * @param chunk An array to contain the chunks
     *
     * @return boolean indicating whether or not the chunk was read
     * successfully. If part of a chunk was read successfully, the return value is <code>true</code>.
     * Note that if only part of a chunk is read, the remaining bytes will contain the original
     * contents of <code>chunk</code>
     *
     */
    public static boolean getChunk(byte[] chunk) {
        
        int totalBytesRead = 0;
        
        try {
            readBytes:
                while (totalBytesRead < chunk.length) {
                int bytesRead = in.read(chunk, totalBytesRead, chunk.length - totalBytesRead);
                
                if (bytesRead == -1) {
                    // eof
                    break readBytes;
                }
                
                totalBytesRead += bytesRead;
                
                }
        }
        catch (IOException e) {
            // note: an exception is not thrown here for EOF
            // in case of EOF, bytesRead is -1
            throw new LoggedException(e);
        }
        
        // return true if we read at least one byte
        if (totalBytesRead > 0) {
            return true;
        }
        
        return false;
    }
    
    
    /**
     * Skip bytes from the input stream.
     *
     * @param bytesToSkip the number of bytes to skip.
     * @return <code>true</code> if the bytes were skipped successfully, <code>false</code> if EOF was found
     * before enough bytes could be skipped.
     *
     */
    public static boolean skip(int bytesToSkip) {
        
        if (bytesToSkip == 0) {
            return true;
        }
        if (skipChunkSize == 0) {
            throw new LoggedException("Internal error: skip chunk size is zero!");
        }
        
        // the skip method in InputStream simply reads and discards one
        // byte at a time, so we won't use that
        int chunksToSkip = bytesToSkip / skipChunkSize;
        
        int remainderSkip = bytesToSkip % skipChunkSize;
        
        // total bytes we read from the stream in this method
        int totalBytesRead = 0;
        
        try {
            
            
            for (int i = 0; i < chunksToSkip; i++) {
                
                int bytesRead = 0;
                
                while (bytesRead < skipChunkSize) {
                    int bytesThisRead = in.read(skipChunk, bytesRead, skipChunkSize - bytesRead);
                    
                    if (bytesThisRead == -1) {
                        // EOF
                        return false;
                    }
                    
                    bytesRead += bytesThisRead;
                }
                
                totalBytesRead += bytesRead;
            }
            
            
            int bytesRead = 0;
            
            while (bytesRead < remainderSkip) {
                int bytesThisRead = in.read(skipChunk, bytesRead, remainderSkip - bytesRead);
                
                if (bytesThisRead == -1) {
                    // EOF
                    return false;
                }
                
                bytesRead += bytesThisRead;
            }
            
            totalBytesRead += bytesRead;
            
            if (totalBytesRead == bytesToSkip) {
                // success
                return true;
            }
            else {
                throw new
                LoggedException("Needed to skip " + bytesToSkip + " bytes, but skipped " + totalBytesRead);
            }
            
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }
    }
    
    
    /**
     * Outputs the chunk to the standard output.
     *
     * @param chunk The chunk to output.
     */
    public static void output(byte[] chunk) {
        try {
            out.write(chunk);
        }
        catch (IOException e) {
            throw new LoggedException(e);
        }
    }
    
    
    /**
     * Reverses the order of the bytes in the array.
     *
     * @param The chunk to reverse.
     */
    public static void arrFlip(byte[] chunk) {
        for (int i = 0; i < chunk.length / 2; i++) {
            byte temp = chunk[i];
            chunk[i] = chunk[chunk.length - i - 1];
            chunk[chunk.length - i - 1] = temp;
        }
    }
    
}


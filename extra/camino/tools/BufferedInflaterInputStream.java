package tools;

import misc.*;

import java.io.*;
import java.util.zip.*;


/**
 * This class is useful to reduce the decompression overhead that comes from reading
 * from a <code>InflaterInputStream</code> a few bytes at a time. This class buffers on the
 * decompressed side. The input stream passed to this class upon construction is responsible
 * for buffering on the disk access side. 
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class BufferedInflaterInputStream extends FilterInputStream {

    private int bufferSize;

    // When we reach EOF, this buffer shrinks to the number of remaining bytes
    private byte[] buffer;

    private boolean streamEOF = false;

    private int bufferIndex = 0;


    /**
     * Create a stream from a compressed input stream. 
     * 
     * @param infIn the underlying decompression stream that is to be buffered.
     * @param bufferSize the amount of decompressed data to buffer.
     */
    public BufferedInflaterInputStream(InflaterInputStream infIn, int bufferSize) {
	super(infIn);

	this.bufferSize = bufferSize;

	buffer = new byte[bufferSize];

	fillBuffer();
    }

  

    public int read() throws IOException {

	int retVal = -1;

	if (bufferIndex < bufferSize) {

	    retVal = buffer[bufferIndex] & 0x00ff;
	    
	    bufferIndex++;
		
	    if (bufferIndex == bufferSize && !streamEOF) {
		fillBuffer();
	    } 
	}
	
	return retVal;
	
    }

    
    
    public int read(byte[] b) throws IOException {
	
	return read(b, 0, b.length);
    }


    public int read(byte[] b, int start, int length) throws IOException {
	
	if (bufferIndex == bufferSize) {
	    // only happens at EOF
	    return -1;
	}

	int bytesRead = 0;	    

	if (length <= bufferSize - bufferIndex) {
	    System.arraycopy(buffer, bufferIndex, b, start, length);
	    bufferIndex += length;

	    if (bufferIndex == bufferSize && !streamEOF) {
		fillBuffer();
	    }

	    bytesRead = length;
	}
	else {
	    
	    System.arraycopy(buffer, bufferIndex, b, start, bufferSize - bufferIndex);

	    bytesRead += bufferSize - bufferIndex;

	    bufferIndex += bytesRead; // ie, bufferIndex now equals bufferSize

	    if (!streamEOF) {
		fillBuffer();
	    }
	    else {
		return bytesRead;
	    }
	    
	    while (bytesRead < length && bufferIndex < bufferSize) {
		
		int bytesToRead = 
		    (bytesRead + bufferSize) < length ? bufferSize : (length - bytesRead);

		System.arraycopy(buffer, bufferIndex, b, start + bytesRead, bytesToRead);

		bytesRead += bytesToRead;

		bufferIndex += bytesToRead;

		if (bufferIndex == bufferSize && !streamEOF) {
		    fillBuffer();
		}
		else {
		    return bytesRead;
		}

	    }
	
	}
	
	return bytesRead;    	
	
    }


    /**
     * Inefficiently skips bytes by reading them, but faster than the superclass method, 
     * which seems to take forever.
     */
    public long skip(long bytesToSkip) throws IOException {

        long remainder = bytesToSkip;

        while (remainder > 0) {

            int toSkip = 0;

            if (remainder > bufferSize) {
                toSkip = bufferSize;
            }
            else {
                // safe because we only get here if skip <= bufferSize, an int
                toSkip = (int)remainder;
            }

            int skipped = read(new byte[toSkip], 0, toSkip);

            if (skipped < 0) {
                throw new IOException("Tried to skip past EOF");
            }

            remainder -= skipped;
        }

        return bytesToSkip;
        
    }


    private void fillBuffer() {

	if (streamEOF) {
	    throw new LoggedException("Tried to read past EOF");
	}
	
	int totalBytesRead = 0;

	try {
	    readBytes:
	    while (totalBytesRead < bufferSize) {
		int bytesRead = in.read(buffer, totalBytesRead, bufferSize - totalBytesRead);

		if (bytesRead == -1) {
		    // eof
		    break readBytes;
		}

		totalBytesRead += bytesRead;

	    }

	    // if at EOF, shrink buffer. If we got no bytes (already reached EOF)
	    // then we do not shrink the array or reset the index
	    if (totalBytesRead < bufferSize) {

		if (totalBytesRead > 0) {

		    byte[] tmp = new byte[totalBytesRead];
		    
		    System.arraycopy(buffer, 0, tmp, 0, totalBytesRead);
		    
		    buffer = tmp;
		    
		    bufferIndex = 0;

		    bufferSize = totalBytesRead;
		    
		}

		// if we got 0 bytes, leave the buffer and index alone (nothing more to read)
		// and set streamEOF
		
		streamEOF = true;
		
	    }
	    else {
		// filled the buffer
		bufferIndex = 0;
	    }


	}
	catch (IOException e) {
	    // note: an exception is not thrown here for EOF
	    // in case of EOF, bytesRead is -1
	    throw new LoggedException(e);
	}

	
    }


    

}

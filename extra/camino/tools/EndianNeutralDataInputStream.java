package tools;

import java.io.*;

public class EndianNeutralDataInputStream implements DataInput {

    
    private final DataInput reader;

    private final InputStream in;

    public EndianNeutralDataInputStream(InputStream in, boolean intelByteOrder) {

	this.in = in;

	if (intelByteOrder) {
	    reader = new tools.LEFilterInputStream(in);
	}
	else {
	    reader = new DataInputStream(in);
	}

    }


    public void close() throws IOException {
	in.close();
    }


    /**
   * Reads a <code>boolean</code> from the underlying input stream by 
   * reading a single byte. If the byte is zero, false is returned.
   * If the byte is positive, true is returned. 
   *
   * @return      b   the <code>boolean</code> value read.
   * @exception  EOFException  if the end of the underlying input stream
   *              has been reached
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public boolean readBoolean() throws IOException {
      return reader.readBoolean();
  }

  /**
   * Reads a signed <code>byte</code> from the underlying input stream
   * with value between -128 and 127
   *
   * @return     the <code>byte</code> value read.
   * @exception  EOFException  if the end of the underlying input stream
   *              has been reached
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public byte readByte() throws IOException {
      return reader.readByte();
  }

  /**
   * Reads an unsigned <code>byte</code> from the underlying 
   * input stream with value between 0 and 255
   *
   * @return     the <code>byte</code> value read.
   * @exception  EOFException  if the end of the underlying input 
   *              stream has been reached
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public int readUnsignedByte() throws IOException {
      return reader.readUnsignedByte();
  }

  /**
   * Reads a two byte signed <code>short</code> from the underlying 
   * input stream. 
   *
   * @return     the <code>short</code> read.
   * @exception  EOFException  if the end of the underlying input stream 
   *              has been reached
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public short readShort() throws IOException {
      return reader.readShort();
  }

  /**
   * Reads a two byte unsigned <code>short</code> from the underlying 
   * input stream in little endian order, low byte first. 
   *
   * @return     the int value of the unsigned short read.
   * @exception  EOFException  if the end of the underlying input stream 
   *              has been reached
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public int readUnsignedShort() throws IOException {
      return reader.readUnsignedShort();
  }

  /**
   * Reads a two byte Unicode <code>char</code> from the underlying 
   * input stream.
   *
   * @return     the int value of the unsigned short read.
   * @exception  EOFException  if the end of the underlying input stream 
   *              has been reached
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public char readChar() throws IOException {
      return reader.readChar();
  }


  /**
   * Reads a four byte signed <code>int</code> from the underlying 
   * input stream. 
   *
   * @return     the <code>int</code> read.
   * @exception  EOFException  if the end of the underlying input stream 
   *              has been reached
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public int readInt() throws IOException {
      return reader.readInt();
  }


  /**
   * Reads a four byte unsigned <code>int</code> from the underlying 
   * input stream. 
   *
   * @return     the <code>int</code> read, as a <code>long</code>.
   * @exception  EOFException  if the end of the underlying input stream 
   *              has been reached
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public long readUnsignedInt() throws IOException {
      int unsigned = readInt();

      return (unsigned & 0xffffffffL);	  
  }



  /**
   * Reads an eight byte signed <code>long</code> from the underlying 
   * input stream.
   *
   * @return     the <code>long</code> read.
   * @exception  EOFException  if the end of the underlying input stream 
   *              has been reached
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public long readLong() throws IOException {
      return reader.readLong();    
  }

  /**
   * Reads a string of no more than 65,535 characters 
   * from the underlying input stream using UTF-8 
   * encoding. This method first reads a two byte short 
   * in <b>big</b> endian order as required by the 
   * UTF-8 specification. This gives the number of bytes in 
   * the UTF-8 encoded version of the string.
   * Next this many bytes are read and decoded as UTF-8
   * encoded characters. 
   *
   * @return     the decoded string
   * @exception  UTFDataFormatException if the string cannot be decoded
   * @exception  IOException  if the underlying stream throws an IOException.
   */
  public String readUTF() throws IOException {
      return reader.readUTF();     
  }

  /**
   *
   * Reads an eight byte <code>double</code> from the underlying 
   * input stream.
   *
   * @return     the <code>double</code> read.
   * @exception  EOFException if end of stream occurs before eight bytes 
   *             have been read.
   * @exception  IOException   if an I/O error occurs.
   */
  public final double readDouble() throws IOException {
      return reader.readDouble();     
  }
  
  /**
   *  
   * Reads a four byte <code>float</code> from the underlying 
   * input stream.
   * @return     the <code>float</code> read.
   * @exception  EOFException if end of stream occurs before four bytes 
   *             have been read.
   * @exception  IOException  if an I/O error occurs.
   */
  public final float readFloat() throws IOException {
      return reader.readFloat();
  }
  
  /**
   * Skip exactly <code>n</code> bytes of input in the underlying 
   * input stream. This method blocks until all the bytes are skipped, 
   * the end of the stream is detected, or an exception is thrown. 
   *
   * @param      n   the number of bytes to skip.
   * @return     the number of bytes skipped, generally n
   * @exception  EOFException  if this input stream reaches the end before
   *               skipping all the bytes.
   * @exception  IOException  if the underlying stream throws an IOException.
   */
   public final int skipBytes(int n) throws IOException {
       return reader.skipBytes(n);
   }

    /**
     * @inheritDoc
     * @see java.io.DataInput#readFully(byte[])
     */
    public final void readFully ( byte b[] ) throws IOException {
        reader.readFully(b);
    }
    
    
    /**
     * @inheritDoc
     * @see java.io.DataInput#readFully(byte[], int, int)
     */
    public final void readFully ( byte b[], int off, int len ) throws IOException {
        reader.readFully( b, off, len );
    }

    
    /**
     * @return a rough approximation of the 8-bit stream as a 16-bit unicode
     *         string
     * @throws IOException
     * @deprecated This method does not properly convert bytes to characters.
     *             Use a Reader instead with a little-endian encoding.
     */
    @Deprecated public final String readLine () throws IOException {
	throw new UnsupportedOperationException("Can't use this deprecated method");
    }








}

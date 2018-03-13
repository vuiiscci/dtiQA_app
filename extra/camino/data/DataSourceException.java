package data;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Exception class used by the DataSource interface.
 * 
 * <dt>Description:
 * 
 * <dd>Standard extension of the <code>LoggedException</code> class. 
 * Instances of this class will be logged on creation.
 * 
 * 
 * @version $Id$
 * @author Danny Alexander
 *  
 */
public class DataSourceException extends misc.LoggedException {

    public DataSourceException() {
        super();
    }

    public DataSourceException(String s) {
        super(s);
    }
}


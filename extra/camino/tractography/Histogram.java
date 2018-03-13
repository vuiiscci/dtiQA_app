package tractography;

import numerics.*;
import misc.LoggedException;
import java.util.Hashtable;
import java.util.Set;

public class Histogram {
    
    /**
     * contents
     */
    private Hashtable <Object, HistBin> histogram;
    private int dimensions;
    private double binSize;
    private boolean useLog;

    /*
     * Creates a new Histogram object that is either 1 or two dimensional.
     * 
     */
    public Histogram(int dims, double binSz, boolean uLog) {

	dimensions = dims;
	binSize = binSz;
	useLog = uLog;

	if(dimensions==1)
	    histogram = new Hashtable<Object, HistBin>();
	else if (dimensions==2)
	    histogram = new Hashtable<Object, HistBin>();
	else
	    throw new LoggedException("Histogram must be either 1 or 2 dimensional.  requested " + dimensions + "D histogram");
    }

    /*
     * Creates a new 1D Histogram object.
     * 
     */
    public Histogram(double binSz, boolean uLog){
	this(1, binSz, uLog);
    }

    /*
     * Adds a fibre-orientation estimate to the histogram.
     * 
     */
    public void add(Vector3D direction, double eig1, double eig2) {
	
	String key = "";

	if(dimensions==1)
	    key = getKey1D(eig1, eig2);
	else if(dimensions ==2)
	    key = getKey2D(eig1, eig2);

	HistBin bin = (HistBin)histogram.get(key);
	//System.err.println(Math.log(((int)(eig1/binSize))*binSize) + " " + Math.log(((int)(eig2/binSize))*binSize));
	if(bin==null){
	    int e1 = (int)(eig1/binSize);
	    int e2 = (int)(eig2/binSize);
	    if(useLog) {
		// binsize will be log if using logs... need to account for this!
		e1 = (int)(Math.log(eig1)/binSize);
		e2 = (int)(Math.log(eig2)/binSize);
		bin = new HistBin(Math.exp(e1*binSize), Math.exp(e2*binSize));
	    }
	    else{
		bin = new HistBin(e1*binSize, e2*binSize);
	    }
	    histogram.put(key, bin);
	}

	bin.addDirToList(direction);
    }

    
    public Set getEntrySet() {
	return histogram.entrySet();
    }


    /*
     * returns a key for a 1D histogram bin.
     * 
     */
    private String getKey1D(double eig1, double eig2) {

	int e1 = (int)(eig1/binSize);
	int e2 = (int)(eig2/binSize);

	int trace = (int)((e1*binSize + e2*binSize)/binSize);

	if(useLog) {
	    e1 = (int)(Math.log(eig1)/binSize);
	    e2 = (int)(Math.log(eig2)/binSize);
	    trace = (int)(Math.log(Math.exp(e2*binSize) + Math.exp(e2*binSize))/binSize);
	}

	return Integer.toString(trace);
    }

    /*
     * returns a key for a 2D histogram bin.
     * 
     */
    private String getKey2D(double eig1, double eig2) {
	
	int e1 = (int)(eig1/binSize);
	int e2 = (int)(eig2/binSize);
	
	if(useLog) {
	    e1 = (int)(Math.log(eig1)/binSize);
	    e2 = (int)(Math.log(eig2)/binSize);
	    //System.err.println(e1 + " " + e2);
	}

	return Integer.toString(e1) + Integer.toString(e2);
    }
}

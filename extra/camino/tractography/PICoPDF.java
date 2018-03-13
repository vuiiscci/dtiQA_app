package tractography;

import misc.*;

import java.util.Random;

/**
 * PDF types: ACG, BINGHAM, WATSON
 */
public enum PICoPDF {
    
    ACG("acg", 3),
    BINGHAM("bingham", 2),
    WATSON("watson", 1);
    

    PICoPDF(String pdfName, int nParams) {
	name = pdfName;
	numParams = nParams;
    }
    
    
    public String toString() {
	return name;
    }
    

    public static PICoPDF getPDF(String s) {

	for (PICoPDF pdf : PICoPDF.values()) {
	    if (s.equals(pdf.name)) {
		return pdf;
	    }
	}
	
	throw new LoggedException("Unsupported PICo PDF type " + s);
    }
    
    public PICoRandomizer getRandomizer(PICoTractographyImage im, Random ran) {

	// Sun advises use of constant-specific methods, ie declare this in braces after each 
	// PDF declaration, but that's less compact
	switch (this) { 
		
	case ACG : return new PICoACGRandomizer(im, ran);
	case BINGHAM : return new PICoBinghamRandomizer(im, ran);
	case WATSON : return new PICoWatsonRandomizer(im, ran);		
	    
	}
	
	// should never happen as long as the switch is correct
	throw new LoggedException("Unknown PICo PDF " + name); 
    }

    /**
     * The name of the PDF.
     */
    public final String name;


    /**
     * Number of parameters for this PDF.
     */
    public final int numParams;
    
}

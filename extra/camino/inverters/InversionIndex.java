package inverters;

import misc.*;

import java.util.Random;

/**
 * Codes for model-based inversions.
 *
 */
public enum InversionIndex {
   
    // call it LDT rather than DT to avoid confusion with DT class
    LDT("ldt", 1, 1),
    LDT_ALIAS("dt", 1, 1),
    ALGDT("algdt", 0, 1),
    NLDT_POS("nldt_pos", 2, 1),
    NLDT("nldt", 4, 1), 
    LDT_WTD("ldt_wtd", 7, 1), 
    ADC("adc", -1, 0),
    RESTORE("restore", -2, 1),
    BALL_STICK("ball_stick", -3, 0),
    CYLCYL("cylcyl", 10, 2),
    CYLCYL_EQ("cylcyl_eq", 20, 2),
    POSPOS("pospos", 30, 2),
    POSPOS_EQ("pospos_eq", 40, 2),
    POSCYL("poscyl", 50, 2),
    POSCYL_EQ("poscyl_eq", 60, 2),
    DTDT("dtdt", 70, 2), // uncoonstrained two tensor
    CYLCYLCYL("cylcylcyl", 210, 3),
    CYLCYLCYL_EQ("cylcylcyl_eq", 220, 3),
    POSPOSPOS("pospospos", 230, 3),
    POSPOSPOS_EQ("pospospos_eq", 240, 3),
    POSPOSCYL("posposcyl", 250, 3),
    POSPOSCYL_EQ("posposcyl_eq", 260, 3),
    POSCYLCYL("poscylcyl", 270, 3),
    POSCYLCYL_EQ("poscylcyl_eq", 280, 3),
    DTDTDT("dtdtdt", 290, 2), // uncoonstrained three tensor
    MFR("mfr", 1000, 0); // multi fibre reconstruction index
    

    /**
     * @param invName a text string indicating the name of the inversion, exactly as entered
     * on the command line.
     *
     * @param numericalIndex the numerical index for the inversion. This is for legacy 
     * compatibility, use numDTs to find out how many tensors are returned by the inversion.
     *
     * @param ndts the number of diffusion tensors calculated by the inversion.
     */
    InversionIndex(String invName, int numericalIndex, int ndts) {
	name = invName;
	numIndex = numericalIndex;
	numDTs = ndts;
    }
    
    
    public String toString() {
	return name;
    }
    

    /**
     * Gets the indices from one or more strings. This method should be passed all inversion
     * codes that are supplied  on the command line, so that it can distinguish legacy inversion 
     * indices from the text codes.
     * 
     * @param codeStrings contains one or more strings. If there is one string, and the string
     * represents an integer index, the method will attempt to parse two indices. For example,
     * <code>getIndices("32")</code> returns <code>[POSPOS, NLDT_POS]</code>. If multiple 
     * Strings are in the array, a single index is returned per string, for example  
     * <code>getIndices({"1", "1", "22"})</code> returns <code>[LDT, LDT, CYLCYL_EQ]</code>. 
     *
     */
    public static InversionIndex[] getIndices(String[] codeStrings) {


	InversionIndex[] codes = null;

	if (codeStrings.length == 1) {
	    // might be a number
	    try {
		int numericalIndex = Integer.parseInt(codeStrings[0]);

		if (numericalIndex > 10) {
		    codes = new InversionIndex[2];

		    codes[0] = getIndex(numericalIndex - (numericalIndex % 10));
		    codes[1] = getIndex(numericalIndex % 10);
		}
		else {
		    
		    codes = new InversionIndex[1];

		    codes[0] = getIndex(numericalIndex);
		}
		
		return codes;
	    }
	    catch (NumberFormatException e) {
		// carry on
	    }
	    
	}
	
	codes = new InversionIndex[codeStrings.length];
	
	findCode: 
	for (int i = 0; i < codes.length; i++) {
	    
	    for (InversionIndex index : InversionIndex.values()) {
		if (codeStrings[i].equals(index.name)) {
		    codes[i] = index;
		    continue findCode;
		}

	    }  
	    
	    // if we get to here, code is either wrong or a number
	    try {
		codes[i] = getIndex(Integer.parseInt(codeStrings[i]));
	    }
	    catch (NumberFormatException e) {
		throw new LoggedException("Unknown inversion code " + codeStrings[i]);
	    }

	}

	return codes;
	
    }
    

    /**
     * @return a code for a single inversion, given a numerical index. Note that
     * while numerical indices can contain two inversions, eg 21, this method will
     * only return the first one. Thus <code>getIndex(1)</code> returns <code>LDT</code> 
     * and <code>getIndex(21)</code> returns <code>CYLCYL</code>.
     *
     */
    public static InversionIndex getIndex(int ind) {

	InversionIndex code = null;
	
	int numericalIndex = ind;
	
	if (numericalIndex > 10 && numericalIndex % 10 != 0) {
	    numericalIndex = numericalIndex - (numericalIndex % 10);
	}
	
	for (InversionIndex index : InversionIndex.values()) {
	    if (index.numIndex == numericalIndex) {
		return index;
	    }
	}  
	
	throw new LoggedException("Unrecognized tensor inversion code " + ind);

    }
    
    /**
     * The name of the index.
     */
    public final String name;


    /**
     * Numerical index for this inversion. The value 1000 is reserved for NO index. The number
     * corresponds to the old integer inversion codes, eg linear DT inversion is 1. 
     * Two tensor inversions are [1-6]0, three tensor inversions are 2[1-8]0.
     */
    public final int numIndex;


    /**
     * Number of diffusion tensors returned by the inversion. May be zero. Useful to identify
     * which inversion to instantiate.
     *
     */
    public final int numDTs;


   
}

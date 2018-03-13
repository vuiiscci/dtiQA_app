package apps;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;
import tractography.*;


/**
 *
 * Like shredder, but used to sample periodically from a set of streamlines
 *
 * @author Philip Cook
 * @version $Id$
 *
 * 
 *
 */
public class TractShredder {



    public static void main(String[] args) {
	

	OutputManager.outputDataType = "float";
	OutputManager.outputFile = null;

	CL_Initializer.inputModel = "raw";
	
	if (args.length != 3) {
            System.err.println("Usage: tractshredder <offset> <number of tracts> <skip>\n" + 
			       "Program processes I/O input on stdin / stdout.\n");
            System.exit(0);
	}

        int offset = Integer.parseInt(args[0]);
        int bunchSize = Integer.parseInt(args[1]);
        int space = Integer.parseInt(args[2]);
	
	
	OutputManager om = new OutputManager();
	    
	TractSource source = new TractSource(null);
	
	// initial offset
	for (int i = 0; i < offset; i++) {
	    source.nextTract();
	}
	    
	// if we run out of tracts part way through a bunch, output what we have so far

	// if we run out midway through a space, then stop

	// if we get part of a tract, an exception is thrown

	moreTracts: 
	while (source.more()) {

	    for (int i = 0; i < bunchSize; i++) {
		Tract t = source.nextTract();
		    
		om.output(t.toArray());

		if (!source.more()) {
		    break moreTracts;
		}
	    }
	    for (int i = 0; i < space; i++) {
		source.nextTract();
		    
		if (!source.more()) {
		    break moreTracts;
		}

	    }
		
	}
	    
	om.close();
    }




}

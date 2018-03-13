package apps;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;
import tractography.*;


/**
 * Counts the number of streamlines in a file.
 *
 * @author Philip Cook
 * @version $Id$
 *
 */
public class TractCounter {



    public static void main(String[] args) {
	
	CL_Initializer.inputModel = "raw";

       	CL_Initializer.CL_init(args);
        CL_Initializer.checkParsing(args);
	    
	TractSource source = new TractSource(CL_Initializer.inputFile);
	
        int counter = 0;

	while (source.more()) {

            Tract t = source.nextTract();
            counter++;
            
	}
	    
        System.out.println(counter);
    }




}

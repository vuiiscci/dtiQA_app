    package imaging;

import misc.*;

import java.text.*;
import java.util.Vector;
import java.util.Scanner;
import java.util.logging.Logger;



/**
 *
 * Twice-refocused scheme with rectangular gradient pulses.
 *
 *
 * @author Matt Hall, Philip Cook
 * @version $Id$
 *
 */
public final class RectGradTRSE_Scheme extends TRSE_Scheme implements SimulableScheme {


    /**
     * String that identifies this scheme version in the scheme file.
     */
    public static final String VERSION = "TRSE";


    /** logging object */
    private static final Logger logger= Logger.getLogger("camino.imaging.RectGradTRSE_Scheme");

    /** times of onset of gradient blocks */
    private final double[] blockStart = new double[4];
    
    /** times gradient blocks end */
    private final double[] blockEnd = new double[4];
    
    protected RectGradTRSE_Scheme(double[][] gDir, double[] modG, double[] del1, double[] t_del1, 
		       double[] del2, double[] t_del2, double[] del3, double[] t_del4, double[] te) {
	    
    	super(gDir, modG, del1, t_del1, del2, t_del2, del3, t_del4, te);
   
   }
	
    
    protected RectGradTRSE_Scheme(double[][] gDir, RectGradTRSE_Scheme source) {
    	super(gDir, source);
    }

    /**
     * Copy constructor for copying scheme.
     */
    protected RectGradTRSE_Scheme(RectGradTRSE_Scheme source) {
        super(source);
    }

    /**
     *
     * @return the b-value for this sequence.
     *
     */
    public double getB_Value(int i) {

		// Definition of b from Danny
	
		double G = getModG(i);
	
		double del1 = getDel1(i);
		double del2 = getDel2(i);
		double del3 = getDel3(i);
		double del4 = getDel4(i);
	
		double t1 = getT_Del1(i);
		double t2 = getT_Del2(i);
		double t3 = getT_Del3(i);
		double t4 = getT_Del4(i);
	
	 	double b = GAMMA * GAMMA * G * G * (del1*del1*(-t1 + t3 - del1/3.0 + del2 - del3) + 
						2.0*del1*del2*(t3 - t2) + del2*del2*(-t2 + t3 - del2/3.0 + del3) +
						2.0*del1*del3*(t2 - t3 + del3) + del2*del3*(2.0*t2 - 2.0*t3 + del3) +
						del3*del3*(t3 - t2 - del3));
	
		return b;

    }

  
    /**
     * Flip a component of the gradient vectors.
     */
    private DW_Scheme flip(int comp) {
		double[][] flippedG_Dir = new double[numMeas][3];
	
		for (int i = 0; i < numMeas; i++) {
		    if (!zero(i)) {
			flippedG_Dir[i] = getG_Dir(i);
			
			flippedG_Dir[i][comp] = -flippedG_Dir[i][comp];
		    }
		}
	
		return new RectGradTRSE_Scheme(flippedG_Dir, this);

    }


    /**
     * Negates the X component of the gradient directions.
     */
    public DW_Scheme flipX() {
    	return flip(0);
    }


    /**
     * Negates the Y component of the gradient directions.
     */
    public DW_Scheme flipY() {
    	return flip(1);
    }


    /**
     * Negates the Z component of the gradient directions.
     */
    public DW_Scheme flipZ() {
    	return flip(2);
    }


    public DW_Scheme getSubsetScheme(int[] indices) {
	
    	int subset = indices.length;

    	double[][] gDirSub = new double[subset][3];
        double[] modG_Sub = new double[subset];
        double[] del1Sub = new double[subset];
        double[] tDel1Sub= new double[subset];
        double[] del2Sub= new double[subset];
        double[] tDel2Sub= new double[subset];
        double[] del3Sub= new double[subset];
        double[] tDel4Sub= new double[subset];
        double[] teSub = new double[subset];

    	for (int i = 0; i < subset; i++) {
	    
	    double[] gDir = getG_Dir(indices[i]);

    	    gDirSub[i][0] = gDir[0];
    	    gDirSub[i][1] = gDir[1];
    	    gDirSub[i][2] = gDir[2];
    	    
    	    modG_Sub[i] = getModG(indices[i]);
    	    
    	    tDel1Sub[i] = getT_Del1(indices[i]);
    	    del1Sub[i] = getDel1(indices[i]);
    	    tDel2Sub[i] = getT_Del2(indices[i]);
    	    del2Sub[i] = getDel2(indices[i]);
    	    del3Sub[i] = getDel3(indices[i]);
    	    tDel4Sub[i] = getT_Del4(indices[i]);
    	    
    	    teSub[i] = getTE(indices[i]);
    	    
    	}

    	return new RectGradTRSE_Scheme(gDirSub, modG_Sub, del1Sub, tDel1Sub, del2Sub, 
				       tDel2Sub, del3Sub, tDel4Sub, teSub);
    	
    }
    
    
    public DW_Scheme gradOrder(int[] order) {
	//expects a 3 element array consisting of 0 (x dir), 1 (y dir) and 2 (z dir)

	double[][] swappedG_Dir = new double[numMeas][3];
	
		for(int i = 0; i < numMeas; i++) {
		    double[] dir = getG_Dir(i);
	
		    swappedG_Dir[i][0] = dir[order[0]];
		    swappedG_Dir[i][1] = dir[order[1]];
		    swappedG_Dir[i][2] = dir[order[2]];
	
		}
	
		return new RectGradTRSE_Scheme(swappedG_Dir, this);
    } 
  
    
    public double[] getGradImpulse(int i, double t, double tLast) {
        
        final double[] zero = new double[3];
        final double[] netShift = new double[3];
        
        // Deal with some easy cases first
        if (tLast > t) {
            throw new LoggedException("tLast (" + tLast + ") must be < t (" + t + ")");
        }
        
        if (tLast == t) {
            return zero;
        }

        if (zero(i)) {
            return zero;
        }

        blockStart[0] = getT_Del1(i);
        blockStart[1] = getT_Del2(i);
        blockStart[2] = getT_Del3(i);
        blockStart[3] = getT_Del4(i);

        blockEnd[0] = blockStart[0] + getDel1(i);
        blockEnd[1] = blockStart[1] + getDel2(i);
        blockEnd[2] = blockStart[2] + getDel3(i);
        blockEnd[3] = blockStart[3] + getDel4(i);

    	if(tLast > blockEnd[3]) {
    	    return zero;
    	}

        double sgn=1.0;

        double modG = getModG(i);
            
        for(int b=0; b<4; b++) {

            if(b>1){
                sgn=-1.0;
            }

            if (t < blockStart[b]) {
                // dt ends before block b, so we're done
                return netShift;
            }
            
            if (tLast > blockEnd[b]) {
                // tLast after block - do nothing on this iteration
                continue;
            }
            

            if (tLast < blockStart[b]) {
                // tLast before block
                double dt = t < blockEnd[b] ? (t - blockStart[b]) : (blockEnd[b] - blockStart[b]);
                
                //netShift+=sgn * modG * dt;
                for(int j=0; j<3; j++){
                    netShift[j] += sgn*modG*gDir[i][j]*dt;
                }
            }
            else { 
                // tLast >= blockStart[b] && tLast <= blockEnd[b]
 
                // condition is implied otherwise we would have continued (if tLast > blockEnd[b])
                // or gone into if (tLast < blockStart[b])

                double dt = t < blockEnd[b] ? (t - tLast) : (blockEnd[b] - tLast);

                //netShift+= sgn * modG * dt;
                for(int j=0; j<3; j++){
                    netShift[j] += sgn*modG*gDir[i][j]*dt;
                }
            }
                       
        }

        return netShift;
    	
    }


    public double getDuration(){
    	
    	double max=0.0;
    	
    	for(int i=0; i<numMeas; i++){
    		
    		double duration = getTE(i);
    		
    		if(duration>max){
    			max=duration;
    		}
    		
    	}
    	
    	return max;
    	
    }

    /**
     *  
     * Reads a RectGradTRSE_Scheme from a String representation.
     *
     *
     * @param lines the contents of the scheme file, parsed by line, with comments and version 
     * information removed.<BR>
     * <BR> The format is <BR>
     *  g1 g2 g3 modG t_del1 del1 t_del2 del2 del3 t_del4 TE
     *  ...
     *  ...
     *  <BR>
     *  where t_delN is the time of the start of gradient pulse N.
     * 
     *
     */
    protected static RectGradTRSE_Scheme readScheme(Vector<String> lines) {
        // Create the pulse-sequence parameter arrays.
	int numMeas = lines.size();

	if (numMeas == 0) {
	    throw new LoggedException("No measurements in scheme file");
	}

        double[][] gDir = new double[numMeas][3];
        double[] modG = new double[numMeas];
        
        double[] del1 = new double[numMeas];
        double[] t_del1 = new double[numMeas];

        double[] del2 = new double[numMeas];
        double[] t_del2 = new double[numMeas];
        
        double[] del3 = new double[numMeas];
        double[] t_del3 = new double[numMeas];
        
	double[] t_del4 = new double[numMeas];
        
        
        double[] TE = new double[numMeas];
      
        // Parse the lines into the arrays.
        try {

            for(int i=0; i<numMeas; i++) {
            
                String[] tokens = lines.elementAt(i).trim().split("\\s+");
                
                if (tokens.length != 11) {
                    throw new LoggedException("Incorrect number of values in measurement " + i);
                }
                
                
                gDir[i][0] = Double.parseDouble(tokens[0]);
                gDir[i][1] = Double.parseDouble(tokens[1]);
                gDir[i][2] = Double.parseDouble(tokens[2]);
                modG[i] = Double.parseDouble(tokens[3]);
                
                t_del1[i] = Double.parseDouble(tokens[4]);
                del1[i] = Double.parseDouble(tokens[5]);
                
                t_del2[i] = Double.parseDouble(tokens[6]);
                del2[i] = Double.parseDouble(tokens[7]);
                
                del3[i]= Double.parseDouble(tokens[8]);
                
                t_del4[i] = Double.parseDouble(tokens[9]);
                
                TE[i] = Double.parseDouble(tokens[10]);
                
                if (modG[i] == 0.0) {
                    // measurement is zero
                    
                    if (gDir[i][0] != 0.0 || gDir[i][1] != 0.0 || gDir[i][2] != 0.0) {
                        logger.info("Zero diffusion weighting in measurement " + i + 
                                    ", setting gradient direction to zero");

                        gDir[i][0] = 0.0;
                        gDir[i][1] = 0.0;
                        gDir[i][2] = 0.0;
		    
                    }
                }
                
            }
        }
        catch (NumberFormatException e) {
            throw new LoggedException(e);
        }
	
	return new RectGradTRSE_Scheme(gDir, modG, del1, t_del1, del2, t_del2, del3, t_del4, TE);
	
    }

 


    public String toString() {

	DecimalFormat df = new DecimalFormat("   0.000000;  -0.000000");
	
	StringBuffer sb = new StringBuffer();


	sb.append("#Twice-refocused scheme with rectangular gradient pulses\n");
	sb.append("#g_x\tg_y\tg_z\t|g|\tt_del1\tdel1\tt_del2\tdel2\tdel3\tt_del4\tTE\n");
	sb.append("VERSION: " + VERSION);
	sb.append("\n");


	for (int i = 0; i < numMeas; i++) {
	    
	    double[] gDir = getG_Dir(i);
	    
	    sb.append(df.format(gDir[0]));
	    sb.append(df.format(gDir[1]));
	    sb.append(df.format(gDir[2]));
	    sb.append(df.format(getModG(i)));
	    sb.append(df.format(getT_Del1(i)));
	    sb.append(df.format(getDel1(i)));
	    sb.append(df.format(getT_Del2(i)));
	    sb.append(df.format(getDel2(i)));
	    sb.append(df.format(getDel3(i)));
	    sb.append(df.format(getT_Del4(i)));
	    sb.append(df.format(getTE(i)));
	    sb.append("\n");
	}
	
	return sb.toString();

    }

}

package imaging;

import java.util.Vector;
import java.util.logging.Logger;

import misc.LoggedException;

/**
 * This behaves like a PGSE sequence, except for the fact that gradients are
 * parabolic. The gradient direction determines the orientation of the parabola
 * which is centred in the middle of the substrate. There's also a need to subtract 
 * the initial phase shift from subsequent values to ensure good convergence.
 * 
 * There's no a lot of difference between this and RectGradSteTanScheme,  
 * 
 * @author matt (matt.hall@ucl.ac.uk)
 *
 */
public final class RectQuadraticGradSteTanScheme extends RectGradSteTanScheme {

    /**
     * String that identifies this scheme version in the scheme file.
     */
    public static final String VERSION = "QUAD_GRAD_PGSE";

    /**
     * logging object
     */
    private static final Logger logger= Logger.getLogger("imaging.RectQuadraticGradSteTanScheme");
	
    /** 
     * enumerator for contrast mechanism: magnitude or phase
     * 
     */
    private enum Contrast{MAG, PHASE};
    
    /**
     * array of labels for phase or magnitude on each entry in the scheme
     */
    private final Contrast[] contrast;
    
    
	protected RectQuadraticGradSteTanScheme(double[][] initG_Dir,
			double[] initModG, double[] initBigDel, double[] initSmallDel,
			double[] initTE, Contrast[] contrast) {
		super(initG_Dir, initModG, initBigDel, initSmallDel, initTE);
		
		this.contrast= contrast;
	}

		
    /**
     * Copy constructor for flipping / gradient component swapping methods.
     */
    protected RectQuadraticGradSteTanScheme(double[][] initG_Dir, RectQuadraticGradSteTanScheme source) {
        super(initG_Dir, source);
        
        this.contrast= source.contrast;
    }

    
    /**
     * checks the contrast mechanism on the measurement.
     * 
     * @return true of magnitude, false if phase
     */
    public boolean isMagnitude(int i){
    	if(contrast[i]==Contrast.MAG){
    		return true;
    	}
    	
    	return false;
    }
    
    
    /**
     * method that returns the gradient strength multiplied by the duration of the
     * pulse in the update. This is an alternative version of the gradient impulse
     * approach, which is not (quite!) compatible with nonlinear gradients.
     */
    public final double getGradientWeight(int i, double t, double tLast){
    	
    	
        final double[] zero= new double[3];
        final double[] gradImpulse= new double[3];
        
        double sgn= 1.0;

        double te = getTE(i);

        if (!(t >= 0.0 && t <= te)) {
            return 0.0;
        }

        if (zero(i)) {
            return 0.0;
        }

        double DELTA = getDELTA(i);
        double delta = getDelta(i);


        // time before and after the gradient blocks
        double pad = (te - DELTA - delta) / 2.0;


        double firstBlockStart= pad;
        double firstBlockEnd= pad+delta;
        double secondBlockStart= pad+DELTA;
        double secondBlockEnd= pad+DELTA+delta;
        
        
        // dt is between pad and first block 
        if((t>=firstBlockStart)&&(t<firstBlockEnd)){
            if(tLast < firstBlockStart) {
                double dt = t - firstBlockStart;
                
                double modG= getModG(i);
                                
                return sgn*modG*dt;
            }
        }


        // if between blocks
        if((t>=firstBlockEnd)&&(t<secondBlockStart)){
            if (tLast<firstBlockEnd) {
                // the block ended between this call and the last one
                // so need to calculate the partial contribution
                double dt= firstBlockEnd-tLast;

                double modG= getModG(i);
                
                return sgn*modG*dt;

            }

            return 0.0;
        }

        // flip the sign of the gradient if in the second block 

        // WARNING - we assume that we never cross blocks, ie there is never a part contribution
        // from different blocks

        if(t>=secondBlockStart){
            sgn=-1.0;
        }

        // if after second block
        if(t>=secondBlockEnd) {
            if(tLast<secondBlockEnd) {
                // the black ended between this call and the last one
                // so need to calculate the partial contribution
                double dt= secondBlockEnd-tLast;

                double modG= getModG(i);
                
                return sgn*modG*dt;

            }


            return 0.0;
        }


        if((t>=secondBlockStart)&&(tLast<secondBlockStart)){
            // the block ended between this call and the last one
            // so need to calculate the partial contribution
            double dt= t-secondBlockStart;

            double modG= getModG(i);
            
            return sgn*modG*dt;

        }

        // if before first block or after second
        if (t < pad || t > te - pad) {
            return 0.0;
        }

        // return gradient from entire timestep    	
        double modG= getModG(i);
        
        return sgn*modG*(t-tLast);

    }
    
    /**
     * Reads a RectGradSteTanScheme from a String representation.
     *
     *
     * @param lines the contents of the scheme file, parsed by line, with comments and version 
     * information removed.
     * 
     */
    protected static RectQuadraticGradSteTanScheme readScheme(Vector<String> lines) {
        // Create the pulse-sequence parameter arrays.
        int numMeas = lines.size();

        if (numMeas == 0) {
            throw new LoggedException("No measurements in scheme file");
        }

        double[][] gDir = new double[numMeas][3];
        double[] modG = new double[numMeas];
        double[] bigDel = new double[numMeas];
        double[] smallDel = new double[numMeas];
        double[] TE = new double[numMeas];
        Contrast[] contrast= new Contrast[numMeas];

        // Parse the lines into the arrays.

        try {

            for(int i=0; i<numMeas; i++) {
                
                String[] tokens = lines.elementAt(i).trim().split("\\s+");
                
                if (tokens.length != 8) {
                    throw new LoggedException("Incorrect number of values in measurement " + i);
                }
                
                
                gDir[i][0] = Double.parseDouble(tokens[0]);
                gDir[i][1] = Double.parseDouble(tokens[1]);
                gDir[i][2] = Double.parseDouble(tokens[2]);
                modG[i] = Double.parseDouble(tokens[3]);
                bigDel[i] = Double.parseDouble(tokens[4]);
                smallDel[i] = Double.parseDouble(tokens[5]);
                TE[i] = Double.parseDouble(tokens[6]);
                if(tokens[7].equalsIgnoreCase("M")){
                	contrast[i]=Contrast.MAG;
                }
                else if(tokens[7].equalsIgnoreCase("P")){
                	contrast[i]=Contrast.PHASE;
                }
                else{
                	throw new LoggedException("Unrecognised contrast specifier '"
                								+tokens[7]+"'. Please use M (magnitude) or P (phase)");
                }
                
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
        

        return new RectQuadraticGradSteTanScheme(gDir, modG, bigDel, smallDel, TE, contrast);

    }

    
    
    
}

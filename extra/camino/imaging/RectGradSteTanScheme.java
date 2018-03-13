package imaging;

import misc.*;

import java.text.*;
import java.util.logging.Logger;
import java.util.Scanner;
import java.util.Vector;


/**
 * 
 * Stejskal-Tanner scheme with rectangular gradient blocks.
 *
 * @author Philip Cook
 * @version $Id$
 *  
 */
public class RectGradSteTanScheme extends StejskalTannerScheme{


    /**
     * String that identifies this scheme version in the scheme file.
     */
    public static final String VERSION = "STEJSKALTANNER";


    /**
     * logging object
     */
    private static final Logger logger = Logger.getLogger("camino.imaging.RectGradSteTanScheme");



    /**
     * Construct a scheme from the imaging parameters.
     *
     * @param initG_Dir gradient directions.
     * @param initModG the average gradient strengths.
     * @param initBigDel gradient separations.
     * @param initSmallDel gradient durations.
     * @param initTE echo times.
     *
     */
    protected RectGradSteTanScheme(double[][] initG_Dir, double[] initModG, double[] initBigDel, 
            double[] initSmallDel, double[] initTE) {

        super(initG_Dir, initModG, initBigDel, initSmallDel, initTE);

    }

    /**
     * Copy constructor for flipping / gradient component swapping methods.
     */
    protected RectGradSteTanScheme(double[][] initG_Dir, RectGradSteTanScheme source) {
        super(initG_Dir, source);
    }
    
    /**
     * Copy constructor for copying scheme.
     */
    protected RectGradSteTanScheme(RectGradSteTanScheme source) {
        super(source);
    }


    /**
     * Gets the b-value, assuming an "effective diffusion time" of (DELTA - delta / 3).
     *
     * @return |q|^2 * (DELTA - delta / 3).
     */
    public double getB_Value(int i) {

        double modQ = getModQ(i);

        return modQ * modQ * (getDELTA(i) - getDelta(i) / 3.0);
    }

    /**
     *
     */
    public double[] getB_Values() {
        return null;
    }
        /**
     *
     */



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

        return new RectGradSteTanScheme(flippedG_Dir, this);

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
        double[] bigDelSub = new double[subset];
        double[] smallDelSub = new double[subset];
        double[] teSub = new double[subset];

        for (int i = 0; i < subset; i++) {
            gDirSub[i] = getG_Dir(indices[i]);
            modG_Sub[i] = getModG(indices[i]);
            bigDelSub[i] = getDELTA(indices[i]);
            smallDelSub[i] = getDelta(indices[i]);
            teSub[i] = getTE(indices[i]);

        }

        return new RectGradSteTanScheme(gDirSub, modG_Sub, bigDelSub, smallDelSub, teSub);

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

        return new RectGradSteTanScheme(swappedG_Dir, this);
    } 


    /**
     * Gets the gradient strength at a specific time point during a measurement. The gradient 
     * block is assumed to take place in the middle of the echo time, ie the time between the 
     * start of the echo time and the first diffusion-weighting gradient is the same as the time 
     * between the end of the second diffusion-weighting gradient and the end of the echo time.
     *  
     * @return gradient strength averaged over a duration (tLast - t) during a measurement i.  
     *
     */
    public double[] getGradImpulse(int i, double t, double tLast) {

        final double[] zero= new double[3];
        final double[] gradImpulse= new double[3];
        
        double sgn= 1.0;

        double te = getTE(i);

        if (!(t >= 0.0 && t <= te)) {
            return zero;
        }

        if (zero(i)) {
            return zero;
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

                /*double rawGrad= (getModG(i))*dt;

                return  sgn*rawGrad;*/
                
                double modG= getModG(i);
                
                for(int j=0; j< gradImpulse.length; j++){
                    gradImpulse[j]= sgn*modG*gDir[i][j]*dt;
                }
                return gradImpulse;
            }
        }


        // if between blocks
        if((t>=firstBlockEnd)&&(t<secondBlockStart)){
            if (tLast<firstBlockEnd) {
                // the block ended between this call and the last one
                // so need to calculate the partial contribution
                double dt= firstBlockEnd-tLast;

                /*double rawGrad= (getModG(i))*dt;

                return  sgn*rawGrad;*/
                
                double modG= getModG(i);
                
                for(int j=0; j< gradImpulse.length; j++){
                    gradImpulse[j]= sgn*modG*gDir[i][j]*dt;
                }
                
                return gradImpulse;

            }

            return zero;
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
                /*double rawGrad= getModG(i)*dt;

                return sgn*(rawGrad);*/
                double modG= getModG(i);
                
                for(int j=0; j< gradImpulse.length; j++){
                    gradImpulse[j]= sgn*modG*gDir[i][j]*dt;
                }

                return gradImpulse;

            }


            return zero;
        }


        if((t>=secondBlockStart)&&(tLast<secondBlockStart)){
            // the block ended between this call and the last one
            // so need to calculate the partial contribution
            double dt= t-secondBlockStart;
            /*double rawGrad= getModG(i)*dt;

            return sgn*(rawGrad);*/
            double modG= getModG(i);
            
            for(int j=0; j< gradImpulse.length; j++){
                gradImpulse[j]= sgn*modG*gDir[i][j]*dt;
            }

            return gradImpulse;

        }

        // if before first block or after second
        if (t < pad || t > te - pad) {
            return zero;
        }

        // return gradient from entire timestep    	
        /*double retVal= getModG(i)*(t-tLast);

        return  sgn*(retVal);*/
        double modG= getModG(i);
        
        for(int j=0; j< gradImpulse.length; j++){
            gradImpulse[j]= sgn*modG*gDir[i][j]*(t-tLast);
        }

        return gradImpulse;


    }

    /**
     * @return the scan duration in seconds
     */
    public double getDuration(){

        double max=0.0;

        for(int i=0; i<numMeas; i++){

            double duration= getTE(i);

            if(duration>max){
                max= duration;
            }
        }

        return max;
    }


    /**
     * Reads a RectGradSteTanScheme from a String representation.
     *
     *
     * @param lines the contents of the scheme file, parsed by line, with comments and version 
     * information removed.
     * 
     */
    protected static RectGradSteTanScheme readScheme(Vector<String> lines) {
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

        // Parse the lines into the arrays.

        try {

            for(int i=0; i<numMeas; i++) {
                
                String[] tokens = lines.elementAt(i).trim().split("\\s+");
                
                if (tokens.length != 7) {
                    throw new LoggedException("Incorrect number of values in measurement " + i);
                }
                
                
                gDir[i][0] = Double.parseDouble(tokens[0]);
                gDir[i][1] = Double.parseDouble(tokens[1]);
                gDir[i][2] = Double.parseDouble(tokens[2]);
                modG[i] = Double.parseDouble(tokens[3]);
                bigDel[i] = Double.parseDouble(tokens[4]);
                smallDel[i] = Double.parseDouble(tokens[5]);
                TE[i] = Double.parseDouble(tokens[6]);
                
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
        

        return new RectGradSteTanScheme(gDir, modG, bigDel, smallDel, TE);

    }


    /**
     * Returns the String representation of this object.
     *
     */
    public String toString() {
        DecimalFormat df = new DecimalFormat("   0.000000;  -0.000000");

        StringBuffer sb = new StringBuffer();


        sb.append("#Stejskal-Tanner scheme with rectangular gradient pulses\n");
        sb.append("#g_x\tg_y\tg_z\t|g|\tDELTA\tdelta\tTE\n");
        sb.append("VERSION: " + VERSION);
        sb.append("\n");

        for (int i = 0; i < numMeas; i++) {

            double[] gDir = getG_Dir(i);

            sb.append(df.format(gDir[0]));
            sb.append(df.format(gDir[1]));
            sb.append(df.format(gDir[2]));
            sb.append(df.format(getModG(i)));
            sb.append(df.format(getDELTA(i)));
            sb.append(df.format(getDelta(i)));
            sb.append(df.format(getTE(i)));
            sb.append("\n");
        }

        return sb.toString();

    }






}

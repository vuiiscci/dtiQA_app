    package imaging;

import misc.*;

import java.text.*;
import java.util.logging.Logger;
import java.util.Scanner;
import java.util.Vector;

/**
 * DPFG scheme with rectangular blocks
 * 
 * @author gmorgan
 * @version $Id$
 *
 */

public final class RectGradDPFG_Scheme extends DPFG_Scheme implements SimulableScheme {
	
    /**
     * String that identifies this scheme version in the scheme file.
     */
    public static final String VERSION = "DPFG";


    /**
     * logging object
     */
    private static final Logger logger = Logger.getLogger("camino.imaging.RectDPFG_Scheme");



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
    protected RectGradDPFG_Scheme(double[][] initG_Dir1, double[] initModG1, double[] initBigDel1, 
    		double[] initSmallDel1, double[] initTM, double[][] initG_Dir2, double[] initModG2, double[] initBigDel2,
    		double[] initSmallDel2, double[] initTE) {

        super(initG_Dir1, initModG1, initBigDel1, initSmallDel1, initTM, initG_Dir2, 
       	     initModG2, initBigDel2, initSmallDel2, initTE);

    }

    /**
     * Copy constructor for flipping / gradient component swapping methods.
     */
    protected RectGradDPFG_Scheme(double[][] initG_Dir1, double[][] initG_Dir2, RectGradDPFG_Scheme source) {
        super(initG_Dir1, initG_Dir2, source);
    }

    /**
     * Copy constructor for copying scheme.
     */
    protected RectGradDPFG_Scheme(RectGradDPFG_Scheme source) {
        super(source);
    }
    /**
     * Haven't implemented bvalue yet (but shouldn't matter for MC simulation)
     * Gets the b-value, assuming an "effective diffusion time" of (DELTA - delta / 3).
     *
     * @param i - the measurement in the scheme file
     * @return |q|^2 * (DELTA - delta / 3).
     */
    public double getB_Value(int i) {
    	
    	/* For a double PGSE scheme file, I'm pretty certain that the total b value will just
    	 * be the sum of the two individual b values - will check this though!
    	 */
    	
    	// final b value to be returned
        double b= 0.0;
        
        // b values for individual PGSE blocks
        double b1 = 0.0;
        double b2 = 0.0;
        
        // get parameters for first PGSE block
        double modG1 = getModG1(i);
        double del1 = getSmallDel1(i);
        double DEL1 = getBigDel1(i);
        
        // get parameters for second PGSE block
        double modG2 = getModG2(i);
        double del2 = getSmallDel2(i);
        double DEL2 = getBigDel2(i);
        
        // calculate b values components, up to the gyromagnetic constant, for both b1 and b2
        b1 = modG1*modG1*del1*del1*(DEL1-(del1/3.0));
        b2 = modG2*modG2*del2*del2*(DEL2-(del2/3.0));
        
        // calculate total b
        b = DW_Scheme.GAMMA*DW_Scheme.GAMMA*(b1+b2);
    
        return b;

    	
    }


    /**
     * Flip a component of the gradient vectors.
     */
    private DW_Scheme flip(int comp) {
        double[][] flippedG_Dir1 = new double[numMeas][3];
        double[][] flippedG_Dir2 = new double[numMeas][3];

        for (int i = 0; i < numMeas; i++) {
            if (!zero(i)) {
                flippedG_Dir1[i] = getG_Dir(i);
                flippedG_Dir2[i] = getG2(i);

                flippedG_Dir1[i][comp] = -flippedG_Dir1[i][comp];
                flippedG_Dir2[i][comp] = -flippedG_Dir2[i][comp];
            }
        }

        return new RectGradDPFG_Scheme(flippedG_Dir1, flippedG_Dir2, this);

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

        double[][] g1DirSub = new double[subset][3];
        double[] modG1Sub = new double[subset];
        double[] bigDel1Sub = new double[subset];
        double[] smallDel1Sub = new double[subset];
        double[] tmSub = new double[subset];
        double[][] g2DirSub = new double[subset][3];
        double[] modG2Sub = new double[subset];
        double[] bigDel2Sub = new double[subset];
        double[] smallDel2Sub = new double[subset];
        double[] teSub = new double[subset];

        for (int i = 0; i < subset; i++) {
            g1DirSub[i] = getG_Dir(indices[i]);
            modG1Sub[i] = getModG1(indices[i]);
            bigDel1Sub[i] = getBigDel1(indices[i]);
            smallDel1Sub[i] = getSmallDel1(indices[i]);
            tmSub[i] = getTM(indices[i]);
            g2DirSub[i] = getG_Dir(indices[i]);
            modG2Sub[i] = getModG2(indices[i]);
            bigDel2Sub[i] = getBigDel2(indices[i]);
            smallDel2Sub[i] = getSmallDel2(indices[i]);
            teSub[i] = getTE(indices[i]);

        }

        return new RectGradDPFG_Scheme(g1DirSub, modG1Sub, bigDel1Sub, smallDel1Sub, tmSub, g2DirSub, 
          	     modG2Sub, bigDel2Sub, smallDel2Sub, teSub);

    }


    public DW_Scheme gradOrder(int[] order) {
        //expects a 3 element array consisting of 0 (x dir), 1 (y dir) and 2 (z dir)

        double[][] swappedG1_Dir = new double[numMeas][3];
        double[][] swappedG2_Dir = new double[numMeas][3];

        for(int i = 0; i < numMeas; i++) {
            double[] dir1 = getG_Dir(i);
            double[] dir2 = getG2(i);

            swappedG1_Dir[i][0] = dir1[order[0]];
            swappedG1_Dir[i][1] = dir1[order[1]];
            swappedG1_Dir[i][2] = dir1[order[2]];
            
            swappedG2_Dir[i][0] = dir2[order[0]];
            swappedG2_Dir[i][1] = dir2[order[1]];
            swappedG2_Dir[i][2] = dir2[order[2]];

        }

        return new RectGradDPFG_Scheme(swappedG1_Dir, swappedG2_Dir, this);
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
        double dt;
        
        double sgn= 1.0;

        double te = getTE(i);

        if (!(t >= 0.0 && t <= te)) {
            return zero;
        }

        if (zero(i)) {
            return zero;
        }

        double DELTA1 = getBigDel1(i);
        double delta1 = getSmallDel1(i);
        double DELTA2 = getBigDel2(i);
        double delta2 = getSmallDel2(i);
        double TM = getTM(i);
        double[] gDir2 = getG2(i);


        // time before and after the gradient blocks
        double pad = (te - (DELTA1 + delta1 +TM + DELTA2 + delta2)) / 2.0;


        double firstBlockStart= pad;
        double firstBlockEnd= pad+delta1;
        double secondBlockStart= pad+DELTA1;
        double secondBlockEnd= pad+DELTA1+delta1;
        double thirdBlockStart = pad+DELTA1+delta1+TM;
        double thirdBlockEnd = pad+DELTA1+delta1+TM+delta2;
        double fourthBlockStart = pad+DELTA1+delta1+TM+DELTA2;
        double fourthBlockEnd = pad+DELTA1+delta1+TM+DELTA2+delta2;

        
        
        // dt is between pad and first block 
        if((t>=firstBlockStart)&&(t<firstBlockEnd)){
            if(tLast < firstBlockStart) {
                dt = t - firstBlockStart;
            }
            else{
            	dt = t - tLast;
            }

            /*double rawGrad= (getModG(i))*dt;

                return  sgn*rawGrad;*/

            double modG= getModG1(i);

            for(int j=0; j< gradImpulse.length; j++){
            	gradImpulse[j]= sgn*modG*gDir[i][j]*dt;
            }
            return gradImpulse;
            
        }


        // if between blocks
        if((t>=firstBlockEnd)&&(t<secondBlockStart)){
            if (tLast<firstBlockEnd) {
                // the block ended between this call and the last one
                // so need to calculate the partial contribution
                dt= firstBlockEnd-tLast;

                /*double rawGrad= (getModG(i))*dt;

                return  sgn*rawGrad;*/
                
                double modG= getModG1(i);
                
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

        if(t>=secondBlockStart && t<fourthBlockStart){
            sgn=-1.0;
        }

        // if after second block
        if((t>=secondBlockStart) && (t<secondBlockEnd)) {
        	if(tLast<secondBlockStart) {
        		// the black ended between this call and the last one
        		// so need to calculate the partial contribution
        		dt= t - secondBlockStart;
        	}
        	else{
        		dt = t - tLast;
        	}

        	/*double rawGrad= getModG(i)*dt;

                return sgn*(rawGrad);*/
        	double modG= getModG1(i);

        	for(int j=0; j< gradImpulse.length; j++){
        		gradImpulse[j]= sgn*modG*gDir[i][j]*dt;
        	}

        	return gradImpulse;
        }


        if((t>=secondBlockEnd)&&(t<thirdBlockStart)){
        	if (tLast<secondBlockEnd) {
                // the block ended between this call and the last one
                // so need to calculate the partial contribution
                dt= secondBlockEnd-tLast;

                /*double rawGrad= (getModG(i))*dt;

                return  sgn*rawGrad;*/
                
                double modG= getModG1(i);
                
                for(int j=0; j< gradImpulse.length; j++){
                    gradImpulse[j]= sgn*modG*gDir[i][j]*dt;
                }
                
                return gradImpulse;

            }

            return zero;

        }
        
     // dt is between pad and first block 
        if((t>=thirdBlockStart)&&(t<thirdBlockEnd)){
            if(tLast < thirdBlockStart) {
                dt = t - thirdBlockStart;
            }
            else{
            	dt = t - tLast;
            }

            /*double rawGrad= (getModG(i))*dt;

                return  sgn*rawGrad;*/

            double modG= getModG2(i);

            for(int j=0; j< gradImpulse.length; j++){
            	gradImpulse[j]= sgn*modG*gDir2[j]*dt;
            }
            return gradImpulse;
            
        }


        // if between blocks
        if((t>=thirdBlockEnd)&&(t<fourthBlockStart)){
            if (tLast<thirdBlockEnd) {
                // the block ended between this call and the last one
                // so need to calculate the partial contribution
                dt= thirdBlockEnd-tLast;

                /*double rawGrad= (getModG(i))*dt;

                return  sgn*rawGrad;*/
                
                double modG= getModG2(i);
                
                for(int j=0; j< gradImpulse.length; j++){
                    gradImpulse[j]= sgn*modG*gDir2[j]*dt;
                }
                
                return gradImpulse;

            }

            return zero;
        }


        // if after second block
        if((t>=fourthBlockStart) && (t<fourthBlockEnd)) {
        	if(tLast<fourthBlockStart) {
        		// the black ended between this call and the last one
        		// so need to calculate the partial contribution
        		dt= t - fourthBlockStart;
        	}
        	else{
        		dt = t - tLast;
        	}

        	/*double rawGrad= getModG(i)*dt;

                return sgn*(rawGrad);*/
        	double modG= getModG2(i);

        	for(int j=0; j< gradImpulse.length; j++){
        		gradImpulse[j]= sgn*modG*gDir2[j]*dt;
        	}

        	return gradImpulse;
        }


        if(t>=fourthBlockEnd){
        	if (tLast<fourthBlockEnd) {
                // the block ended between this call and the last one
                // so need to calculate the partial contribution
                dt= fourthBlockEnd-tLast;

                /*double rawGrad= (getModG(i))*dt;

                return  sgn*rawGrad;*/
                
                double modG= getModG2(i);
                
                for(int j=0; j< gradImpulse.length; j++){
                    gradImpulse[j]= sgn*modG*gDir2[j]*dt;
                }
                
                return gradImpulse;

            }

            return zero;

        }

        // if before first block or after fourth
        if (t < pad || t > te - pad) {
            return zero;
        }
        
        return zero;


/*        // return gradient from entire timestep    	
        double retVal= getModG(i)*(t-tLast);

        return  sgn*(retVal);
        double modG= getModG(i);
        
        for(int j=0; j< gradImpulse.length; j++){
            gradImpulse[j]= sgn*modG*gDir[i][j]*(t-tLast);
        }

        return gradImpulse;*/


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
     * Reads a RectGradDPFG_Scheme from a String representation.
     *
     *
     * @param lines the contents of the scheme file, parsed by line, with comments and version 
     * information removed.
     * 
     */
    protected static RectGradDPFG_Scheme readScheme(Vector<String> lines) {
        // Create the pulse-sequence parameter arrays.
        int numMeas = lines.size();

        if (numMeas == 0) {
            throw new LoggedException("No measurements in scheme file");
        }

        double[][] g1Dir = new double[numMeas][3];
        double[] modG1 = new double[numMeas];
        double[] bigDel1 = new double[numMeas];
        double[] smallDel1 = new double[numMeas];
        double[] TM = new double[numMeas];
        double[][] g2Dir = new double[numMeas][3];
        double[] modG2 = new double[numMeas];
        double[] bigDel2 = new double[numMeas];
        double[] smallDel2 = new double[numMeas];
        double[] TE = new double[numMeas];

        // Parse the lines into the arrays.

        try {

            for(int i=0; i<numMeas; i++) {
                
                String[] tokens = lines.elementAt(i).trim().split("\\s+");
                
                if (tokens.length != 14) {
                    throw new LoggedException("Incorrect number of values in measurement " + i);
                }
                
                
                g1Dir[i][0] = Double.parseDouble(tokens[0]);
                g1Dir[i][1] = Double.parseDouble(tokens[1]);
                g1Dir[i][2] = Double.parseDouble(tokens[2]);
                modG1[i] = Double.parseDouble(tokens[3]);
                bigDel1[i] = Double.parseDouble(tokens[4]);
                smallDel1[i] = Double.parseDouble(tokens[5]);
                TM[i] = Double.parseDouble(tokens[6]);
                g2Dir[i][0] = Double.parseDouble(tokens[7]);
                g2Dir[i][1] = Double.parseDouble(tokens[8]);
                g2Dir[i][2] = Double.parseDouble(tokens[9]);
                modG2[i] = Double.parseDouble(tokens[10]);
                bigDel2[i] = Double.parseDouble(tokens[11]);
                smallDel2[i] = Double.parseDouble(tokens[12]);
                TE[i] = Double.parseDouble(tokens[13]);
                
                if (modG1[i] == 0.0 && modG2[i] == 0.0) {
                    // measurement is zero
                    
                    if (g1Dir[i][0] != 0.0 || g1Dir[i][1] != 0.0 || g1Dir[i][2] != 0.0 || g2Dir[i][0] !=0.0 || g2Dir[i][1] != 0.0 || g2Dir[i][2] !=0.0) {
                        logger.info("Zero diffusion weighting in measurement " + i + 
                    ", setting gradient direction to zero");
                        
                        g1Dir[i][0] = 0.0;
                        g1Dir[i][1] = 0.0;
                        g1Dir[i][2] = 0.0;
                        g2Dir[i][0] = 0.0;
                        g2Dir[i][1] = 0.0;
                        g2Dir[i][2] = 0.0;
                        
                    }
                    
                }
                
            }

        }
        catch (NumberFormatException e) {
            throw new LoggedException(e);
        }
        

        return new RectGradDPFG_Scheme(g1Dir, modG1, bigDel1, smallDel1, TM, g2Dir, 
         	     modG2, bigDel2, smallDel2, TE);

    }


    /**
     * Returns the String representation of this object.
     *
     */
    public String toString() {
        DecimalFormat df = new DecimalFormat("   0.000000;  -0.000000");

        StringBuffer sb = new StringBuffer();


        sb.append("#dPFG scheme with rectangular gradient pulses\n");
        sb.append("#g1_x\tg1_y\tg1_z\t|g1|\tDELTA1\tdelta1\tTM\tg2_x\tg2_y\tg2_z\t|g2|\tDELTA1\tdelta2\tTE\n");
        sb.append("VERSION: " + VERSION);
        sb.append("\n");

        for (int i = 0; i < numMeas; i++) {

            double[] g1Dir = getG_Dir(i);
            double[] g2Dir = getG2(i);

            sb.append(df.format(g1Dir[0]));
            sb.append(df.format(g1Dir[1]));
            sb.append(df.format(g1Dir[2]));
            sb.append(df.format(getModG1(i)));
            sb.append(df.format(getBigDel1(i)));
            sb.append(df.format(getSmallDel1(i)));
            sb.append(df.format(getTM(i)));
            sb.append(df.format(g2Dir[0]));
            sb.append(df.format(g2Dir[1]));
            sb.append(df.format(g2Dir[2]));
            sb.append(df.format(getModG2(i)));
            sb.append(df.format(getBigDel2(i)));
            sb.append(df.format(getSmallDel2(i)));
            sb.append(df.format(getTE(i)));
            sb.append("\n");
        }

        return sb.toString();

    }

}

package data;

import misc.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Simple model for the diffusion displacement density function .
 * 
 * <dt>Description:
 * 
 * <dd>Models the displacement by a mixture of Gaussian densities.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id: GaussianMixture.java 58 2006-10-26 11:38:18Z ucacmgh $
 *  
 */
public class BallStick implements ModelPDF {

	/**
	 * The orientation of the stick
	 */
	private double[] orientation;
	/**
	 * diffusivity 
	 */
	private double diffusivity;
	/**
	 * volume fraction of the stick
	 */
	private double volfrac;
	/**
	 * allows us to implement at and ft
	 */
	private GaussianMixture gm;
	
	private DT [] dt;

	/**
	 * constructor for specified diffusivity, volume fraction and
	 * orientation
	 */
	public BallStick(double diff, double vf, double[] ori){
		orientation = ori;
		volfrac = vf;
		diffusivity = diff;
		dt = new DT [2];
	}
	
	public double at(double[] x, double tau) {

    	if (gm == null){
    		makeGaussianMix();
    	}
    	
    	return gm.at(x, tau);
    }
	
    public double ftAt(double[] q, double tau) {
	  
    	if (gm == null){
    		makeGaussianMix();
    	}
        return gm.ftAt(q, tau);
    }


    public double ftAtB_Vec(double[] g, double b) {
	  
    	if (gm == null){
    		makeGaussianMix();
    	}
        return gm.ftAt(g, b);
    }

    
    private void makeGaussianMix() {
    	//remember ordering of DT values is dxx, dxy, dxz, dyy, dyz, dzz
    	//DT[] dt = new DT[2];
    //	dt[0] = new DT(1.0,2.0,3.0,4.0,5.0,6.0);
    //	dt[1] = new DT(4.0,2.0,3.0,4.0,5.0,9.0);
    	dt[0] = new DT(diffusivity, 0.0, 0.0, diffusivity, 0.0, diffusivity);
    	dt[1] = new DT(diffusivity*orientation[0]*orientation[0], diffusivity*orientation[0]*orientation[1], diffusivity*orientation[0]*orientation[2], diffusivity*orientation[1]*orientation[1], diffusivity*orientation[1]*orientation[2], diffusivity*orientation[2]*orientation[2]);
    	double [] mixtures = new double [2];
    //	mixtures[1] = 0.6;
    //	mixtures[2] = 0.4;
    	mixtures[0] = 1-volfrac;
    	mixtures[1] = volfrac;
    	gm = new GaussianMixture(dt, mixtures);
    }


    /**
     * Returns a list of principal directions.
     * 
     * @return Array of principal directions.
     */
    public double[][] getPDs() {
        double[][] dirs = new double[dt.length][3]; //should this be dt.length - 1? 
        for (int i = 0; i < dt.length; i++) {
            double[] dir = dt[i].getPD();
            for (int j = 0; j < 3; j++) {
                dirs[i][j] = dir[j];
            }
        }

        return dirs;
    }
    
    public double[] getInformationArray(){
    	double[] infoarray = new double[5];
    	infoarray[0] = diffusivity;
    	infoarray[1] = volfrac;
    	infoarray[2] = orientation[0];
    	infoarray[3] = orientation[1];
    	infoarray[4] = orientation[2];
    	return infoarray;
    }


}

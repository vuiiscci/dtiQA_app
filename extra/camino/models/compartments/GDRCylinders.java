package models.compartments;

import numerics.GammaFunctions;
import numerics.GammaRandom;
import numerics.RealMatrix;
import tools.CL_Initializer;
import misc.LoggedException;
import models.ParametricModel;
import imaging.DW_Scheme;
import imaging.StejskalTannerScheme;

public class GDRCylinders extends ListDistributionCompartment {

	/** gyromagnetic ratio */
	private final double GAMMA = DW_Scheme.GAMMA;

	/** gamma shape parameter */
	private static double gamma_k;
	
	/** gamma scale parameter */
	private static double gamma_beta;

	/** orientation of cylinders, theta */
	private static double theta;
	
	/** orientation of cylinders, phi */
	private static double phi;
	
	/** diffusivity */
	private static double d= CL_Initializer.DIFF_CONST; 
	
	/** number of bins in gamma distrtibution */
	private final int numBins;
	

	
	
	/** 
	 * constructor with nothing
	 */
	public GDRCylinders(){
		super(CompartmentType.GAMMADISTRIBRADIICYLINDERS.numParams, new CylinderGPD());
		
		this.numBins=5;
	}
	
	
	
	/** 
	 * constructor with specified number of bins
	 * 
	 * @param bins number of bins
	 */
	public GDRCylinders(double k, double beta, double theta, double phi, double d, int numBins){
		// construct params list and specify CylinderGPC compartment
		// in superclass.
		super(getGammaDistrnListParams(numBins), "CylinderGPC", 5);
		
		GDRCylinders.gamma_k= k;
		GDRCylinders.gamma_beta= beta;
		GDRCylinders.theta= theta;
		GDRCylinders.phi= phi;
		GDRCylinders.d=d;
		
		this.numBins=numBins;
	}
	
	/**
	 * constructor with default number of bins
	 */
	public GDRCylinders(double k, double beta, double theta, double phi, double d){
		this(k, beta, theta, phi, d, 5);
	}
	
	/**
	 * translate gamma parameters and number of bins into 
	 * a list of cylinder parameters and weights. This finds a range of values
	 * for the distribution but searching downwards and upwards from the 
	 * mode (or mean if k<=1). It then divides this range into bins,
	 * calculates a weight from the gamma pdf and assembles the list
	 * of cylinder params for the list distrn to evaluate.
	 * 
	 * @param numBins number of bins in gamma distrn
	 * 
	 * @return array of parameters arrays for the 
	 */
	private static double[][] getGammaDistrnListParams(int numBins){
		
		
		
		
		double lower;		  // lower limit to find
		double upper;         // upper limit to find
		
		
		
		
//		 search downwards from the mode to find the
//		 point at which the Cumulative distribution fn
//		 falls below 1/num bins.
//		while(gammaCDF(gamma_k, gamma_beta, current)>1.0/numBins){
//			
//			current-= searchInc;
//			
//			if(current<=0.0){
//				current=0.0;
//				break;
//			}
//			
//		}
//		lower= current;
		
		//find upper limit 1-1.0/numBins using brents algorithm
		lower=findGammaCDFCrossing(0, gamma_beta*gamma_k, 1.0/numBins, 1e-20);
		
	
		
		
		
		
		// reset search start pos
//		   current= searchStart;
		
		// search upwards from the start point to find the 
		// point at which the cumulative distribution fn
		// climbs about 1-1/num bins
//		int it_count=0;		
//		while(gammaCDF(gamma_k, gamma_beta, current)<(1.0-1.0/numBins)){			
//			it_count++;
//			current+=  searchInc;
//		}
		//upper= current;
		
		//find upper limit 1-1.0/numBins using brents algorithm
		upper=findGammaCDFCrossing(lower, numBins*gamma_beta*gamma_k, (1-1.0/numBins), 1e-20);			
		
		
		
		
		
		
		
	
		double range= (upper-lower);
		double binWidth= range/numBins;
		
		
		double[][] params= new double[numBins][CompartmentType.CYLINDERGPD.numParams+1];
		
				for(int i=0; i<numBins; i++){
			
			double r= lower + (i+0.5)*binWidth;
			double bottom= lower + i*binWidth;
			double top= lower + (i+1)*binWidth;
			
//			double f= gammaPDF(gamma_k, gamma_beta, bottom, top)*binWidth;
		
			//f is the cumulative PDF of gamma between bottom and top (normalized to the bining range of the CDF)   
			double f= (gammaCDF(gamma_k, gamma_beta, top)-gammaCDF(gamma_k, gamma_beta, bottom))/(1-2*(1.0/numBins));
			
			params[i][0]= f;				// weight
			params[i][1]= d;				// diffusivity
			params[i][2]= theta; 			// theta
			params[i][3]= phi; 				// phi
			params[i][4]= r;				// radius			
			
			
		}
		
		return params;
		
	}
	
	
	
	/**
	 * returns the cumulative gamma distribution between the given limits.
	 * this is evaluated as
	 *     F(xmin, xmax; k, theta) = \frac{1}{\Gamma(k)}(\gamma(k, xmax) - \gamma(k, xmin))
	 *  where
	 *    \Gamma(k) is the gamma function and
	 *    \gamma(k, x) is the incomplete gamma function
	 *    
	 * @see http://en.wikipedia.org/wiki/Gamma_distribution
	 * 
	 * @param k gamma shape param
	 * @param theta gamma scale param
	 * @param xmin lower limit of parameter range
	 * @param xmax upper limit of parameter range
	 * 
	 * @return cumulative gamma dist at xmax minus same at xmin
	 */
	private static double gammaPDF(double k, double theta, double xmin, double xmax){
	    
	    double Gammak= GammaFunctions.gammln(k);
	    
	    double gammaXmin= GammaFunctions.gammp(k, xmin/theta);
	    
	    double gammaXmax= GammaFunctions.gammp(k, xmax/theta);
	    
	    return (gammaXmax - gammaXmin)/Math.exp(Gammak)/(xmax-xmin);
		
	}
	
	
	/**
	 * calculates the cumulative gamma function up to the value given
	 * 
	 * @param k gamma shape param
	 * @param theta gamma scale param
	 * @param x top end upper limit of integral
	 * 
	 * @return gamma(x/theta)/Gamma(k)
	 * gamma(k, z)= incomplete gamma fn
	 * Gamma(k)= gamma function
	 * 
	 * 
	 */
	private static double gammaCDF(double k, double thet, double x){
	    
	    double gammaX= GammaFunctions.gammp(k, x/thet);
	    return gammaX;
	}

	
	
	
	/**
	 * generates signals from this compartment for each line
	 * in the scheme given.
	 * 
	 * @param scheme scan specifics
	 */
	public RealMatrix getSignals(double[] params, DW_Scheme rawScheme) {

		RealMatrix signals = new RealMatrix(rawScheme.numMeasurements(),1);
		
		//Initialise class variables from model parameters
		GDRCylinders.gamma_k= params[0];
		GDRCylinders.gamma_beta= params[1];
		GDRCylinders.d=params[2];
		GDRCylinders.theta= params[3];
		GDRCylinders.phi= params[4];
		

		StejskalTannerScheme scheme;

		try {
			scheme = (StejskalTannerScheme) rawScheme;
		} catch (ClassCastException cce) {
			throw new LoggedException(
					"scheme object passed to cylinder compartment is not a StejskalTanner sequence");
		}

        double[][] subParams= getGammaDistrnListParams(numBins);

		for (int m=0; m < rawScheme.numMeasurements();m++)
		{
			signals.setEntry(m, 0 , this.getWeightedCylinderSignal(subParams, scheme, m));
		}
		
        return signals;	    
	}
	
	
	
	public double getSignal(double[] params, DW_Scheme rawScheme, int i){
        
		//Initialise class variables from model parameters
		GDRCylinders.gamma_k= params[0];
		GDRCylinders.gamma_beta= params[1];
		GDRCylinders.d=params[2];
		GDRCylinders.theta= params[3];
		GDRCylinders.phi= params[4];
		

		StejskalTannerScheme scheme;

		try {
			scheme = (StejskalTannerScheme) rawScheme;
		} catch (ClassCastException cce) {
			throw new LoggedException(
					"scheme object passed to cylinder compartment is not a StejskalTanner sequence");
		}

		double[][] subParams= getGammaDistrnListParams(numBins);
        		
		
		return getWeightedCylinderSignal(subParams, scheme, i);
	}
	
	private double getWeightedCylinderSignal(double[][] subParams, StejskalTannerScheme scheme, int i){
		
		
		double signal = 0.0;
		final double[] cylParams= new double[subParams[0].length];
        
		
		for(int ci=0; ci<subParams.length; ci++){
			
        	// construct params array for the cylinder
        	for(int cj=1; cj<subParams[0].length; cj++){
        		cylParams[cj-1]= subParams[ci][cj];
        	}
        	
        	// add weighted signal to accumulated signal
        	double contrib= model.getSignal(cylParams, scheme, i);
        	
        	contrib=contrib*subParams[ci][0];
        	
        	signal+=contrib;
        	
        }
		
		
        
        return signal; 
	}
	
	
	//Using Brent root finding to determine cdfs
	private static double findGammaCDFCrossing(double startx, double stopx, double offset, double convergence)
	{
		double fstartx=gammaCDF(gamma_k, gamma_beta, startx) - offset;
		double fstopx=gammaCDF(gamma_k, gamma_beta, stopx) - offset;
		double delta=Math.abs(stopx-startx);
		
		
		if (fstartx * fstopx > 0)
		{
			//try to fix start and stop points
			
			//if lower bound is too high, start from 0 again
			if (fstartx>0)
			{
				fstartx=gammaCDF(gamma_k, gamma_beta, 0) - offset;
			}
			
			//if upper bound is too low, ignore kappa
			else if (fstopx<0)
			{
				fstopx=gammaCDF(gamma_k, gamma_beta, stopx/gamma_k) - offset;
			}
			
			else //there is something wrong with this gamma distribution
			{
				throw new LoggedException("fa="+fstartx+" and fb="+fstopx+"same sign! Can't find root for GDR.\n+" +
					"Offending parameters were: k="+gamma_k+" beta="+ gamma_beta);
			}
		}
		
		
		
		double root=startx;
		double froot=fstartx;
		
		boolean mflag=true;
		
		
		double s = 0;
		double de = 0;
		
		while (!(delta<convergence || fstartx==0 || fstopx==0)){
			
			 //System.out.println("a="+startx+" fa="+fstartx+" b="+stopx+" fb="+fstopx);
			
			 if (fstartx != froot && fstopx != froot)
			 {
				 //inverse interpolation
				 s=startx*fstopx*froot/((fstartx-fstopx)*(fstartx-froot));
				 s+=stopx*fstartx*froot/((fstopx-fstartx)*(fstopx-froot));
				 s+=root*fstartx*fstopx/((froot-fstartx)*(froot-stopx));
			 }
			 else
			 {
				 //secant method
				 s=stopx-fstopx*(stopx-startx)/(fstopx-fstartx);
			 }
			 
			 boolean condition1=!((s>=(3*startx+stopx)/4) && (s<=stopx) || (s<=(3*startx+stopx)/4) && (s>stopx));
			 boolean condition2=mflag && (Math.abs(s-stopx) >= Math.abs(stopx-root)/2);
			 boolean condition3=!mflag && (Math.abs(s-stopx) >= Math.abs(root-de)/2);
			 boolean condition4=mflag && (Math.abs(stopx-root) < delta);
			 boolean condition5=!mflag && (Math.abs(root-delta) < delta);//??
				
			 //bisection
			 if (condition1 || condition2 ||  condition3 || condition4 || condition5)
			 {
				 s=(startx+stopx)/2;
				 mflag=true;
			 }
			 else
			 {
				 mflag=false;
			 }
			 double fs=gammaCDF(gamma_k, gamma_beta, s) - offset;
			 
			 //update start,stop,root,d
			 
			 de=root;
			 root=stopx;
			 froot=fstopx;
			 
			 if ((fstartx * fs) < 0)
			 {
				 stopx=s;
				 fstopx=fs;
				 
			 }
			 else
			 {
				 startx = s;
				 fstartx = fs;
			 }
			 
			 if (Math.abs(fstartx) < Math.abs(fstopx))
			 {
				 //swap startx and stopx
				 double tmp=stopx;
				 double ftmp=fstopx;
				 stopx=startx;
				 fstopx=fstartx;
				 startx=tmp;
				 fstartx=ftmp;
			 }
			 
			 
			 //calculate new iteration fs
			 //fstartx=gammaCDF(gamma_k, gamma_beta, startx) - offset;
			 //fstopx=gammaCDF(gamma_k, gamma_beta, stopx) - offset;
			 delta=Math.abs(stopx-startx);
		}
		return s;		
	}

}
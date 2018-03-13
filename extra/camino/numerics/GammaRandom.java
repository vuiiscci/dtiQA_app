package numerics;

import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Logger;


/**
 * generates gamma-distributed random numbers using a
 * mersienne twister and a rejection-sampling technique
 *
 *
 * @author matt (m.hall@cs.ucl.ac.uk)
 *
 */
public class GammaRandom {

	/** logging object */
	private final Logger logger= Logger.getLogger(this.getClass().getName());
	
	/** twister, for uniform distributed random numbers */
	private final MTRandom rng;
	
	/** shape parameter of distribution */
	private final double k;
	
	/** inverse scale parameter */
	private final double beta;
	
	/** integer part of shape parameter */
	private final int kFloor;
	
	/** non-integer part of k */
	private final double delta;
	
	/** rejection sampling parameter for Best's algorithm */
	private final double b;
	
	/** rejection sampling parameter for Best's algorithm */	
	private final double c;
	
	/**
	 * constructor. instantiates rng and sets distribution parameters.
	 * 
	 * @param seed rng seed
	 * @param k distribution shape parameter
	 * @param beta distribution inverse-scale parameter (beta=1/theta) 
	 */
	public GammaRandom(long seed, double k, double beta){
		this.rng= new MTRandom(seed);
		this.k=k;
		this.kFloor=(int)Math.floor(k);
		this.delta=k-kFloor;
		this.beta=beta;
		
		this.b= k-1;
		this.c= 3*k-0.75;	
	}
	
	
	/**
	 * returns a sample drawn from the gamma distribution given by
	 * k and beta. This uses Best's algorithm. 
	 * 
	 * @see Devroye, Non-uniform random variate generation, Springer-Verlag, New York (1986) Chapter 9, pg 410 
	 * 
	 * This is a rejection-sampling method based on cubic transforms. It has uniform performance across all 
	 * parameter values and shows superior convergence to the method used previously. 
	 * 
	 * @return gamma distributed random number
	 */
	public double nextGamma(){
		
		while(true){
			
			double U=rng.nextDouble();
			double V=rng.nextDouble();
			
			double W=U*(1-U);
			double Y=Math.sqrt(c/W)*(U-0.5);
			double X=b+Y;
			
			if(X>=0.0){
				
				double Z= 64*(W*W*W)*(V*V);
				
				if(Z<=1.0-2.0*Y*Y/X){
					return beta*X;
				}
				else{
					if(Math.log(Z)<=2*(b*Math.log(X/b)-Y)){
						return beta*X;
					}
				}
			}
		}
	}
	
	
	/**
	 * generates a uniformly distrubuted random number from
	 * the interval [0,1). This is basically the same as the 
	 * rng.nextDouble() method, but rejects any zeros that 
	 * crop up.
	 * 
	 * @return value uniformly distributed on [0,1)
	 */
	private double getUniformNonZeroValue(){
		double U=0.0;
		
		while(U==0.0){
			U=rng.nextDouble();
		}
		
		return U;
	}
	
	
	/**
	 * returns the weight for the specified bin in the specified distribution,
	 * weighted by bin width.
	 * 
	 * @param k shape parameter
	 * @param beta scale parameter
	 * @param xmin lower end of bin
	 * @param xmax upper end of bin
	 * 
	 * @return integral over range 
	 */
	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		int samples=100000;
		
		int bins=1000;
		double[] distrn= new double[bins];
		
		double min=Double.MAX_VALUE;
		double max=0.0;
		
		//double binsize=(max-min)/bins;
				
		
		GammaRandom grng= new GammaRandom(34756037485l, 2.331206, 6.909243E-7);
		
		System.err.print("generating and binning "+samples+" numbers... ");
		
		for(int i=0; i<samples; i++){
			double num= grng.nextGamma();
			
			if(num<min){
				min=num;
			}
			if(num>max){
				max=num;
			}
		}

		double binsize=(max-min)/bins;
		for(int i=0; i<samples; i++){
						
			double num= grng.nextGamma();
			
			int bin=Math.min((int)(num/binsize), bins-1);
			
			if(bin==-1){
				System.err.println("i="+i+" bin is -1. that shouldn't happen...");
			}
			
			distrn[bin]++;
		}

		
		
		
		System.err.println("done");
		
		System.err.print("outputting distributiion... ");
		
		FileWriter writer;
		
		try{
			writer= new FileWriter("gammaDist.csv");
			
			for(int i=0; i<bins; i++){
				writer.write((i*binsize)+","+distrn[i]+"\n");
			}
			
			writer.flush();
			writer.close();
		}
		catch(IOException ioe){
			throw new RuntimeException(ioe);
		}
		
		System.err.println("done");

	}

}

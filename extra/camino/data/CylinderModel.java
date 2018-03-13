package data;

import java.util.logging.Logger;

import misc.DT;
import misc.LoggedException;

import data.BallStick;
import data.ModelPDF;

import tools.CL_Initializer;
import imaging.DW_Scheme;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd> Simple model for the diffusion displacement density function in white matter with two
 * compartments. The model assumes parallel cylindrical axon cells with equal
 * radii and impermeable walls.
 * 
 * <dt>Description: Generates synthetic data analytically.
 * 
 * <dd>Input is a scheme file and some imaging parameters and output is a set
 * of measurements.
 * 
 * </dl>
 * 
 * @author Eleftheria Panagiotaki email:panagio@cs.ucl.ac.uk
 */
public class CylinderModel implements ModelPDF {
	
	/** gyromagnetic ratio */
	double GAMMA = DW_Scheme.GAMMA;
	/** Volume fraction */
	private final double f ;
	/** parallel diffusivity */
	private final double dparal;
	/** perpendicular diffusivity */
	private final double dperp ;
	/** Cylinder radius */
	private final double R ;
	/** fibre direction n */
	private final double[] n;

	
	/**
	 * constructor for volume fraction, specified diffusivities, axon radius and fibre
	 * direction.
	 */
	public CylinderModel(double f, double dparal, double dperp,double R,double[] n){
		this.f=f;
		this.dparal=dparal;
		this.dperp=dperp;
		this.R=R;
		this.n=new double [3];
		for(int i=0 ;i<3; i++){
			this.n[i] = n[i];
		}
		
		
	}
	

	public double signalAt(double[] G, double DELTA,double delta) {

		double modG=0.0;
		
		for(int i=0; i<G.length; i++){
			modG += G[i]*G[i];
		}
		
		modG= Math.sqrt(modG);



		/** 60 first roots from the equation j'1(am*x)=0 */
		double[] am = { 1.84118307861360, 5.33144196877749, 8.53631578218074,
				11.7060038949077, 14.8635881488839, 18.0155278304879,
				21.1643671187891, 24.3113254834588, 27.4570501848623,
				30.6019229722078, 33.7461812269726, 36.8899866873805,
				40.0334439409610, 43.1766274212415, 46.3195966792621,
				49.4623908440429, 52.6050411092602, 55.7475709551533,
				58.8900018651876, 62.0323477967829, 65.1746202084584,
				68.3168306640438, 71.4589869258787, 74.6010956133729,
				77.7431620631416, 80.8851921057280, 84.0271895462953,
				87.1691575709855, 90.3110993488875, 93.4530179063458,
				96.5949155953313, 99.7367932203820, 102.878653768715,
				106.020498619541, 109.162329055405, 112.304145672561,
				115.445950418834, 118.587744574512, 121.729527118091,
				124.871300497614, 128.013065217171, 131.154821965250,
				134.296570328107, 137.438311926144, 140.580047659913,
				143.721775748727, 146.863498476739, 150.005215971725,
				153.146928691331, 156.288635801966, 159.430338769213,
				162.572038308643, 165.713732347338, 168.855423073845,
				171.997111729391, 175.138794734935, 178.280475036977,
				181.422152668422, 184.563828222242, 187.705499575101 };

		double[] am1 = new double[am.length];
		for (int i1 = 0; i1 < am.length; i1++) {
			am1[i1] = am[i1] / R;
			// System.err.println("am1["+i+"] = "+ am1[i]);
		}

		/** tau = diffusion time */
		double tau = DELTA - (delta / 3);

		/** modulo of n */
		double modn = 0.0;
		for (int i1 = 0; i1 < n.length; i1++) {
			modn += n[i1] * n[i1];
		}
		modn = Math.sqrt(modn);

		/** Diffusion tensor, cylindrically symmetric */
		// constructing the compartments of the tensor
		double dcomp1 = (dparal - dperp);
		double[] dnnT = new double[6];

		dnnT[0] = dcomp1 * n[0] * n[0];
		dnnT[1] = dcomp1 * n[0] * n[1];
		dnnT[2] = dcomp1 * n[0] * n[2];
		dnnT[3] = dcomp1 * n[1] * n[1];
		dnnT[4] = dcomp1 * n[1] * n[2];
		dnnT[5] = dcomp1 * n[2] * n[2];

		double[] e1 = new double[] { 1.0, 0.0, 0.0 };

		double[] dperpI = new double[6];

		dperpI[0] = dperp * e1[0];
		dperpI[1] = dperp * e1[1];
		dperpI[2] = dperp * e1[1];
		dperpI[3] = dperp * e1[0];
		dperpI[4] = dperp * e1[1];
		dperpI[5] = dperp * e1[0];

		double[] Dcomparts = new double[6];

		Dcomparts[0] = dnnT[0] + dperpI[0];
		Dcomparts[1] = dnnT[1] + dperpI[1];
		Dcomparts[2] = dnnT[2] + dperpI[2];
		Dcomparts[3] = dnnT[3] + dperpI[3];
		Dcomparts[4] = dnnT[4] + dperpI[4];
		Dcomparts[5] = dnnT[5] + dperpI[5];

		DT d = new DT(Dcomparts[0], Dcomparts[1], Dcomparts[2], Dcomparts[3],
				Dcomparts[4], Dcomparts[5]);

		/** calculating the sum for the perpendicular intra-cellular signal */
		double sum = 0;
		for (int i1 = 0; i1 < am1.length; i1++) {
			// d*am^2
			double dam = dparal * am1[i1] * am1[i1];
			// -d*am^2*delta
			double e11 = -dam * delta;
			// -d*am^2*DELTA
			double e2 = -dam * DELTA;
			// -d*am^2*(DELTA-delta)
			double dif = DELTA - delta;
			double e3 = -dam * dif;
			// -d*am^2*(DELTA+delta)
			double plus = DELTA + delta;
			double e4 = -dam * plus;
			// nominator of the fraction
			double nom = 2 * dam * delta - 2 + 2 * Math.exp(e11) + 2
					* Math.exp(e2) - Math.exp(e3) - Math.exp(e4);
			// denominator
			double denom = dam * dam * am1[i1] * am1[i1]
					* (R * R * am1[i1] * am1[i1] - 1);
			// the sum of the fraction
			sum += (nom / denom);

		}

		/**Normalised signal*/
		double normSignal = 0.0;
		/** theta is the angle between n and G */
		double theta = 0.0;
		double dotprodGn = 0.0;
		/** q the wavenumber */
		double[] q = new double[3];
		double qDq = 0.0;

		dotprodGn = n[0] * G[0] + n[1] * G[1] + n[2] * G[2];
		theta = Math.acos(dotprodGn / modG * modn);

		q[0] = GAMMA * delta * G[0];
		q[1] = GAMMA * delta * G[1];
		q[2] = GAMMA * delta * G[2];

		/** Contract the diffusion tensor by a vector */
		qDq = d.contractBy(q);

		/** Hindered diffusion in the extra-cellular component */
		double extracelSignal = Math.exp(-tau * qDq);

		/** Intra-cellular parallel signal */
		double intraPar1 = (GAMMA * delta) * modG * Math.cos(theta);
		double intraPar2 = intraPar1 * intraPar1 * dparal;
		double intraParallel = Math.exp(-tau * intraPar2);

		/** -2g^2*Gperp^2 */
		double gg = -2 * GAMMA * GAMMA * modG * modG * Math.sin(theta)
				* Math.sin(theta);

	
		double intraPerpend = Math.exp(gg * sum);


		double intracelSignal = intraParallel * intraPerpend;

		normSignal = f * intracelSignal + (1 - f) * extracelSignal;
		

		return normSignal;

	}
    public double ftAt(double[] q, double tau) {                                          
  	 
    	throw new LoggedException("CylinderModel not defined given only q vector and tau. requires more detailed scan parameters");
    	
    }
	
	public double ftAtB_Vec(double[] g, double b){
	    
	    throw new LoggedException("CylinderModel not defined given only g and b.");
	}

	public double at(double[] x, double tau) {                                                  

		throw new LoggedException("CylinderModel not defined given only x vector and tau. requires more detailed scan parameters");
    }
	

	
    public double[][] getPDs() {

    	double[][] pds= new double[1][3];
    	
    	for(int i=0; i<3; i++){
    		pds[0][i]=n[i];
    	}
    	
        return pds;
    }
	
}




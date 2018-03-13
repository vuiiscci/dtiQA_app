package numerics;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Contains implementation of the various Bessel functions.
 * 
 * <dt>Description:
 * 
 * <dd>Adapted from NRC++ (2002) code.
 * 
 * </dl>
 * 
 * @author Andy Sweet
 *  
 */

public class BesselFunctions {
	
	// accuracy of besselIn
	private static final double ACC = 40.0;
	// used to prevent overflow in besselIn
	private static final double BIGNO = 1.0e10;
	private static final double BIGNI = 1.0e-10;
	
	/**
	 * Zeroth-order Bessel function of the first kind. See NRC++ (2002), pg 237.
	 *
	 * @param x value to evaluate function at
	 */
	public static double besselJ0(double x) {
		
		double ax,z,xx,y,ans,ans1,ans2;
		
		// direct rational approximation
		if ((ax=Math.abs(x)) < 8.0) {
			y=x*x;
			
			ans1 = 57568490574.0 + y*(-13362590354.0 + y*(651619640.7 + y*(-11214424.18 + y*(77392.33017 + y*(-184.9052456)))));
				
			ans2 = 57568490411.0 + y*(1029532985.0 + y*(9494680.718 + y*(59272.64853 + y*(267.8532712 + y*1.0))));
			
			ans = ans1/ans2;
		} else {
		// fitting function
			z = 8.0/ax;
			y = z*z;
			xx = ax -0.785398164;
			
			ans1 = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4 + y*(-0.2073370639e-5 + y*0.2093887211e-6)));
			
			ans2 = -0.1562499995e-1 + y*(0.1430488765e-3 + y*(-0.6911147651e-5 + y*(0.7621095161e-6 - y*(0.934945152e-7))));
			
			ans = Math.sqrt(0.636619772/ax)*(Math.cos(xx)*ans1-z*Math.sin(xx)*ans2);
		}
		return ans;
	}
	
	/**
	 * First-order Bessel function of the first kind. See NRC++ (2002), pg 238.
	 *
	 * @param x value to evaluate function at
	 */
	public static double besselJ1(double x) {
		
		double ax,z,xx,y,ans,ans1,ans2;
		ax = Math.abs(x);
		
		// direct rational approximation
		if (ax < 8.0) {
			y = x*x;
			
			ans1 = x*(72362614232.0 + y*(-7895059235.0 + y*(242396853.1 + y*(-2972611.439 + y*(15704.48260 + y*(-30.16036606))))));
				
			ans2 = 144725228442.0 + y*(2300535178.0 + y*(18583304.74 + y*(99447.43394 + y*(376.9991397 + y*1.0))));
			
			ans = ans1/ans2;
		} else {
		// fitting function
			z = 8.0/ax;
			y = z*z;
			xx = ax - 2.356194491;
			
			ans1 = 1.0 + y*(-0.183105e-2 + y*(-0.3516396496e-4 + y*(-0.2457520174e-5 + y*(-0.240337019e-6))));
			
			ans2 = 0.04687499995 + y*(-0.2002690873e-3 + y*(0.8449199096e-5 + y*(-0.88228987e-6 - y*(0.105787412e-6))));
			
			ans = Math.sqrt(0.636619772/ax)*(Math.cos(xx)*ans1-z*Math.sin(xx)*ans2);
			
			if(x < 0.0) ans = -ans;
		}
		return ans;
	}
	
	/**
	 * Zeroth-order modified Bessel function of the first kind. See NRC++ (2002), pg 241-242.
	 *
	 * @param x value to evaluate function at
	 */
	public static double besselI0(double x) {
		double ax, ans, y;
		ax = Math.abs(x);
		
		// polynomial fit
		if (ax < 3.75) {
			y = x/3.75;
			y *= y;
			
			ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
			
		} else {
			y = 3.75/ax;
			
			ans = (Math.exp(ax)/Math.sqrt(ax))*0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2 + y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1 + y*(0.392377e-2))))))));
		}
		return ans;
	} 
	
	/**
	 * First-order modified Bessel function of the first kind. See NRC++ (2002), pg 241-242.
	 *
	 * @param x value to evaluate function at
	 */
	public static double besselI1(double x) {
		double ax, ans, y;
		ax = Math.abs(x);
		
		// polynomial fit
		if (ax < 3.75) {
			y = x/3.75;
			y *= y;
			
			ans = ax*(0.5 + y*(0.87890594 + y*(0.51498869 + y*(0.15084934 + y*(0.2658733e-1 + y*(0.301532e-2 + y*0.32411e-2))))));
			
		} else {
			y = 3.75/ax;
			
			ans = 0.2282967e-1 + y*(-0.2895312e-1 + y*(0.1787654e-1 - y*0.420059e-2));
			
			ans = 0.39894228 + y*(-0.3988024e-1 + y*(-0.362018e-2 + y*(0.163801e-2 + y*(-0.1031555e-1 + y*ans))));
			
			ans *= Math.exp(ax)/Math.sqrt(ax);
		}
		return ans;
	}
	
	/**
	 * nth-order modified Bessel function of the first kind. See NRC++ (2002), pg 242-243??
	 *
	 * @param x value to evaluate function at
	 */
	public static double besselIn(int n, double x) {
		
		int j;
		double bi,bim,bip,tox,ans;
		
		// functions are defined for n=0,1 already
		if (n == 0) {
			return besselI0(x);
		} else if (n == 1) { 
			return besselI1(x);
		} 
		
		// no point computing for 0
		if (x == 0.0) { 
			return  0.0;	
		} else {
			tox = 2.0/Math.abs(x);
			bip = ans = 0.0;
			bi = 1.0;
			
			// downward sum of series approximation
			for (j=2*(n + (int)Math.sqrt(ACC*n)); j > 0; j--) {
				bim = bip+j*tox*bi;
				bip = bi;
				bi = bim;
				// renormalize to prevent overflows
				if (Math.abs(bi) > BIGNO) {
					ans *= BIGNI;
					bi *= BIGNI;
					bip *= BIGNI;
				}
				
				if (j == n)	ans = bip;
			}
			
			// normalize with bessI0 and flip sign for negative x
			ans *= besselI0(x)/bi;
			if(x < 0.0) ans = -ans;
			return ans;
		}
	}
	
	// public static void main(String[] args) {
	// 	
	// 	// try some values of different order and compare them with MATLAB
	// 	double[] testVals = {0.1, 0.5, 1.7, 23.5, 100.1, 2000}; 
	// 	for (int i=0; i<testVals.length; i++) {
	// 		System.out.println("BesselI0("+ testVals[i] + "): " + BesselFunctions.besselI0(testVals[i]));
	// 	} 
	// 	for (int i=0; i<testVals.length; i++) {
	// 		System.out.println("BesselI1("+ testVals[i] + "): " + BesselFunctions.besselI1(testVals[i]));
	// 	}
	// 	for (int i=0; i<testVals.length; i++) {
	// 		System.out.println("BesselI2("+ testVals[i] + "): " + BesselFunctions.besselIn(2, testVals[i]));
	// 	}
	// }
}
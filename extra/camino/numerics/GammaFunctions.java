package numerics;

import misc.LoggedException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Contains implementation of the gamma function.
 * 
 * <dt>Description:
 * 
 * <dd>Adapted from NRC code.
 * 
 * </dl>
 * 
 * @author Danny Alexander $Id: GammaFunctions.java,v 1.3 2005/08/18 11:12:21
 *         ucacmgh Exp $
 *  
 */
public class GammaFunctions {

	/**
	 * k proxy
	 */
	
	/**
	 * gammln(k) proxy
	 */
	
    /**
     * Euler's constant. (NRC p. 223)
     */
    public static final double EULERGAMMA = 0.5772156649;

    /**
     * max number of iterations in series evaluation
     */
    private static final int ITMAX=100;
    
    /**
     * small number constant used in gamma series evaluation
     */
    private static final double EPS=3E-7;
    
    /**
     * small number constant used in continued fraction gamma evaluation
     */
    private static final double FPMIN= 1E-30;
    
    
    
    /**
     * The gamma function. See NRC p. 214.
     * 
     * @param xx
     *            Argument at which to evaluate the gamma function.
     */
    public static double gammln(double xx) {
        double x, y, tmp, ser;
        double[] cof = { 76.18009172947146, -86.50532032941677, 24.01409824083091,
                -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
        int j;

        y = x = xx;
        tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        ser = 1.000000000190015;
        for (j = 0; j <= 5; j++)
            ser += cof[j] / ++y;
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }
    
    /**
     * Returns the incomplete gamma function P(a; x).
     * see NRC p. 218.
     * @param a
     * @param x
     * @return
     */
    public static double gammp(double a, double x){
        
        if(x<0.0 || a<= 0.0){
            throw new LoggedException("Invalid args in incomplete gamma function");
        }
        
        if(x<(a+1.0)){
            return gser(a, x);
        }
        else{
            return 1.0-gcf(a, x);
        }
    }
    
    
    
    
    /**
     * Returns the incomplete gamma function P(a; x) evaluated by its
     * series representation as gamser.
     */
    private static double gser(double a, double x){
        int n;
        double sum,del,ap;
        double gln=gammln(a);
        if (x <= 0.0) {
            if (x < 0.0) 
                throw new LoggedException("x less than 0 in routine gser");
            
            return 0.0;
        } 
        else{
            ap=a;
            del=sum=1.0/a;
            for (n=1;n<=ITMAX;n++) {
                ++ap;
                del *= x/ap;
                sum += del;
                if (Math.abs(del) < Math.abs(sum)*EPS) {
                    return sum*Math.exp(-x+a*Math.log(x)-(gln));
                }
            }
            throw new LoggedException("a too large, ITMAX too small in routine gser");
        }
    }

    
    /**
     * Returns the incomplete gamma function Q(a; x) evaluated by its continued 
     * fraction representation.
     * 
     * @param a
     * @param x
     * @return
     */
    private static double gcf(double a, double x){
        int i;
        double an,b,c,d,del,h;
        double gln=gammln(a);

        //Set up for evaluating continued fraction by modified Lentz's method (x5.2) with b0 = 0.
        b=x+1.0-a; 
        c=1.0/FPMIN;
        d=1.0/b;
        h=d;
        for (i=1;i<=ITMAX;i++) { //Iterate to convergence.
            an = -i*(i-a);
            b += 2.0;
            d=an*d+b;
            if (Math.abs(d) < FPMIN) 
                d=FPMIN;

            c=b+an/c;

            if (Math.abs(c) < FPMIN) 
                c=FPMIN;

            d=1.0/d;
            del=d*c;
            h *= del;
            if (Math.abs(del-1.0) < EPS)
                break;
        }
        if (i > ITMAX)
            throw new LoggedException("a too large, ITMAX too small in gcf");

        return Math.exp(-x+a*Math.log(x)-(gln))*h; //Put factors in front.
    }

}
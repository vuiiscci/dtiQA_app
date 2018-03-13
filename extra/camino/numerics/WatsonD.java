package numerics;

import numerics.WatsonDistribution;


/**
 * <dl>
 * <dt> Description: Evaluates D_3(kappa), eqn 10.3.32 in Mardia and Jupp "Directional Statistics" (2000)
 * </dl>
 *
 * @version $Id$
 * @author  Philip Cook
 * @see sphereDistFit.WatsonFitter
 *
 * 
 */
public class WatsonD extends SVAnalyticalFunction {  

    public WatsonD() {
	super(true);
    }

   
    /**
     * @return the value of the function for the given parameter.
     */ 
    public double evaluate(double kappa) {

	if (kappa < 700.0) {
	    return 
		WatsonDistribution.hyper1F1(1.5, 2.5, kappa, 1.0e-10) 
		/ 
		( 3.0 * WatsonDistribution.hyper1F1(0.5, 1.5, kappa, 1.0e-10) );
	}
	else {
	    return Math.exp( WatsonDistribution.logHyper1F1(1.5, 2.5, kappa, 1.0e-10) 
			     - Math.log(3.0) - WatsonDistribution.logHyper1F1(0.5, 1.5, kappa, 1.0e-10) );
	}
    }
    

    /**
     * @return the derivative of this function. 
     *
     */
    public SVAnalyticalFunction differentiated() {

	return new SVAnalyticalFunction(false) {
		
		public double evaluate(double kappa) {
		    
		    if (kappa < 700.0) {
			double a = WatsonDistribution.hyper1F1(0.5, 1.5, kappa, 1.0e-10);
			double b = WatsonDistribution.hyper1F1(1.5, 2.5, kappa, 1.0e-10);
			double c = WatsonDistribution.hyper1F1(2.5, 3.5, kappa, 1.0e-10);
			
			// do the multiplication to avoid huge numbers
			return (3.0 / 15.0) * (c / a) -  (1.0 / 9.0) * (b / a) * (b / a);
		    }
		    else {
			double a = WatsonDistribution.logHyper1F1(0.5, 1.5, kappa, 1.0e-10);
			double b = WatsonDistribution.logHyper1F1(1.5, 2.5, kappa, 1.0e-10);
			double c = WatsonDistribution.logHyper1F1(2.5, 3.5, kappa, 1.0e-10);

			//			return Math.log(3.0 / 15.0) + c - Math.log(1.0 / 9.0) - 2.0 * b + a;
			return (3.0 / 15.0) * Math.exp(c - a) - (1.0 / 9.0) * Math.exp(2.0 * (b - a));

		    }

		}
		
		public String toString() {
		    return "\\frac{3}{p^2 + 2p} " +
			"\\frac{M(\\frac{5}{2}, \\frac{p}{2} + 2, \\kappa)}" + 
			"{M(\\frac1{2}, \\frac{p}{2}, \\kappa)} " +
			"- \\left( \\frac1{p} \\frac{M(\\frac{3}{2}, \\frac{p}{2} + 1, \\kappa)}" +
			"{M(\\frac{1}{2}, \\frac{p}{2}, \\kappa)} \\right)^2";

		}

		
	    };
	
    }


    /**
     *
     * @return a LaTeX string describing this function.
     */
    public String toString() {
	return "D_3(\\kappa) = \\frac{M(\\frac{3}{2}, \\frac{p}{2} + 1, \\kappa)}{pM(\\frac1{2}, \\frac{p}{2}, \\kappa)}";
    }
    
    
}

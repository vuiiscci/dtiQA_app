package numerics;


import java.text.DecimalFormat;

/**
 * 
 * Translation of John Kent's code that implements his paper "Asymptotic Expansions for the Bingham Distribution" 
 * (Applied Statistics 36(2):139-144, 1987). According to John Kent, the C code from which this class is 
 * translated was in itself a translation of the Fortran code that he originally wrote.
 * 
 * @version $Id$
 * @author John T Kent
 * @author Philip Cook
 *
 */
public class BinghamFitter extends SphericalDistributionFitter {
 

    private static final int c__9 = 9;
    private static final int c__1 = 1;
    private static final int c__3 = 3;
    private static final int c__5 = 5;
    private static final double c_b28 = 8.;

    private BinghamFitter() {
    }

    
    /*    VERSION 1.1  4 NOVEMBER 1985                                      BIN000
	  10*/
    /*    COPYRIGHT J.T. KENT                                               BIN000
	  20*/
    /*    THIS IS A DRIVER PROGRAM TO CALCULATE THE MOMENTS OF THE          BIN000
	  30*/
    /*    BINGHAM DISTRIBUTION FROM PARAMETERS AK1, AK2,                    BIN000
	  40*/
    /*    OR                                                                BIN000
	  50*/
    /*    TO CALCULATE THE PARAMETERS CORRESPONDING TO GIVEN VALUES OF      BIN000
	  60*/
    /*    THE MOMENTS.                                                      BIN000
	  70*/
    /*    THIS PROGRAM IS DESIGNED TO RUN INTERACTIVELY.                    BIN000
	  80*/
    /*                                                                      BIN000
									    90*/

    /**
     * Get a Bingham distribution from a scatter matrix of samples.
     */
    public static BinghamDistribution getBinghamDistribution(EigenSystem3D tBarEig) throws ConvergenceException {
	return getBinghamDistribution(tBarEig, new java.util.Random());
    }


    /**
     * Get a Bingham distribution from a scatter matrix of samples, and also specify the random number generator used by the distribution.
     */
    public static BinghamDistribution getBinghamDistribution(EigenSystem3D tBarEig, java.util.Random ran) throws ConvergenceException {

	double t2 = tBarEig.eigenvalues[1];
	double t3 = tBarEig.eigenvalues[2];
	
	/* L1: */
	double[] akfc = bngpar(t2, t3);

	
	return new BinghamDistribution(tBarEig.eigenvectors, akfc[0], akfc[1],  akfc[2], ran);
    }

    


    /**
     *
     * Gets the Bingham parameters <code>k1, k2</code> (concentration parameters) and <code>bingC</code> (normalization constant).
     * <p> Specify the eigenvalues in descending order. t2 > t1 (numbering scheme of Watson).
     *
     * @return <code>{k1, k2, bingc}</code>.
     */
    public static double[] bngpar(double t2, double t1) throws ConvergenceException
    {
	double[] hi = new double[4];	/* was [2][2] */;

	/* System generated locals */
	double d__1 = 0.0;
	double d__2 = 0.0;

	double fc = 0.0;

	/* Local variables */
	int mode = 0;
	double rerr = 0.0;
	int i = 0;
	int[] level = new int[1]; // do this because level gets assigned in another method
	double fd1, fd2, bk1, bk2;
	double[] der = new double[2];

	double[] hes = new double[4];	/* was [2][2] */;

	double ak1 = 0.0;
	double ak2 = 0.0;
	
	/*    THIS ROUTINE COMPUTES THE VALUES OF THE PARAMETERS AK1            BI
	      N00880*/
	/*      AND AK2 WHICH CORRESPOND TO THE GIVEN MOMENTS T1 AND T2,        BI
		N00890*/
	/*      THE EXPECTED VALUES OF X1 AND  X2  FROM A BINGHAM DISTRIBUTION. BI
		N00900*/
	/*                                                                      BI
										N00910*/
	/*    INPUT:   T1,T2  MOMENTS FROM A BINGHAM DISTRIBUTION,              BI
	      N00920*/
	/*                    OR THE SMALLEST 2 EIGENVALUES OF THE SAMPLE       BI
			      N00930*/
	/*                    MOMENT OF INERTIA MATRIX FROM A SAMPLE OF SIZE N  BI
			      N00940*/
	/*                      0 < T1 <= T2,   T1 + 2*T2 <= 1                  BI
				N00950*/
	/*    OUTPUT:  AK1, AK2  VALUES OF THE CONCENTRATION PARAMETERS         BI
	      N00960*/
	/*                       CORRESPONDING TO T1 AND T2, OR                 BI
				 N00970*/
	/*                       THE MAXIMUM LIKELIHOOD ESTIMATES OF AK1 AND    BI
				 N00980*/
	/*                       AK2 FROM A SAMPLE.                             BI
				 N00990*/
	/*             HI(2,2)   N TIMES THE ASYMPTOTIC VARIANCE MATRIX OF THE  BI
		       N01000*/
	/*                       MLES OF AK1 AND AK2  FROM A SAMPLE             BI
				 N01010*/
	/*             IFAIL = 0  IF ALL OK                                     BI
		       N01020*/
	/*                   = 1 IF T1 AND T2 LIE OUTSIDE REQUIRED DOMAIN       BI
			     N01030*/
	/*                   = 2 IF ITERATIONS DO NOT CONVERGE                  BI
			     N01040*/
	/*                        THIS SHOULD NOT HAPPEN!                       BI
				  N01050*/
	/*                                                                      BI

										N01060*/
	// 	/* Parameter adjustments */
	// 	hi -= 3;

	

	// catch precision errors
	double eigEps = 1E-15;
	
	// if t1 and t2 are almost zero, we have very high concentration.
	// use a small delta to give approximate result without causing the 
	// algorithm to break down
	if (t1 < 0.0 && t1 > -eigEps) {
	    t1 = eigEps;
	    
	    if (t2 < t1) {
		t2 = t1;
	    }
	}

	/* Function Body */
	if (t1 <= t2 && t1 > 0.0 && t1 + t2 * 2.0 <= 1.0) {
	} 
	else {

	    // something wrong
	    throw new misc.LoggedException("Can't use t1 == " + t1 + " and t2 == " + t2 + " to calculate A");
	}
	
	if (t2 <= 0.04) {
	    mode = 3;
	    ak1 = -1.0 / (t1 * 2.0);
	    ak2 = -1.0 / (t2 * 2.0);
	} else if (t1 <= 0.04 && t2 / t1 >= 3.0) {
	    mode = 2;
	    if (0.5 - t2 <= 0.35) {
		ak2 = (1.0 - t2 * 2.0) * -2.0;
	    } else {
		ak2 = -1.0 / (0.5 - t2);
	    }
	    ak1 = -1.0 / (t1 * 2.0);
	} else {
	    mode = 1;
	    ak1 = 0.0;
	    ak2 = 0.0;
	}
	for (i = 1; i <= 20; ++i) {
	    bk1 = ak1;
	    bk2 = ak2;
	    level[0] = 8;
	    //	    System.out.println("level == " + level[0]);
	    fc = bingc(ak1, ak2, der, hes, mode, level);
	    fd1 = der[0] - t1;
	    fd2 = der[1] - t2;
	    
	    inv(hes, hi);

	    ak1 = bk1 - hi[0] * fd1 - hi[2] * fd2;
	    ak2 = bk2 - hi[1] * fd1 - hi[3] * fd2;

	    if (ak2 > 0.0) {
		ak2 = 0.0;
	    }
	    if (ak1 > ak2) {
		ak1 = ak2;
	    }
	    d__1 = ak1 - bk1;
	    d__2 = ak2 - bk2;
	    // 	    rerr = (d__1 = *ak1 - bk1, fabs(d__1)) / (1 - *ak1) + (d__2 = *ak2 - 
	    //								   bk2, fabs(d__2)) / (1 - *ak2);
	    
 	    rerr = Math.abs(d__1) / (1 - ak1) + Math.abs(d__2) / (1 - ak2);

	    if (rerr < 1e-7) {
		return new double[] {ak1, ak2, fc / (4.0 * Math.PI)};
	    }
	    
	}

	// couldn't converge
	throw new ConvergenceException("bngpar couldn't converge");

    } /* bngpar_ */
    

    /**
     * Automatically chooses which type of numerical method to use.
     * @param ak1 kappa1 < 0
     * @param ak2 0 > kappa2 > kappa1 
     */
   public static double bingc(double ak1, double ak2) throws ConvergenceException {

       if (ak2 > 0.0 || ak1 > ak2) {
	   throw new 
	       IllegalArgumentException("Illegal concentration parameters {ak1, ak2} = {" + 
					ak1 + ", " + ak2 + "}");
       }

	double[] der = new double[2];

	double[] hes = new double[4];	/* was [2][2] */;
	
	int[] level = new int[1]; level[0] = 8;

	int mode = 0;

	// from Kent's paper
	if (ak2 < -8.5) {
	    mode = 3;
	}
	else {
	    mode = 1;

	    if (ak1 <= -10.0 && ak1 / ak2 >= 2.0) {
		mode = 2;
	    }

	}

	
	return bingc(ak1, ak2, der, hes, mode, level) / (4.0 * Math.PI);
   }



    /**
     * Get the normalization constant and its derivatives. 
     * Automatically chooses which type of numerical method to use.
     * @param ak1 kappa1 < 0
     * @param ak2 kappa2 > kappa1
     * @param der {dC/dk_1 dC/dk_2}
     * @param hes {d^2C/dk_1^2, d^2C/dk_1dk_2, d^2C/dk_2dk_1, d^2C/dk_2^2} 
     */
   public static double bingc(double ak1, double ak2, double[] der, double[] hes) throws ConvergenceException {

	
	int[] level = new int[1]; level[0] = 8;

	int mode = 0;

	// from Kent's paper
	if (ak2 < -8.5) {
	    mode = 3;
	}
	else {
	    mode = 1;

	    if (ak1 <= -10.0 && ak1 / ak2 >= 2.0) {
		mode = 2;
	    }

	}


	double bingc = bingc(ak1, ak2, der, hes, mode, level);

	der[0] = der[0] / (4.0 * Math.PI);
	der[1] = der[1] / (4.0 * Math.PI);

	hes[0] = hes[0] / (4.0 * Math.PI);
	hes[1] = hes[1] / (4.0 * Math.PI);
	hes[2] = hes[2] / (4.0 * Math.PI);
	hes[3] = hes[3] / (4.0 * Math.PI);

	return bingc / (4.0 * Math.PI);
   }


    private static double bingc(double ak1, double ak2, double[] der, 
		 double[] hes, int mode, int[] level) throws ConvergenceException {
	/* Initialized data */
	
	final double c1 = 15.74960995;
	final double pi2 = 6.283185308;
	int init = 0;
	
	/* System generated locals */
	int i__1, i__2, i__3;
	double ret_val = 0.0;
	
	/* Local variables */
	double dera, derb, fact, derd, derl, trma, trmb, term, fact0, term0;

	double[] a = new double[10];
	double[] b = new double[8];
	double[] d = new double[8];
	double[] e = new double[9];

	int i, j, k, n;
	double hesaa, hesab;
	int ifail;
	double hesbb;
	double[] cfold = new double[9];
	double y, termb, termd, terml, b2;
	double[] bb = new double[9]; 
	double ea,be;
	double[] dd = new double[9];
	double[] cf = new double[9];
	double eb, de, al;

	double add, eps;
	int imx = 0;
	int jmx = 0;
	double sum, add0, bes0;


	/*    THIS FUNCTION CALCULATES THE NORMALIZATION CONSTANT FOR THE       BI
	      N00050*/
	/*      BINGHAM DISTRIBUTION AND ITS DERIVATIVES WITH RESPECT TO        BI
		N00060*/
	/*       THE PARAMETERS AK1 AND AK2.                                    BI
		 N00070*/
	/*                                                                      BI
										N00080*/
	/*    INPUT:  AK1 <= AK2 <= 0   CONCENTRATION PARAMETERS                BI
	      N00090*/
	/*            MODE = 1  USE TAYLOR SERIES EXPANSION                     BI
		      N00100*/
	/*                   2  NORMAL-VON MISES ASYMPTOTIC SERIES              BI
			     N00110*/
	/*                   3  BIVARIATE NORMAL ASYMPTOTIC SERIES              BI
			     N00120*/
	/*             LEVEL >= 0  NUMBER OF TERMS USED IN ASYMPTOTIC SERIES    BI
		       N00130*/
	/*                         (MODES 2 AND 3)                              BI
				   N00140*/
	/*     OUTPUT: BINGC   NORMALIZATION CONSTANT                           BI
	       N00150*/
	/*             DER(2)  DERIVATIVES OF BINGC WRT  AK1 AND AK2            BI
		       N00160*/
	/*             HES(2,2) SECOND DERIVATIVES                              BI
		       N00170*/
	/*                      BASED ON POWER SERIES(MODE 1) OR THE DOMINANT   BI
				N00180*/
	/*                      TERM IN THE ASYMPTOTIC EXPANSION(MODE = 2 OR 3) BI
				N00190*/
	/*             LEVEL    NUMBER OF TERMS USED IN EACH SUM (MODE 1)       BI
		       N00200*/
	/*             IND = IFAIL = 1 IF SERIES DO NOT CONVERGE                BI
		       N00210*/
	/* Parameter adjustments */
	//	hes -= 3;
	// 	--der;

	/* Function Body */
	//    *ind = 0;
	eps = 1e-6;
	/*    LEVUP = MAXIMUM NO OF LEVELS                                      BI
	      N00310*/
	if (ak1 > ak2 || ak2 > 0.0) {
	    ifail = 1;
	    return ret_val;
	}
	//     if (init == 1) {
	// 	goto L6;
	//     }
	if (init != 1) {
	    
	    init = 1;
	    e[0] = 1.0;
	    for (i = 1; i <= 8; ++i) {
		/* L5: */
		e[i] = e[i - 1] * (i * 2.0 - 1.0);
	    }
	    /*    THE E(I) REPRESENT THE 2*I'TH CENTRAL MOMENT OF THE STANDARD      BI
		  N00430*/
	    /*    NORMAL DISTRIBUTION                                               BI
		  N00440*/
	    
	}


	// L6:
	/* L100: */
	//     if (*mode != 1) {
	// 	goto L200;
	//     }
    
	if (mode == 1) {
	
	    /*                                                                      BI
										    N00480*/
	    /*    NOW WE FIND BOUNDS ON THE NUMBER OF TERMS NEEDED FOR GIVEN        BI
		  N00490*/
	    /*     ACCURACY IN TERMS OF THE EXP FUNCTION                            BI
		   N00500*/
	    /*                                                                      BI
										    N00510*/
	    //	*ind = 0;
	    al = -(ak1);
	    be = -(ak1) + ak2;
	    ea = Math.exp(al);
	    trma = al;
	    for (i = 2; i <= 100; ++i) {
		trma = trma * al / i;
		if ((double) i > al && trma / ea < eps) {
		    imx = i;
		    break; // goto L62;
		}
	    }
	    if (i == 101) {
		throw new ConvergenceException("L61: Normalization constant couldn't converge");
	    }
		/* L61: */
	   
	    // 	    *ind = 1;
	    //	L62:
	    eb = Math.exp(be);
	    trmb = be;
	    for (i = 2; i <= 100; ++i) {
		trmb = trmb * be / i;
		if ((double) i > be && trmb / eb < eps) {
		    jmx = i;
		    break; // goto L64;
		}
		/* L63: */
	    }
	    if (i == 101) {
		throw new ConvergenceException("L63: Normalization constant couldn't converge");
	    }
	    //	*ind = 1;
	    //	L64:
	    term0 = pi2 * 2;
	    ret_val = 0.0;
	    dera = 0.0;
	    derb = 0.0;
	    hesaa = 0.0;
	    hesab = 0.0;
	    hesbb = 0.0;
	    if (Math.abs(al) < 1e-15) {
		al = 1e-15;
	    }
	    if (Math.abs(be) < 1e-15) {
		be = 1e-15;
	    }
	    /*    BINGC = EXP(AK1)*SUM( TERM(I,J) )  OVER I,J >=0                   BI
		  N00860*/
	    /*     TERM0 = TERM(I,0)                                                BI
		   N00870*/
	    /*                                                                      BI
										    N00880*/
	    i__1 = imx;
	    for (i = 0; i <= i__1; ++i) {
		if (i >= 1) {
		    term0 = term0 * (i - 0.5) / (i + 0.5) / i * al;
		}
		term = term0;
		sum = term;
		dera += term * i / al;
		hesaa += term * i * (i - 1) / al / al;
		i__2 = jmx;
		for (j = 1; j <= i__2; ++j) {
		    term = term * (j - .5) / ((i + j + .5) * j) * be;
		    dera += term * i / al;
		    derb += term * j / be;
		    hesaa += term * i * (i - 1) / al / al;
		    hesab += term * i * j / al / be;
		    hesbb += term * j * (j - 1) / be / be;
		    sum += term;
		    if (term < 1e-16 && j >= 2) {
			break; //    goto L34;
		    }
		    /* L32: */
		}
		// L34:
		ret_val += sum;
		if (sum < 1e-16 && i >= 2) {
		    break; // goto L33;
		}
		/* L31: */
	    }
	    //	L33:
	    level[0] = (imx>jmx)?imx:jmx;
	    der[0] = 1 - (dera + derb) / ret_val;
	    der[1] = derb / ret_val;
	    b2 = ret_val * ret_val;
	    hesaa = hesaa / ret_val - dera * dera / b2;
	    hesab = hesab / ret_val - dera * derb / b2;
	    hesbb = hesbb / ret_val - derb * derb / b2;
	    hes[0] = hesaa + hesab * 2 + hesbb;
	    hes[2] = -hesab - hesbb;
	    hes[1] = hes[2];
	    hes[3] = hesbb;
	    ret_val *= Math.exp(ak1);

	    // exception already thrown if needed
	    // 	if (*ind == 1) {
	    // 	    cout << " NO CONVERGENCE IN BINGC; MODE 1\n";
	    // 	    exit(0);
	    // 	}
	    return ret_val;

	} // end if mode == 1
    
	// L200:
	//    if (*mode != 2) {
	//        goto L300;
	//    }
	else if (mode == 2) {

	    /*    A(K)=  BESRAT   (BE) = I (BE)/I   (BE)                            BI
		  N01260*/
	    /*                 K-1        K      K-1                                BI
			       N01270*/
	    /*    THEN SET                                                          BI
		  N01280*/
	    /*     A(K) = K'TH MOMENT FOR VON MISES DSN                             BI
		   N01290*/
	    /*          = I (BE)/I (BE) = OLD  A(1)...A(K)                          BI
			N01300*/
	    /*             K      0                                                 BI
			   N01310*/
	    /*                                                                      BI
										    N01320*/
	    al = ak2 - ak1 * 2.0;
	    be = -(ak2) / 2.0;
	    if (al <= 0.001) {
		ret_val = c1;
		der[0] = 3.0;
		der[1] = 1.0;
		hes[0] = 1.0;
		hes[2] = 0.0;
		hes[1] = 0.0;
		hes[3] = 1.0;
		return ret_val;
	    }
	    y = be * be / 4.0;
	    bes0 = bstrs0(y);
	    /*           = I (BE) * 1 * EXP(-BE)                                    BI
			 N01470*/
	    /*              0                                                       BI
			    N01480*/
	    /*                                                                      BI
										    N01490*/
	    a[9] = besrat(c_b28, be, eps);
	    for (i = 8; i >= 1; --i) {
		/* L1: */
		a[i] = be / (be * a[i + 1] + (i << 1));
	    }
	    for (i = 2; i <= 9; ++i) {
		if (a[i - 1] < 1e-15) {
		    a[i - 1] = 0.;
		}
		/* L2: */
		a[i] = a[i - 1] * a[i];
	    }
	    a[0] = 1.;
	    bb[0] = a[1];
	    for (i = 0; i <= 8; ++i) {
		cfold[i] = 0.0;
		/* L50: */
		cf[i] = 0.0;
	    }
	    cf[0] = 1.0;
	    /* L56: */
	    term = 1.0;
	    fact0 = 1.0;
	    termb = bb[0];
	    terml = 0.0;
	    // 	if (*level == 0) {
	    // 	    goto L8;
	    // 	}
	    if (level[0] != 0) {

		i__1 = level[0];
		for (k = 1; k <= i__1; ++k) {
		    n = k;
		    i__2 = n - 1;
		    for (i = 0; i <= i__2; ++i) {
			/* L52: */
			cfold[i] = cf[i];
		    }
		    if (n == 1) {
			cf[0] = 0.0;
			cf[1] = 1.0;
		    }
		    if (n == 2) {
			cf[0] = 0.5;
			cf[1] = 0.0;
			cf[2] = 0.5;
		    }
		    if (n > 2) {
			cf[0] = cfold[1] / 2.0;
			cf[1] = cfold[0] + cfold[2] / 2.0;
			cf[n - 1] = cfold[n - 2] / 2.0;
			cf[n] = cfold[n - 1] / 2.0;
		    }
		    // 		if (n <= 3) {
		    // 		    goto L54;
		    // 		}
		    if (n > 3) {
		    
			i__2 = n - 2;
			for (i = 2; i <= i__2; ++i) {
			    /* L53: */
			    cf[i] = (cfold[i - 1] + cfold[i + 1]) / 2.0;
			}
		    }
		    //		L54:
		    b[n - 1] = 0.0;
		    bb[n] = 0.0;
		    i__2 = n;
		    for (i = 0; i <= i__2; ++i) {
			if (i == 0) {
			    bb[n] += cf[0] * a[1];
			} else {
			    bb[n] += cf[i] * (a[i - 1] + a[i + 1]) / 2.0;
			}
			/* L55: */
			b[n - 1] += cf[i] * a[i];
		    }
		    if (k == 1) {
			fact0 = fact0 * (-1 / al) / k;
		    } else {
			fact0 = fact0 * (-be / al) / k;
		    }
		    fact = fact0 * be;
		    termb += e[k] * (b[k - 1] * fact0 * k + bb[k] * fact);
		    terml -= k * e[k] * b[k - 1] * fact / al;
		    add = e[k] * b[k - 1] * fact;
		    term += add;
		    if (Math.abs(add) / term < 1e-15) {
			break; // goto L8;
		    }
		    /* L3: */
		}
	    } // if level != 0

 	
	    // 	L8:
	    /* L4: */
	    ret_val = c1 / Math.sqrt(al) * bes0 * term;
	    derl = -0.5 / al + terml / term;
	    derb = termb / term - 1;
	    der[0] = derl * -2;
	    der[1] = derl - derb * 0.5;
	    hes[0] = 1.0 / (ak1 * 2.0 * ak1);
	    hes[2] = 0.0;
	    hes[1] = hes[2];
	    hes[3] = (a[2] * 0.5 + 0.5 - a[1] * a[1]) / 4.0;
	    return ret_val;

	

	} // end if mode == 2
	else if (mode == 3) {
	    // 	 L300:
	    // 	if (*mode != 3) {
	    // 	goto L400;
	    //     }

	    al = ak2 * -2.0;
	    if (al == 0.0) {
		ret_val = 1.0;
		der[0] = 1.0;
		der[1] = 1.0;
		hes[0] = 1.0;
		hes[2] = 0.0;
		hes[1] = 0.0;
		hes[3] = 1.0;
		return ret_val;
	    }
	    de = ak2 / ak1;
	    for (i = 0; i <= 8; ++i) {
		cf[i] = 0.0;
		/* L41: */
		cfold[i] = 0.0;
	    }
	    cf[0] = 1.0;
	    /* L47: */
	    term = 1.0;
	    fact = 1.0;
	    add = 1.0;
	    termd = 0.0;
	    terml = 0.0;
	    // 	if (*level == 0) {
	    // 	    goto L9;
	    // 	}
	    if (level[0] != 0) {
 	    
		i__1 = level[0];
		for (k = 1; k <= i__1; ++k) {
		    n = k;
		    i__2 = n - 1;
		    for (i = 0; i <= i__2; ++i) {
			/* L43: */
			cfold[i] = cf[i];
		    }
		    cf[0] = cfold[0];
		    cf[n] = cfold[n - 1];
		    // 		if (n < 2) {
		    // 		    goto L44;
		    // 		}
		    if (n >= 2) {
			i__2 = n - 1;
			for (i = 1; i <= i__2; ++i) {
			    /* L45: */
			    cf[i] = cfold[i - 1] + cfold[i];
			}
			/*    THE CF REPRESENT THE BINOMIAL COEFFICIENTS                      
			      BIN02560*/
		    }

		    //		L44:
		    d[n - 1] = 0.0;
		    dd[n] = 0.0;
		    i__2 = n;
		    for (i = 0; i <= i__2; ++i) {
			if (i > 0) {
			    i__3 = i - 1;
			    dd[n] += i * cf[i] * e[i] * e[n - i] * Math.pow(de, i__3);
			}
			/* L46: */
			d[n - 1] += cf[i] * e[i] * e[n - i] * Math.pow(de, i);
		    }
		    fact = -fact * (-0.5 - k + 1) / k / al;
		    termd += fact * dd[k];
		    terml -= fact * d[k - 1] * k / al;
		    add0 = add;
		    add = fact * d[k - 1];
		    term += add;
		    if (Math.abs(add) / term < 1e-10) {
			break; // goto L9;
		    }
		    /* L12: */
		}
	    
	    } // if level != 0
	
	    //	L9:
	    /* L11: */
	    ret_val = 2 / al * Math.sqrt(de) * pi2 * term;
	    derl = -1 / al + terml / term;
	    derd = 0.5 / de + termd / term;
	    der[0] = derd * (-(ak2)) / (ak1 * ak1);
	    der[1] = derl * -2 + derd / ak1;
	    hes[0] = 1.0 / (ak1 * 2.0 * ak1);
	    hes[2] = 0.0;
	    hes[1] = 0.0;
	    hes[3] = 1.0 / (ak2 * 2.0 * ak2);
	    /*    DER(1)=DERL                                                       BI
		  N02810*/
	    /*    DER(2)=DERD                                                       BI
		  N02820*/
	    /*    THIS GIVES DERIVATIVES WRT LAMBDA AND SPIKE                       BI
		  N02830*/
	    /*                                                                      BI
										    N02840*/
	    return ret_val;

	} // end if mode == 3

	else { 
	    //	L400:
	    System.out.println("INVALID MODE IN BINGC");
	    return ret_val;
	}

    } /* bingc_ */	



    private static double besrat(double v, double x, double eps)
    {
	/* Format strings */
	//	String fmt_3 = "(\002 FAULTY RATIO\002,4(2x,f12.5))";
	
	/* System generated locals */
	int i__1;
	
	/* Local variables */
	int m, n;
	double r, a1, u0, u1, y0, x2, y1, er, cx;
	double[] d1 = new double[20];
	double rat = 0.0;
	
	/*    THIS ROUTINE CALCULATES THE BESSEL FUNCTION RATIO                 BI
	      N02920*/
	/*                                                                      BI
										N02930*/
	/*           I   (X)/I (X)  = RAT, SAY                                  BI
		     N02940*/
	/*            V+1     V                                                 BI
		      N02950*/
	/*                                                                      BI
										N02960*/
	/*      WHERE I(.) IS THE MODIFIED BESSEL FUNCTION,                     BI
		N02970*/
	/*      BY THE METHOD OF AMOS(1974)                                     BI
		N02980*/
	/*    INPUT:  V  ORDER OF THE BESSEL FUNCTION                           BI
	      N02990*/
	/*            X  ARGUMENT OF THE BESSEL FUNCTION                        BI
		      N03000*/
	/*            EPS  ACCURACY DESIRED                                     BI
		      N03010*/
	/*    OUTPUT: RAT   BESSEL FUNCTION RATIO                               BI
	      N03020*/
	/*                                                                      BI
										N03030*/
	x2 = x * x;
	cx = (v + 1.0) * 2.0;
	u0 = v + 1.5;
	y0 = u0 * u0;
	d1[0] = x / (u0 - 1. + Math.sqrt(y0 + x2));
	for (n = 1; n <= 19; ++n) {
	    u0 = v + (double) n + 1.5;
	    y0 = u0 * u0;
	    a1 = x / (u0 - 1.0 + Math.sqrt(y0 + x2));
	    i__1 = n;
	    for (m = 1; m <= i__1; ++m) {
		r = a1 / d1[m - 1];
		d1[m - 1] = a1;
		u1 = v + (double) n - (double) m + 1.0;
		y1 = u1 * u1;
		/* L2: */
		a1 = x / (u1 + Math.sqrt(y1 + x2 * r));
	    }
	    d1[n] = a1;
	    rat = a1;
	    er = x / (cx + x * d1[n - 1]);
	    er -= a1;
	    er = Math.abs(er) / a1;
	    /* L1: */
	    if (er < eps) {
		return rat;
	    }
	}
	return rat;
    } /* besrat_ */

/*VERSION 1 OF LBJK11 FORTRAN AT 13:41:51 ON FRIDAY 04/03/81            BIN033
00*/



    private static double bstrs0(double y)
    {
	/* Format strings */
	//	String fmt_2 = "(\002 BSTRS0 NO CONVERGENCE, BSTRS0=\002,e12.5\ ,\002 COEF=\002,e12.5,\002 Y=\002,e12.5);
	//	String fmt_8 = "(\002 BSTRS0  ASYMPTOTICS NOT CONVERGED\002)";
	
	/* System generated locals */
	double ret_val;
	
	/* Local variables */
	double coef;
	int i, k;
	double t, u, x, z, b1, f1;

	final double pi2 = Math.atan(1.0) * 8.0;
	
	/*                                                                      BI
										N03350*/
	/*   INPUT:                                                             BI
	     N03360*/
	/*        Y = X*X/4   WHERE X IS THE ARGUMENT OF THE BESSEL FUNCTION    BI
		  N03370*/
	/*                                                                      BI
										N03380*/
	/*   OUTPUT:                                                            BI
	     N03390*/
	/*                                                                      BI
										N03400*/
	/*        BSTRS0 = I (X) * EXP(-X)                                      BI
		  N03410*/
	/*                  0                                                   BI
			    N03420*/
	/*                                                                      BI
										N03430*/
	/*                WHERE  X = SQRT(4*Y)                                  BI
			  N03440*/
	/*                                                                      BI
										N03450*/
	/*          SO THAT   0   LT  BSTRS0  LE  1                             BI
		    N03460*/
	/*                                                                      BI
										N03470*/
	/*                                                                      BI
										N03480*/
	x = Math.sqrt(y * 4.0);
	//     if (*y >= 100.) {
	// 	goto L3;
	//     }
	if (y < 100.0) {
	    ret_val = Math.exp(-x);
	    coef = ret_val;
	    for (k = 1; k <= 200; ++k) {
		coef = coef * y / (k * k);
		b1 = ret_val;
		ret_val += coef;
		if (b1 == ret_val) {
		    return ret_val;
		}
		/* L1: */
	    }
	    return ret_val;
	}
	
	//L3:
	
	ret_val = 1.0 / Math.sqrt(pi2 * x);
	coef = ret_val;
	u = 0.0;
	z = x * 8;
	for (i = 1; i <= 40; ++i) {
	    t = i * 2.0 - 1.0;
	    coef = -coef * (u - t * t) / i / z;
	    f1 = ret_val;
	    ret_val += coef;
	    if (f1 == ret_val) {
		return ret_val;
	    }
	/* L7: */
	}
	return ret_val;
    } /* bstrs0_ */
    








    private static int inv(double[] a, double[] b)
    {
	double det;

	/* Parameter adjustments */
	//	b -= 3;
	//	a -= 3;
	// 
	
	/* Function Body */
	det = a[0] * a[3] - a[2] * a[1];
	b[0] = a[3] / det;
	b[2] = -a[1] / det;
	b[1] = -a[2] / det;
	b[3] = a[0] / det;
	return 0;
    } /* inv_ */



    public static void main(String[] args) {

	System.out.println("Testing Bingham distribution\n");

	double t2 = 0.0;
	double t3 = 0.0;
	
	double ak1 = 0.0;
	double ak2 = 0.0;
	double afc = 0.0;
	
	int i = 0;
	int j = 0;

	System.out.print("t2\t\tt3\t\tak1\t\tak2\n");
	
	t2 = 0.0;
	t3 = 0.0;


	DecimalFormat df = new DecimalFormat("0.000000");
	
	for (i = 10; i > 0; i--) {
	    for (j = 10; j >= i; j--) {
		t2 = 0.33333332 / (double)i;
		t3 = 0.33333332 / (double)j;
		
		if (t3 <= t2) {
		   System.out.print(df.format(t2) + "\t" + df.format(t3));
		    
		    //getbingpars(t3, t2, ak1, ak2, afc); // oddly, program asks for them in reverse order
		    //printf("\t%f\t%f\t%f\n", ak1, ak2, afc);
		    
		   try {
		       double[] bngpars = bngpar(t2, t3);
		       
		       System.out.print("\t" + df.format(bngpars[0]) + "\t" + df.format(bngpars[1]) + "\t" + df.format(bngpars[2] * 4.0 * Math.PI) + "\n");
		   }
		   catch (ConvergenceException e) {

		   }
		    
		}
	    }
	}
	
	System.out.println("\n");

	for (i = 10; i > 0; i--) {
	    for (j = 10; j >= i; j--) {
		t2 = 0.33333332 / (double)(2*i);
		t3 = 0.33333332 / (double)(2*j);
		
		if (t3 <= t2) {
		    
		    //getbingpars(t3, t2, ak1, ak2, afc); // oddly, program asks for them in reverse order
		    //printf("\t%f\t%f\t%f\n", ak1, ak2, afc);
		    double[] der = new double[2];
		    
		   try {

		       double[] bngpars = bngpar(t2, t3);
		       double bingc = bingc(bngpars[0], bngpars[1], der, new double[4]);

		       System.out.print(df.format(bngpars[0]) + "\t" + df.format(bngpars[1]) + "\t" + df.format(der[0] * 4.0 * Math.PI) + "\t" + df.format(der[1] * 4.0 * Math.PI) + "\n");

		   }
		   catch (ConvergenceException e) {

		   }   

		}
	    }
	} 
	
    }
    

}

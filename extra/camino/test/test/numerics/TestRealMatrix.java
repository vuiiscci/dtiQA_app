package numerics;

import junit.framework.*;
import junit.extensions.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>RealMatrix.java</code>.
 * <BR>
 * </dl>
 *
 * @version $Id: TestRealMatrix.java,v 1.3 2006/02/24 16:36:28 kseunari Exp $
 * @author  Daniel Alexander
 * @see numerics.RealMatrix
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestRealMatrix extends TestCase {


    public TestRealMatrix(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	// any class variables you declare should be initialized here. This is called before each test
    }
    
    protected void tearDown() {
	// does the opposite of setup. Should be used to free memory, I/O resources, etc
    }
    
    public static Test suite() {
	return new TestSuite(TestRealMatrix.class);
    }


    public void testSVD() {

	// Testing RealMatrix
	RealMatrix rm = new RealMatrix(3,4);
	rm.setEntry(0, 0, 0.9501);
	rm.setEntry(0, 1, 0.4860);
	rm.setEntry(0, 2, 0.4565);
	rm.setEntry(0, 3, 0.4447);
	rm.setEntry(1, 0, 0.2311);
	rm.setEntry(1, 1, 0.8913);
	rm.setEntry(1, 2, 0.0185);
	rm.setEntry(1, 3, 0.6154);
	rm.setEntry(2, 0, 0.6068);
	rm.setEntry(2, 1, 0.7621);
	rm.setEntry(2, 2, 0.8214);
	rm.setEntry(2, 3, 0.7919);
	try {
	    RealMatrix[] svd = rm.svd();

	    assertEquals(svd[0].entries[0][0], -0.553350, 0.000001);
	    assertEquals(svd[0].entries[1][0], -0.451129, 0.000001);
	    assertEquals(svd[0].entries[2][1], 0.114099, 0.000001);

	    assertEquals(svd[1].entries[0][0], 2.103685, 0.000001);
	    assertEquals(svd[1].entries[1][1], 0.667638, 0.000001);
	    assertEquals(svd[1].entries[2][2], 0.389216, 0.000001);

	    assertEquals(svd[2].entries[0][1], 0.578624, 0.000001);
	    assertEquals(svd[2].entries[1][2], -0.131603, 0.000001);
	    assertEquals(svd[2].entries[3][3], -0.784575, 0.000001);
	} catch(Exception e) {
	    fail("SVD failed");
	}


    }


   public void testInverse() {

	/**
	   >> rand(3,3)
	   
	   ans =
	   
	   0.95012928514718   0.48598246870930   0.45646766516834
	   0.23113851357429   0.89129896614890   0.01850364324822
	   0.60684258354179   0.76209683302739   0.82140716429525
	   
	   >> inv(ans)
	   
	   ans =
	   
	   1.67404446726824  -0.11964439680952  -0.92759516260175
	   -0.41647243581381   1.17375817197623   0.20499869641804
	   -0.85035677245799  -1.00061468473081   1.71251901466682

	*/

	RealMatrix square = new RealMatrix(3,3);

	square.entries = new double[][] {{0.95012928514718, 0.48598246870930, 0.45646766516834},
			  {0.23113851357429, 0.89129896614890, 0.01850364324822},
			  {0.60684258354179, 0.76209683302739, 0.82140716429525}};

	double[][] inv = {{1.67404446726824,-0.11964439680952,-0.92759516260175},
			  {-0.41647243581381, 1.17375817197623, 0.20499869641804},
			  {-0.85035677245799,-1.00061468473081, 1.71251901466682}};


	GenTestMethods.assertArraysEqual(inv, square.inverse().entries, 0.00000001);
    }
  

    public void testProduct() {

	/*

	>> a = randn(3,4)
	
	a =
	
	0.38033725171391  -0.04822078914531   1.09500373878749   0.89563847121175
	-1.00911552434078   0.00004319184163  -1.87399025764096   0.73095733842945
	-0.01951066953029  -0.31785945124769   0.42818327304516   0.57785734633080
	
	>> b = randn(4,2)
	
	b =
	
	0.04031403161844  -0.37746895552236
	0.67708918759730  -0.29588711000356
	0.56890020520072  -1.47513450585526
	-0.25564541563196  -0.23400404765603
	
	>> a*b
	
	ans =
	
	0.37666513550642  -1.95415842183436
	-1.29363160513411   2.97423771989152
	-0.12013878094357  -0.66543369244201
	
	>> 
	
	*/

	double[][] ae = {{0.38033725171391, -0.04822078914531, 1.09500373878749, 0.89563847121175},
			{-1.00911552434078, 0.00004319184163,-1.87399025764096, 0.73095733842945},
			{-0.01951066953029,-0.31785945124769, 0.42818327304516, 0.57785734633080}};


	
	double[][] be = {{0.04031403161844, -0.37746895552236},
			{0.67708918759730, -0.29588711000356},
			{0.56890020520072, -1.47513450585526},
			{-0.25564541563196, -0.23400404765603}};


	double[][] prodE = {{0.37666513550642,  -1.95415842183436},
			   {-1.29363160513411,  2.97423771989152},
			   {-0.12013878094357, -0.66543369244201}};

	RealMatrix a = new RealMatrix(3,4);
	a.entries = ae;

	RealMatrix b = new RealMatrix(4,2);
	b.entries = be;
	
	RealMatrix prod = a.product(b);

	assertEquals(3,prod.r);
	assertEquals(2,prod.c);

	GenTestMethods.assertArraysEqual(prodE, prod.entries, 0.0000001);
	
    }


    
    public void testJacobi() {
	
	
	RealMatrix a = new RealMatrix(3,3);

	a.entries[0][0] = 0.81316649730376;
	a.entries[0][1] = 0.89364953091353;
	a.entries[0][2] = 0.05789130478427;
	a.entries[1][0] = a.entries[0][1];
	a.entries[1][1] = 0.41027020699095;
	a.entries[1][2] = 0.00986130066092;
	a.entries[2][0] = a.entries[0][2];
	a.entries[2][1] = a.entries[1][2];
	a.entries[2][2] = 0.19872174266149;

	RealMatrix[] j = a.jacobi();

	assertEquals(-0.30596402890815, j[0].entries[1][1], 1E-6);
	assertEquals(0.19834720084370, j[0].entries[2][2], 1E-6);
	assertEquals(1.52977527502065, j[0].entries[0][0], 1E-6);

	assertEquals(-0.62477669321024, j[1].entries[0][1], 1E-6);
	assertEquals(0.77876021516487, j[1].entries[1][1], 1E-6);
	assertEquals(0.05645007438119, j[1].entries[2][1], 1E-6);

	assertEquals(0.00515904076713, j[1].entries[0][2], 1E-6);
	assertEquals(-0.06817852319658, j[1].entries[1][2], 1E-6);
	assertEquals(0.99765979836470, j[1].entries[2][2], 1E-6);

	assertEquals(0.78078644194149, j[1].entries[0][0], 1E-6);
	assertEquals(0.62360581800613, j[1].entries[1][0], 1E-6);
	assertEquals(0.03857869657286, j[1].entries[2][0], 1E-6);


    }
     
    public void testPseudoInverse() {
	
	/**
	   >> rand(3,3)
	   
	   ans =
	   
	   0.95012928514718   0.48598246870930   0.45646766516834
	   0.23113851357429   0.89129896614890   0.01850364324822
	   0.60684258354179   0.76209683302739   0.82140716429525
	   
	   >> pinv(a)

	   ans =

	    1.67404446726824  -0.11964439680953  -0.92759516260175
	   -0.41647243581382   1.17375817197623   0.20499869641805
	   -0.85035677245799  -1.00061468473080   1.71251901466683
	*/

	RealMatrix square = new RealMatrix(3,3);

	square.entries = new double[][] {{0.95012928514718, 0.48598246870930, 0.45646766516834},
			  {0.23113851357429, 0.89129896614890, 0.01850364324822},
			  {0.60684258354179, 0.76209683302739, 0.82140716429525}};

	double[][] inv = {{1.67404446726824,-0.11964439680952,-0.92759516260175},
			  {-0.41647243581381, 1.17375817197623, 0.20499869641804},
			  {-0.85035677245799,-1.00061468473081, 1.71251901466682}};


	GenTestMethods.assertArraysEqual(inv, square.pseudoInv().entries, 0.00000001);
    }
	 

    public void testAdd() {

	/*

	>> a = rand(2)
	
	a =
	
	0.7922    0.6557
	0.9595    0.0357
	
	
	>> b = rand(2)

	b =
	
	0.8491    0.6787
	0.9340    0.7577

	>> a + b
	
	ans =
	
	1.6413    1.3345
	1.8935    0.7935
	
	*/

	RealMatrix a = new RealMatrix(new double[][] {{ 0.7922, 0.6557}, {0.9595, 0.0357}});

	RealMatrix b = new RealMatrix(new double[][] {{ 0.8491, 0.6787}, {0.9340, 0.7577}});

	RealMatrix c = a.add(b);

	assertEquals(1.6413, c.entries[0][0], 1E-4);
	assertEquals(1.3345, c.entries[0][1], 1E-4);
	assertEquals(1.8935, c.entries[1][0], 1E-4);
	assertEquals(0.7935, c.entries[1][1], 1E-4);

	
    }
 
	
    public void testSub() {
	
	/*
	  a, b as above


	  >> a - b
	  
	  ans =
	  
	  -0.0569   -0.0230
	  0.0255   -0.7220
	*/


	RealMatrix a = new RealMatrix(new double[][] {{ 0.7922, 0.6557}, {0.9595, 0.0357}});

	RealMatrix b = new RealMatrix(new double[][] {{ 0.8491, 0.6787}, {0.9340, 0.7577}});

	RealMatrix c = a.sub(b);

	assertEquals(-0.0569, c.entries[0][0], 1E-4);
	assertEquals(-0.0230, c.entries[0][1], 1E-4);
	assertEquals(0.0255, c.entries[1][0], 1E-4);
	assertEquals(-0.7220, c.entries[1][1], 1E-4);

	
    }
  
}

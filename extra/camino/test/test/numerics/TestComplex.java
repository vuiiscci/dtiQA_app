package numerics;

import junit.framework.*;
import junit.extensions.*;

/**
 * <dl>
 * <dt>Purpose: Automated tests for <code>Complex.java</code>.
 * <BR><BR>
 *
 * <dt>Description:
 * <dd> This class is used to perform tests on <code>Complex</code> with JUnit 3.8.
 * 
 * </dl>
 *
 * @version $Id: TestComplex.java,v 1.1 2005/09/26 10:36:47 ucacpco Exp $
 * @author  Philip Cook
 * @see numerics.Complex
 * @see <A HREF="http://www.junit.org">Junit Homepage</A>
 *
 *
 */
public class TestComplex extends TestCase {

    private Complex c1 = null;

    public TestComplex(String name) {
	super(name);
    }

    public static void main (String[] args) {
	junit.textui.TestRunner.run(suite());
    }


    /** 
     * Initialise global resources for tests.
     */
    protected void setUp() {
	c1 = new Complex(10.0, -4.0);
    }
    
    public static Test suite() {
	return new TestSuite(TestComplex.class);
    }



    // compare with output from Matlab
    public void testInverse() {
	try {
	    Complex inv = c1.inverse();
	    
	    assertEquals(0.08620689655172, inv.real(), 0.0000001);
	    assertEquals(0.03448275862069, inv.imag(), 0.0000001);

	}
	catch (ComplexNumberException e) {
	    fail("Unexpected exception " + e);
	}

    }
  
    
    public void testSqrt() {
	Complex sqrt = c1.sqrt();
	
	assertEquals(3.22260217947151, sqrt.real(), 0.0000001);
	assertEquals(-0.62061647346369, sqrt.imag(), 0.0000001);
	
    }

    
    /**
     * Test sin of complex number. If this works, sinh and cosh of real numbers also works.
     *
     */
    public void testSin() {

	Complex sin = c1.sin();
	
	assertEquals(-14.85625516387525, sin.real(), 0.0000001);
	assertEquals(22.89819255096376, sin.imag(), 0.0000001);

    }
    
    /**
     * Test cos of complex number. If this works, sinh and cosh of real numbers also works.
     *
     */
    public void testCos() {

	Complex cos = c1.cos();
	
	assertEquals(-22.91356068209214, cos.real(), 0.0000001);
	assertEquals(-14.84629106966036, cos.imag(), 0.0000001);

    }

    public void testTanhReal() {
	assertEquals(0.83697920564276, Complex.tanh(1.211), 0.00000001);
    }


    /**
     * Tests sinh of Complex number.
     */
    public void testSinh() {
	Complex sinh = c1.sinh();
	
	assertEquals(-7.198729413635292E03, sinh.real(), 0.0000001);
	assertEquals(8.334842155341617E03, sinh.imag(), 0.0000001);
    }


    /**
     * Tests cosh of Complex number.
     */
    public void testCosh() {
	Complex cosh = c1.cosh();
	
	assertEquals(-7.198729443310667E03, cosh.real(), 0.0000001);
	assertEquals(8.334842120982836E03, cosh.imag(), 0.0000001);
    }



    /**
     * Tests tanh of Complex number.
     */
    public void testTanh() {
	Complex tanh = c1.tanh();
	
	assertEquals(1.00000000059980, tanh.real(), 0.0000001);
	assertEquals(-0.00000000407844, tanh.imag(), 0.0000001);
    }
    

   

  
}

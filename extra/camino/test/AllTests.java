import junit.framework.*;

import apps.*;
import misc.*;
import data.*;
import numerics.*;
import imaging.*;
import inverters.*;
import tools.*;
import tractography.*;
import sphfunc.*;
import models.compartments.*;


public class AllTests extends TestCase {

    private static boolean testApps = true;
    private static boolean testData = true;
    private static boolean testImaging = true;
    private static boolean testInverters = true;
    private static boolean testMisc = true;
    private static boolean testNumerics = true;
    private static boolean testOptimizers = true;
    private static boolean testTools = true;
    private static boolean testSphFunc = true;
    private static boolean testTractography = true;
    private static boolean testGreatCircleIntegrals = true;
    private static boolean testModels = true;


    public static Test suite() {

        TestSuite suite = new TestSuite();

	// apps 
	if (testApps) {
		suite.addTest(TestSphFuncPICoCalibrationData.suite());
		suite.addTest(TestRGB_ScalarImage.suite());

	}
	

	// data
	if (testData) {
	    suite.addTest(TestDataSynthesizer.suite());
	    suite.addTest(TestGaussianMixture.suite());
	    suite.addTest(TestVoxelOrderDataSource.suite());
	    suite.addTest(TestStandardTestFunctions.suite());
	}
	
	// inverters
	if (testInverters) {
	    suite.addTest(TestInverters.suite());
	    suite.addTest(TestEvenSphHarmFitter.suite());
	    suite.addTest(TestModelIndex.suite());
	}

	// misc
	if (testMisc) {
	    suite.addTest(TestDT.suite());
	    suite.addTest(TestOrderedAcqSubsetMinimizer.suite());
	    suite.addTest(TestOrderedAcqWeightedMinimizer.suite());
	    suite.addTest(TestOrderedAcqSingleSubsetMinimizer.suite());
	    suite.addTest(TestScalarImage.suite());
	    suite.addTest(TestDynamicScalarImage.suite());
	    suite.addTest(TestSparseVectorImage.suite());
	    
	}


	// numerics 
	if (testNumerics) {
	    suite.addTest(TestACG_Distribution.suite());
	    suite.addTest(TestACG_Fitter.suite());
	    suite.addTest(TestBinghamDistribution.suite());
	    suite.addTest(TestBinghamFitter.suite());
	    suite.addTest(TestComplex.suite());
	    suite.addTest(TestEigenSystem3D.suite());
	    suite.addTest(GenTestMethods.suite());
	    suite.addTest(TestRealMatrix.suite());
	    suite.addTest(TestRotations.suite());
	    suite.addTest(TestSphericalHarmonics.suite());
	    suite.addTest(TestSphericalDistributionFitter.suite());
            suite.addTest(TestSymmetricMatrix.suite());
	    suite.addTest(TestTwoFibreACGFitter.suite());
	    suite.addTest(TestTwoFibreWatsonFitter.suite());
	    suite.addTest(TestTwoFibreFixedPropWatsonFitter.suite());
	    suite.addTest(TestTwoFibreBipolarWatsonFitter.suite());
	    suite.addTest(TestTwoFibreBinghamFitter.suite());
	    suite.addTest(TestVector3D.suite());
	    suite.addTest(TestWatsonDistribution.suite());
	    suite.addTest(TestWatsonFitter.suite());
	}
	
	// imaging
	if (testImaging) {
	    suite.addTest(TestAnalyzeHeader.suite());
	    suite.addTest(TestMetaImageHeader.suite());
	    suite.addTest(TestNifti1Dataset.suite());
	    suite.addTest(TestQ_VectorScheme.suite());
	    suite.addTest(TestB_VectorScheme.suite());
	    suite.addTest(TestRectGradSteTanScheme.suite());
	    suite.addTest(TestRectGradTRSE_Scheme.suite());
	}

	// tools
	if (testTools) {
	    suite.addTest(TestArrayOps.suite());
	}

        // sphfunc
        if (testSphFunc) {
            suite.addTest(TestEvenSHS.suite());
            suite.addTest(TestRBF_Sum.suite());
            suite.addTest(TestMaxEntProfile.suite());
            suite.addTest(TestTuchRBF.suite());
            suite.addTest(TestImagSH.suite());
            suite.addTest(TestRealSH.suite());
        }

	// tractography
	if (testTractography) {
 	    suite.addTest(TestConnectionProbabilityImage.suite());
	    suite.addTest(TestConnectivitySegmentedImage.suite());
	    suite.addTest(TestCubicVoxelRegion.suite());
	    suite.addTest(TestDT_LinearInterpolator.suite());
 	    suite.addTest(TestDT_LookupTableGenerator.suite());
	    suite.addTest(TestDT_NC_Interpolator.suite());
	    suite.addTest(TestDT_NN_Interpolator.suite());
            suite.addTest(TestDT_TractographyImage.suite());
            suite.addTest(TestEightNeighbourInterpolator.suite());
	    suite.addTest(TestEulerFibreTracker.suite());
	    suite.addTest(TestFACT_FibreTracker.suite());
 	    suite.addTest(TestLookupTable.suite());
	    suite.addTest(TestNearestNeighbourInterpolator.suite());
	    suite.addTest(TestNeighbourChoiceInterpolator.suite());
	    suite.addTest(TestOneTensorLUTGenerator.suite());
            suite.addTest(TestPICoACGRandomizer.suite());
	    suite.addTest(TestPICoBinghamRandomizer.suite());
            suite.addTest(TestPICoTractographyImage.suite());	    
	    suite.addTest(TestPICoWatsonRandomizer.suite());
	    suite.addTest(TestPointListROI.suite());
	    suite.addTest(TestRK4FibreTracker.suite());
            suite.addTest(TestSF_TractographyImage.suite());	    
 	    suite.addTest(TestStreamlineROI_Filter.suite());
 	    suite.addTest(TestStreamlineConnectivityGraph.suite());
 	    suite.addTest(TestTargetCP_Image.suite());
 	    suite.addTest(TestTract.suite());
 	    suite.addTest(TestTractCollection.suite());
 	    suite.addTest(TestTractSource.suite());
	    suite.addTest(TestTractStatisticFilter.suite());
	    suite.addTest(TestTractStatisticImage.suite());
 	    suite.addTest(TestTwoTensorLUTGenerator.suite());
	    suite.addTest(TestVectorLinearInterpolator.suite());
	    suite.addTest(TestVoxelList.suite());
	}
	
	if(testGreatCircleIntegrals){

	    suite.addTest(TestGreatCircleIntegrals.suite());
	}

	if(testModels){
	    suite.addTest(TestCylinderGPD.suite());
	    suite.addTest(TestBall.suite());
            suite.addTest(TestZeppelin.suite());
            suite.addTest(TestStick.suite());
            suite.addTest(TestTensor.suite());
            suite.addTest(TestSphereGPD.suite());
            suite.addTest(TestAstrocylinders.suite());
            suite.addTest(TestAstrosticks.suite());
            suite.addTest(TestGDRCylinders.suite());
            suite.addTest(TestDot.suite());
            
	}

	return suite;
    }

    // runs all tests in suite
    public static void main(String args[]) {
        
	// default is to run all tests
	
	if (args.length > 0) {
            
            testApps = false;
	    testData = false;
	    testImaging = false;
	    testInverters = false;
	    testMisc = false;
	    testNumerics = false;
	    testOptimizers = false;
            testSphFunc = false;
	    testTools = false;
	    testTractography = false;
	    testGreatCircleIntegrals = false;
	    testModels = false;
            
	    for (int i = 0; i < args.length; i++) {
		if (args[i].equals("-apps")) {
		    testApps = true;
		}
		if (args[i].equals("-data")) {
		    testData = true;
		}
		if (args[i].equals("-imaging")) {
		    testImaging = true;
		}
		if (args[i].equals("-inverters")) {
		    testInverters = true;
		}
		if (args[i].equals("-misc")) {
		    testMisc = true;
		}
		if (args[i].equals("-numerics")) {
		    testNumerics = true;
		}
		if (args[i].equals("-optimizers")) {
		    testNumerics = true;
		}
		if (args[i].equals("-tools")) {
		    testTools = true;
		}
		if (args[i].equals("-sphfunc")) {
		    testSphFunc = true;
		}
		if (args[i].equals("-tractography")) {
		    testTractography = true;
		}
		if(args[i].equals("-greatcircles")) {
		    testGreatCircleIntegrals = true;
		}
		if(args[i].equals("-models")) {
		    testModels = true;
		}

	    }

	}
	
	
	
        junit.textui.TestRunner.run(suite());
    }
}

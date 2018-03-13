package tools;

import java.util.Random;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.logging.Logger;
import java.io.*;

import apps.SphFuncPD_Stats;
import numerics.*;
import misc.*;
import models.compartments.CompartmentFactory;
import models.compartments.CompartmentModel;
import models.compartments.CompartmentType;
import imaging.*;
import data.*;
import simulation.DiffusionSimulation;
import simulation.SimulationParams;
import simulation.dynamics.StepGeneratorFactory;
import simulation.dynamics.StepGeneratorFactory.StepType;
import simulation.geometry.elements.CylinderFactory.CylType;
import simulation.geometry.substrates.ParallelCylinderSubstrate;
import simulation.geometry.substrates.SubstrateFactory.SubstrateType;
import sphfunc.*;
import inverters.*;
import mesd.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Parses command line arguments common to many of the Camino applications
 * and sets up the common objects.
 * 
 * <dt>Description:
 * 
 * <dd>Several Camino applications have common command line options. This class
 * reads them all, parses them and stored them as publicly accessible structure
 * fields.
 * 
 * The class also creates objects that are common to the Camino applications
 * using the command line options. For example, it sets up the input data
 * stream, the output stream, the DW_Scheme object.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 * 
 */
public class CL_Initializer {

	/**
	 * logging object
	 */
	private static Logger logger = Logger
			.getLogger("camino.tools.CL_Initializer");

	/**
	 * The number of values per voxel for the ball and stick model.
	 */
	public static final int BALLSTICKVALSPERVOXEL = 7;

	/**
	 * The number of values per voxel in single tensor data.
	 */
	public static final int ONETENVALSPERVOXEL = DT_Inversion.ITEMSPERVOX;

	/**
	 * The number of values per voxel in two-tensor data.
	 */
	public static final int TWOTENVALSPERVOXEL = TwoTensorInversion.ITEMSPERVOX;

	/**
	 * The number of values per voxel in three-tensor data.
	 */
	public static final int THREETENVALSPERVOXEL = ThreeTensorInversion.ITEMSPERVOX;

	/**
	 * Standard test function index.
	 */
	public static int testFunction = -1;

	/**
	 * Indicates whether or not the user specified a non-standard test function.
	 */
	public static boolean nonStandardTestFunction = false;

	/**
	 * Indicates whether to synthesise data using simulation
	 */
	public static boolean brownianSimulation = false;

	/**
	 * Indicates whether to synthesise data using a compartment model
	 */
	public static boolean compartmentModel = false;

	/**
	 * indicates if we are using our compartment model for fitting or not
	 */
	public static boolean compartmentModelForFitting = false;

	/**
	 * type of tissue model for fitting
	 */
	public static String fitModel = null;
	/**
	 * Start point for compartment model fitting from commandline
	 */
	public static boolean fixedStartPoint = false;

	/**
	 * type of noise model for fitting
	 */
	public static String noiseModel = "GAUSSIAN";

	/**
	 * type of fitting algorithm
	 */
	public static String fitAlgorithm = null;

	/**
	 * diffusion simulation parameters object
	 */
	public static SimulationParams simParams = null;

	/**
	 * simulation diffusivity
	 */
	public static double DIFF_CONST = 2.02E-9; // m^2 s^{-1}

	/**
	 * volume fraction
	 */
	public static double volFrac = 1;

	/**
	 * Dimensions of voxel grid
	 */
	public static int[] dataDims = { 0, 0, 0 };

	/**
	 * Voxel size in millimeters
	 */
	public static double[] voxelDims = { 0, 0, 0 };

	/**
	 * distance of centre of circle from centre of voxel grid
	 */
	public static double centreDist = 0;

	/**
	 * If a non-standard test function was specified, the parameters are read
	 * from the command line into this array.
	 */
	public static double[] testFunctionData = {};

	/**
	 * Index specifying the inversion.
	 */
	public static ModelIndex[] inversionIndices = new ModelIndex[] { ModelIndex.LDT };

	/**
	 * Index of the random rotation index to apply to the test function.
	 */
	public static int rotationIndex = 0;

	/**
	 * The number of voxels to process.
	 */
	public static int numVoxels = 0;

	/**
	 * Alternative setting for lambda1 in the standard test functions. This can
	 * be used to change the anisotropy. If left negative, defaults are used.
	 */
	public static double lambda1 = -1;

	/**
	 * Specifies the angle of rotation for StandardTestFunctions.dt2. If left
	 * zero, dt2 is not rotated.
	 */
	public static double dt2rotangle = 0;

	/**
	 * Specifies the mixing parameter for StandardTestFunctions.dt2. If left
	 * negative, the mixing parameter stays at 0.5.
	 */
	public static double dt2mix = -1;

	/**
	 * Alternative setting for scaling the standard test functions. If left
	 * negative the default scale is used.
	 */
	public static double scale = -1;

	// The parameters of a standard sequence with fixed b and
	// points from an electrostatic point set.

	/**
	 * number of elements in a voxel
	 */
	public static int numElements = -1;

	/**
	 * The number of zero measurements.
	 */
	public static int M = -1;

	/**
	 * The number of non-zero measurements.
	 */
	public static int N = -1;

	/**
	 * Fixed b-value for schemes created at run time.
	 */
	public static double bValue = 0.0;

	/**
	 * The signal to noise ratio in synthetic data.
	 */
	public static double SNR = -1;

	/**
	 * The number of bootstraps used in bootstrapping. The bootstrapper is only
	 * used if bootstrap > 0. If bootstrap == 1, the wild bootstrapper is used.
	 */
	public static int bootstrap = -1;

	/**
	 * Model used in wild bootstrapping.
	 * 
	 */
	public static String wildBootstrapModel = null;

	/**
	 * List of files containing repeated measurements of the same data for use
	 * in bootstrapping.
	 * 
	 */
	public static String[] bsDataFiles = null;

	/**
	 * Random number generator seed.
	 */
	public static int seed = 36558013;

	/**
	 * Name of file containing voxel classification. Used in classified model
	 * fit.
	 */
	public static String voxelClassMap = null;

	/**
	 * Maximum number of components in classified model fitting.
	 */
	public static int maxTensorComponents = 2;

	/**
	 * shape parameter for gamma distributed random numbers
	 */
	public static double gamma_k = 1.84037;

	/**
	 * scale parameter for gamma distribution
	 */
	public static double gamma_beta = 7.8E-7;

    /**
     * Intrinsic diffusivity in the MMWMD models.  If set
     * greater than zero it replaces the default values in
     * specific classes.
     */
    public static double mmwmddiff = -1;

	/**
	 * Array of indices of the models in classified model fitting. The default
	 * expects the output of VoxelClassify and does linear DT fitting where SH
	 * order is 0 or 2 and full two-tensor fitting where SH order is 4.
	 */
	public static ModelIndex[][] classifiedModelIndices = new ModelIndex[][] {
			{ ModelIndex.LDT }, { ModelIndex.LDT }, { ModelIndex.LDT },
			{ ModelIndex.POSPOS, ModelIndex.LDT } };

	/**
	 * Name of the input data file.
	 */
	public static String inputFile = null;

	/**
	 * Name of the background mask file.
	 */
	public static String bgMaskFile = null;
	
	/**
	 * Name of the gradient adjust file.
	 */
	public static String gradAdjFile = null;
	
	/**
	 * Name of the imaging scheme file.
	 */
	public static String schemeFile = null;

	/**
	 * The expected standard deviation of the noise in the data.
	 */
	public static double sigma = -1;

	/**
	 * Specifies the type of noise to add.
	 */
	public static String noiseType = "rician";

	/**
	 * Name of the outlier map for RESTORE.
	 */
	public static String outlierMapFile = null;

	/**
	 * Name of the noise variance map file for inversions that compute it.
	 * 
	 */
	public static String noiseVarianceMapFile = null;

	/**
	 * Name of the residual vector map file for inversions that compute it.
	 * 
	 */
	public static String residualVectorMapFile = null;

	/**
	 * Name of the matrix file for linear reconstructions.
	 */
	public static String matrixFile;

	/**
	 * Specifies whether or not the data should be normalized before linear
	 * reconstruction.
	 */
	public static boolean lrNormalize = false;

	/**
	 * Specifies whether or not the linear reconstruction should use the log
	 * measurements.
	 */
	public static boolean lrLog = false;

	/**
	 * The type of the input data
	 */
	public static String inputDataType = "float";

	/**
	 * The type of the background mask
	 */
	public static String bgMaskType = "short";

	/**
	 * The type of the gradient adjust 
	 */
	public static String gradAdjType = "float";

	/**
	 * The type of model represented by the input data. Normally this is null,
	 * but synthetic data can be generated from model parameters on the input.
	 */
	public static String inputModel = null;

	/**
	 * list of names for each compartment in compartment model
	 */
	public static String[] compartmentNames = null;

	/**
	 * list of parsed parameters for compartmentModel
	 */
	public static double[] compParams = null;

	/**
	 * The type of basis function used for linear reconstructions. The default
	 * setting is rbf.
	 */
	public static int basisType = BasisSumFactory.TUCH_RBF;

	/**
	 * Index of the ISCodes point set to use for computing spherical integrals.
	 */
	public static int pointSetInd = 3;

	/**
	 * Flag indicating that the IS codes point set index was set.
	 */
	public static boolean pointSetIndSet = false;

	/**
	 * Number of principal directions in the input or output stream.
	 */
	public static int numPDsIO = 3;

	/**
	 * The maximum order in the fitted spherical harmonic series.
	 */
	public static int maxOrder = 4;

	/**
	 * The f-test thresholds. If they are left at these default values, the
	 * program does not do the thresholding, but outputs the f-statistic
	 * probabilities instead.
	 */
	public static double f1 = -1;

	public static double f2 = -1;

	public static double f3 = -1;

	/**
	 * Specification of the maximum entropy deconvolution filter. The default
	 * filter is spike deconvolution with bd = 1.0.
	 */
	public static double[] kernelParams = { 0.0, 1.0 };

	/**
	 * Specifies the set of reconstruction directions for maximum entropy
	 * deconvolution if it is different from the direction set specified in the
	 * scheme file.
	 */
	public static int mePointSet = -1;

	/**
	 * Threshold on the q=0 measurements to detect background.
	 */
	public static double BACKGROUNDTHRESHOLD = 0.0;

	/**
	 * Threshold on the q=0 measurements to detect CSF.
	 */
	public static double CSFTHRESHOLD = -1.0;

	/**
	 * The parameters of the imaging sequence corresponding to the input data.
	 */
	public static DW_Scheme imPars = null;

	/**
	 * The test function to use for simulations.
	 */
	public static ModelPDF p = null;

	/**
	 * The data source to process.
	 */
	public static DataSource data = null;

	/**
	 * The background mask.
	 */
	public static DataSource bgMask = null;

	/**
	 * The gradient adjust.
	 */
	public static DataSource gradAdj = null;

	/**
	 * Flags to indicate parsing of arguments
	 */
	private static boolean[] parsed = null;

	/**
	 * Name of a file containing an affine image transformation.
	 */
	public static String transformFile = null;

	/**
	 * Names of files containing an image warp.
	 */
	public static String transformFileX = null;
	public static String transformFileY = null;
	public static String transformFileZ = null;


	/**
	 * Name of the reorientation strategy to use.
	 */
	public static String reorientation = "ppd";

	/**
	 * Burn-in period of Markov Chain Monte Carlo algorithm .
	 */
	public static int burnInTime = 1000;
	
	/**
	 * Number of samples produced by Markov Chain Monte Carlo algorithm .
	 */
	public static int samples = 100;
	
	/**
	 * Number of steps between samples in Markov Chain Monte Carlo algorithm .
	 */
	public static int interval = 100;
	
	/**
	 * Name of a file to output MCMC statistics.
	 */
	public static String mcmcStatsFile = null;


    /**
     * Header file defining the input / output physical space
     *
     */
    public static String headerTemplateFile;

    /**
     * Header defining the image space, data type, and other options for writing image data.
     *
     */
    public static ImageHeader headerTemplate;


	/**
	 * Empty constructor.
	 */
	private CL_Initializer() {

	}

	/**
	 * The constructor requires only the array of command line arguments.
	 * 
	 * @param args
	 *            The array of command line arguments to a DiffusionInversions
	 *            application.
	 */
	public static void CL_init(String[] args) {

		parsed = new boolean[args.length];

		for (int i = 0; i < parsed.length; i++) {
			parsed[i] = false;
		}

		boolean numPDsSet = false;
		boolean maxTensorComponentsSet = false;

		// Loop over the command line arguments and
		// store the required values in the appropriate variables.
		for (int i = 0; i < args.length; i++) {

			if (args[i].equals("-argfile"))
				try {

					// Read extra arguments from argfile and augment args array.
					args = augmentArgs(args, args[i + 1]);

					// Also need to augment the parsed array.
					boolean[] parsedNew = new boolean[args.length];
					for (int j = 0; j < i; j++) {
						parsedNew[j] = parsed[j];
					}
					for (int j = i; j < parsedNew.length; j++) {
						parsedNew[j] = false;
					}

					parsed = parsedNew;

					// Mark the argfile arguments as parsed.
					markAsParsed(i);
					markAsParsed(i + 1);

				} catch (Exception e) {
					throw new LoggedException(e.toString()
							+ "\nAttempting to continue...");
				}

			if (args[i].equals("-testfunc")) {
				// parse the arguments
				testFunction = Integer.parseInt(args[i + 1]);

				// flag and value both processed
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-gaussmix")) {
				nonStandardTestFunction = true;
				int components = Integer.parseInt(args[i + 1]);
				testFunctionData = new double[1 + components * 7];
				testFunctionData[0] = (double) components;

				// flag parsed
				markAsParsed(i);
				// components flag parsed
				markAsParsed(i + 1);

				for (int j = 0; j < components * 7; j++)
					try {
						testFunctionData[j + 1] = Double.parseDouble(args[i + j
								+ 2]);
						parsed[i + j + 2] = true;
					} catch (Exception e) {

						throw new LoggedException(
								"Enter {Dxx, Dxy, Dxz, Dyy, Dyz, Dzz, mix} for each component diffusion tensor.");
					}
			}
			/*
			 * if(args[i].equals("-delta")){ throw new LoggedException("-delta
			 * option disabled. Needs fixing.");
			 * //SchemeV0.setDelta(Double.parseDouble(args[i+1]));
			 * //markAsParsed(i, 2); }
			 */
			if (args[i].equals("-rotation")) {
				rotationIndex = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-inversion")) {
				inversionIndices = ModelIndex
						.getIndices(new String[] { args[i + 1] });
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equalsIgnoreCase("-fitmodel")) {
				fitModel = args[i + 1];
				compartmentModel = true;
				compartmentModelForFitting = true;
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-noiseModel")) {
				noiseModel = args[i + 1];
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-fitalgorithm")) {
				fitAlgorithm = args[i + 1];
				markAsParsed(i, 2);
			}
			if (args[i].equals("-model")) {
				/*
				 * if(args[i+1].equalsIgnoreCase("compartment")){ // parse the
				 * compartment model. this option is used by // modelfit. but
				 * not datasynth initCompartmentModel(args, i+1);
				 * 
				 * markAsParsed(i, 2);
				 *  } else{
				 */
				// use the old inversion indices path
				String[] tmp = new String[1000];
				int t = 0;

				while (t + i + 1 < args.length
						&& !args[i + t + 1].startsWith("-")) {
					tmp[t] = args[i + t + 1];
					t++;
				}

				String[] inversionStrings = new String[t];

				System.arraycopy(tmp, 0, inversionStrings, 0, t);

				inversionIndices = ModelIndex.getIndices(inversionStrings);

				markAsParsed(i, t + 1);
				// }

			}
			if (args[i].equals("-bootstrap")) {
				bootstrap = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-seed")) {
				seed = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-voxels")) {
				numVoxels = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-bmx")) {
				logger
						.warning("This option has been deprecated and is no longer functional. Use -schemefile.");
				markAsParsed(i);
				markAsParsed(i + 1);

			}
			if (args[i].equals("-voxelelements")) {
				numElements = Integer.parseInt(args[i + 1]);
				markAsParsed(i, 2);
			}
			if (args[i].equals("-fixedmodq")) {
				M = Integer.parseInt(args[i + 1]);
				N = Integer.parseInt(args[i + 2]);

				double modQ = Double.parseDouble(args[i + 3]);
				double diffusionTime = Double.parseDouble(args[i + 4]);

				bValue = modQ * modQ * diffusionTime;

				markAsParsed(i);
				markAsParsed(i + 1);
				markAsParsed(i + 2);
				markAsParsed(i + 3);
				markAsParsed(i + 4);
			}
			if (args[i].equals("-fixedbvalue")) {
				M = Integer.parseInt(args[i + 1]);
				N = Integer.parseInt(args[i + 2]);

				bValue = Double.parseDouble(args[i + 3]);

				markAsParsed(i);
				markAsParsed(i + 1);
				markAsParsed(i + 2);
				markAsParsed(i + 3);
			}
			if (args[i].equals("-qscale")) {
				logger
						.warning("This option has been deprecated and is no longer functional.");
				markAsParsed(i, 2);
			}
			if (args[i].equals("-tau")) {
				logger
						.warning("This option has been deprecated and is no longer functional.");
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-tauscale")) {
				logger
						.warning("This option has been deprecated and is no longer functional.");
				markAsParsed(i, 2);
			}
			if (args[i].equals("-lambda1")) {
				lambda1 = Double.parseDouble(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-scale")) {
				scale = Double.parseDouble(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-dt2rotangle")) {
				dt2rotangle = Double.parseDouble(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-dt2mix")) {
				dt2mix = Double.parseDouble(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-snr")) {
				SNR = Double.parseDouble(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-outputfile")) {
				OutputManager.outputFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-gzip")) {
				OutputManager.gzipOut = true;
				markAsParsed(i);
			}
			if (args[i].equals("-inputfile")) {
				inputFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);

			}
			if (args[i].equals("-bsdatafiles")) {
				String[] tmp = new String[1000];
				int t = 0;
				while (t + i + 1 < args.length
						&& !args[t + i + 1].startsWith("-")) {
					tmp[t] = args[i + t + 1];
					t++;
				}
				bsDataFiles = new String[t];

				bootstrap = t;

				System.arraycopy(tmp, 0, bsDataFiles, 0, t);
				markAsParsed(i, t + 1);

			}
			if (args[i].equals("-wildbsmodel")) {
				wildBootstrapModel = args[i + 1];
				bootstrap = 1;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-bgmask") || args[i].equals("-brainmask")) {

                            // (PAC) the rest of the world calls it a brain mask, we should
                            // start doing the same

				bgMaskFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);

			}
			if (args[i].equals("-gradadj")) {

				gradAdjFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);

			}
			if (args[i].equals("-schemefile")) {
				schemeFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-voxclassmap")) {
				voxelClassMap = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-outliermap")) {
				outlierMapFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-noisemap")) {
				noiseVarianceMapFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-residualmap")) {
				residualVectorMapFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-matrix")) {
				matrixFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-normalize")) {
				lrNormalize = true;
				markAsParsed(i);
			}
			if (args[i].equals("-log")) {
				lrLog = true;
				markAsParsed(i);
			}
			if (args[i].equals("-sigma")) {
				sigma = Double.parseDouble(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-noisetype")) {
				// better to do lower case conversion here than at every call to
				// nextVoxel() in the DataSynthesizer
				noiseType = args[i + 1].toLowerCase();
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-inputdatatype")) {
				inputDataType = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-outputdatatype")) {
				OutputManager.outputDataType = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-maskdatatype")) {
				bgMaskType = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-gradadjdatatype")) {
				gradAdjType = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-inputmodel")) {
				inputModel = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-synthmodel")) {
				// parse the compartment model. this option is used by
				// datasynth but not modelfit
				initCompartmentModel(args, i + 1);

				compartmentModel = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-startpoint")) {
				// parse the compartment model. this option is used by
				// modelfit
				initCompartmentModel(args, i + 1);

				fixedStartPoint = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-basistype")) {
				if (args[i + 1].equals("rbf")) {
					basisType = BasisSumFactory.TUCH_RBF;
				} else if (args[i + 1].equals("sh")) {
					basisType = BasisSumFactory.SPHERICAL_HARMONICS;
				} else {
					throw new LoggedException("Unrecognized basis type: "
							+ args[i + 1]);
				}
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-order")) {
				maxOrder = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-fastmesd")) {
				MESD_Inversion.doNumIntTests = false;
				markAsParsed(i);
			}
			if (args[i].equals("-basepointset")) {
				MESD_Inversion.basePointSet = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-pointset")) {
				pointSetInd = Integer.parseInt(args[i + 1]);
				pointSetIndSet = true;
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-maxcomponents")) {
				maxTensorComponents = Integer.parseInt(args[i + 1]);
				maxTensorComponentsSet = true;
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-classifiedmodels")) {
				int numModels = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);

				String[] inversionStrings = new String[numModels];
				for (int j = 0; j < numModels; j++) {
					inversionStrings[j] = args[i + 2 + j];
					markAsParsed(i + 2 + j);
				}

				classifiedModelIndices = ModelIndex
						.getClassifiedModelIndices(inversionStrings);

			}
			if (args[i].equals("-rbfpointset")) {
				RBF_Sum.setPoints(SphericalPoints.getElecPointSet(Integer
						.parseInt(args[i + 1])));
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-rbfsigma")) {
				TuchRBF_Sum.setSigma(Double.parseDouble(args[i + 1]));
				;
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-filter")) {
				markAsParsed(i);
				if (args[i + 1].equals("SPIKE")) {
					kernelParams = new double[2];
					kernelParams[0] = 0.0;
					kernelParams[1] = Double.parseDouble(args[i + 2]);
					markAsParsed(i + 1);
					markAsParsed(i + 2);
				} else if (args[i + 1].equals("PAS")) {
					kernelParams = new double[2];
					kernelParams[0] = 1.0;
					kernelParams[1] = Double.parseDouble(args[i + 2]);
					markAsParsed(i + 1);
					markAsParsed(i + 2);
				} else if (args[i + 1].equals("TENSOR")) {
					kernelParams[0] = 2.0;
					kernelParams[1] = Double.parseDouble(args[i + 2]);
					markAsParsed(i + 1);
					markAsParsed(i + 2);
				} else {
					throw new LoggedException("Unrecognized filter type: "
							+ args[i + 1]);
				}
			}
			if (args[i].equals("-mepointset")) {
				mePointSet = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-ftest")) {
				f1 = Double.parseDouble(args[i + 1]);
				f2 = Double.parseDouble(args[i + 2]);
				f3 = Double.parseDouble(args[i + 3]);
				markAsParsed(i);
				markAsParsed(i + 1);
				markAsParsed(i + 2);
				markAsParsed(i + 3);
			}
			if (args[i].equals("-bgthresh")) {
				BACKGROUNDTHRESHOLD = Double.parseDouble(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-csfthresh")) {
				CSFTHRESHOLD = Double.parseDouble(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-numpds")) {
				numPDsIO = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
				numPDsSet = true;
			}
			if (args[i].equals("-searchradius")) {
				SphericalFunction.setSearchRadius(Double
						.parseDouble(args[i + 1]));
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-paramscale")) {
				SphericalFunction.SFPD_PARAMSCALE = Double
						.parseDouble(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-v")) {
				DB.setVerbosity(Integer.parseInt(args[i + 1]));
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-trans")) {
				transformFile = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-xtrans")) {
				transformFileX = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-ytrans")) {
				transformFileY = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-ztrans")) {
				transformFileZ = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-reorientation")) {
				reorientation = args[i + 1];
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-diffusivity")) {
				DIFF_CONST = Double.parseDouble(args[i + 1]);
				markAsParsed(i, 2);
			}
			if (args[i].equals("-mmwmddiff")) {
				mmwmddiff = Double.parseDouble(args[i + 1]);
				markAsParsed(i, 2);
			}
			if (args[i].equals("-gamma")) {
				gamma_k = Double.parseDouble(args[i + 1]);
				gamma_beta = Double.parseDouble(args[i + 2]);
				markAsParsed(i, 3);
			}
			if (args[i].equalsIgnoreCase("-drawcrosssection")) {
				SimulationParams.sim_drawCrossSection = true;
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			// simulation options (all set brownianSimulation to true)
			if (args[i].equals("-walkers")) {
				SimulationParams.sim_N_walkers = Integer.parseInt(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-tmax")) {
				SimulationParams.sim_tmax = Integer.parseInt(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-p")) {
				SimulationParams.sim_p = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-pstick")) {
				SimulationParams.sim_p_stick = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-punstick")) {
				SimulationParams.sim_p_unstick = Double
						.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-voxelsizefrac")) {
				SimulationParams.sim_voxelSizeFrac = Double
						.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-crossangle")) {
				SimulationParams.sim_cAngle = Double.parseDouble(args[i + 1]);
				markAsParsed(i, 2);
			}
			if (args[i].equals("-D1")) {
				SimulationParams.sim_cyl_D1 = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-D2")) {
				SimulationParams.sim_cyl_D2 = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-compartmentsignal")) {
				if (args[i + 1].equalsIgnoreCase("EXTRAONLY")) {
					SimulationParams.sim_compartmentSignal = SimulationParams.EXTRAONLY;
					brownianSimulation = true;
					markAsParsed(i, 2);
				}
				if (args[i + 1].equalsIgnoreCase("INTRAONLY")) {
					SimulationParams.sim_compartmentSignal = SimulationParams.INTRAONLY;
					brownianSimulation = true;
					markAsParsed(i, 2);
				}
				if (args[i + 1].equalsIgnoreCase("ALL")) {
					SimulationParams.sim_compartmentSignal = SimulationParams.ALLCOMPS;
					brownianSimulation = true;
					markAsParsed(i, 2);
				}
			}
			if (args[i].equals("-initial")) {
				if (args[i + 1].equalsIgnoreCase("uniform")) {
					SimulationParams.sim_initial = SimulationParams.UNIFORM;
					markAsParsed(i + 1);
				} else if ((args[i + 1].equalsIgnoreCase("spike"))
						|| (args[i + 1].equalsIgnoreCase("delta"))) {
					SimulationParams.sim_initial = SimulationParams.SPIKE;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("special")) {
					SimulationParams.sim_initial = SimulationParams.SPECIAL;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("intra")) {
					SimulationParams.sim_initial = SimulationParams.INTRACELLULAR;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("extra")) {
					SimulationParams.sim_initial = SimulationParams.EXTRACELLULAR;
					markAsParsed(i + 1);
				} else {
					logger.warning("simulation initial conditions '"
							+ args[i + 1]
							+ "' default value (UNIFORM) will be used");
				}
				brownianSimulation = true;
				markAsParsed(i);
			}
			if (args[i].equalsIgnoreCase("-duration")) {
				SimulationParams.duration = Double.parseDouble(args[i + 1]);
				SimulationParams.trajectories = true;
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-trajfile")) {
				SimulationParams.trajFile = args[i + 1];
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-separateruns")) {
				SimulationParams.sim_separate_runs = true;
				brownianSimulation = true;
				markAsParsed(i);
			}
			if (args[i].equalsIgnoreCase("-facets")) {
				SimulationParams.sim_num_facets = Integer.parseInt(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if ((args[i].equals("-geometry")) || (args[i].equals("-substrate"))) {
				if (args[i + 1].equalsIgnoreCase("cell-striped")) {
					SimulationParams.sim_geomType = SubstrateType.CELL_STRIPED;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("cell-perc")) {
					SimulationParams.sim_geomType = SubstrateType.CELL_PERC;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("cylinder")) {
					SimulationParams.sim_geomType = SubstrateType.CYL_1_FIXED;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("inflammation")) {
					SimulationParams.sim_geomType = SubstrateType.CYL_1_INFLAM;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("facetted")) {
					SimulationParams.sim_geomType = SubstrateType.CYL_1_FACET;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("ply")) {
					SimulationParams.sim_geomType = SubstrateType.TRI_PLY_MESH;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("percolation")) {
					SimulationParams.sim_geomType = SubstrateType.CYL_1_PERC;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("myelinated")) {
					SimulationParams.sim_geomType = SubstrateType.CYL_1_MYELIN;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("crossing")) {
					SimulationParams.sim_geomType = SubstrateType.CYL_2_FIXED;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("distribcylinder")) {
					SimulationParams.sim_geomType = SubstrateType.CYL_1_DISTRIB;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("stickycylinder")) {
					SimulationParams.sim_geomType = SubstrateType.CYL_1_STICKY;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("regsphere")){
					SimulationParams.sim_geomType = SubstrateType.SPH_1_pC;
					markAsParsed(i+1);
				} else if (args[i + 1].equalsIgnoreCase("empty")) {
					SimulationParams.sim_geomType = SubstrateType.EMPTY;
					markAsParsed(i + 1);
				} else {
					logger.warning("unrecognised geometry type '" + args[i + 1]
							+ "' will use default ("
							+ SimulationParams.sim_geomType + ")");
				}
				brownianSimulation = true;
				markAsParsed(i);
			}
			if (args[i].equalsIgnoreCase("-substrateinfo")) {
				SimulationParams.substrateInfo = true;
				markAsParsed(i);
			}
			if (args[i].equals("-packing")) {
				if (args[i + 1].equalsIgnoreCase("square")) {
					SimulationParams.sim_cyl_pack = ParallelCylinderSubstrate.SQUARE;
					markAsParsed(i + 1);
				} else if (args[i + 1].equalsIgnoreCase("hex")) {
					SimulationParams.sim_cyl_pack = ParallelCylinderSubstrate.HEX;
					markAsParsed(i + 1);
				} else {
					logger
							.warning("cylinder packing '"
									+ args[i + 1]
									+ "' not supported. will default to hexagonal packing");
				}
				brownianSimulation = true;
				markAsParsed(i);
			}
			if (args[i].equals("-numcylinders")) {
				SimulationParams.sim_cyl_dist_size = Integer
						.parseInt(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-cylinderminrad")) {
				SimulationParams.sim_cyl_min_r = Double
						.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-cylindermaxrad")) {
				SimulationParams.sim_cyl_max_r = Double
						.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-cylinderrad")||args[i].equalsIgnoreCase("-radius")) {
				SimulationParams.sim_r = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-cylinderinnerrad")) {
				SimulationParams.sim_cyl_r1 = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-cylindersep")||args[i].equalsIgnoreCase("-separation")) {
				SimulationParams.sim_R = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if(args[i].equals("-meshsep")){
				double[] meshSep= new double[DiffusionSimulation.D];
				
				for(int j=0; j<meshSep.length; j++){
					meshSep[j]= Double.parseDouble(args[i+j+1]);
				}
				
				SimulationParams.sim_mesh_sep= meshSep;
				brownianSimulation= true;
				markAsParsed(i, 4);
				
			}
			if(args[i].equals("-cylindertype")){
				// these are the only ones anyone actually uses...
				if(args[i+1].equalsIgnoreCase("basic")){
					SimulationParams.cylinderType= CylType.BASIC;
				}
				if(args[i+1].equalsIgnoreCase("nested")){
					SimulationParams.cylinderType= CylType.NESTED;
				}
				
				brownianSimulation= true;
				markAsParsed(i, 2);
			}
			if(args[i].equals("-nestdepth")){
				SimulationParams.cyl_nest_depth_max=Integer.parseInt(args[i+1]);
				
				brownianSimulation=true;
				markAsParsed(i,2);
			}
			if (args[i].equals("-increments")) {
				SimulationParams.sim_inflamm_increments = Integer
						.parseInt(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if(args[i].equals("-cylfile")){
				SimulationParams.sim_cylFile= new String(args[i+1]);
				brownianSimulation= true;
				markAsParsed(i,  2);
			}
			if (args[i].equalsIgnoreCase("-onlyrun")) {
				SimulationParams.sim_onlyRun = Integer.parseInt(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-voxelsize")) {
				SimulationParams.sim_voxelSize = Double
						.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-plyfile")) {
				SimulationParams.sim_plyfile = new String(args[i + 1]);
				String suffix = SimulationParams.sim_plyfile
						.substring(SimulationParams.sim_plyfile.length() - 4);
				if (!suffix.equalsIgnoreCase(".ply")) {
					SimulationParams.sim_plyfile += ".ply";
				}
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-steptype")) {
				if (args[i + 1].equals("fixedlength")) {
					SimulationParams.sim_stepType = StepType.FIXEDLENGTH;
					markAsParsed(i + 1);
				} else {
					logger.warning("step type '" + args[i + 1]
							+ "' unknown. default "
							+ "(fixedlegnth) generator will be used");
				}
				brownianSimulation = true;
				markAsParsed(i);
			}
			if (args[i].equals("-substratesize")) {
				SimulationParams.sim_L = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-latticesize")) {
				SimulationParams.sim_L = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-cellsize")) {
				SimulationParams.sim_l = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-stripethickness")) {
				SimulationParams.sim_stripethickness = Integer
						.parseInt(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equals("-pperc")) {
				SimulationParams.sim_p_perc = Double.parseDouble(args[i + 1]);
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			if (args[i].equalsIgnoreCase("-statsfile")) {
				SimulationParams.sim_statsfile = args[i + 1];
				brownianSimulation = true;
				markAsParsed(i, 2);
			}
			// end of simulation options
			if (args[i].equals("-volfrac")) {
				volFrac = Double.parseDouble(args[i + 1]);
				markAsParsed(i, 2);
			}
			if (args[i].equals("-gridsize") || args[i].equals("-datadims")) {
				dataDims[0] = Integer.parseInt(args[i + 1]);
				dataDims[1] = Integer.parseInt(args[i + 2]);
				dataDims[2] = Integer.parseInt(args[i + 3]);
				markAsParsed(i, 4);
			}
			if (args[i].equals("-voxeldims")) {
				voxelDims[0] = Double.parseDouble(args[i + 1]);
				voxelDims[1] = Double.parseDouble(args[i + 2]);
				voxelDims[2] = Double.parseDouble(args[i + 3]);
				markAsParsed(i, 4);

			}
			if (args[i].equals("-centredist")) {
				centreDist = Double.parseDouble(args[i + 1]);
				markAsParsed(i, 2);
			}
			if (args[i].equals("-burnin")) {
				burnInTime = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-samples")) {
				samples = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equals("-interval")) {
				interval = Integer.parseInt(args[i + 1]);
				markAsParsed(i);
				markAsParsed(i + 1);
			}
			if (args[i].equalsIgnoreCase("-mcmcstatsfile")){
				mcmcStatsFile = args[i+1];
				markAsParsed(i, 2);
			}			
                        if (args[i].equals("-header")) {

                            headerTemplateFile = args[i+1];

                            markAsParsed(i,2);
                        }

		
		}
                
                // if someone has done -outputfile something.gz, assume -gzip
                if (OutputManager.outputFile != null && OutputManager.outputFile.endsWith(".gz")) {
                    OutputManager.gzipOut = true;
                }

                
                // set numPDsIO
		if (!numPDsSet && inputModel != null) {

			if (inputModel.endsWith("dt")) {
				numPDsIO = 1;
			} else if (inputModel.equals("twotensor")) {
				numPDsIO = 2;
			} else if (inputModel.equals("threetensor")) {
				numPDsIO = 3;
			} else if (inputModel.equals("pico")) {
				numPDsIO = 1;
			} else if (inputModel.equals("ballstick")) {
				numPDsIO = 1;
			} else if (inputModel.startsWith("bayesdirac")) {
				numPDsIO = 1;
			} else if (inputModel.endsWith("multitensor")) {
				numPDsIO = maxTensorComponents;
			}
		}
		else if (numPDsSet && inputModel != null) {

			boolean wrongNumberOfPDs = false;

			if (inputModel.equals("dt")) {
				if (numPDsIO != 1) {
					wrongNumberOfPDs = true;
				}
			} else if (inputModel.equals("twotensor")) {
				if (numPDsIO != 2) {
					wrongNumberOfPDs = true;
				}
			} else if (inputModel.equals("threetensor")) {
				if (numPDsIO != 3) {
					wrongNumberOfPDs = true;
				}
			} else if (inputModel.equals("ballstick")) {
				if (numPDsIO != 1) {
					wrongNumberOfPDs = true;
				}
			}

			if (wrongNumberOfPDs) {
				throw new LoggedException("Wrong number of PDs (" + numPDsIO
						+ ") specified for input model " + inputModel);
			}

		}
                
		// now do maxTensorComponents
		if (!maxTensorComponentsSet && inputModel != null) {

			// (PAC) Better to tell people to use multitensor with
			// -maxcomponents, rather than
			// twotensor or threetensor. track and dteig no longer recommend
			// twotensor and threetensor.

			if (inputModel.endsWith("dt")) {
				maxTensorComponents = 1;
			} else if (inputModel.equals("twotensor")) {
				maxTensorComponents = 2;
			} else if (inputModel.equals("threetensor")) {
				maxTensorComponents = 3;
			}
		} else if (maxTensorComponentsSet && inputModel != null) {

			boolean wrongComponents = false;

			if (inputModel.equals("dt")) {
				if (maxTensorComponents != 1) {
					wrongComponents = true;
				}
			} else if (inputModel.equals("twotensor")) {
				if (maxTensorComponents != 2) {
					wrongComponents = true;
				}
			} else if (inputModel.equals("threetensor")) {
				if (maxTensorComponents != 3) {
					wrongComponents = true;
				}
			}

			if (wrongComponents) {
				throw new LoggedException("Wrong number of tensor components ("
						+ maxTensorComponents + ") specified for input model "
						+ inputModel);
			}

		}
		/*
		 * if (brownianSimulation) { throw new LoggedException("Simulation is
		 * broken :-( "); }
		 */

	}

	/**
	 * Checks the arguments and logs any args that have not been marked as
	 * parsed.
	 */
	public static final void checkParsing(String[] args) {

		if (parsed == null) {
			RuntimeException e = new RuntimeException(
					"attempt to check command line argument parsing made"
							+ " before parsing has been attempted!");

			StackTraceElement[] stackTrace = e.getStackTrace();
			String stString = new String();

			// log the exceptions message
			logger.severe(e.toString());

			// log the stack trace
			for (int ii = 0; ii < stackTrace.length; ii++) {
				stString += stackTrace[ii] + "\n";
			}
			logger.severe(stString);

			throw e;
		}

		for (int i = 0; i < parsed.length; i++) {
			if (!parsed[i]) {
				logger.warning("WARNING: couldn't parse arg " + i + ": '"
						+ args[i] + "'");
			}
		}
	}

	/**
	 * Augments the command-line arguments array with arguments read from a
	 * config file.
	 * 
	 * @param args
	 *            The current args array
	 * 
	 * @param filename
	 *            The name of the config file.
	 * 
	 * @return The augmented args array.
	 */
	public static String[] augmentArgs(String[] args, String filename)
			throws FileNotFoundException {

		// Read the file.
		Scanner s = new Scanner(new File(filename));
		Vector<String> newArgs = new Vector<String>();
		while (s.hasNext()) {
			newArgs.add(s.next());
		}
		s.close();

		// Create the augmented args array.
		String[] augArgs = new String[args.length + newArgs.size()];
		for (int i = 0; i < args.length; i++) {
			augArgs[i] = args[i];
		}
		for (int i = 0; i < newArgs.size(); i++) {
			augArgs[args.length + i] = newArgs.elementAt(i);
		}

		return augArgs;
	}

	/**
	 * Initializes the DW_Scheme object containing the acqusition information.
	 */
	public static void initImagingScheme() {

		// Sort out the parameters of the imaging sequence.
		if (M >= 0) {
			imPars = B_VectorScheme.elecPointSetScheme(M, N, bValue);
		} else if (schemeFile != null) {
			imPars = DW_Scheme.readScheme(schemeFile);
		} else {
			if (!brownianSimulation) {
				throw new LoggedException(
						"Tried to initialize imaging scheme, but no scheme parameters were specified");
			}
		}
	}

	/**
	 * Initializes the reconstruction direction set for maximum entropy
	 * deconvolution.
	 */
	public static void initMaxEnt() {

		if (mePointSet > 0) {
			// Set the reconstruction directions from the specified
			// electrostatic point set.
			MaxEntProfile.setReconDirs(SphericalPoints
					.getElecPointSet(mePointSet));

		} else {
			// Set the central points in the maximum entry profile to the
			// set of gradient directions in the schemefile.
			// (PAC) Now setting them to actual gradient directions
			MaxEntProfile.setReconDirs(imPars.getNonZeroG_Dirs());
		}

	}

	/**
	 * Specifies a specific data source object.
	 * 
	 * @param newDataSource
	 *            The data source.
	 */
	public static void setDataSource(DataSource newDataSource) {
		data = newDataSource;
	}

	/**
	 * Specifies a specific data source object for the background mask.
	 * 
	 * @param newMaskSource
	 *            The data source.
	 */
	public static void setMaskSource(DataSource newMaskSource) {
		bgMask = newMaskSource;
	}

	/**
	 * Specifies a specific data source object for the gradient adjust.
	 * 
	 * @param newGradAdjSource
	 *            The data source.
	 */
	public static void setGradAdjSource(DataSource newGradAdjSource) {
		gradAdj = newGradAdjSource;
	}

	/**
	 * Initializes an external data source for spherical function parameters.
	 */
	public static void initSphFuncDataSource() {
		if (inputModel == null) {

                    // If no model was specified for the input data, the input data
                    // is diffusion-weighted MRI measurements.
                    
                    data = ExternalDataSource.getDataSource(inputFile, imPars.numMeasurements(), inputDataType);
		} else if (inputModel.equals("dt")) {
                    data = ExternalDataSource.getDataSource(inputFile, ONETENVALSPERVOXEL, inputDataType);
		} else if (inputModel.equals("ballstick")) {
			data = ExternalDataSource.getDataSource(inputFile, BALLSTICKVALSPERVOXEL, inputDataType);
		} else if (inputModel.equals("sh")) {
			data = ExternalDataSource.getDataSource(inputFile, SphericalHarmonics.evenFuncsUpTo(maxOrder) + 2, inputDataType);
		} else if (inputModel.equals("rbf")) {
			data = ExternalDataSource.getDataSource(inputFile, RBF_Sum.numPoints() + 2, inputDataType);
		} else if (inputModel.equals("maxent")) {
			data = ExternalDataSource.getDataSource(inputFile, MaxEntProfile.numLambdas() + 2, inputDataType);
		} else {
			RuntimeException e = new RuntimeException(
					"Unknown input model type: " + inputModel);

			StackTraceElement[] stackTrace = e.getStackTrace();
			String stString = new String();

			// log the exceptions message
			logger.severe(e.toString());

			// log the stack trace
			for (int ii = 0; ii < stackTrace.length; ii++) {
				stString += stackTrace[ii] + "\n";
			}
			logger.severe(stString);

			throw e;
		}

		initMaskSource();
        initGradAdjSource();
	}

	/**
	 * Initializes an external data source for tensor data
	 */
	public static void initTensorDataSource() {

            if (inputModel.equals("dt")) {
                data = ExternalDataSource.getDataSource(inputFile, ONETENVALSPERVOXEL, inputDataType);
            } else if (inputModel.equals("twotensor")) {
                data = ExternalDataSource.getDataSource(inputFile, TWOTENVALSPERVOXEL, inputDataType);
            } else if (inputModel.equals("threetensor")) {
                data = ExternalDataSource.getDataSource(inputFile, THREETENVALSPERVOXEL, inputDataType);
            } else if (inputModel.equals("multitensor")) {
                data = ExternalDataSource.getDataSource(inputFile, 3 + 7 * maxTensorComponents, inputDataType);
            } else {
                throw new LoggedException(
                                          "Cannot process non-tensor data from model: " + inputModel);
            }

            initMaskSource();
            initGradAdjSource();
	}


	/**
	 * Initializes the data source for the background mask.
	 */
	public static void initMaskSource() {
            if (bgMaskFile != null) {
                bgMask = ExternalDataSource.getDataSource(bgMaskFile, 1, bgMaskType);
            }
	}

	/**
	 * Initializes the data source for the gradient adjust.
	 */
	public static void initGradAdjSource() {
            if (gradAdjFile != null) {
                gradAdj = ExternalDataSource.getDataSource(gradAdjFile, 9, gradAdjType);
            }
	}


	/**
	 * parse the command line arguments for compartment models into a concrete
	 * instance. Argument should be formatted as follows: -model compartment
	 * <num compartments> <comp1 name> <comp1 params> <comp2 name> <comp2
	 * params> ...
	 * 
	 * the name of the compartment should be one of those from
	 * data.compartments.CompartmentType (@see data.compartment.CompartmentType)
	 * and is not case sensitive. The number of parameters will vary from
	 * compartment to compartment. params should be specified in the order in
	 * which they are expected by the compartment.
	 * 
	 * This method will mark as parsed everything in the model specification
	 * after the "compartment" argument, provided it is formatted correctly.
	 * Exceptions will be thrown if: o volume fractions sum to greater than 1 o
	 * a negative volume fraction is specified o the wrong number of paramters
	 * is specified for a compartment
	 * 
	 * @param args
	 *            the list of command line args
	 * @param i
	 *            the index in the arg list of the "compartment" switch
	 * 
	 * @return a new compartment model object
	 * 
	 */
	private static void initCompartmentModel(String[] args, int i) {

		int numComps = Integer.parseInt(args[i + 1]);

		markAsParsed(i + 1);

		String[] comps = new String[numComps];
		CompartmentType[] cTypes = new CompartmentType[numComps];

		double[][] params = new double[numComps][];

		int counter = 2;

		double totalVolFrac = 0.0;

		// space for vol fracs and number of compartments
		int totalNumParams = numComps + 1;

		for (int j = 0; j < numComps; j++) {
			comps[j] = args[i + counter];
			CompartmentType compType = CompartmentType
					.getCompartmentType(comps[j]);
			cTypes[j] = compType;
			int numParams = compType.numParams;
			params[j] = new double[numParams + 1];
			markAsParsed(i + counter);
			counter++;

			// count space for parameters
			totalNumParams += numParams;

			if (j < numComps - 1) {
				params[j][0] = Double.parseDouble(args[i + counter]);
				if (params[j][0] < 0.0) {
					throw new LoggedException(
							"negative volume fraction specified in compartment model initialisation");
				}
				totalVolFrac += params[j][0];
				markAsParsed(i + counter);
				counter++;
			} else {
				if (totalVolFrac < 1.0) {
					params[j][0] = 1.0 - totalVolFrac;
				} else {
					throw new LoggedException(
							"specified volume fractions in compartment model sum to greater than 1.");
				}
			}

			for (int k = 0; k < numParams; k++) {
				params[j][k + 1] = Double.parseDouble(args[i + counter]);
				markAsParsed(i + counter);
				counter++;
			}
		}

		// let's see how we've done...
		logger.info("Compartment model specified");
		logger.info("There are " + numComps + " compartments:");
		String compString = new String();
		for (int j = 0; j < comps.length; j++) {
			compString += ("\n\t" + cTypes[j] + ": vol frac= " + params[j][0] + " params = ");
			for (int k = 1; k < params[j].length; k++) {
				compString += (" " + params[j][k]);
			}
		}
		logger.info(compString + "\n");

		// linearise the params array
		double[] allParams = new double[totalNumParams];

		// start the counter after the number of comps and the vol fracs
		int paramIndex = numComps + 1;

		// first entry is number of compartments
		allParams[0] = 1;// not the number of compartments but the unweighted
							// signal
		for (int j = 0; j < params.length; j++) {
			// volume fraction goes at the start of the linearised array
			allParams[j + 1] = params[j][0];
			for (int k = 1; k < params[j].length; k++) {
				// copy the params sequentially into the space after the vol
				// fracs
				allParams[paramIndex] = params[j][k];
				paramIndex++;
			}
		}

		String linCompString = new String("linearised params array is:");
		for (int j = 0; j < allParams.length; j++) {
			linCompString += (allParams[j] + " ");
		}
		logger.info(linCompString);

		compartmentNames = comps;
		compParams = allParams;

		// return new CompartmentModel(comps, allParams);

	}



    public static void initTestFunction() {

        // Set up the test function if specified.
        if (testFunction >= 0) {
            
            // Set up the specified options.
            if (rotationIndex > 0) {
                StandardTestFunctions.setTransformation(Rotations.randomRotMat(new Random(rotationIndex)));
            }
            if (scale > 0) {
                StandardTestFunctions.setScale(scale);
            }
            if (dt2rotangle != 0.0) {
                StandardTestFunctions.setDT2RotationAngle(dt2rotangle);
            }
            if (dt2mix > 0) {
                StandardTestFunctions.setMix3(1.0 - dt2mix);
            }
            if (lambda1 > 0) {
                StandardTestFunctions.setLambda1(lambda1);
            }

            p = StandardTestFunctions.getFunction(testFunction);
        }
        
        // A non-standard test function may be specified on the
        // command line.
        if (nonStandardTestFunction) {
            
            // The test function is assumed to be a mixture of Gaussian
            // densities.

            // The first element of the array specified the number of
            // components in the model.
            int components = (int) testFunctionData[0];

            // Construct the arrays of diffusion tensors and mixing
            // parameters
            DT[] dts = new DT[components];
            double[] mixpars = new double[components];
            for (int i = 0; i < components; i++) {
                dts[i] = new DT(testFunctionData[i * 7 + 1],
                                testFunctionData[i * 7 + 2], testFunctionData[i * 7 + 3],
                                testFunctionData[i * 7 + 4], testFunctionData[i * 7 + 5],
                                testFunctionData[i * 7 + 6]);
                mixpars[i] = testFunctionData[i * 7 + 7];

                if (rotationIndex > 0) {
                    dts[i] = dts[i].transform(Rotations.randomRotMat(new Random(rotationIndex)));
                }

            }

            p = new GaussianMixture(dts, mixpars);

            
        }


    }


    /**
     * Creates the source of data from options specified on the command line.
     *
     */
    public static void initDataSynthesizerFromTestFunction() {

        MTRandom ran = new MTRandom(seed);
        
        if (bootstrap < 1) {

            data = new DataSynthesizer(p, imPars, SNR, numVoxels, ran);
        } 
        else if (bootstrap == 1) {
            
            DataSource input = new DataSynthesizer(p, imPars, SNR, 1, ran);
                   
            if (wildBootstrapModel.equals("dt")) {
                data = new WildBS_DT_DataSynth(input, imPars, numVoxels, seed);
            } else {
                throw new LoggedException("Wild bootstrap is not supported for model: " + wildBootstrapModel);
            }
        }
        else {
            data = new SyntheticBootstrapper(p, imPars, SNR, bootstrap, numVoxels, seed);
        }

    }


    public static void initDataSynthesizerFromInputModelData() {

        // if bootstrap then make a bootstrapper using whatever was
        // specified as input
        if (bootstrap > 0) {
            
            DataSynthFromInput input = null;
            
            // If repetition bootstrapping, we set snr = -1 in the data
            // source.
            // The bootstrap object will duplicate each voxel's noise-free
            // data
            // and add noise itself.
            double inputSNR = -1.0;
            
            // If wild bootstrapping, we need one voxel of noisy data from
            // the input source
            if (wildBootstrapModel != null) {
                inputSNR = SNR;
            }
            
            if (inputModel.equals("dt")) {
                input = new DataSynthFromDT_Input(inputFile, inputDataType, imPars, inputSNR, seed);
            } else if (inputModel.equals("twotensor")) {
                input = new DataSynthFromTwoTensorInput(inputFile, inputDataType, imPars, inputSNR, seed);
            } else if (inputModel.equals("threetensor")) {
                input = new DataSynthFromThreeTensorInput(inputFile, inputDataType, imPars, inputSNR, seed);
            } else if (inputModel.equals("multitensor")) {
                input = new DataSynthFromMultiTensorInput(inputFile, inputDataType, imPars, maxTensorComponents, inputSNR, seed);
            } else if (inputModel.equals("ballstick")) {
                input = new DataSynthFromBallStickInput(inputFile, inputDataType, imPars, inputSNR, seed);
            } else {
                throw new LoggedException("Unrecognized input model: " + inputModel);
            }
            
            if (wildBootstrapModel != null) {
                if (wildBootstrapModel.equals("dt")) {
                    data = new WildBS_DT_DataSynth(input, imPars, numVoxels, seed);
                }
                else {
                    throw new LoggedException("Wild bootstrap is not supported for model: " + wildBootstrapModel);
                }
            }
            else {
                data = new BootstrapDataSynthFromInput(input, bootstrap,
                                                       imPars, numVoxels, SNR, seed);
            }
            
        } 
        else {
            // The data is the parameters of models that must be used
            // to generate synthetic data.
            if (inputModel.equals("dt")) {
                data = new DataSynthFromDT_Input(inputFile, inputDataType, imPars, SNR, seed);
            } else if (inputModel.equals("twotensor")) {
                data = new DataSynthFromTwoTensorInput(inputFile, inputDataType, imPars, SNR, seed);
            } else if (inputModel.equals("threetensor")) {
                data = new DataSynthFromThreeTensorInput(inputFile, inputDataType, imPars, SNR, seed);
            } else if (inputModel.equals("multitensor")) {
                data = new DataSynthFromMultiTensorInput(inputFile, inputDataType, imPars, maxTensorComponents, SNR, seed);
            } else if (inputModel.equals("ballstick")) {
                data = new DataSynthFromBallStickInput(inputFile, inputDataType, imPars, SNR, seed);
            } else if (inputModel.equals("sh") || inputModel.equals("rbf")) {
                // data = new DataSynthFromESHS_Input(inputFile,
                // inputDataType,
                // imPars, SNR, seed);
            } else {
                throw new LoggedException("Unrecognized input model: "
                                          + inputModel);
            }
            
        }

    }


    public static void initDataSynthesizerFromSimulation() {

        // if we've set up one or more diffusion simulation
        // parameters, then initialise the data source as a
        // difusion simulation
        double[] stepParams = null;

        SimulableScheme simScheme = null;

        
        if (SimulationParams.trajectories) {
            simParams = new SimulationParams(
                                             SimulationParams.sim_N_walkers,
                                             SimulationParams.sim_tmax, SimulationParams.sim_p,
                                             SimulationParams.sim_initial,
                                             SimulationParams.sim_geomType,
                                             SimulationParams.sim_stepType,
                                             SimulationParams.sim_voxelSize,
                                             SimulationParams.duration);
            
        } 
        else {
            
            simParams = new SimulationParams(
                                             SimulationParams.sim_N_walkers,
                                             SimulationParams.sim_tmax, SimulationParams.sim_p,
                                             SimulationParams.sim_initial,
                                             SimulationParams.sim_geomType,
                                             SimulationParams.sim_stepType,
                                             SimulationParams.sim_voxelSize,
                                             (SimulableScheme) imPars);

            
            try {
                simScheme = (SimulableScheme) imPars;
            } 
            catch (ClassCastException cce) {
            	
            	logger.severe("The scheme file provided doesn't contain enough information to generate measurements from simulation");
            	logger.severe("Simulation-based data synthesis requires gradent strength, pulse durations and timings to be specified");
            	logger.severe("Compatible scheme file formats include STEJSKALTANNER, TRSE, and genwave");
            	
              	throw new LoggedException("Scheme file not compatible with simulation");
            }
        

        }
        
        stepParams = StepGeneratorFactory.getStepParamsArray(SimulationParams.sim_stepType, simParams);
        
        simParams.setStepParams(stepParams);
        
        if (bootstrap > 0) {
            
            int tmpVoxels = numVoxels;
            
            
            // Simulation automatically reads CL_Initializer.numVoxels, so we have to change it here
            //
            // We want bootstrap voxels of simulated data, that we will resample via the bootstrap method
            numVoxels = bootstrap;
            
            DiffusionSimulation sim = new DiffusionSimulation(simParams, simScheme);
            
            // restore numVoxels to the number we actually want
            numVoxels = tmpVoxels;
            
            if (wildBootstrapModel != null) {
                if (wildBootstrapModel.equals("dt")) {
                    data = new WildBS_DT_DataSynth(sim, imPars, numVoxels, seed);
                } 
                else {
                    throw new LoggedException("Wild bootstrap is not supported for model: " + wildBootstrapModel);
                }
            }
            else { 
                data = new BootstrapDataSynthFromSimulationInput(sim, bootstrap, imPars.numMeasurements(), numVoxels, seed);
            }
        } 
        else {
            if (SimulationParams.trajectories) {
                data = new DiffusionSimulation(simParams);
            } 
            else {
                
                // if we've got this far, the scheme object is simulable
                data = new DiffusionSimulation(simParams, simScheme);
            }
        }
        
        
    }


    /**
     * Initializes the synthetic data source. Call initImagingScheme before this
     * method.
     */
    public static void initDataSynthesizer() {

        initTestFunction();
            
        // Create the data source.
        if (p != null) {
            initDataSynthesizerFromTestFunction();
        } 
        else if (inputModel != null) {
            
            // no p, but have an input model
            initDataSynthesizerFromInputModelData();

        } 
        else if (compartmentModel && !compartmentModelForFitting) {
            // compartment model has been specified for data synthesis
            data = new CompartmentModel(compartmentNames, compParams, imPars);
        }
        else if (brownianSimulation) {
            initDataSynthesizerFromSimulation();
        }
        else if (bootstrap > 0) {
                // bootstrap with no input model and no p
                // this implies bootstrapping of real data
                if (SNR > -1) {
                    logger
                        .warning("Cannot set SNR when bootstrapping raw data."
                                 + "Program will continue but will not add noise");
                }

                if (wildBootstrapModel != null) {

                    DataSource dwiData = ExternalDataSource.getDataSource(inputFile, imPars.numMeasurements(), inputDataType);

                    if (wildBootstrapModel.equals("dt")) {
                        data = new WildBS_DT_DataSynth(dwiData, imPars, numVoxels, seed);
                    } else {
                        throw new LoggedException("Wild bootstrap is not supported for model: " + wildBootstrapModel);
                    }
                } else {
                    data = new BootstrapDataSynthesizer(bsDataFiles, imPars.numMeasurements(), numVoxels, inputDataType, seed);
                }
                
        }
        else {
            data = ExternalDataSource.getDataSource(inputFile, imPars.numMeasurements(), inputDataType); 
            initMaskSource();
            initGradAdjSource();
        }
        
        
        
    }
  



	/**
	 * Mark a single argument <code>i</code> as parsed.
	 */
	public static final void markAsParsed(int i) {
		markAsParsed(i, 1);
	}

	/**
	 * Mark sequential arguments as parsed, from <code>i</code> to
	 * <code>i + numArgs - 1</code>.
	 * 
	 * @param i
	 *            the index of the first argument to mark as parsed.
	 * @param numArgs
	 *            the number of arguments to mark as parsed.
	 */
	public static final void markAsParsed(int i, int numArgs) {
		for (int a = 0; a < numArgs; a++) {
			parsed[i + a] = true;
		}
	}

	/**
	 * routine that throws an exeption, to test exception handling.
	 * 
	 * @throws Exception
	 */
	public static void exceptionTest() throws Exception {

		throw new Exception("test exception thrown");
	}

	public static void testCompartmetnParsing() {

		// line of args
		System.err.println("well-formatted, 3-comp model:");
		String argline = new String(
				"-model compartment 3 Ball 0.01 1.0 STICK 0.02 2.0 0.2 0.2 cYlindergpd 0.3 1.3 0.3 0.03");

		// convert to a list of strings
		StringTokenizer tok = new StringTokenizer(argline);

		int numArgs = tok.countTokens();

		String[] args = new String[numArgs];

		for (int count = 0; tok.hasMoreTokens(); count++) {
			args[count] = tok.nextToken();
		}

		parsed = new boolean[args.length];
		for (int i = 0; i < parsed.length; i++) {
			parsed[i] = false;
		}

		initCompartmentModel(args, 1);

		CompartmentModel compModel = new CompartmentModel(compartmentNames,
				compParams);

		checkParsing(args);

		// line of args
		System.err.println("\n\n3-comps, negative vol frac:");
		argline = new String(
				"-model compartment 3 Ball 0.09 1.0 tensor -0.02 1.0 0.0 0.0 1.0 0.0 1.0 cYlindergpd 0.1 1.1 0.3 0.4");

		// convert to a list of strings
		tok = new StringTokenizer(argline);

		numArgs = tok.countTokens();

		args = new String[numArgs];

		for (int count = 0; tok.hasMoreTokens(); count++) {
			args[count] = tok.nextToken();
		}

		parsed = new boolean[args.length];
		for (int i = 0; i < parsed.length; i++) {
			parsed[i] = false;
		}

		initCompartmentModel(args, 1);

		try {
			compModel = new CompartmentModel(compartmentNames, compParams);
		} catch (LoggedException le) {
			System.err.println("exception caught:");
			le.printStackTrace();
		}

		checkParsing(args);

		System.err.println("\n\n3-comps, vol fracs > 1:");
		argline = new String(
				"-model compartment 3 Ball 0.09 1.0 tensor 0.99 1.0 0.0 0.0 1.0 0.0 1.0 cYlindergpd 0.1 1.1 0.3 0.4");

		// convert to a list of strings
		tok = new StringTokenizer(argline);

		numArgs = tok.countTokens();

		args = new String[numArgs];

		for (int count = 0; tok.hasMoreTokens(); count++) {
			args[count] = tok.nextToken();
		}

		parsed = new boolean[args.length];
		for (int i = 0; i < parsed.length; i++) {
			parsed[i] = false;
		}

		initCompartmentModel(args, 1);

		try {
			compModel = new CompartmentModel(compartmentNames, compParams);
		} catch (LoggedException le) {
			System.err.println("exception caught:");
			le.printStackTrace();
		}

		checkParsing(args);

		// line of args
		System.err.println("\n\n3-comps, too few parameters");
		argline = new String(
				"-model compartment 3 Ball 0.09 1.0 tensor 0.02 1.0 0.0 0.0 cYlindergpd 0.1 1.1 0.4");

		// convert to a list of strings
		tok = new StringTokenizer(argline);

		numArgs = tok.countTokens();

		args = new String[numArgs];

		for (int count = 0; tok.hasMoreTokens(); count++) {
			args[count] = tok.nextToken();
		}

		parsed = new boolean[args.length];
		for (int i = 0; i < parsed.length; i++) {
			parsed[i] = false;
		}

		initCompartmentModel(args, 1);

		try {
			compModel = new CompartmentModel(compartmentNames, compParams);
		} catch (NumberFormatException nfe) {
			System.err.println("exception caught:");
			nfe.printStackTrace();
		}

		checkParsing(args);

	}


    


	public static void main(String[] args) {

		testCompartmetnParsing();
	}


    
    /**
     * Sets up various options relating to the definition of physical space of input data, and the header 
     * for output data.
     *
     */
    public static void initInputSpaceAndHeaderOptions() {
        
        ImageHeader inputHeader;

        // Use input file as header template if it is an image and the user has not specified -header
        if (headerTemplateFile == null && inputFile != null && ImageHeader.imageExists(inputFile)) {
            headerTemplateFile = inputFile;
        }

        
        if (headerTemplateFile != null) {
            try {
                inputHeader = ImageHeader.readHeader(headerTemplateFile);
            }
            catch (Exception e) {
                throw new LoggedException(e + "\n\n    Can't use " + headerTemplateFile + 
                                          " to define input space.");
            }
        
            if (dataDims[0] > 0) {
                
                // check for consistency
                int[] hdrDataDims = inputHeader.getDataDims();
                
                for (int i = 0; i < 3; i++) {
                    if (hdrDataDims[i] != dataDims[i]) {
                        throw new LoggedException("Command line data dimensions inconsistent with input header");
                    }
                }

                double[] hdrVoxelDims = inputHeader.getVoxelDims();
                
                for (int i = 0; i < 3; i++) {
                    if ( Math.abs(hdrVoxelDims[i] -  dataDims[i]) > 1E-6 ) {
                        throw new LoggedException("Command line voxel dimensions inconsistent with input header");
                    }
                }
                
                
            }
            else {
                
                dataDims = inputHeader.getDataDims();
                voxelDims = inputHeader.getVoxelDims();

                for (int i = 0; i < 3; i++) {
                    if (voxelDims[i] < 0.0) {
                        logger.warning("Negative voxel dimensions detected, using absolute values, physical space may be misaligned");
                        voxelDims[i] = Math.abs(voxelDims[i]);
                    }
                }
            
            }
            
            headerTemplate = inputHeader;
            
            headerTemplate.setDataType(OutputManager.outputDataType);
            
            // if file is .nii.gz, gzip nii output
            if (headerTemplateFile.endsWith(".gz")) {
                headerTemplate.setGzip(true);
            }
            else {
                headerTemplate.setGzip(OutputManager.gzipOut);
            }
            
        }
        else {

            logger.info("No definition of physical space provided, will attempt to proceed with specified data " + 
                        "and voxel dimensions");

            // no header file - default nii, identity transform
            Nifti1Dataset nds = new Nifti1Dataset();

            if (dataDims[0] == 0) {
                throw new LoggedException("Data dimensions are zero, cannot initialize image space");
            }

            if (voxelDims[0] == 0.0) {
              logger.warning("No header or voxel dimensions specified, assuming identity transform to physical space");
              voxelDims = new double[] {1.0, 1.0, 1.0};
            }

            nds.setDims(4, dataDims[0], dataDims[1], dataDims[2], 1, 0, 0, 0);
            nds.setPixDims(1.0f, (float)voxelDims[0], (float)voxelDims[1], (float)voxelDims[2], 0.0f, 1.0f, 1.0f, 1.0f);

            nds.setQuaternion((short)1, (short)1, new float[] {0.0f, 0.0f, 0.0f}, new float[] {0.0f, 0.0f, 0.0f});

            nds.setDataType(OutputManager.outputDataType);
            
            nds.setFilename("caminoHeaderTemplate", true, OutputManager.gzipOut);
            
            headerTemplate = nds;
        }

    }
    


}

package apps;

import java.io.*;
import java.util.logging.*;

import optimizers.MinimizerException;
import optimizers.MultiRunMinimizer;

import tools.*;
import imaging.*;
import inverters.*;
import data.*;

import fitters.*;

import misc.LoggedException;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Diffusion data inversion application using model fitting.
 * 
 * <dt>Description:
 * 
 * <dd>The program fits simple models, specified by the user, to input data and
 * outputs the model parameters in each voxel.
 * 
 * </dl>
 * 
 * @author Danny Alexander
 * @version $Id$
 * 
 */
public class ModelFit extends Executable {

	/**
	 * Logging object
	 */
	protected static Logger logger = Logger.getLogger("camino.apps.ModelFit");

	// The fitter
	private Fitter fitter;
	// Or the inverter (only one is used; inverters will become
	// deprecated, but code remains for backward compatibility).
	protected DiffusionInversion inv;

	public ModelFit(String[] args) {
		super(args);
	}

	public void initDefaultVals() {
	}

	public void initOptions(String[] args) {

		// Parse the command line arguments
		CL_Initializer.CL_init(args);
		CL_Initializer.checkParsing(args);
		CL_Initializer.initImagingScheme();
		CL_Initializer.initDataSynthesizer();
	}

	public void initVariables() {

		if (CL_Initializer.compartmentModel) {

			FitModel fm = FitModel.getFitModel(CL_Initializer.fitModel);
			NoiseModel nm = NoiseModel.getNoiseModel(CL_Initializer.noiseModel);
			FitAlgorithm fa = FitAlgorithm
					.getFitAlgorithm(CL_Initializer.fitAlgorithm);

			fitter = getFitter(fm, nm, fa);

		} else {

			// Choose the inversion to run.
			inv = getIndexedInversion(CL_Initializer.inversionIndices,
					CL_Initializer.imPars);
		}
	}

	public void execute(OutputManager om) {

            // Loop over the data
            int voxelNumber = 0;
            while (CL_Initializer.data.more())
                try {
                    
                    double[] nextVoxel = CL_Initializer.data.nextVoxel();
                    double[] nextGradAdj = {};
                    if (CL_Initializer.gradAdj!= null) {
                        nextGradAdj = CL_Initializer.gradAdj.nextVoxel();
                    }

                    // Fit or output background default.
                    double backgroundB0 = CL_Initializer.imPars
                        .geoMeanZeroMeas(nextVoxel);
                    boolean bg = isBG(backgroundB0);
                    if (bg) {
                        int ipv = 0;
                        int vps = 0;
                        if (CL_Initializer.compartmentModel) {
                            ipv = fitter.getNumValuesPerRun();
                            vps = fitter.getNumValuesPerSolution();
                        }
                        else {
                            ipv = inv.itemsPerVoxel();
                            vps = ipv;
                        }

                        double[] voxOut = new double[ipv];
                        
                        for(int i=0; i<ipv; i=i+vps) {

                        	// Set the exitcode to -1 to indicate background.
                        	voxOut[i] = -1;
                        
                        	if (backgroundB0 > 0.0) {
                        		voxOut[i+1] = Math.log(backgroundB0);
                        	} else {
                        		voxOut[i+1] = 0.0;
                        	}
                        }
	
                        om.output(voxOut);
                        
                        // The inverter may have other operations to perform
                        // in background voxels.
                        if (!CL_Initializer.compartmentModel)
                            inv.background();
                        
                    } else try {
                        // Fit the model and output the result.
                        if (CL_Initializer.compartmentModel) {
                            double[][] fit;
                            if (CL_Initializer.gradAdj!=null) {
                                fit = fitter.fit(nextVoxel, nextGradAdj);
                            }
                            else {
                                fit = fitter.fit(nextVoxel);
                            }
                            for (int i = 0; i < fit.length; i++) {
                                om.output(fit[i]);
                            }
                        }
                        
                        else {

                            double[] fittedData;
                            if (CL_Initializer.gradAdj!=null) {
                                fittedData = inv.invert(nextVoxel, nextGradAdj);
                            }
                            else {
                                fittedData = inv.invert(nextVoxel);
                            }

                            om.output(fittedData);

                        }

                    } catch (MinimizerException e) {

                        throw new LoggedException(e);

                    }

                    logger.fine("Completed voxel: " + voxelNumber + " " + bg);

                    voxelNumber++;
                    
                } catch (DataSourceException e) {
                    throw new LoggedException("The data file does not contain a whole number of voxels. Check the scheme file. Got Exception " + e);
                }
            
            // Output statistics compiled by the inverter.
            if (!CL_Initializer.compartmentModel)
                inv.close();
            
            // Flush output
            om.close();
            
        }



	/**
	 * assembles the combination of tissue model, noise model and fitting
	 * algorithm into a fitter object that fits the specified model to the data
	 * using the specified algorithm and noise.
	 * 
	 * 
	 * @param fm
	 *            model to fit
	 * @param nm
	 *            noise model to adopt
	 * @param fa
	 *            fitting algorithm to use
	 * 
	 * @return a fitter, if the combination is implemented.
	 */
	public static final Fitter getFitter(FitModel fm, NoiseModel nm,
			FitAlgorithm fa) {

            //if (nm == NoiseModel.GAUSSIAN) {
			if (fm == FitModel.BALLSTICK) {
				if (fa == FitAlgorithm.LM) {
					return new BallStickLM_Fitter(CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallStickMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				} else if (fa == FitAlgorithm.MCMC) {
					return new BallStickMCMC_Fitter(
							CL_Initializer.imPars);
				}
			} else if (fm == FitModel.BIZEPPELIN) {
				if (fa == FitAlgorithm.LM) {
					return new BiZeppelinLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BiZeppelinMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}  else if (fm == FitModel.ZEPPELINGDRCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinGDRCylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinGDRCylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}else if (fm == FitModel.ZEPPELINGDRCYLINDERSDOT) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinGDRCylindersDotLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinGDRCylindersDotMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.ZEPPELINGDRCYLINDERSSPHERE) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinGDRCylindersSphereLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinGDRCylindersSphereMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.ZEPPELINGDRCYLINDERSASTROSTICKS) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinGDRCylindersAstrosticksLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinGDRCylindersAstrosticksMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.ZEPPELINGDRCYLINDERSASTROCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinGDRCylindersAstrocylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinGDRCylindersAstrocylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			
			else if (fm == FitModel.BALLCYLINDER) {
				if (fa == FitAlgorithm.LM) {
					return new BallCylinderLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallCylinderMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINSTICK) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinStickLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinStickMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINSTICKDIRECT) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinStickLM_DirectFitter(CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
				    return new ZeppelinStickMultiRunLM_DirectFitter(CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINSTICKTORT) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinStickTortLM_Fitter(CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
				    return new ZeppelinStickTortMultiRunLM_Fitter(CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINCYLINDER) {
				if (fa == FitAlgorithm.LM) {
				    return new ZeppelinCylinderLM_Fitter(CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinCylinderMultiRunLM_Fitter(CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINCYLINDERDIRECT) {
				if (fa == FitAlgorithm.LM) {
				    return new ZeppelinCylinderLM_DirectFitter(CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinCylinderMultiRunLM_DirectFitter(CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINCYLINDERTORT) {
				if (fa == FitAlgorithm.LM) {
				    return new ZeppelinCylinderTortLM_Fitter(CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinCylinderTortMultiRunLM_Fitter(CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.TENSORSTICK) {
				if (fa == FitAlgorithm.LM) {
					return new TensorStickLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorStickMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.TENSORCYLINDER) {
				if (fa == FitAlgorithm.LM) {
					return new TensorCylinderLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorCylinderMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}

			else if (fm == FitModel.BALLSTICKDOT) {
				if (fa == FitAlgorithm.LM) {
					return new BallStickDotLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallStickDotMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.BALLSTICKASTROSTICKS) {
				if (fa == FitAlgorithm.LM) {
					return new BallStickAstrosticksLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallStickAstrosticksMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.BALLSTICKASTROCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new BallStickAstrocylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallStickAstrocylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}

			else if (fm == FitModel.ZEPPELINSTICKDOT) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinStickDotLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinStickDotMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINSTICKSPHERE) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinStickSphereLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinStickSphereMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.ZEPPELINSTICKASTROSTICKS) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinStickAstrosticksLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinStickAstrosticksMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINSTICKASTROCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinStickAstrocylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinStickAstrocylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.TENSORSTICKDOT) {
				if (fa == FitAlgorithm.LM) {
					return new TensorStickDotLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorStickDotMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} 
			else if (fm == FitModel.TENSORSTICKSPHERE) {
				if (fa == FitAlgorithm.LM) {
					return new TensorStickSphereLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorStickSphereMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			
			else if (fm == FitModel.TENSORSTICKASTROSTICKS) {
				if (fa == FitAlgorithm.LM) {
					return new TensorStickAstrosticksLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorStickAstrosticksMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.TENSORSTICKASTROCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new TensorStickAstrocylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorStickAstrocylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.BALLCYLINDERDOT) {
				if (fa == FitAlgorithm.LM) {
					return new BallCylinderDotLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallCylinderDotMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINCYLINDERDOT) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinCylinderDotLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinCylinderDotMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}else if (fa == FitAlgorithm.MCMC) {
					return new ZeppelinCylinderDotMCMC_GaussianFitter(CL_Initializer.imPars);
				}
			} else if (fm == FitModel.ZEPPELINCYLINDERDOTDIRECT) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinCylinderDotLM_DirectFitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
				    return new ZeppelinCylinderDotMultiRunLM_DirectFitter(CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINCYLINDERDOTCSFDIRECT) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinCylinderDotCSF_LM_DirectFitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
				    return new ZeppelinCylinderDotCSF_MultiRunLM_DirectFitter(CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.TENSORCYLINDERDOT) {
				if (fa == FitAlgorithm.LM) {
					return new TensorCylinderDotLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorCylinderDotMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.BALLCYLINDERASTROSTICKS) {
				if (fa == FitAlgorithm.LM) {
					return new BallCylinderAstrosticksLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallCylinderAstrosticksMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.BALLCYLINDERASTROCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new BallCylinderAstrocylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallCylinderAstrocylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} 
			else if (fm == FitModel.BALLCYLINDERSPHERE) {
				if (fa == FitAlgorithm.LM) {
					return new BallCylinderSphereLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallCylinderSphereMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}else if (fm == FitModel.TENSORCYLINDERASTROSTICKS) {
				if (fa == FitAlgorithm.LM) {
					return new TensorCylinderAstrosticksLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorCylinderAstrosticksMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.TENSORCYLINDERASTROCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new TensorCylinderAstrocylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorCylinderAstrocylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}   else if (fm == FitModel.TENSORCYLINDERSPHERE) {
				if (fa == FitAlgorithm.LM) {
					return new TensorCylinderSphereLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorCylinderSphereMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINCYLINDERASTROSTICKS) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinCylinderAstrosticksLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinCylinderAstrosticksMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.ZEPPELINCYLINDERASTROCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new ZeppelinCylinderAstrocylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new ZeppelinCylinderAstrocylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			} else if (fm == FitModel.MMWMDBASIC) {
				if (fa == FitAlgorithm.LM) {
					return new MMWMD_BasicLM_DirectFitter(CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new MMWMD_BasicMultiRunLM_DirectFitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				} else if (fa == FitAlgorithm.MCMC) {
					return new MMWMD_BasicMCMC_Fitter(CL_Initializer.imPars);
				}

			} else if (fm == FitModel.MMWMDINVIVO) {
				if (fa == FitAlgorithm.LM) {
					return new MMWMD_InVivoLM_DirectFitter(CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new MMWMD_InVivoMultiRunLM_DirectFitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				} else if (fa == FitAlgorithm.MCMC) {
					return new MMWMD_InVivoMCMC_Fitter(CL_Initializer.imPars);
				}

			} else if (fm == FitModel.MMWMDFIXED) {
				if (fa == FitAlgorithm.LM) {
                                    return new MMWMD_FixedLM_DirectFitter(CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
                                    return new MMWMD_FixedMultiRunLM_DirectFitter(CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				} else if (fa == FitAlgorithm.MCMC) {
					return new MMWMD_FixedMCMC_Fitter(CL_Initializer.imPars);
				}

			}

			 else if (fm == FitModel.ZEPPELINCYLINDERSPHERE) {
					if (fa == FitAlgorithm.LM) {
						return new ZeppelinCylinderSphereLM_Fitter(
								CL_Initializer.imPars);
					} else if (fa == FitAlgorithm.MULTIRUNLM) {
						return new ZeppelinCylinderSphereMultiRunLM_Fitter(
								CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
					}
				}
			else if (fm == FitModel.BALLGDRCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new BallGDRCylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallGDRCylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.BALLGDRCYLINDERSDOT) {
				if (fa == FitAlgorithm.LM) {
					return new BallGDRCylindersDotLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallGDRCylindersDotMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}

			else if (fm == FitModel.BALLGDRCYLINDERSSPHERE) {
				if (fa == FitAlgorithm.LM) {
					return new BallGDRCylindersSphereLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallGDRCylindersSphereMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.BALLGDRCYLINDERSASTROSTICKS) {
				if (fa == FitAlgorithm.LM) {
					return new BallGDRCylindersAstrosticksLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallGDRCylindersAstrosticksMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.BALLGDRCYLINDERSASTROCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new BallGDRCylindersAstrocylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallGDRCylindersAstrocylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.TENSORGDRCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new TensorGDRCylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorGDRCylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}

			else if (fm == FitModel.TENSORGDRCYLINDERSDOT) {
				if (fa == FitAlgorithm.LM) {
					return new TensorGDRCylindersDotLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorGDRCylindersDotMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			
			else if (fm == FitModel.TENSORGDRCYLINDERSSPHERE) {
				if (fa == FitAlgorithm.LM) {
					return new TensorGDRCylindersSphereLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorGDRCylindersSphereMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.TENSORGDRCYLINDERSASTROSTICKS) {
				if (fa == FitAlgorithm.LM) {
					return new TensorGDRCylindersAstrosticksLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorGDRCylindersAstrosticksMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.TENSORGDRCYLINDERSASTROCYLINDERS) {
				if (fa == FitAlgorithm.LM) {
					return new TensorGDRCylindersAstrocylindersLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new TensorGDRCylindersAstrocylindersMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.BALLSTICKSPHERE) {
				if (fa == FitAlgorithm.LM) {
					return new BallStickSphereLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new BallStickSphereMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
			else if (fm == FitModel.VERDICTCOLORECTAL) {
				if (fa == FitAlgorithm.LM) {
					return new VERDICTcolorectalLM_Fitter(
							CL_Initializer.imPars);
				} else if (fa == FitAlgorithm.MULTIRUNLM) {
					return new VERDICTcolorectalMultiRunLM_Fitter(
							CL_Initializer.imPars, CL_Initializer.samples, CL_Initializer.seed);
				}
			}
              
                        //}

		

		throw new LoggedException("no fitter available for combination of "
				+ fm + ", " + nm + " and " + fa);

	}

	/**
	 * Creates and returns the indexed diffusion inversion.
	 * 
	 * @param indices
	 *            The indices for the inversion(s). Either {index} or {two
	 *            tensor index, DT index}, or {three tensor index, DT index}.
	 * 
	 * @return The inverter.
	 */
	public static DiffusionInversion getIndexedInversion(ModelIndex[] indices,
			DW_Scheme imPars) {

		// test for non-DT inversions first, then go through the 1,2,3 tensor
		// possibilities

		if (indices[0] == ModelIndex.BALL_STICK) {
			return new BallStickInversion(imPars);
		}
		if (indices[0] == ModelIndex.RESTORE) {

			// Use the RESTORE algorithm
			if (CL_Initializer.sigma == -1) {
				logger.severe("Noise level must be specified (-sigma <std>) for RESTORE.");
				System.exit(1);
			}
			return new RestoreDT_Inversion(imPars, CL_Initializer.sigma);
		}
		if (indices[0] == ModelIndex.ADC) {

			// ADC only inversion.
			return new LinearADC_Inversion(imPars);
		}
		if (indices[0].numDTs == 1) {

			// Single DT inversion.
			return DT_Inversion.getIndexedDT_Inversion(indices[0], imPars);
		}
		if (indices[0].numDTs == 2) {

			if (indices.length == 1) {
				return new TwoTensorInversion(imPars, indices[0],
						ModelIndex.LDT);
			} else {
				return new TwoTensorInversion(imPars, indices[0], indices[1]);
			}
		}
		if (indices[0].numDTs == 3) {

			if (indices.length == 1) {
				return new ThreeTensorInversion(imPars, indices[0],
						ModelIndex.LDT);
			} else {
				return new ThreeTensorInversion(imPars, indices[0], indices[1]);
			}
		}

		return null;
	}

	/**
	 * Check whether the next voxel is background
	 * 
	 * @return boolean indicating background status.
	 */
	public static boolean isBG(double backgroundB0) {

		boolean background = false;
		if (CL_Initializer.BACKGROUNDTHRESHOLD > 0.0) {
			background = backgroundB0 < CL_Initializer.BACKGROUNDTHRESHOLD;
		}
		if (CL_Initializer.bgMask != null)
			try {
				background = (CL_Initializer.bgMask.nextVoxel()[0] == 0.0);
			} catch (DataSourceException e) {
				throw new LoggedException("Error reading background mask." + e);
			}

		return background;
	}

}

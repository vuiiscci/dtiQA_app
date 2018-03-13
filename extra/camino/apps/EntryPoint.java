package apps;

import java.util.logging.Logger;

import data.OutputManager;
//import data.*;

import misc.LoggedException;

import apps.AddNoise;

/**
 * Unified entrypoint for all camino commands,
 * incorporating generic main method structure
 * and command checking.
 * 
 * @author matt, shahrum (s.nedjati-gilani@cs.ucl.ac.uk)
 */
public class EntryPoint {

	private static final Logger logger= Logger.getLogger("apps.EntryPoint");

	/**
	 * parses the first agument from the commandline and
	 * identifies which command we want to run.
	 *
	 * @param execname name of executable from commandline
	 *
	 */
	private static Executable getExecutable(String execname, String[] newArgs){
	
		//Executable b; 
		// identify command to run
		if(execname.equalsIgnoreCase("AddNoise")){
			return new AddNoise(newArgs);
		}
		else if(execname.equalsIgnoreCase("AverageDWI")){
			return new AverageDWI(newArgs);
		}
/*		else if(execname.equalsIgnoreCase("BallSticksOrientationViewer")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("TrajToSignals")){
			return new Scan(newArgs);
		}
		else if(execname.equalsIgnoreCase("BinaryToText")){
			return new BinaryToText(newArgs);
		}
		else if(execname.equalsIgnoreCase("ChunkStats")){
			return new ChunkStats(newArgs);
        }
		else if(execname.equalsIgnoreCase("ClassifiedModelFit")){
			return new ClassifiedModelFit(newArgs);
        }
		/*
		else if(execname.equalsIgnoreCase("CombineTwoFibreLUTs")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("ConnectivityMatrix")){
			return new ConnectivityMatrix(newArgs);
        }
		else if(execname.equalsIgnoreCase("ConsistencyFraction")){
			return new ConsistencyFraction(newArgs);
        }
                /*
		else if(execname.equalsIgnoreCase("CountSeeds")){
			return new blah(newArgs);
		}
		else if(execname.equalsIgnoreCase("CP_Stats")){
			return new blah(newArgs);
		}
		else if(execname.equalsIgnoreCase("DataStats")){
			return new blah(newArgs);
		}
		else if(execname.equalsIgnoreCase("DeconvToCamenoFormat")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("DT_EigenSystem")){
			return new DT_EigenSystem(newArgs);
		}
		else if(execname.equalsIgnoreCase("DT_FitMatrix")){
			return new DT_FitMatrix(newArgs);
		}
		else if(execname.equalsIgnoreCase("DT_ShapeStatistics")){
			return new DT_ShapeStatistics(newArgs);
		}
		else if(execname.equalsIgnoreCase("DT_ToCamino")){
			return new DT_ToCamino(newArgs);
		}
		else if(execname.equalsIgnoreCase("DT_ToImage")){
			return new DT_ToImage(newArgs);
		}
		else if(execname.equalsIgnoreCase("ElectrostaticSubsets")){
			return new ElectrostaticSubsets(newArgs);
		}
		else if(execname.equalsIgnoreCase("EstimateSNR")){
			return new EstimateSNR(newArgs);
		}/*
		else if(execname.equalsIgnoreCase("F_TestThresholdSelector")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("FanningGrid")){
			return new FanningGrid(newArgs);
		}
                else if(execname.equalsIgnoreCase("FixNonSPD_Tensors")){
                        return new FixNonSPD_Tensors(newArgs);
                }
		else if(execname.equalsIgnoreCase("FracAnis")){
			return new FracAnis(newArgs);
		}
		else if(execname.equalsIgnoreCase("FSL_ToScheme")){
			return new FSL_ToScheme(newArgs);
		}/*
		else if(execname.equalsIgnoreCase("GatherStats")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("GenerateDTLUT")){
			return new GenerateDTLUT(newArgs);
		}/*
		else if(execname.equalsIgnoreCase("GenerateMFR_LUT")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("ImageMath")){
			return new ImageMath(newArgs);
		}
		else if(execname.equalsIgnoreCase("ImageToVoxel")){
			return new ImageToVoxel(newArgs);
		}/*
		else if(execname.equalsIgnoreCase("ImportSH_Coeffs")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("InversionStats")){
			return new InversionStats(newArgs);
		}
		else if(execname.equalsIgnoreCase("LinearRecon")){
			return new LinearRecon(newArgs);
		}
		else if(execname.equalsIgnoreCase("MeanDiff")){
			return new MeanDiff(newArgs);
		}/*
		else if(execname.equalsIgnoreCase("MESD")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("ModelFit")){
			return new ModelFit(newArgs);
                }
                else if(execname.equalsIgnoreCase("MultiFibreReconStats")){
			return new MultiFibreReconStats(newArgs);
		}/*
		else if(execname.equalsIgnoreCase("NavPanel")){
			return new blah(newArgs);
                        }*/
		else if(execname.equalsIgnoreCase("OrderElectrostaticPoints")){
			return new OrderElectrostaticPoints(newArgs);
		}
		else if(execname.equalsIgnoreCase("OrientationBiasMap")){
			return new OrientationBiasMap(newArgs);
		}
                /*
		else if(execname.equalsIgnoreCase("PD_OrientationViewer")){
			return new blah(newArgs);
		}
		else if(execname.equalsIgnoreCase("PICoPDF_OrientationViewer")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("PICoPDFs")){
			return new PICoPDFs(newArgs);
		}
		else if(execname.equalsIgnoreCase("PointSetToScheme")){
			return new PointSetToScheme(newArgs);
		}/*
		else if(execname.equalsIgnoreCase("PointSetToScheme2")){
			return new blah(newArgs);
		}
*/
		else if(execname.equalsIgnoreCase("ProcessStreamlines")){
                    return new ProcessStreamlines(newArgs);
		}         
        else if(execname.equalsIgnoreCase("QBallMX")){
			return new QBallMX(newArgs);
		}
		/*else if(execname.equalsIgnoreCase("Reorient")){
			return new Reorient(newArgs);
		}*/
/*		else if(execname.equalsIgnoreCase("RGB_ScalarImage")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("ScannerToVoxel")){
			return new ScannerToVoxel(newArgs);
		}/*
		else if(execname.equalsIgnoreCase("SchemePanel")){
			return new blah(newArgs);
		}
		else if(execname.equalsIgnoreCase("SchemeToFSL")){
			return new blah(newArgs);
		}
		else if(execname.equalsIgnoreCase("SequenceStats")){
			return new blah(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("SelectShells")){
			return new SelectShells(newArgs);
		}
		else if(execname.equalsIgnoreCase("Shredder")){
			return new Shredder(newArgs);
		}
		else if(execname.equalsIgnoreCase("SphFuncAnisotropy")){
			return new SphFuncAnisotropy(newArgs);
		}
		else if(execname.equalsIgnoreCase("SphFuncBitMap")){
			return new SphFuncBitMap(newArgs);
		}	
		else if(execname.equalsIgnoreCase("SphHarmFitter")){
			return new SphHarmFitter(newArgs);
		}		
		else if(execname.equalsIgnoreCase("SphFuncSkewness")){
			return new SphFuncSkewness(newArgs);
		}	
		else if(execname.equalsIgnoreCase("SphPDF_Fit")){
			return new SphPDF_Fit(newArgs);			
		}			
		else if(execname.equalsIgnoreCase("SphFuncKurtosis")){
			return new SphFuncKurtosis(newArgs);
		}			
		else if(execname.equalsIgnoreCase("Split4D_NiiImage")){
			return new Split4D_NiiImage(newArgs);
		}
		else if(execname.equalsIgnoreCase("StreamlineTractography")){
			return new StreamlineTractography(newArgs);
		}
		else if(execname.equalsIgnoreCase("SyntheticData")){
			return new SyntheticData(newArgs);
		}				
		else if(execname.equalsIgnoreCase("VoxelClassify")){
			return new VoxelClassify(newArgs);
		}
		else if(execname.equalsIgnoreCase("SubsetScheme")){
			return new SubsetScheme(newArgs);
		}
		else if(execname.equalsIgnoreCase("SphFuncPD_Stats")){
			return new SphFuncPD_Stats(newArgs);
		}
		/*		else if(execname.equalsIgnoreCase("SyntheticData")){
			return new SyntheticData(newArgs);
		}*/
		else if(execname.equalsIgnoreCase("TraceD")){
			return new TraceD(newArgs);
		}			
		else if(execname.equalsIgnoreCase("VoxelToScanner")){
			return new VoxelToScanner(newArgs);
		}
		else if(execname.equalsIgnoreCase("VoxelToImage")){
			return new VoxelToImage(newArgs);
		}
		else if(execname.equalsIgnoreCase("VoxelwiseImageStats")){
			return new VoxelwiseImageStats(newArgs);
		}
		/*else if(execname.equalsIgnoreCase("WriteZeros")){
			return new WriteZeros(newArgs);
		}*/
		else{
			throw new LoggedException("Unrecognised command "+execname);
		}
		//return b;
		
	}


	/**
	 * unified entrypoint for all camino commands from
	 * commandline.
	 *
	 * @param args array of arguments from commandline
	 */
	public static void main(String[] args){
	
		String[] newArgs = null;
		
		// construct arguments array
		//if(args.length>1){
		newArgs= new String[args.length-1];

		if(args.length>1){	
			for(int i=1; i<args.length; i++){
				newArgs[i-1] = args[i];
			}
		}
			
		Executable exec = getExecutable(args[0], newArgs);
		OutputManager om = new OutputManager();
		exec.execute(om);
	}



}

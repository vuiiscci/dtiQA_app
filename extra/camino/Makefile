# Make file for the camino applications.
# 
# $Id: Makefile,v 1.31 2006/08/10 14:33:12 ucacpco Exp $

# Default javac call
# JAVAC := javac -classpath .

# In a cluster environment, it might be necessary to limit memory usage
JAVAC := javac -J-Xmx500M -classpath .


all: 				AddNoise AnalyzeHeader BinaryToText ClassifiedModelFit ConnectivityMatrix ConsistencyFraction CountSeeds DataStats DeconvToCamenoFormat DT_EigenSystem DT_FitMatrix DT_ShapeStatistics DT_ToCamino DT_ToImage ElectrostaticSubsets EntryPoint EstimateSNR FixNonSPD_Tensors FSL_ToScheme F_TestThresholdSelector FanningGrid FracAnis GatherStats MeanDiff MESD ModelFit ImageMath ImageToVoxel InversionStats LinearRecon MultiFibreReconStats OrderElectrostaticPoints OrientationBiasMap PICoApps PointSetToScheme ProcessStreamlines QBallMX Reorient RGB_ScalarImage ScannerToVoxel SchemeToFSL SequenceStats Shredder Simulation Split4D_Nii SphFuncAnisotropy SphFuncSkewness SphFuncKurtosis SphFuncBitMap SphHarmFitter SphFuncPD_Stats SphPDF_Fit StreamlineTractography SubsetScheme SyntheticData TargetProbsToASCII TensorOrientationViewer ThresholdB0 TraceD TractCounter TractShredder TractStatistics VoxelClassify VoxelMean VoxelToImage VoxelToScanner VoxelwiseImageStats VTK_Streamlines WriteZeros MatlabArrayDataSource jar 

# Builds everything not in the test directory
allclasses::			
				$(JAVAC) `find . -name "*.java" | grep -v \./test`			

AddNoise:			apps/AddNoise.class
apps/AddNoise.class:		apps/AddNoise.java
				$(JAVAC) apps/AddNoise.java

AnalyzeHeader:			imaging/AnalyzeHeader.class
imaging/AnalyzeHeader.class:	imaging/AnalyzeHeader.java
				$(JAVAC) imaging/AnalyzeHeader.java

BinaryToText:			apps/BinaryToText.class

apps/BinaryToText.class: 	apps/BinaryToText.java
				$(JAVAC) apps/BinaryToText.java

ClassifiedModelFit:		apps/ClassifiedModelFit.class

apps/ClassifiedModelFit.class: apps/ClassifiedModelFit.java
				$(JAVAC) apps/ClassifiedModelFit.java

ConnectivityMatrix:		apps/ConnectivityMatrix.class

apps/ConnectivityMatrix.class: apps/ConnectivityMatrix.java
				$(JAVAC) apps/ConnectivityMatrix.java

ConsistencyFraction:		apps/ConsistencyFraction.class

apps/ConsistencyFraction.class: apps/ConsistencyFraction.java
				$(JAVAC) apps/ConsistencyFraction.java

CountSeeds:			apps/CountSeeds.class

apps/CountSeeds.class:		apps/CountSeeds.java
				$(JAVAC) apps/CountSeeds.java


DataStats:			apps/DataStats.class

apps/DataStats.class:		apps/DataStats.java
				$(JAVAC) apps/DataStats.java	


DeconvToCamenoFormat:		apps/DeconvToCamenoFormat.class

apps/DeconvToCamenoFormat.class: apps/DeconvToCamenoFormat.java
				$(JAVAC) apps/DeconvToCamenoFormat.java

DT_EigenSystem:			apps/DT_EigenSystem.class

apps/DT_EigenSystem.class: 	apps/DT_EigenSystem.java
				$(JAVAC) apps/DT_EigenSystem.java

DT_FitMatrix:			apps/DT_FitMatrix.class

apps/DT_FitMatrix.class: 	apps/DT_FitMatrix.java
				$(JAVAC) apps/DT_FitMatrix.java

DT_ShapeStatistics:		apps/DT_ShapeStatistics.class

apps/DT_ShapeStatistics.class: 	apps/DT_ShapeStatistics.java
				$(JAVAC) apps/DT_ShapeStatistics.java

DT_ToCamino:			apps/DT_ToCamino.class

apps/DT_ToCamino.class:		apps/DT_ToCamino.java
				$(JAVAC) apps/DT_ToCamino.java

DT_ToImage:			apps/DT_ToImage.class

apps/DT_ToImage.class:		apps/DT_ToImage.java
				$(JAVAC) apps/DT_ToImage.java


ElectrostaticSubsets:		apps/ElectrostaticSubsets.class

apps/ElectrostaticSubsets.class:	apps/ElectrostaticSubsets.java
				$(JAVAC) apps/ElectrostaticSubsets.java

EntryPoint:			apps/EntryPoint.class
apps/EntryPoint.class:		apps/EntryPoint.java
				$(JAVAC) apps/EntryPoint.java

EstimateSNR:			apps/EstimateSNR.class
apps/EstimateSNR.class:		apps/EstimateSNR.java
				$(JAVAC) apps/EstimateSNR.java

F_TestThresholdSelector:	apps/F_TestThresholdSelector.class

apps/F_TestThresholdSelector.class: apps/F_TestThresholdSelector.java
				$(JAVAC) apps/F_TestThresholdSelector.java

FanningGrid:    apps/FanningGrid.class

apps/FanningGrid.class:		apps/FanningGrid.java
				$(JAVAC) apps/FanningGrid.java

FixNonSPD_Tensors:              apps/FixNonSPD_Tensors.java
				$(JAVAC) apps/FixNonSPD_Tensors.java

FracAnis:			apps/FracAnis.class

apps/FracAnis.class: 		apps/FracAnis.java
				$(JAVAC) apps/FracAnis.java

FSL_ToScheme:			apps/FSL_ToScheme.class

apps/FSL_ToScheme.class: 	apps/FSL_ToScheme.java
				$(JAVAC) apps/FSL_ToScheme.java

GatherStats:			apps/GatherStats.class

apps/GatherStats.class:		apps/GatherStats.java
				$(JAVAC) apps/GatherStats.java

ImageMath:			apps/ImageMath.class

apps/ImageMath.class:		apps/ImageMath.java
				$(JAVAC) apps/ImageMath.java

ImageToVoxel:			apps/ImageToVoxel.class

apps/ImageToVoxel.class:	apps/ImageToVoxel.java
				$(JAVAC) apps/ImageToVoxel.java

MeanDiff:			apps/MeanDiff.class

apps/MeanDiff.class:		apps/MeanDiff.java	
				$(JAVAC) apps/MeanDiff.java

MESD:				apps/MESD.class

apps/MESD.class:		apps/MESD.java	
				$(JAVAC) apps/MESD.java

ModelFit:			apps/ModelFit.class

apps/ModelFit.class:		apps/ModelFit.java	
				$(JAVAC) apps/ModelFit.java

InversionStats:			apps/InversionStats.class

apps/InversionStats.class:	apps/InversionStats.java
				$(JAVAC) apps/InversionStats.java

LinearRecon:			apps/LinearRecon.class

apps/LinearRecon.class:		apps/LinearRecon.java
				$(JAVAC) apps/LinearRecon.java

MultiFibreReconStats:		apps/MultiFibreReconStats.class

apps/MultiFibreReconStats.class:	apps/MultiFibreReconStats.java
				$(JAVAC) apps/MultiFibreReconStats.java

OrderElectrostaticPoints:		apps/OrderElectrostaticPoints.class

apps/OrderElectrostaticPoints.class:	apps/OrderElectrostaticPoints.java
				$(JAVAC) apps/OrderElectrostaticPoints.java

OrientationBiasMap:		apps/OrientationBiasMap.class

apps/OrientationBiasMap.class:	apps/OrientationBiasMap.java
				$(JAVAC) apps/OrientationBiasMap.java

PICoApps:			apps/GenerateDTLUT.class apps/CombineTwoFibreLUTs.class apps/PICoPDFs.class apps/GenerateMFR_LUT.class apps/SphFuncPICoCalibrationData.class

apps/GenerateDTLUT.class:	apps/GenerateDTLUT.java
				$(JAVAC) apps/GenerateDTLUT.java

apps/SphFuncPICoCalibrationData.class:	apps/SphFuncPICoCalibrationData.java
				$(JAVAC) apps/SphFuncPICoCalibrationData.java

apps/GenerateMFR_LUT.class:	apps/GenerateMFR_LUT.java
				$(JAVAC) apps/GenerateMFR_LUT.java


apps/CombineTwoFibreLUTs.class:	apps/CombineTwoFibreLUTs.java
				$(JAVAC) apps/CombineTwoFibreLUTs.java

apps/PICoPDFs.class:		apps/PICoPDFs.java
				$(JAVAC) apps/PICoPDFs.java

PointSetToScheme:		apps/PointSetToScheme.class

apps/PointSetToScheme.class:	apps/PointSetToScheme.java
				$(JAVAC) apps/PointSetToScheme.java

ProcessStreamlines:		apps/ProcessStreamlines.class

apps/ProcessStreamlines.class:	apps/ProcessStreamlines.java
				$(JAVAC) apps/ProcessStreamlines.java

QBallMX:			apps/QBallMX.class

apps/QBallMX.class:		apps/QBallMX.java
				$(JAVAC) apps/QBallMX.java

Reorient:			apps/Reorient.class

apps/Reorient.class:		apps/Reorient.java
				$(JAVAC) apps/Reorient.java

RGB_ScalarImage:		apps/RGB_ScalarImage.class

apps/RGB_ScalarImage.class:	apps/RGB_ScalarImage.java
				$(JAVAC) apps/RGB_ScalarImage.java


ScannerToVoxel:			apps/ScannerToVoxel.class

apps/ScannerToVoxel.class:	apps/ScannerToVoxel.java
				$(JAVAC) apps/ScannerToVoxel.java

SchemeToFSL:			apps/SchemeToFSL.class

apps/SchemeToFSL.class:		apps/SchemeToFSL.java
				$(JAVAC) apps/SchemeToFSL.java


SequenceStats:			apps/SequenceStats.class

apps/SequenceStats.class:	apps/SequenceStats.java
				$(JAVAC) apps/SequenceStats.java

Shredder:			apps/Shredder.class
apps/Shredder.class:		apps/Shredder.java
				$(JAVAC) apps/Shredder.java

Simulation:			apps/SyntheticData.class
				$(JAVAC) simulation/*.java

Split4D_Nii:			apps/Split4D_NiiImage.class
apps/Split4D_NiiImage.class:	apps/Split4D_NiiImage.java
				$(JAVAC) apps/Split4D_NiiImage.java


SphFuncAnisotropy:		apps/SphFuncAnisotropy.class

apps/SphFuncAnisotropy.class:	apps/SphFuncAnisotropy.java
				$(JAVAC) apps/SphFuncAnisotropy.java

SphFuncSkewness:		apps/SphFuncSkewness.class

apps/SphFuncSkewness.class:	apps/SphFuncSkewness.java
				$(JAVAC) apps/SphFuncSkewness.java

SphFuncKurtosis:		apps/SphFuncKurtosis.class

apps/SphFuncKurtosis.class:	apps/SphFuncKurtosis.java
				$(JAVAC) apps/SphFuncKurtosis.java

SphFuncBitMap:			apps/SphFuncBitMap.class

apps/SphFuncBitMap.class:	apps/SphFuncBitMap.java
				$(JAVAC) apps/SphFuncBitMap.java

SphHarmFitter:			apps/SphHarmFitter.class

apps/SphHarmFitter.class:	apps/SphHarmFitter.java
				$(JAVAC) apps/SphHarmFitter.java

SphFuncPD_Stats:		apps/SphFuncPD_Stats.class

apps/SphFuncPD_Stats.class:	apps/SphFuncPD_Stats.java
				$(JAVAC) apps/SphFuncPD_Stats.java

SphPDF_Fit:			apps/SphPDF_Fit.class

apps/SphPDF_Fit.class:		apps/SphPDF_Fit.java
				$(JAVAC) apps/SphPDF_Fit.java


StreamlineTractography:		apps/StreamlineTractography.class
apps/StreamlineTractography.class:	apps/StreamlineTractography.java
					$(JAVAC) apps/StreamlineTractography.java


SubsetScheme:		        apps/SubsetScheme.class

apps/SubsetScheme.class:	apps/SubsetScheme.java
				$(JAVAC) apps/SubsetScheme.java


SyntheticData:		        apps/SyntheticData.class

apps/SyntheticData.class:	apps/SyntheticData.java
				$(JAVAC) apps/SyntheticData.java

TargetProbsToASCII:		apps/TargetProbsToASCII.class

apps/TargetProbsToASCII.class: apps/TargetProbsToASCII.java
				$(JAVAC) apps/TargetProbsToASCII.java


TensorOrientationViewer:	apps/PD_OrientationViewer.class

apps/PD_OrientationViewer.class: apps/PD_OrientationViewer.java
				$(JAVAC) apps/PD_OrientationViewer.java

ThresholdB0:			apps/ThresholdB0.class

apps/ThresholdB0.class: 	apps/ThresholdB0.java
				$(JAVAC) apps/ThresholdB0.java

TractCounter:			apps/TractCounter.class

apps/TractCounter.class: 	apps/TractCounter.java
				$(JAVAC) apps/TractCounter.java

TractShredder:			apps/TractShredder.class

apps/TractShredder.class: 	apps/TractShredder.java
				$(JAVAC) apps/TractShredder.java

TractStatistics:		apps/TractStatistics.class

apps/TractStatistics.class: 		apps/TractStatistics.java
				$(JAVAC) apps/TractStatistics.java


TraceD:				apps/TraceD.class

apps/TraceD.class: 		apps/TraceD.java
				$(JAVAC) apps/TraceD.java

VoxelClassify:			apps/VoxelClassify.class

apps/VoxelClassify.class:	apps/VoxelClassify.java
				$(JAVAC) apps/VoxelClassify.java

VoxelMean:			apps/VoxelMean.class

apps/VoxelMean.class:		apps/VoxelMean.java
				$(JAVAC) apps/VoxelMean.java

VoxelToImage:			apps/VoxelToImage.class

apps/VoxelToImage.class:	apps/VoxelToImage.java
				$(JAVAC) apps/VoxelToImage.java

VoxelToScanner:			apps/VoxelToScanner.class

apps/VoxelToScanner.class:	apps/VoxelToScanner.java
				$(JAVAC) apps/VoxelToScanner.java

VoxelwiseImageStats:		apps/VoxelwiseImageStats.class

apps/VoxelwiseImageStats.class:	apps/VoxelwiseImageStats.java
				$(JAVAC) apps/VoxelwiseImageStats.java

VTK_Streamlines:		apps/VTK_Streamlines.class

apps/VTK_Streamlines.class:	apps/VTK_Streamlines.java
				$(JAVAC) apps/VTK_Streamlines.java

WriteZeros:			apps/WriteZeros.class

apps/WriteZeros.class:		apps/WriteZeros.java
				$(JAVAC) apps/WriteZeros.java


MatlabArrayDataSource:		data/MatlabArrayDataSource.class

data/MatlabArrayDataSource.class:		data/MatlabArrayDataSource.java
				$(JAVAC) data/MatlabArrayDataSource.java


docs:
				javadoc -d ./javadoc Jama apps data imaging inverters linearHARDI linearReconstruction mesd misc numerics optimizers simulation sphfunc tools tractography			
manhtml:			
				doc/makeManHTML.pl

jar:
				jar cf camino.jar `find . -name "*\.class"`


cleandocs:
				rm -rf javadoc

clean:
				rm -f `find . -name "*~"` `find . -name "*\.class"`

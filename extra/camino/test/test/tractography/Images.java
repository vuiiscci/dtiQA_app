package tractography;

import junit.framework.*;
import junit.extensions.*;

import data.*;
import misc.DT;
import numerics.*;

import java.io.*;
import java.util.Random;

/**
 * Provides test images of various types and some methods to get and test values obtained from them.
 *
 *
 */
public class Images {

    private static final String crossingOneDT_File = "./test/tractography/crossingTest.oneDT.Bdouble";
    private static final String crossingTwoDT_File = "./test/tractography/crossingTest.twoDT.Bdouble";
  
    private static final String crossingTwoDTPICoWatsonFile = 
	"./test/tractography/crossingTestPICoWatson.twoDT.WatsonPDF.Bdouble";

    private static final String crossingTwoDTPICoBinghamFile = 
	"./test/tractography/crossingTestPICoBingham.twoDT.BinghamPDF.Bdouble";

    private static final String crossingTwoDTPICoACG_File = 
	"./test/tractography/crossingTestPICoACG.twoDT.ACGPDF.Bdouble";


    private static final String crossingOneSF_File = "./test/tractography/crossingTest.onePD.Bdouble";
    private static final String crossingTwoSF_File = "./test/tractography/crossingTest.twoPD.Bdouble";
 

    private static final String linearOneDT_File = "./test/tractography/linearTest.oneDT.Bdouble";

    private static final String snakeOneDT_File = "./test/tractography/snakeTest.oneDT.Bdouble";

    private static final String cubeOneDT_File = "./test/tractography/cubeTest.oneDT.Bdouble";

    private static final String twoTensorCubeTwoDT_File = "./test/tractography/twoTensorCubeTest.twoDT.Bdouble";


    protected static final double xVoxelDim = 2.0;
    protected static final double yVoxelDim = 3.0;
    protected static final double zVoxelDim = 4.0;
    

    private static final long seed = 73820l;


    /**
     * Just a cube. Each tensor has a different orientation, so you can see which PD you are getting.
     */
    private static DT_TractographyImage cube;


    /*
      slice 0:
      
      |3 | 2|
      |-----|
      |0 | 1|
      
      slice 1:
       
      |7 | 6|
      |-----|
      |4 | 5|
      
      voxel numbers defined in counter clockwise order. All tensors prolate with FA == 0.6.

      voxel 0 has e1 == x
      voxel (1-3) is 0 rotated about z by (1-3) * 30 degrees.

      voxel 4 has e1 == xz
      voxel (5-7) is 0 rotated about z by (1-3) * 30 degrees.

      
    */



    /**
     * Linear anisotropic region with isotropic padding
     */
    private static DT_TractographyImage linear;


    /**
     * Orthogonal fibre crossing - tensor version. 
     */
    private static DT_TractographyImage crossing;


    /**
     * Orthogonal fibre crossing - tensor version. Includes Watson PICo parameters. 
     */
    private static PICoTractographyImage crossingWatson;


    /**
     * Orthogonal fibre crossing - tensor version. Includes Bingham PICo parameters. 
     */
    private static PICoTractographyImage crossingBingham;


    /**
     * Orthogonal fibre crossing - tensor version. Includes ACG PICo parameters.
     */
    private static PICoTractographyImage crossingACG;


    /**
     * Orthogonal fibre crossing - vector (sfpeaks output) version. 
     *
     */
    private static SF_TractographyImage crossingSF;


    /**
     * snake
     */
    private static DT_TractographyImage snake;

    //                 _ _ _ _ _  image is used to test that trackers follow curved paths
    //        _ _ _ _ /



    /**
     * Cube with identical two-tensor model in each voxel.
     */
    private static DT_TractographyImage twoTensorCube;


    // initialize images
    static {


	File f = new File(crossingOneDT_File);
	
	if (!f.exists()) {
	    createImageFiles();
	}
	
	crossing = 
	    DT_TractographyImage.getTractographyImage(crossingTwoDT_File,
                                                         "double", 2, 
                                                         null, 0.0,
                                                         new int[] {3, 3, 2},
                                                         new double[] {xVoxelDim, yVoxelDim, zVoxelDim});

	
	crossingSF = 
	    SF_TractographyImage.getTractographyImage(crossingTwoSF_File,
                                                         "double", 2, 
                                                         null, 0.0, 
                                                         new int[] {3, 3, 2},
                                                         new double[] {xVoxelDim, yVoxelDim, zVoxelDim});
	
	
	crossingBingham = 
	    PICoTractographyImage.getTractographyImage(crossingTwoDTPICoBinghamFile,
                                                           "double", 2, PICoPDF.BINGHAM,
                                                           null, 0.0, 
                                                           new int[] {3, 3, 2},
                                                           new double[] {xVoxelDim, yVoxelDim, zVoxelDim}, 
							   new Random(seed));
	
	crossingACG = 
	    PICoTractographyImage.getTractographyImage(crossingTwoDTPICoACG_File,
                                                           "double", 2, PICoPDF.ACG,
                                                           null, 0.0, 
                                                           new int[] {3, 3, 2},
                                                           new double[] {xVoxelDim, yVoxelDim, zVoxelDim}, 
							   new Random(seed));
        
	crossingWatson =
	    PICoTractographyImage.getTractographyImage(crossingTwoDTPICoWatsonFile,
                                                           "double", 2, PICoPDF.WATSON,
                                                           null, 0.0, 
                                                           new int[] {3, 3, 2},
                                                           new double[] {xVoxelDim, yVoxelDim, zVoxelDim}, 
							   new Random(seed)); 

	
	linear = 
	    DT_TractographyImage.getTractographyImage(linearOneDT_File,
                                                         "double", 1,
                                                         null, 0.0,
                                                         new int[] {5, 1, 1},
                                                         new double[] {xVoxelDim, yVoxelDim, zVoxelDim});
	
	cube = 
	    DT_TractographyImage.getTractographyImage(cubeOneDT_File,
                                                         "double", 1,
                                                         null, 0,
                                                         new int[] {2, 2, 2},
                                                         new double[] {xVoxelDim, yVoxelDim, zVoxelDim});

	snake = 
	    DT_TractographyImage.getTractographyImage(snakeOneDT_File,
                                                         "double", 1,
                                                         null, 0,
                                                         new int[] {5, 2, 1},
                                                         new double[] {xVoxelDim, yVoxelDim, zVoxelDim});

	twoTensorCube =
	    DT_TractographyImage.getTractographyImage(twoTensorCubeTwoDT_File,
                                                      "double", 2,
                                                      null, 0,
                                                      new int[] {2, 2, 2},
                                                      new double[] {xVoxelDim, yVoxelDim, zVoxelDim}); 
    }
    
    
    // define images and write to disk. This only needs to be done once 
    private static void createImageFiles() {
	
	DT[][][][] crossingData = new DT[3][3][2][];

	// simple line
	DT[][][][] linearData = new DT[5][1][1][];

	DT[][][][] snakeData = new DT[5][2][1][];

		
	// FA == 0.5
	DT leftRight = new DT(1143E-12, 0.0, 0.0, 479E-12, 0.0, 479E-12);
	DT upDown = new DT(479E-12, 0.0, 0.0, 1143E-12, 0.0, 479E-12);
	
	// actually, FA == 0.1
	DT isotropic = new DT(781E-12, 0.0, 0.0, 660E-12, 0.0, 660E-12);
	
	DT[] hor = new DT[] {leftRight};
	DT[] ver = new DT[] {upDown};
	DT[] cross = new DT[] {leftRight, upDown};

	DT[] background = new DT[] {isotropic};

	linearData[0][0][0] = background;
	linearData[1][0][0] = hor;
	linearData[2][0][0] = hor;
	linearData[3][0][0] = hor;
	linearData[4][0][0] = background;

	int[][][] linearClassification = new int[5][1][1];
	linearClassification[0][0][0] = 0;
	linearClassification[1][0][0] = 2;
	linearClassification[2][0][0] = 2;
	linearClassification[3][0][0] = 2;
	linearClassification[4][0][0] = 0;


	snakeData[0][0][0] = hor;
	snakeData[1][0][0] = hor;
	snakeData[2][0][0] = new DT[] {leftRight.transform
				       (Rotations.getRotMat(new Vector3D(0.0,0.0,1.0), Math.PI * 0.4))};
	snakeData[3][0][0] = background;
	snakeData[4][0][0] = background;

	snakeData[0][1][0] = background;
	snakeData[1][1][0] = background;
	snakeData[2][1][0] = hor;
	snakeData[3][1][0] = hor;
	snakeData[4][1][0] =  new DT[] {leftRight.transform
				       (Rotations.getRotMat(new Vector3D(0.0,1.0,0.0), Math.PI * 0.2))};

	int[][][] snakeClassification = new int[5][2][1];

	snakeClassification[0][0][0] = 2;
	snakeClassification[1][0][0] = 2;
	snakeClassification[2][0][0] = 2;
	snakeClassification[3][0][0] = 0;
	snakeClassification[4][0][0] = 0;
	snakeClassification[0][1][0] = 0;
	snakeClassification[1][1][0] = 0;
	snakeClassification[2][1][0] = 2;
	snakeClassification[3][1][0] = 2;
	snakeClassification[4][1][0] = 2;
	

	crossingData[0][0][0] = background;
	crossingData[1][0][0] = ver;
	crossingData[2][0][0] = background;

	crossingData[0][1][0] = hor;
	crossingData[1][1][0] = cross;
	crossingData[2][1][0] = hor;
	crossingData[0][2][0] = hor; // classified as -1, so should always come out isotropic
	crossingData[1][2][0] = ver;
	crossingData[2][2][0] = background;

	int[][][] crossingClassification = new int[3][3][2];

	crossingClassification[0][0][0] = 0;
	crossingClassification[1][0][0] = 2;
	crossingClassification[2][0][0] = 0;

	crossingClassification[0][1][0] = 2;
	crossingClassification[1][1][0] = 4;
	crossingClassification[2][1][0] = 2;

	crossingClassification[0][2][0] = -1;
	crossingClassification[1][2][0] = 2;
	crossingClassification[2][2][0] = -1;


	for (int i = 0; i < 3; i++) {
	    for (int j = 0; j < 3; j++) {
		crossingData[i][j][1] = crossingData[i][j][0];
		crossingClassification[i][j][1] = crossingClassification[i][j][0];
	    }
	}


	// for cube
	DT x = new DT(1256E-12, 0.0, 0.0, 422E-12, 0.0, 422E-12);
	DT xz = x.transform(Rotations.getRotMat(Rotations.Y_AXIS, Math.PI / 4.0));

	DT[][][][] cubeData = new DT[2][2][2][];
	
	cubeData[0][0][0] = new DT[] {x};
	cubeData[1][0][0] = new DT[] {x.transform
				      (Rotations.getRotMat(Rotations.Z_AXIS, Math.PI / 6.0))};
	cubeData[1][1][0] = new DT[] {x.transform
				      (Rotations.getRotMat(Rotations.Z_AXIS, 2.0 * Math.PI / 6.0))};

	cubeData[0][1][0] = new DT[] {x.transform
				      (Rotations.getRotMat(Rotations.Z_AXIS, 3.0 * Math.PI / 6.0))};


	cubeData[0][0][1] = new DT[] {xz};
	cubeData[1][0][1] =new DT[] {xz.transform
				     (Rotations.getRotMat(Rotations.Z_AXIS, Math.PI / 6.0))};

	cubeData[1][1][1] = new DT[] {xz.transform
				      (Rotations.getRotMat(Rotations.Z_AXIS, 2.0 * Math.PI / 6.0))};

	cubeData[0][1][1] = new DT[] {xz.transform
				      (Rotations.getRotMat(Rotations.Z_AXIS, 3.0 * Math.PI / 6.0))};


	int[][][] cubeClassification = new int[2][2][2];

	cubeClassification[0][0][0] = 2;
	cubeClassification[0][0][1] = 2;
	cubeClassification[0][1][0] = 2;
	cubeClassification[0][1][1] = 2;
	cubeClassification[1][0][0] = 2;
	cubeClassification[1][0][1] = 2;
	cubeClassification[1][1][0] = 2;
	cubeClassification[1][1][1] = 2;


	DT[][][][] twoTensorCubeData = new DT[2][2][2][];

	// simple

	double fa = 0.1;
	
	for (int i = 0; i < 2; i++) {
	    for (int j = 0; j < 2; j++) {
		for (int k = 0; k < 2; k++) {
		    double t = 2100E-12;
		    
		    double l1 = (-3.0 * t + 2 * fa * fa * t - 2.0 * 
				 Math.sqrt(3.0 * fa * fa * t * t - 2 * fa * fa * fa * fa * t * t))
			/(3.0 * (-3.0 + 2.0 * fa * fa));
		    
		    double l2 = (t - l1) / 2.0;
		    double l3 = l2;
		    
		    DT d1 = new DT(l1, 0.0, 0.0, l2, 0.0, l3);
		    DT d2 = new DT(l2, 0.0, 0.0, l1, 0.0, l3);
		    
		    twoTensorCubeData[i][j][k] = new DT[] {d1, d2};
		    
		    fa += 0.1;
		}
	    }
	}

	int[][][] twoTensorCubeClassification = new int[2][2][2];

	twoTensorCubeClassification[0][0][0] = 4;
	twoTensorCubeClassification[0][0][1] = 4;
	twoTensorCubeClassification[0][1][0] = 4;
	twoTensorCubeClassification[0][1][1] = 4;
	twoTensorCubeClassification[1][0][0] = 4;
	twoTensorCubeClassification[1][0][1] = 4;
	twoTensorCubeClassification[1][1][0] = 4;
	twoTensorCubeClassification[1][1][1] = 4;


	try {

	    // write crossing for test
	    FileOutputStream fout = new FileOutputStream(crossingTwoDT_File);

	    DataOutputStream twoDTout = new DataOutputStream(new BufferedOutputStream(fout, 1024*10));
	    

	    fout = new FileOutputStream(crossingTwoDTPICoWatsonFile);
	    
	    DataOutputStream twoDTPICoWatsonout = 
		new DataOutputStream(new BufferedOutputStream(fout, 1024*10));

	    fout = new FileOutputStream(crossingTwoDTPICoBinghamFile);
	    
	    DataOutputStream twoDTPICoBinghamout = 
		new DataOutputStream(new BufferedOutputStream(fout, 1024*10));
	    
	    fout = new FileOutputStream(crossingTwoDTPICoACG_File);
	    
	    DataOutputStream twoDTPICoACGout = 
		new DataOutputStream(new BufferedOutputStream(fout, 1024*10));
	    
	    fout = new FileOutputStream(crossingTwoSF_File);
	    
	    DataOutputStream twoSFout = 
		new DataOutputStream(new BufferedOutputStream(fout, 1024*10));

	    
	    double[] lrComps = leftRight.getComponents();
	    double[] udComps = upDown.getComponents();
	    
	    DT oblate = new DT((lrComps[0] + udComps[0]) / 2.0, 
			       (lrComps[1] + udComps[1]) / 2.0, 
			       (lrComps[2] + udComps[2]) / 2.0, 
			       (lrComps[3] + udComps[3]) / 2.0, 
			       (lrComps[4] + udComps[4]) / 2.0, 
			       (lrComps[5] + udComps[5]) / 2.0);
			       
	    for (int k = 0; k < 2; k++) {
		for (int j = 0; j < 3; j++) {
		    for (int i = 0; i < 3; i++) {

                        // exit
                        if (crossingClassification[i][j][k] == -1) {
                            twoDTout.writeDouble(-1.0);    
                            twoSFout.writeDouble(-1.0);    
                        }
                        else {
                            twoDTout.writeDouble(0.0);    
                            twoSFout.writeDouble(0.0);    
                        }


                        // ln A(0)
			twoDTout.writeDouble(1.0);
                        twoSFout.writeDouble(1.0);    
                        

			if (i == 1 && j == 1) {
                            
                            // numPDs mix
                            twoDTout.writeDouble(2.0);
                            twoDTout.writeDouble(0.5);

                            twoDTPICoWatsonout.writeDouble(2.0);

                            twoDTPICoBinghamout.writeDouble(2.0);

                            twoDTPICoACGout.writeDouble(2.0);

                            // dt file
			    double[] comps = crossingData[i][j][k][0].getComponents();
			    for (int c = 0; c < 6; c++) {
				twoDTout.writeDouble(comps[c]);
			    }

                            twoDTout.writeDouble(0.5);

			    comps = crossingData[i][j][k][1].getComponents();
			    for (int c = 0; c < 6; c++) {
				twoDTout.writeDouble(comps[c]);
			    }

                            // write others
                            twoSFout.writeDouble(2.0);    // number of peaks

                            // flag for  consistency  with  repeated  run
                            // mean(f). 
                            // std(f).  
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    


                            double[][] sEig1 = crossingData[i][j][k][0].sortedEigenSystem();
                            double[][] sEig2 = crossingData[i][j][k][1].sortedEigenSystem();
                            
                            twoSFout.writeDouble(sEig1[1][0]);
                            twoSFout.writeDouble(sEig1[2][0]);
                            twoSFout.writeDouble(sEig1[3][0]);  

                            // f
                            twoSFout.writeDouble(1.0);    
                            
                            // H
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    

                            twoSFout.writeDouble(sEig2[1][0]);
                            twoSFout.writeDouble(sEig2[2][0]);
                            twoSFout.writeDouble(sEig2[3][0]);  
                                
                            // f
                            twoSFout.writeDouble(1.0);    
                            
                            // H
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    

                            Vector3D[][] evecs = new Vector3D[2][3];
                            
                            // first eigenvector is the e1 of each tensor
                            evecs[0][0] = new Vector3D(sEig1[1][0], sEig1[2][0], sEig1[3][0]);
                            evecs[1][0] = new Vector3D(sEig2[1][0], sEig2[2][0], sEig2[3][0]);
                            
                           
                            // third eigenvector is cross product of the two firsts, ie normal to plane of crossing
                            evecs[0][2] = evecs[0][0].cross(evecs[1][0]).normalized();
                            evecs[1][2] = evecs[0][2];
                            
                            // second is the first rotated about the third by 90 degrees
                            evecs[0][1] = Rotations.rotateVector(evecs[0][0], evecs[0][2], Math.PI / 2.0);
                            evecs[1][1] = Rotations.rotateVector(evecs[1][0], evecs[1][2], Math.PI / 2.0);
                            

                            for (int n = 0; n < 2; n++) {

                                twoDTPICoWatsonout.writeDouble(0.5);
                                twoDTPICoBinghamout.writeDouble(0.5);
                                twoDTPICoACGout.writeDouble(0.5);

                                for (int m = 0; m < 3; m++) {
                                    twoDTPICoWatsonout.writeDouble(evecs[n][m].x);
                                    twoDTPICoWatsonout.writeDouble(evecs[n][m].y);
                                    twoDTPICoWatsonout.writeDouble(evecs[n][m].z);
                                    
                                    twoDTPICoBinghamout.writeDouble(evecs[n][m].x);
                                    twoDTPICoBinghamout.writeDouble(evecs[n][m].y);
                                    twoDTPICoBinghamout.writeDouble(evecs[n][m].z);
                                    
                                    twoDTPICoACGout.writeDouble(evecs[n][m].x);
                                    twoDTPICoACGout.writeDouble(evecs[n][m].y);
                                    twoDTPICoACGout.writeDouble(evecs[n][m].z);
                                }
                                
                                twoDTPICoWatsonout.writeDouble(50.0);
                                
                                twoDTPICoBinghamout.writeDouble(-100.0);
                                twoDTPICoBinghamout.writeDouble(-50.0);
                                
                                twoDTPICoACGout.writeDouble(2.942083193342585);
                                twoDTPICoACGout.writeDouble(0.08377478758707287);
                                twoDTPICoACGout.writeDouble(0.008377478758707287);
                                
                            }
                        
                            
                        
                        }
			else {

                            // numPDs mix
                            twoDTout.writeDouble(1.0);
                            twoDTout.writeDouble(1.0);
                            
                             if (crossingClassification[i][j][k] == -1) { 
                                 twoDTPICoWatsonout.writeDouble(0.0);
                                 twoDTPICoBinghamout.writeDouble(0.0);
                                 twoDTPICoACGout.writeDouble(0.0);
                             }
                             else {
                                 
                                 twoDTPICoWatsonout.writeDouble(1.0);
                                 twoDTPICoBinghamout.writeDouble(1.0);
                                 twoDTPICoACGout.writeDouble(1.0);
                             }

                            twoDTPICoWatsonout.writeDouble(1.0);

                            twoDTPICoBinghamout.writeDouble(1.0);

                            twoDTPICoACGout.writeDouble(1.0);

                            // number of peaks
                            twoSFout.writeDouble(1.0);  
                            // flag for  consistency  with  repeated  run
                            // mean(f). 
                            // std(f).  
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    


                            double[] comps = crossingData[i][j][k][0].getComponents();
			    
                            for (int c = 0; c < 6; c++) {
                                twoDTout.writeDouble(comps[c]);
                            }
                            
                            twoDTout.writeDouble(0.0);
                            
                            for (int c = 0; c < 6; c++) {
                                twoDTout.writeDouble(0.0);
                            }
                            

			    double[][] sEig = crossingData[i][j][k][0].sortedEigenSystem();

                            twoSFout.writeDouble(sEig[1][0]);
			    twoSFout.writeDouble(sEig[2][0]);
			    twoSFout.writeDouble(sEig[3][0]);  

			    for (int m = 0; m < 3; m++) {
                                twoDTPICoWatsonout.writeDouble(sEig[1][m]);
                                twoDTPICoWatsonout.writeDouble(sEig[2][m]);
                                twoDTPICoWatsonout.writeDouble(sEig[3][m]);

                                twoDTPICoBinghamout.writeDouble(sEig[1][m]);
                                twoDTPICoBinghamout.writeDouble(sEig[2][m]);
                                twoDTPICoBinghamout.writeDouble(sEig[3][m]);

                                twoDTPICoACGout.writeDouble(sEig[1][m]);
                                twoDTPICoACGout.writeDouble(sEig[2][m]);
                                twoDTPICoACGout.writeDouble(sEig[3][m]);
                            }
                            
                            twoDTPICoWatsonout.writeDouble(100.0);
                            
                            twoDTPICoBinghamout.writeDouble(-200.0);
                            twoDTPICoBinghamout.writeDouble(-100.0);
                            
                            twoDTPICoACGout.writeDouble(2.942083193342585);
                            twoDTPICoACGout.writeDouble(0.08377478758707287);
                            twoDTPICoACGout.writeDouble(0.008377478758707287);
                            
                            // f
                            twoSFout.writeDouble(1.0);    
                            
                            // H
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    
                            twoSFout.writeDouble(1.0);    
                            
                            
			    for (int n = 0; n < 11; n++) {
                                twoDTPICoWatsonout.writeDouble(0.0);
                            }

			    for (int n = 0; n < 12; n++) {
                                twoDTPICoBinghamout.writeDouble(0.0);
                            }

			    for (int n = 0; n < 13; n++) {
                                twoDTPICoACGout.writeDouble(0.0);
                            }
                            
			    for (int n = 0; n < 8; n++) {
                                twoSFout.writeDouble(0.0);
                            }
                            
                            

			}
			
		    }
		}
	    }


	    twoDTout.close();

	    twoSFout.close();

	    twoDTPICoWatsonout.close();
	    twoDTPICoBinghamout.close();
	    twoDTPICoACGout.close();


	    // linear
 
	    fout = new FileOutputStream(linearOneDT_File);
	    DataOutputStream oneDTout = new DataOutputStream(new BufferedOutputStream(fout, 1024*10));

	    for (int k = 0; k < 1; k++) {
		for (int j = 0; j < 1; j++) {
                    for (int i = 0; i < 5; i++) {
			// vcout.writeInt(linearClassification[i][j][k]);

                        if (linearClassification[i][j][k] == -1) {
                            oneDTout.writeDouble(-1.0);    
                        }
                        else {
                            oneDTout.writeDouble(0.0);    
                        }

                        oneDTout.writeDouble(1.0);
			
			double[] comps = linearData[i][j][k][0].getComponents();
			for (int c = 0; c < 6; c++) {
			    oneDTout.writeDouble(comps[c]);
			}
		    }
		}
	    }

	    oneDTout.close();



	    // snake

	    fout = new FileOutputStream(snakeOneDT_File);
	    oneDTout = new DataOutputStream(new BufferedOutputStream(fout, 1024*10));

	    for (int k = 0; k < 1; k++) {
		for (int j = 0; j < 2; j++) {
		    for (int i = 0; i < 5; i++) {

                        if (snakeClassification[i][j][k] == -1) {
                            oneDTout.writeDouble(-1.0);    
                        }
                        else {
                            oneDTout.writeDouble(0.0);    
                        }

			oneDTout.writeDouble(1.0);
			
			double[] comps = snakeData[i][j][k][0].getComponents();
			for (int c = 0; c < 6; c++) {
			    oneDTout.writeDouble(comps[c]);
			}
		    }
		}
	    }

	    oneDTout.close();


	    // cube

	    fout = new FileOutputStream(cubeOneDT_File);
	    oneDTout = new DataOutputStream(new BufferedOutputStream(fout, 1024*10));

	    for (int k = 0; k < 2; k++) {
		for (int j = 0; j < 2; j++) {
		    for (int i = 0; i < 2; i++) {
			
                        if (cubeClassification[i][j][k] == -1) {
                            oneDTout.writeDouble(-1.0);    
                        }
                        else {
                            oneDTout.writeDouble(0.0);    
                        }

			oneDTout.writeDouble(1.0);
			
			double[] comps = cubeData[i][j][k][0].getComponents();
			for (int c = 0; c < 6; c++) {
			    oneDTout.writeDouble(comps[c]);
			}
		    }
		}
	    }

	    oneDTout.close();

	    // twoTensorCube

	    fout = new FileOutputStream(twoTensorCubeTwoDT_File);
	    twoDTout = new DataOutputStream(new BufferedOutputStream(fout, 1024*10));

	    for (int k = 0; k < 2; k++) {
		for (int j = 0; j < 2; j++) {
		    for (int i = 0; i < 2; i++) {
                        
                        if (twoTensorCubeClassification[i][j][k] == -1) {
                            twoDTout.writeDouble(-1.0);    
                        }
                        else {
                            twoDTout.writeDouble(0.0);    
                        }

			twoDTout.writeDouble(1.0);
			twoDTout.writeDouble(2.0);
			
                        twoDTout.writeDouble(0.5);
			    
			double[] comps = twoTensorCubeData[i][j][k][0].getComponents();
			for (int c = 0; c < 6; c++) {
			    twoDTout.writeDouble(comps[c]);
			}
			
                        twoDTout.writeDouble(0.5);
			
			comps = twoTensorCubeData[i][j][k][1].getComponents();
			for (int c = 0; c < 6; c++) {
			    twoDTout.writeDouble(comps[c]);
			}

			
		    }
		}
	    }

	    twoDTout.close();


	    
	    
	}
	catch(IOException e) {
	    throw new RuntimeException("While trying to write images: " + e);
	}


    }





    /**
     *
     * Straight line of anisotropic tensors with isotropic (FA==0.1) voxels at each end.
     *
     */
    public static DT_TractographyImage getLinear() {

	return new DT_TractographyImage(linear);
    }


    /**
     *
     * Orthogonal fibre crossing.
     *
     */
    public static DT_TractographyImage getCrossing() {

	return new DT_TractographyImage(crossing);
    }

    
    /**
     *
     * Orthogonal fibre crossing with PICo parameters.
     *
     */
    public static PICoTractographyImage getCrossingWatson() {

	return new PICoTractographyImage(crossingWatson, new Random(seed));
    }

    /**
     *
     * Orthogonal fibre crossing with PICo parameters.
     *
     */
    public static PICoTractographyImage getCrossingBingham() {

	return new PICoTractographyImage(crossingBingham, new Random(seed));
    }


    /**
     *
     * Orthogonal fibre crossing with PICo parameters.
     *
     */
    public static PICoTractographyImage getCrossingACG() {

	return new PICoTractographyImage(crossingACG, new Random(seed));
    }


    /**
     *
     * Orthogonal fibre crossing, SFs only.
     *
     */
    public static SF_TractographyImage getCrossingSF() {

	return new SF_TractographyImage(crossingSF);
    }


    public static DT_TractographyImage getSnake() {

	return new DT_TractographyImage(snake);
    }


    /**
     *
     * 2x2x2 block, each with unique PD.
     *
     */
    public static DT_TractographyImage getCube() {

	return new DT_TractographyImage(cube);
    }


    /**
     *
     * 2x2x2 block, each with a fibre crossing of a tensor lying along the x axis and another 
     * along the y axis, varying FA.
     *
     */
    public static DT_TractographyImage getTwoTensorCube() {

	return new DT_TractographyImage(twoTensorCube);
    }


    /**
     * Get a point at the centre of a voxel for use in testing
     */
    public static Point3D getPointAtVoxelCentre(int i, int j, int k) {
        
        return new Point3D((i + 0.5) * xVoxelDim, (j + 0.5) * yVoxelDim, (k + 0.5) * zVoxelDim);

    }



    /**
     *
     * Convert a voxel index to mm. Follows the Camino convention where (0,0,0) is the corner of the voxel, 
     * and (0.5, 0.5, 0.5) is the centre.
     *
     * 
     */
    public static Point3D getPointAtVoxel(double x, double y, double z) {
        
        return new Point3D(x * xVoxelDim, y * yVoxelDim, z * zVoxelDim);

    }



    /**
     * @return true if all tensor components are equal within the specified delta
     */
    public static void assertTensorsEqual(TestCase test, DT d1, DT d2, double delta) {

	
	double[] d1c = d1.getComponents();
	double[] d2c = d2.getComponents();

	for (int i = 0; i < 6; i++) {
	    if ( Math.abs(d1c[i] - d2c[i]) > delta ) {
		test.fail("expected\n" + d1 + " but was\n" + d2);
	    }
	    
	}

    }

 
}


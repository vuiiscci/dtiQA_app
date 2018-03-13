package apps;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;

import java.io.*;

import javax.swing.event.*;
import java.text.*;

import java.util.*;
import java.util.logging.Logger;

import data.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;

/**
 *
 * Simple tool to visualize principal directions of tensor data. Reads output of dteig.
 * 
 * @author Philip Cook, Kiran Seunarine
 * @version $Id$
 *  
 */
public abstract class PD_OrientationViewer extends JPanel {

    /**
     * Logging object
     */
    private static Logger logger = Logger.getLogger("camino.apps.PD_OrientationViewer");

    private JFrame frame;

    private JPanel sliceSelectionPanel;

    private NavPanel navPanel;
    private SchemePanel schemePanel;

    private JList<Integer> sliceList;

    // the image that is being displayed in the main window
    private BufferedImage zoomedSlice;
   
    // used to set slice orientation
    public static final int XY = 100;
    public static final int XZ = 200;
    public static final int YZ = 300;
    
    // schemefile editor variables
    public static final int NO_FLIP = 0;
    public static final int FLIP_X = 1;
    public static final int FLIP_Y = 2;
    public static final int FLIP_Z = 3;
    
    private int [] flipDirs = {1, 1, 1};
    private int [] currentXYZ_Order = DW_Scheme.gradXYZ;

    private boolean isFlippedDir = false;
    private boolean isSwappedDirs = false;
    protected boolean showVectors = true;

    protected int sliceOrientation = -1;

    protected int currentSliceIndex = 0;

    protected double[][][] scalarVol = null;

    // between 0.0 and 1.0
    protected double[][][] red = null;
    protected double[][][] green = null;
    protected double[][][] blue = null;

    // reference to image's vectors
    protected Vector3D[][][][] vectors;

    // data coordinates, describes size of array
    protected int xSize;

    protected int ySize;

    protected int zSize;

    // number of slices in current stack - used to save slices
    protected int noSlices;

    
    // default number of voxels in main window is preferredZoomedExtent * preferredZoomedExtent
    // change this to zoom in and out
    protected int preferredZoomedExtent = 30;


    // normalized scalar threshold for display of PDs
    protected double normScalarThresh;


    protected final RGB_ScalarImage image;

    protected final MouseWheelListener wheelListener;


    /**
     * @param im the image to display.
     * @param scalarThresh do not display vectors unless the scalar value is above the threshold.
     */
    protected PD_OrientationViewer(RGB_ScalarImage im, double scalarThresh) {

	image = im;

	vectors = im.vectors;
	scalarVol = im.scalarVol;
	
	if (scalarThresh < im.minScalarValue) {
	    normScalarThresh = -1.0;
	}
	else if (scalarThresh > im.maxScalarValue) {
	    normScalarThresh = 1.0;
	}
	else {
	    normScalarThresh = (scalarThresh - im.minScalarValue) / (im.maxScalarValue - im.minScalarValue);   
	}

	xSize = image.xDataDim;
	ySize = image.yDataDim;
	zSize = image.zDataDim;

        sliceOrientation = XY;

        if (preferredZoomedExtent % 2 > 0) {
            preferredZoomedExtent = preferredZoomedExtent - 1;
        }
	if (preferredZoomedExtent < 2) {
	    preferredZoomedExtent = 2;
	}

	vectors = image.vectors;

	// scroll slices with mouse wheel
	wheelListener = new MouseWheelListener() {
	    
		public void mouseWheelMoved(MouseWheelEvent e) {
		    int count = e.getWheelRotation();
		    
		    int newSliceIndex = sliceList.getSelectedIndex() + count;

		    int max = sliceList.getModel().getSize() - 1;

		    if (newSliceIndex < 0) {
			newSliceIndex = 0;
		    }
		    else if (newSliceIndex > max) {
			newSliceIndex = max;
		    }
		    
                    
		    sliceList.setSelectedIndex(newSliceIndex);
		    sliceList.ensureIndexIsVisible(newSliceIndex); 
		    
		}
	    };
	
	addMouseWheelListener(wheelListener);


	// add top panel after wheel listener because the zoom panel also uses it
        navPanel = new NavPanel(this);

	schemePanel = new SchemePanel(this);

	
    }


    /**
     * Gets e1 and scalar value for a point in the zoomed panel.
     *
     */
    protected String getText(int ex, int ey) {
       
        Dimension panelSize = getSize(); // should be the size of the panel

	int[] bounds = navPanel.getZoomedRegion(); // extent of zoomed region

	// pixels per voxel, assuming isotropic display
        double xScale = (panelSize.getWidth() / (bounds[1] - bounds[0] + 1.0));
        double yScale = (panelSize.getHeight() / (bounds[3] - bounds[2] + 1.0));

	if (xScale < yScale) {
            yScale = xScale;
        }
        else {
            xScale = yScale;
        }
	

	// in data coordinates
        int x = 0;
        
        int y = 0;

	x = (int)(ex / xScale);
	y = (int)((panelSize.getHeight() - ey) / yScale);

	// bounds in data coordinates
	int xMinVoxel = bounds[0];

	int xMaxVoxel = bounds[1];

	int yMinVoxel = bounds[2];

	int yMaxVoxel = bounds[3];

	int dataXIndex = 0;
	int dataYIndex = 0;
	int dataZIndex = 0;
	

        switch (sliceOrientation) {

        case PD_OrientationViewer.XY:
	    dataXIndex = xMinVoxel + x;
	    dataYIndex = yMinVoxel + y;
	    dataZIndex = currentSliceIndex;
	    break;

        case PD_OrientationViewer.XZ:
	    dataXIndex = xMinVoxel + x;
	    dataYIndex = currentSliceIndex;
	    dataZIndex = yMinVoxel + y;
	    break;

        case PD_OrientationViewer.YZ:
	    dataXIndex = currentSliceIndex;
	    dataYIndex = xMinVoxel + x;
	    dataZIndex = yMinVoxel + y;
	    break;
                            
        default:
            throw new java.lang.IllegalStateException("Unrecognised plane index "
                                                      + sliceOrientation);
                            
        }

	DecimalFormat df = new DecimalFormat("0.000000");

	if (dataXIndex < 0 || dataXIndex >= xSize || dataYIndex < 0 || dataYIndex >= ySize ||
	    dataZIndex < 0 || dataZIndex >= zSize) {
	    return "\n";
	}
	if (vectors[dataXIndex][dataYIndex][dataZIndex].length == 1) {
	    return dataXIndex + " " + dataYIndex + " " + dataZIndex + 
		"\n" + df.format(vectors[dataXIndex][dataYIndex][dataZIndex][0].x) + 
		" " + df.format(vectors[dataXIndex][dataYIndex][dataZIndex][0].y) +
		" " + df.format(vectors[dataXIndex][dataYIndex][dataZIndex][0].z) + 
		" (" + df.format(scalarVol[dataXIndex][dataYIndex][dataZIndex]) +
		")\n";
	}
	else if (vectors[dataXIndex][dataYIndex][dataZIndex].length >= 2) {
	    return dataXIndex + " " + dataYIndex + " " + dataZIndex + 
		"\n" + df.format(vectors[dataXIndex][dataYIndex][dataZIndex][0].x) + 
		" " + df.format(vectors[dataXIndex][dataYIndex][dataZIndex][0].y) +
		" " + df.format(vectors[dataXIndex][dataYIndex][dataZIndex][0].z) + 
		" (" + df.format(scalarVol[dataXIndex][dataYIndex][dataZIndex]) +
		")\n" + df.format(vectors[dataXIndex][dataYIndex][dataZIndex][1].x) + 
		" " + df.format(vectors[dataXIndex][dataYIndex][dataZIndex][1].y) +
		" " + df.format(vectors[dataXIndex][dataYIndex][dataZIndex][1].z);
	}
	else {
	    return "\n";
	}
	
    }
    

    public void writeSlices() {
	
	String fileRoot = "";

	JFileChooser chooser = new JFileChooser();
	
	javax.swing.filechooser.FileFilter pngFilter = getFileFilter("PNG (*.png)","png");

	chooser.setAcceptAllFileFilterUsed(false);
	chooser.addChoosableFileFilter(pngFilter);

	chooser.setFileFilter(pngFilter);

	try {
	    File currentDir = new File(new File(".").getCanonicalPath());
	    
	    chooser.setCurrentDirectory(currentDir);
	}
	catch (IOException e) {
	    throw new LoggedException(e);
	}
	    
	String filePath = null;
	
	int r = chooser.showSaveDialog(this);

	if(r == JFileChooser.APPROVE_OPTION) {
	    

	    // Assuming *nix systems will use absolute paths (ie, they will start with /), seems to work
            // Careful about using system file separators, because we might need to use forward slash
            // on Windows (Cygwin)
	    if (chooser.getCurrentDirectory().toString().startsWith("/") ) {
		filePath = chooser.getCurrentDirectory() + "/" + chooser.getSelectedFile().getName();
	    }
	    else {
		filePath = chooser.getCurrentDirectory() + "\\" + chooser.getSelectedFile().getName();
	    }
	    
	    try {
                
                // Strip off .png
                if (filePath.endsWith(".png")) {
                    filePath = filePath.substring(0, filePath.length() - 4);
                }

                int oldSliceIndex = currentSliceIndex;

                for (int i = 0; i < noSlices; i++) {
                    
                    currentSliceIndex = i;

                    sliceList.setSelectedIndex(currentSliceIndex);
                    sliceList.ensureIndexIsVisible(currentSliceIndex);

                    navPanel.setSlice();
                    setZoomedSlice();

                    // doesn't seem to update UI
                    revalidate();
                    repaint();

                    BufferedImage image = new BufferedImage(getSize().width, getSize().height, BufferedImage.TYPE_INT_RGB);
                    paint(image.createGraphics());

                    DecimalFormat df = new DecimalFormat("000");

                    File outputFile = new File(filePath + df.format(currentSliceIndex) + ".png");
                    javax.imageio.ImageIO.write(image, "png", outputFile);
                    
                }
 
                // reset slice index
                sliceList.setSelectedIndex(oldSliceIndex);
                sliceList.ensureIndexIsVisible(oldSliceIndex);

	    }
	    catch (IOException e) {
	        logger.warning("Could not write file (exception follows)");
	        LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
	    }

	}
    }

    
    public void writeRGB() {
	
	String fileRoot = "";

	JFileChooser chooser = new JFileChooser();
	
	javax.swing.filechooser.FileFilter metaFilter = getFileFilter("Meta IO (*.mha)","mha");
	javax.swing.filechooser.FileFilter vtkFilter = getFileFilter("VTK (*.vtk)","vtk");
	javax.swing.filechooser.FileFilter niiFilter = getFileFilter("NIfTI (*.nii)","nii");
	javax.swing.filechooser.FileFilter niiGZ_Filter = getFileFilter("NIfTI (*.nii.gz)","nii.gz");

	chooser.setAcceptAllFileFilterUsed(false);
	chooser.addChoosableFileFilter(metaFilter);
	chooser.addChoosableFileFilter(vtkFilter);
	chooser.addChoosableFileFilter(niiFilter);
	chooser.addChoosableFileFilter(niiGZ_Filter);

	chooser.setFileFilter(metaFilter);

	try {
	    File currentDir = new File(new File(".").getCanonicalPath());
	    
	    chooser.setCurrentDirectory(currentDir);
	}
	catch (IOException e) {
	    throw new LoggedException(e);
	}
	    
	String filePath = null;
	
	int r = chooser.showSaveDialog(this);

	if(r == JFileChooser.APPROVE_OPTION) {
	    

	    // Assuming *nix systems will use absolute paths (ie, they will start with /), seems to work
            // Careful about using system file separators, because we might need to use forward slash
            // on Windows (Cygwin)
	    if (chooser.getCurrentDirectory().toString().startsWith("/") ) {
		filePath = chooser.getCurrentDirectory() + "/" + chooser.getSelectedFile().getName();
	    }
	    else {
		filePath = chooser.getCurrentDirectory() + "\\" + chooser.getSelectedFile().getName();
	    }
	    
	    try {

		    
		javax.swing.filechooser.FileFilter formatFilter = chooser.getFileFilter();
		
		if (formatFilter == metaFilter) {
		    if (!(filePath.endsWith(".mha"))) {
			filePath += ".mha";
		    }
		}
		else if (formatFilter == vtkFilter) {
		    if (!(filePath.endsWith(".vtk"))) {
			filePath += ".vtk";
		    }
		}
		else if (formatFilter == niiFilter) {
		    if (!(filePath.endsWith(".nii"))) {
			filePath += ".nii";
		    }
		}
		else if (formatFilter == niiGZ_Filter) {
		    if (!(filePath.endsWith(".nii.gz"))) {
			filePath += ".nii.gz";
		    }
		}
		else {
		    throw new LoggedException("Unknown file filter");
		}
		
		image.writeImage(filePath);
	    }
	    catch (IOException e) {
	        logger.warning("Could not write file (exception follows)");
	        LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
	    }

	}
    }
    
    public void writeScheme() {

	CL_Initializer.initImagingScheme();
	DW_Scheme newScheme = CL_Initializer.imPars;

	//need to update flips/order
	// note: must reorder first!!
        newScheme = newScheme.gradOrder(currentXYZ_Order);
	
        if(flipDirs[0]==-1)
	    newScheme = newScheme.flipX();
	else if(flipDirs[1]==-1)
            newScheme = newScheme.flipY();
	else if(flipDirs[2]==-1)
            newScheme = newScheme.flipZ();
	

	//set up file chooser
	String fileRoot = "";

	JFileChooser chooser = new JFileChooser();
	
	javax.swing.filechooser.FileFilter schemeFilter = getFileFilter("Scheme file (*.scheme)","scheme");

	chooser.setAcceptAllFileFilterUsed(true);
	chooser.addChoosableFileFilter(schemeFilter);

	chooser.setFileFilter(schemeFilter);

	try {
	    File currentDir = new File(new File(".").getCanonicalPath());
	    
	    chooser.setCurrentDirectory(currentDir);
	}
	catch (IOException e) {
	    throw new LoggedException(e);
	}
	    
	String filePath = null;
	
	int r = chooser.showSaveDialog(this);

	if(r == JFileChooser.APPROVE_OPTION) {

            // Cygwin uses / but Windows Java will report \ as separator. Thus we always use /
	    String slashie = "/";
	    
	    filePath = chooser.getCurrentDirectory() + slashie + chooser.getSelectedFile().getName();
	    
	    try {

		javax.swing.filechooser.FileFilter formatFilter = chooser.getFileFilter();

		if(!chooser.getSelectedFile().getName().contains(".")) {
		    filePath += ".scheme";
		}

		FileWriter fw = new FileWriter(filePath);
		BufferedWriter br = new BufferedWriter(fw);

		//write schemefile!!
 		br.write(newScheme.toString());
 		br.close();

		isFlippedDir = false;
		isSwappedDirs = false;

		flipDirs[0]=1;
		flipDirs[1]=1;
		flipDirs[2]=1;

		currentXYZ_Order = DW_Scheme.gradXYZ;

		schemePanel.resetButtons();
		CL_Initializer.schemeFile=filePath;
		CL_Initializer.initImagingScheme();
		logger.info("Saved scheme file (Note that pds/dt eigs have not been updated.  You must recalculate these using the updated schemefile).");
 		}
 	    catch (IOException e) {
 	        logger.warning("Could not write file (exception follows)");
 	        LoggedException.logExceptionWarning(e, Thread.currentThread().getName());
	    }

	}
    }

    // sets up a file filter
    private javax.swing.filechooser.FileFilter getFileFilter(final String name, final String ext) {
	return new javax.swing.filechooser.FileFilter() {  
		public boolean accept(File f) {
		    return f.getName().toLowerCase().endsWith("." + ext) || f.isDirectory();
		}
		
		public String getDescription() {  
		    return name;
		}
	    };
    }



    /**
     * Should be called by instantiating constructors.
     *
     */
    protected void display() {

	addMouseMotionListener(new MouseInputAdapter() {
		
		public void mouseMoved(MouseEvent e) {
		    int x = e.getX();
		    int y = e.getY();

		    PD_OrientationViewer.this.navPanel.setText(PD_OrientationViewer.this.getText(x,y));
		}
	    });


	navPanel.resetZoomPosition();


        frame = new JFrame("PD Orientation Viewer");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	frame.addWindowListener(new WindowListener(){
		public void windowClosing(WindowEvent e){
		    if(isFlippedDir || isSwappedDirs){
			Object [] options = { "Yes", "No"};
			int n = JOptionPane.showOptionDialog(frame, "Schemefile has not been saved.  Save now?", "Save scheme file?",
							     JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);

			if(n==0) // Yes -> save dialog
			    writeScheme();
		    }
		}

		public void windowClosed(WindowEvent e) {}
		
		public void windowIconified(WindowEvent e) {}
				
		public void windowOpened(WindowEvent e) {}
	
		public void windowDeiconified(WindowEvent e) {}

		public void windowActivated(WindowEvent e) {}
		
		public void windowDeactivated(WindowEvent e) {}
		
		public void windowGainedFocus(WindowEvent e) {}
		
		public void windowLostFocus(WindowEvent e) {}
		
		public void windowStateChanged(WindowEvent e) {}
		
	    });

	
        sliceOrientation = XY;

        initSliceSelectionPanel();

        navPanel.setSlice();
        setZoomedSlice();
        
	// contains scheme panel, nav panel
	JPanel topPanel = new JPanel();

	topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));

	// add scheme panel if scheme file is specified
	if (CL_Initializer.schemeFile != null) {
	    topPanel.add(schemePanel);
	}

	topPanel.add(navPanel);

        //Add stuff to the frame
        frame.getContentPane().add(this, BorderLayout.CENTER);
        frame.getContentPane().add(sliceSelectionPanel, BorderLayout.EAST);
        frame.getContentPane().add(topPanel, BorderLayout.NORTH);

        //Display the window.
        frame.pack();
        frame.setSize(700, 700);
        frame.setVisible(true);
    }

    /** Called to update an existing display. Updates the zoom selection and the zoomed slice. */ 
    public void updateDisplay() {
        navPanel.setSlice();
        setZoomedSlice();
    } 


    public static void main(String[] args) {

        CL_Initializer.inputDataType = "double";
        CL_Initializer.maxTensorComponents = 1;

	CL_Initializer.numPDsIO = -1;

	CL_Initializer.inputModel = "dteig";

        CL_Initializer.CL_init(args);
        
        int xS = 0;
        int yS = 0;
        int zS = 0;

	// default voxel dims are 1.0
	// in the future, may scale display for anisotropic voxels
	// currently only used for writing RGB volume
        double xV = 1.0;
        double yV = 1.0;
        double zV = 1.0;

        double scalarThresh = 0.0;

        double minScalar = 0.0;
        double maxScalar = 0.0;

        String scalarFile = null;

	boolean rgb = true;

	String picoPDF = "bingham";	

	// eigenvector index, 0 - 2, where 0 displays e1
	int eigIndex = 0;
	
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-scalarthresh")) { 
                scalarThresh = Double.parseDouble(args[i+1]);
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-scalarrange")) { 
                minScalar = Double.parseDouble(args[i+1]);
                maxScalar = Double.parseDouble(args[i+2]);
                CL_Initializer.markAsParsed(i, 3);
            }
            if (args[i].equals("-scalarfile")) {
                scalarFile = args[i+1];
                CL_Initializer.markAsParsed(i, 2);
            }
            if (args[i].equals("-norgb")) {
                rgb = false;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-e1")) {
                eigIndex = 0;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-e2")) {
                eigIndex = 1;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-e3")) {
                eigIndex = 2;
                CL_Initializer.markAsParsed(i);
            }
            if (args[i].equals("-pdf")) {
                picoPDF = args[i+1];
                CL_Initializer.markAsParsed(i, 2);
            }

        }


        CL_Initializer.checkParsing(args);

        if ( (CL_Initializer.maxTensorComponents > 2 && CL_Initializer.inputModel.equals("dteig")) || CL_Initializer.numPDsIO > 2 ) {
            logger.info("This program displays a maximum of 2 PDs per voxel. Display shows the first two PDs in each voxel.");
        }

	// set num PDs if default not changed
	// max tensor components set separately
	if (CL_Initializer.numPDsIO < 0) {
	    if (CL_Initializer.inputModel.equals("pds")) {
		// default is 3 for sfpeaks
		CL_Initializer.numPDsIO = 3;
	    }
	    else {
		// default is 1 for everything else
		CL_Initializer.numPDsIO = 1;
	    }
	}
        
        if (ImageHeader.imageExists(scalarFile)) {
            // image space has to be the same as scalar space
            CL_Initializer.headerTemplateFile = scalarFile;
        }
        
        CL_Initializer.initInputSpaceAndHeaderOptions();

	xS = CL_Initializer.dataDims[0];
	yS = CL_Initializer.dataDims[1];
	zS = CL_Initializer.dataDims[2];

	xV = CL_Initializer.voxelDims[0];
	yV = CL_Initializer.voxelDims[1];
	zV = CL_Initializer.voxelDims[2];

	
	double[][][] scalarVol = null;

	if (scalarFile != null) {
	    
	    DataSource scalars = null;
	    
	    if (ImageHeader.imageExists(scalarFile)) {
                
		ImageHeader ih = null;
		
		try {
		    ih = ImageHeader.readHeader(scalarFile);
		}
		catch (IOException e) {
		    throw new LoggedException(e);
                    
		}
		// get data dims, voxel dims, from header
		int[] dataDims = new int[] {xS, yS, zS};
		double[] voxelDims = new double[] {xV, yV, zV};

                if (!ih.sameDimensions(dataDims, voxelDims)) {
                    throw new LoggedException("Scalar file has different dimensions to data");
                }

		scalars = ih.getImageDataSource();
	    }
	    else {
		scalars = new VoxelOrderDataSource(scalarFile, 1, CL_Initializer.inputDataType);
	    }

	    scalarVol = new double[xS][yS][zS];
            
	    
	    for (int k = 0; k < zS; k++) { 
		for (int j = 0; j < yS; j++) {
		    for (int i = 0; i < xS; i++) {
			scalarVol[i][j][k] = scalars.nextVoxel()[0];
		    }
		}
	    }
	}


	RGB_ScalarImage image = null;

        if (CL_Initializer.inputModel.equals("dteig")) {

            // Test that all the required information was provided on the
            // command line.
            if (xS <= 0 || yS <= 0 || zS <= 0) {
                throw new misc.LoggedException("Cannot create image without data dimensions");
            }

            DataSource data = ExternalDataSource.getDataSource(CL_Initializer.inputFile, 12 * CL_Initializer.maxTensorComponents, CL_Initializer.inputDataType);

	    image = RGB_ScalarImage.imageFromTensorEigenSys(data, CL_Initializer.headerTemplate, scalarVol,
							    minScalar, maxScalar, eigIndex);

	}
        else if (CL_Initializer.inputModel.equals("pds")) {

            // Test that all the required information was provided on the
            // command line.
            if (xS <= 0 || yS <= 0 || zS <= 0) {
                throw new misc.LoggedException("Cannot create image without data dimensions");
            }

            DataSource data = ExternalDataSource.getDataSource(CL_Initializer.inputFile, 6 + 8 * CL_Initializer.numPDsIO, CL_Initializer.inputDataType);


	    image = RGB_ScalarImage.imageFromSphFuncPDs(data, CL_Initializer.headerTemplate, scalarVol,
							    minScalar, maxScalar);

        }
        else if (CL_Initializer.inputModel.equals("pico")) {

            // Test that all the required information was provided on the
            // command line.
            if (xS <= 0 || yS <= 0 || zS <= 0) {
                throw new misc.LoggedException("Cannot create image without data dimensions");
            }

	    int paramsPerPD = 1;

	    if (picoPDF.equals("bingham")) {
		paramsPerPD = 2;
	    }
	    if (picoPDF.equals("acg")) {
		paramsPerPD = 3;
	    }


            DataSource data = ExternalDataSource.getDataSource(CL_Initializer.inputFile, 
                                                               1 + (10 + paramsPerPD) * CL_Initializer.numPDsIO, 
                                                               CL_Initializer.inputDataType);


	    image = RGB_ScalarImage.imageFromPICoPDFs(data, CL_Initializer.numPDsIO, paramsPerPD,
                                                      CL_Initializer.headerTemplate, scalarVol,
						      minScalar, maxScalar, eigIndex);

        }
        else if (CL_Initializer.inputModel.equals("ballstick")) {

            // Test that all the required information was provided on the
            // command line.
            if (xS <= 0 || yS <= 0 || zS <= 0) {
                throw new misc.LoggedException("Cannot create image without data dimensions");
            }

            DataSource data = ExternalDataSource.getDataSource(CL_Initializer.inputFile, 7, 
                                                               CL_Initializer.inputDataType);


	    image = RGB_ScalarImage.imageFromBallStick(data, CL_Initializer.headerTemplate,
                                                       scalarVol, minScalar, maxScalar);

        }

        else if (CL_Initializer.inputModel.equals("compartment")) {

            // Test that all the required information was provided on the
            // command line.
            if (xS <= 0 || yS <= 0 || zS <= 0) {
                throw new misc.LoggedException("Cannot create image without data dimensions");
            }

	    image = RGB_ScalarImage.imageFromCompartmentModel(CL_Initializer.headerTemplate, scalarVol,
							    minScalar, maxScalar);
        }

        else {
            throw new misc.LoggedException("Unrecognized input model " + CL_Initializer.inputModel);
        }
        

	if (rgb) {
	    image.setRGB_Gamma(1.0);
	}
	else {
	    image.setRGB_Gamma(0.0);
	}
	
	// create viewer of appropriate type
	if (CL_Initializer.inputModel.equals("dteig")) {
	    TensorOrientationViewer viewer = new TensorOrientationViewer(image, scalarThresh);
	}
        else if (CL_Initializer.inputModel.equals("pds")) {
	    SphFuncPeaksOrientationViewer viewer = new SphFuncPeaksOrientationViewer(image, scalarThresh);
	}
        else if (CL_Initializer.inputModel.equals("pico")) {
	    PICoPDF_OrientationViewer viewer = new PICoPDF_OrientationViewer(image, scalarThresh);
	}
        else if (CL_Initializer.inputModel.equals("ballstick")) {
	    BallSticksOrientationViewer viewer = new BallSticksOrientationViewer(image, scalarThresh);
        }
        else if (CL_Initializer.inputModel.equals("compartment")) {
	    CompartmentModelOrientationViewer viewer = new CompartmentModelOrientationViewer(image, scalarThresh);
	}

    }

    /**
     * Sets the slice orientation to XY, XZ, or YZ, then updates the slice list
     * and display.
     *  
     */
    public void setSliceOrientation(int orientation) {
        sliceOrientation = orientation;
        navPanel.resetZoomPosition();
        setSliceListModel(); // causes set slice
	repaint();
    }

    public void swapDirs(int [] swapOrder) {
	if(swapOrder!= DW_Scheme.gradXYZ)
	    isSwappedDirs = true;
	else
	    isSwappedDirs = false;

	if(currentXYZ_Order != DW_Scheme.gradXYZ) // don't need to reorder dirs if they're already in [x y z] order
	    revertGrads(currentXYZ_Order);

	reorderGrads(swapOrder);
	currentXYZ_Order = swapOrder;

	if(flipDirs[0]==-1)
	    schemePanel.updateFlipButtons(FLIP_X);
	else if(flipDirs[1]==-1)
	    schemePanel.updateFlipButtons(FLIP_Y);
	else if(flipDirs[2]==-1)
	    schemePanel.updateFlipButtons(FLIP_Z);

	image.calculateRGB();
	
	setSliceListModel(); // causes set slice
	updateDisplay();
	repaint();
    }

    private void reorderGrads(int [] reOrder) {
	int numVectsPerVox = 0;
	double [] vec = new double [3];
	for(int i=0;i<image.xDataDim;i++){
	    for(int j=0;j<image.yDataDim;j++){
	        for(int k=0;k<image.zDataDim;k++){
		    numVectsPerVox=image.vectors[i][j][k].length;
		    for(int l=0;l<numVectsPerVox;l++){

			vec[0]=image.vectors[i][j][k][l].x;
			vec[1]=image.vectors[i][j][k][l].y;
			vec[2]=image.vectors[i][j][k][l].z;

			image.vectors[i][j][k][l] = new Vector3D(vec[reOrder[0]],vec[reOrder[1]],vec[reOrder[2]]);
		    }
		}
	    }
	}

	int tmpFlipX = flipDirs[0];
	int tmpFlipY = flipDirs[1];
	int tmpFlipZ = flipDirs[2];
	flipDirs[reOrder[0]]=tmpFlipX;
	flipDirs[reOrder[1]]=tmpFlipY;
	flipDirs[reOrder[2]]=tmpFlipZ;
    }

    private void revertGrads(int [] reOrder) {
	int numVectsPerVox = 0;
	double [] vec = new double [3];
	double x=0;
	double y=0;
	double z=0;
	int [] revertOrder = new int[3];
	for(int m=0;m<3;m++) {
	    if(reOrder[m]==0)
	    	revertOrder[0]=m;		    
	    else if(reOrder[m]==1)
		revertOrder[1]=m;
	    else if(reOrder[m]==2)
		revertOrder[2]=m;
	}

	for(int i=0;i<image.xDataDim;i++){
	    for(int j=0;j<image.yDataDim;j++){
	        for(int k=0;k<image.zDataDim;k++){
		    numVectsPerVox=image.vectors[i][j][k].length;
		    for(int l=0;l<numVectsPerVox;l++){

			vec[0]=image.vectors[i][j][k][l].x;
			vec[1]=image.vectors[i][j][k][l].y;
			vec[2]=image.vectors[i][j][k][l].z;
			
			image.vectors[i][j][k][l] = new Vector3D(vec[revertOrder[0]], vec[revertOrder[1]], vec[revertOrder[2]]);
		    }
		}
	    }
	}
	int tmpFlipX = 0;
	int tmpFlipY = 0;
	int tmpFlipZ = 0;

	tmpFlipX=flipDirs[reOrder[0]];
	tmpFlipY=flipDirs[reOrder[1]];
	tmpFlipZ=flipDirs[reOrder[2]];
	flipDirs[0]=tmpFlipX;
	flipDirs[1]=tmpFlipY;
	flipDirs[2]=tmpFlipZ;
    }
    /**
     * flips the x, y or z direction
     *  
     */
    public void setFlip(int flipDirection) {

	int numVectsPerVox = 0;
	double x, y, z;
	
	switch(flipDirection){
        case FLIP_X:
	    flipDirs[0]=-1 * flipDirs[0];
	    isFlippedDir=true;
	    break;
	case FLIP_Y:
	    flipDirs[1]=-1 * flipDirs[1];
	    isFlippedDir =true;
	    break;
	case FLIP_Z:
	    flipDirs[2]=-1 * flipDirs[2];
	    isFlippedDir =true;
	    break;
	case NO_FLIP:
	    isFlippedDir =false;

	    break;
	}

	for(int i=0;i<image.xDataDim;i++){
	    for(int j=0;j<image.yDataDim;j++){
		for(int k=0;k<image.zDataDim;k++){
		    numVectsPerVox=image.vectors[i][j][k].length;
		    for(int l=0;l<numVectsPerVox;l++){
			x = flipDirs[0] * image.vectors[i][j][k][l].x;
			y = flipDirs[1] * image.vectors[i][j][k][l].y;
			z = flipDirs[2] * image.vectors[i][j][k][l].z;
			image.vectors[i][j][k][l] = new Vector3D(x, y, z);
		    }
		}
	    }
	}

	flipDirs[0]=1;
	flipDirs[1]=1;
	flipDirs[2]=1;
	switch(flipDirection){
        case FLIP_X:
	    flipDirs[0]=-1 * flipDirs[0];
	    break;
	case FLIP_Y:
	    flipDirs[1]=-1 * flipDirs[1];
	    break;
	case FLIP_Z:
	    flipDirs[2]=-1 * flipDirs[2];
	    break;
	case NO_FLIP:
	    //	    isNewScheme=false;
	    break;
	}
	image.calculateRGB();
	repaint();
	
    }


    /**
     * Initialise the slice selection panel.
     */
    private void initSliceSelectionPanel() {
        sliceList = new JList<>();
        sliceSelectionPanel = new JPanel(new BorderLayout());

        setSliceListModel();

        // needed first time only
        currentSliceIndex = sliceList.getSelectedIndex();

        sliceList.setBackground(Color.white);
        sliceList.setFixedCellWidth(22);
        sliceList.addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                if (sliceList.getSelectedIndex() > -1) {
                    // calls to setSliceListModel() will cause -1 to be selected
                    // this traps that situation and allows setSliceListModel()
                    // to correct it
                    currentSliceIndex = sliceList.getSelectedIndex();
                    navPanel.setSlice();
                    setZoomedSlice();
                }
            }
        });

        //	sliceList.setToolTipText("slice selection");
        // forces single selection only
        sliceList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        sliceList.setPrototypeCellValue(new Integer(10000));
        JScrollPane scrollPane = new JScrollPane(sliceList);

        sliceSelectionPanel.add(scrollPane);

    }

    /**
     * Sets up the list of slices in the current plane view.
     */
    protected void setSliceListModel() {

        DefaultListModel<Integer> sliceModel = new DefaultListModel<>();
        
        noSlices = 0;
        
        switch (sliceOrientation) {
        case XY:
            noSlices = zSize;
            break;
        case XZ:
            noSlices = ySize;
            break;
        case YZ:
            noSlices = xSize;
            break;

        default:
            throw new java.lang.IllegalStateException("Unrecognised plane index "
                    + sliceOrientation);

        }

        for (int i = 0; i < noSlices; i++) {
            sliceModel.addElement(new Integer(i));
        }

        sliceList.setModel(sliceModel);
        sliceList.setSelectedIndex(noSlices / 2);
        sliceList.ensureIndexIsVisible(noSlices / 2);
        
    }



    /**
     * Called if zoomed section is moved, the slice is changed, or the plane view is switched.
     */
    public void setZoomedSlice() {
        int sl = currentSliceIndex;

	int[] bounds = navPanel.getZoomedRegion();


	// bounds in data coordinates. 
	int xMinVoxel = bounds[0];

	int xMaxVoxel = bounds[1];

	int yMinVoxel = bounds[2];

	int yMaxVoxel = bounds[3];


        int zoomedXSize = xMaxVoxel - xMinVoxel + 1;
        int zoomedYSize = yMaxVoxel - yMinVoxel + 1;


        zoomedSlice = new BufferedImage(zoomedXSize, zoomedYSize, BufferedImage.TYPE_INT_RGB);

        switch (sliceOrientation) {
            
        case XY:
            for (int i = 0; i < zoomedXSize; i++) {
                for (int j = 0; j < zoomedYSize; j++) {
		    
		    // java draws top down
		    zoomedSlice.setRGB(i, zoomedYSize - j - 1, 
				       image.rgbIndex(xMinVoxel + i, yMinVoxel + j, sl)); 
                    
                }
            }

            break;
        case XZ:
            for (int i = 0; i < zoomedXSize; i++) {
                for (int j = 0; j < zoomedYSize; j++) {
		    zoomedSlice.setRGB(i, zoomedYSize - j - 1, 
				       image.rgbIndex(xMinVoxel + i, sl, yMinVoxel + j)); 
                }
            }

            break;
        case YZ:
            for (int i = 0; i < zoomedXSize; i++) {
                for (int j = 0; j < zoomedYSize; j++) {
		    zoomedSlice.setRGB(i, zoomedYSize - j - 1, 
				       image.rgbIndex(sl, xMinVoxel + i, yMinVoxel + j)); 
                }
            }

            break;

        default:
            throw new java.lang.IllegalStateException("Unrecognised plane index " + sliceOrientation);
            
        }
        
        repaint();
        
    }



    /**
     * Draws the image and PDs.
     *  
     */
    public void paintComponent(Graphics g) {

        super.paintComponent(g);

        Dimension panelSize = getSize(); // should be the size of the panel
        // the BorderLayout should size this panel to fit the window

        // resize image so that it is as big as it could be, given the window
        // size
        // but maintain aspect ratio

        // x and y here refer to screen dimensions, nothing to do with slice
        // orientation

	
        double xScale = panelSize.getWidth() / zoomedSlice.getWidth();
        double yScale = panelSize.getHeight() / zoomedSlice.getHeight();

        if (xScale < yScale) {
            yScale = xScale;
        }
        else {
            xScale = yScale;
        }

	// integer scaling allows us to tell where we are in the zoomed region
        Graphics2D g2 = (Graphics2D) g;
        g2.scale(xScale, yScale);
	g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);

        g2.drawImage(zoomedSlice, 0, 0, this);

        g2.scale(1.0 / xScale, 1.0 / yScale);

	double spacing = xScale;

	int[] bounds = navPanel.getZoomedRegion();


	// bounds in data coordinates
	int xMinVoxel = bounds[0];

	int xMaxVoxel = bounds[1];

	int yMinVoxel = bounds[2];

	int yMaxVoxel = bounds[3];

	if (showVectors) {

	    Color pd1Color = Color.red;
	    Color pd2Color = Color.blue;
	    
	    if (image.rgbGamma() > 0.0) {
		pd1Color = Color.white;
		pd2Color = Color.white;
	    }

	    g2.setColor(pd1Color);

	    int xZoomSize = xMaxVoxel - xMinVoxel + 1;
	    int yZoomSize = yMaxVoxel - yMinVoxel + 1;

	    for (int i = 0; i < xZoomSize; i++) {
		for (int j = 0; j < yZoomSize; j++) {
		
		    // screen coordinates
		    int centreX = (int)(spacing / 2 + i * spacing);
		    int centreY = (int)(spacing / 2 + j * spacing);

		    int x1 = 0;
		    int y1 = 0;
		    int x2 = 0;
		    int y2 = 0;
                
		    double magnitude = spacing;
		    Vector3D vec = null;

		    switch (sliceOrientation) {
		    case XY:

			if (image.normScalarVol[i + xMinVoxel][yMaxVoxel - j][currentSliceIndex] > 
			    normScalarThresh) {
			    
			    vec = vectors[i + xMinVoxel][yMaxVoxel - j][currentSliceIndex][0];

			    // sign is different because Java's y axis starts at the top of the screen
			    x1 = centreX - (int)(magnitude * vec.x / 2.0);
			    y1 = centreY + (int)(magnitude * vec.y / 2.0);
                        
			    x2 = centreX + (int)(magnitude * vec.x / 2.0);
			    y2 = centreY - (int)(magnitude * vec.y / 2.0);

			    g2.drawLine(x1, y1, x2, y2);
                      
			    if (vectors[i + xMinVoxel][yMaxVoxel - j][currentSliceIndex].length == 2) {
				vec = vectors[i + xMinVoxel][yMaxVoxel - j][currentSliceIndex][1];

				if (vec.mod() > 0) {
				    g2.setColor(pd2Color);
				    x1 = centreX - (int)(magnitude * vec.x / 2.0);
				    y1 = centreY + (int)(magnitude * vec.y / 2.0);
                                
				    x2 = centreX + (int)(magnitude * vec.x / 2.0);
				    y2 = centreY - (int)(magnitude * vec.y / 2.0);
                                
				    g2.drawLine(x1, y1, x2, y2);
				    g2.setColor(pd1Color); 
				}

			    }
                        
			}
			break;
                    
		    case XZ:


			if (image.normScalarVol[i + xMinVoxel][currentSliceIndex][yMaxVoxel - j] > 
			    normScalarThresh) {

			    vec = vectors[i + xMinVoxel][currentSliceIndex][yMaxVoxel - j][0];
                        
			    // sign is different because Java's y axis starts at the top of the screen
			    x1 = centreX - (int)(magnitude * vec.x / 2.0);
			    y1 = centreY + (int)(magnitude * vec.z / 2.0);
                        
			    x2 = centreX + (int)(magnitude * vec.x / 2.0);
			    y2 = centreY - (int)(magnitude * vec.z / 2.0);

			    g2.drawLine(x1, y1, x2, y2);                        

			    if (vectors[i + xMinVoxel][currentSliceIndex][yMaxVoxel - j].length == 2) {
				vec = vectors[i + xMinVoxel][currentSliceIndex][yMaxVoxel - j][1];

				if (vec.mod() > 0) {
				    g2.setColor(pd2Color);

				    x1 = centreX - (int)(magnitude * vec.x / 2.0);
				    y1 = centreY + (int)(magnitude * vec.z / 2.0);
                                
				    x2 = centreX + (int)(magnitude * vec.x / 2.0);
				    y2 = centreY - (int)(magnitude * vec.z / 2.0);
                                
				    g2.drawLine(x1, y1, x2, y2);
				    g2.setColor(pd1Color); 
				}

			    }

			}
			break;
                    
		    case YZ:
                    
			if (image.normScalarVol[currentSliceIndex][i + xMinVoxel][yMaxVoxel - j] > 
			    normScalarThresh) {
			    
			    vec = vectors[currentSliceIndex][i + xMinVoxel][yMaxVoxel - j][0];
                        
			    // sign is different because Java's y axis starts at the top of the screen
			    x1 = centreX - (int)(magnitude * vec.y / 2.0);
			    y1 = centreY + (int)(magnitude * vec.z / 2.0);
                        
			    x2 = centreX + (int)(magnitude * vec.y / 2.0);
			    y2 = centreY - (int)(magnitude * vec.z / 2.0);

			    g2.drawLine(x1, y1, x2, y2);

			    if (vectors[currentSliceIndex][i + xMinVoxel][yMaxVoxel - j].length == 2) {
				vec = vectors[currentSliceIndex][i + xMinVoxel][yMaxVoxel - j][1];

				if (vec.mod() > 0) {
				    g2.setColor(pd2Color);
                                
				    x1 = centreX - (int)(magnitude * vec.y / 2.0);
				    y1 = centreY + (int)(magnitude * vec.z / 2.0);
                                
				    x2 = centreX + (int)(magnitude * vec.y / 2.0);
				    y2 = centreY - (int)(magnitude * vec.z / 2.0);
                                
				    g2.drawLine(x1, y1, x2, y2);
				    g2.setColor(pd1Color); 
				}

			    }
			}
                    
			break;

		    default:
			throw new java.lang.IllegalStateException("Unrecognised plane index "
								  + sliceOrientation);

		    }

                
		
		}
	    }
        

	}

    }

    
}



/**
 * Panel shows the whole slice with a box outlining the zoomed region.
 *
 */
class ZoomPanel extends JPanel {

    // voxels contained in zoomed region, data coordinates
    private int zoomedXMin = 0;
    private int zoomedYMin = 0;

    private int zoomedXMax = 0;
    private int zoomedYMax = 0;

    // centre of box in data coordinates
    private int xBoxCentre = 0;
    private int yBoxCentre = 0;
    
    // corner of bounding box in screen coordinates
    private int zoomOutlineXMin = 0;
    private int zoomOutlineYMin = 0;
    
    private BufferedImage slice = null;
    
    private int zoomedXExtent = 0;
    private int zoomedYExtent = 0;
    
    private PD_OrientationViewer container = null;

    private double[][][] normScalarVol = null;

    protected int xSize;

    protected int ySize;

    protected int zSize;

    public ZoomPanel(PD_OrientationViewer container) {
        this.container = container;

        xSize = container.xSize;
        ySize = container.ySize;
        zSize = container.zSize;

	resetZoomPosition();

        //        zoomedXMax = preferredZoomedExtent - 1;
        //   zoomedYMax = preferredZoomedExtent - 1;
       
        normScalarVol = container.image.normScalarVol;
        setLayout(new GridLayout(2, 1));
        addMouseListener(mouseListener);
        addMouseMotionListener(mouseListener);
	setPreferredSize(new Dimension(xSize,ySize));
	setMinimumSize(new Dimension(xSize,ySize));

	addMouseWheelListener(container.wheelListener);
    }


    
    // listens for mouse clicks in the image
    private final MouseInputAdapter mouseListener = new MouseInputAdapter() {
            public void mouseClicked(MouseEvent e) {
                setBoxCentre(e);
                container.setZoomedSlice();
            }

            public void mouseDragged(MouseEvent e) {
                setBoxCentre(e);
                container.setZoomedSlice();
            }
        };        

    
    public void resetZoomPosition() {
        
        zoomedXExtent = container.preferredZoomedExtent;
        zoomedYExtent = container.preferredZoomedExtent;

        switch (container.sliceOrientation) {
        case PD_OrientationViewer.XY:
            
            if (xSize < zoomedXExtent) {
                zoomedXExtent = xSize;
            }
            if (ySize < zoomedYExtent) {
                zoomedYExtent = ySize;
            }

        
            zoomedXMin = 0;
            zoomedXMax = zoomedXMin + zoomedXExtent - 1;

	    zoomedYMin = (ySize - zoomedYExtent);
	    zoomedYMax = zoomedYMin + zoomedYExtent - 1;
	    
	    zoomOutlineXMin = zoomedXMin;
	    zoomOutlineYMin = (int)((ySize - zoomedYMax - 1));

            break;
        case PD_OrientationViewer.XZ:

            if (xSize < zoomedXExtent) {
                zoomedXExtent = xSize;
            }
            if (zSize < zoomedYExtent) {
                zoomedYExtent = zSize;
            }

            zoomedXMin = 0;
            zoomedXMax = zoomedXMin + zoomedXExtent - 1;

            zoomedYMin = (zSize - zoomedYExtent);
            zoomedYMax = zoomedYMin + zoomedYExtent - 1;
            
            zoomOutlineXMin = zoomedXMin;
            zoomOutlineYMin = (int)((zSize - zoomedYMax - 1));


            break;
        case PD_OrientationViewer.YZ:

            if (ySize < zoomedXExtent) {
                zoomedXExtent = ySize;
            }
            if (zSize < zoomedYExtent) {
                zoomedYExtent = zSize;
            }

            zoomedXMin = 0;
            zoomedXMax = zoomedXMin + zoomedXExtent - 1;

            zoomedYMin = (zSize - zoomedYExtent);
            zoomedYMax = zoomedYMin + zoomedYExtent - 1;
            
            zoomOutlineXMin = zoomedXMin;
            zoomOutlineYMin = (int)((zSize - zoomedYMax - 1));

            break;
        }

    }
   
    
    private void setBoxCentre(MouseEvent e) {
	setBoxCentre(e.getX(), e.getY());
    }

    protected int getXBoxCentre() {
	return xBoxCentre;
    }

    protected int getYBoxCentre() {
	return yBoxCentre;
    }

    /**
     * Sets the centre and hence the boundaries of the zoom box.
     *
     * @param ex the x coordinate of the box centre.
     * @param ey the y coordinate of the box centre.
     *
     */
    protected void setBoxCentre(int ex, int ey) {        

	xBoxCentre = ex;
	yBoxCentre = ey;

        zoomedXExtent = container.preferredZoomedExtent;
        zoomedYExtent = container.preferredZoomedExtent;

        switch (container.sliceOrientation) {
        case PD_OrientationViewer.XY:
            
            if (xSize < zoomedXExtent) {
                zoomedXExtent = xSize;
            }
            if (ySize < zoomedYExtent) {
                zoomedYExtent = ySize;
            }

            break;
        case PD_OrientationViewer.XZ:

            if (xSize < zoomedXExtent) {
                zoomedXExtent = xSize;
            }
            if (zSize < zoomedYExtent) {
                zoomedYExtent = zSize;
            }

            break;
        case PD_OrientationViewer.YZ:

            if (ySize < zoomedXExtent) {
                zoomedXExtent = ySize;
            }
            if (zSize < zoomedYExtent) {
                zoomedYExtent = zSize;
            }

            break;
        }


        Dimension panelSize = getSize(); // should be the size of the panel
        // the BorderLayout should size this panel to fit the window
                        
        // resize image so that it is as big as it could be, given the window
        // size
        // but maintain aspect ratio

        // x and y scale refer to screen dimensions, nothing to do with slice
        // orientation
                        
        double xScale = panelSize.getWidth() / slice.getWidth();
        double yScale = panelSize.getHeight() / slice.getHeight();


        if (xScale < yScale) {
            yScale = xScale;
        }
        else {
            xScale = yScale;
        }

        // in data coordinates
        int x = (int)Math.round(ex / xScale);
        
        int y = 0;

        if (container.sliceOrientation == PD_OrientationViewer.XY) {
            y = ySize - (int)Math.round(ey / yScale) - 1;
        }
        else {
            y = zSize - (int)Math.round(ey / yScale) - 1;
        }

        // box is described in data coordinates
        zoomedXMin = x - zoomedXExtent / 2 + 1;
        zoomedXMax = x + zoomedXExtent / 2;
        
        zoomedYMin = y - zoomedYExtent / 2 + 1;
        zoomedYMax = y + zoomedYExtent / 2;

        switch (container.sliceOrientation) {
        case PD_OrientationViewer.XY:

            if (zoomedXMin < 0) {
                zoomedXMin = 0;
                zoomedXMax = zoomedXExtent - 1;
            }
            else if (zoomedXMax >= xSize) {
                zoomedXMax = xSize - 1;
                zoomedXMin = zoomedXMax - zoomedXExtent + 1;
            }
            
            if (zoomedYMin < 0) {
                zoomedYMin = 0;
                zoomedYMax = zoomedYExtent - 1;
            }
            else if (zoomedYMax >= ySize) {
                zoomedYMax = ySize - 1;
                zoomedYMin = zoomedYMax - zoomedYExtent + 1;
            }

            zoomOutlineYMin = (int)((ySize - zoomedYMax - 1));                            

            break;
        case PD_OrientationViewer.XZ:

            if (zoomedXMin < 0) {
                zoomedXMin = 0;
                zoomedXMax = zoomedXExtent - 1;
            }
            else if (zoomedXMax >= xSize) {
                zoomedXMax = xSize - 1;
                zoomedXMin = zoomedXMax - zoomedXExtent + 1;
            }

            if (zoomedYMin < 0) {
                zoomedYMin = 0;
                zoomedYMax = zoomedYExtent - 1;
            }
            else if (zoomedYMax >= zSize) {
                zoomedYMax = zSize - 1;
                zoomedYMin = zoomedYMax - zoomedYExtent + 1;
            }

            zoomOutlineYMin = (int)((zSize - zoomedYMax - 1));

            break;
        case PD_OrientationViewer.YZ:

            if (zoomedXMin < 0) {
                zoomedXMin = 0;
                zoomedXMax = zoomedXExtent - 1;
            }
            else if (zoomedXMax >= ySize) {
                zoomedXMax = ySize - 1;
                zoomedXMin = zoomedXMax - zoomedXExtent + 1;
            }

            if (zoomedYMin < 0) {
                zoomedYMin = 0;
                zoomedYMax = zoomedYExtent - 1;
            }
            else if (zoomedYMax >= zSize) {
                zoomedYMax = zSize - 1;
                zoomedYMin = zoomedYMax - zoomedYExtent + 1;
            }

            zoomOutlineYMin = (int)((zSize - zoomedYMax - 1));

            break;
                            
        default:
            throw new java.lang.IllegalStateException("Unrecognised plane index "
                                                      + container.sliceOrientation);
                            
        }
                        
        zoomOutlineXMin = (int)(zoomedXMin);

        repaint();

    }


    protected void recentreBox() {
	setBoxCentre(xBoxCentre, yBoxCentre);
    }


    /**
     * Draws the image.
     *  
     */
    public void paintComponent(Graphics g) {
                        
        super.paintComponent(g);
                        
        Dimension panelSize = getSize(); // should be the size of the panel
        // the BorderLayout should size this panel to fit the window
                        
        // resize image so that it is as big as it could be, given the window
        // size
        // but maintain aspect ratio

        // x and y here refer to screen dimensions, nothing to do with slice
        // orientation
                        
        double xScale = panelSize.getWidth() / slice.getWidth();
        double yScale = panelSize.getHeight() / slice.getHeight();
                        
        if (xScale < yScale) {
            yScale = xScale;
        }
        else {
            xScale = yScale;
        }
                        
        Graphics2D g2 = (Graphics2D) g;
        g2.scale(xScale, yScale);
	g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);                        
        g2.drawImage(slice, 0, 0, this);

        g2.setColor(Color.yellow);

        g2.drawRect(zoomOutlineXMin, zoomOutlineYMin, 
                    zoomedXExtent, zoomedYExtent);

    }


    /**
     * Called if slice is changed or the plane view is switched.
     */
    public void setSlice() {
        int currentSliceIndex = container.currentSliceIndex;
            
        switch (container.sliceOrientation) {
                
        case PD_OrientationViewer.XY:

            slice = new BufferedImage(xSize, ySize, BufferedImage.TYPE_INT_RGB);
                
            for (int i = 0; i < xSize; i++) {
                for (int j = 0; j < ySize; j++) {
                    slice.setRGB(i, ySize - j - 1, 
				 container.image.rgbIndex(i,j,currentSliceIndex)); // java draws top
                    // down
                }
            }
                
            break;
        case PD_OrientationViewer.XZ:


            slice = new BufferedImage(xSize, zSize, BufferedImage.TYPE_INT_RGB);

            for (int i = 0; i < xSize; i++) {
                for (int j = 0; j < zSize; j++) {
                    slice.setRGB(i, zSize - j - 1, 
				 container.image.rgbIndex(i,currentSliceIndex, j));
                }
            }
                
            break;

        case PD_OrientationViewer.YZ:

            slice = new BufferedImage(ySize, zSize, BufferedImage.TYPE_INT_RGB);

            for (int i = 0; i < ySize; i++) {
                for (int j = 0; j < zSize; j++) {
                    slice.setRGB(i, zSize - j - 1, 
				 container.image.rgbIndex(currentSliceIndex, i, j));
                }
            }
                
            break;
                
        default:
            throw new java.lang.IllegalStateException("Unrecognised plane index " + container.sliceOrientation);
                
        }
            
        repaint();
            
    }

    public int[] getZoomedRegion() {
        return new int[] {zoomedXMin, zoomedXMax, zoomedYMin, zoomedYMax}; 
    }

  
    
                    
}




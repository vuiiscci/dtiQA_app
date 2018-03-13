package apps;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.text.DecimalFormat;
import javax.swing.event.*;
import java.io.*;
import java.util.*;

import tools.*;
import inverters.*;
import data.*;

/**
 * <dl>
 * 
 * <dt>Purpose:
 * 
 * <dd>Loads F-test data for interactive threshold selection.
 * 
 * <dt>Description:
 * 
 * <dd>This class provides simple visualization of the output of VoxelClassify.
 * The user can see the voxel classification resulting from different F-test
 * thresholds.
 * <dd>
 * <dd>The input data for this class is the output of VoxelClassify when the
 * -ftest option is not used.
 * 
 * </dl>
 * 
 * @see apps.VoxelClassify
 * @author Philip Cook
 * @version $Id: F_TestThresholdSelector.java,v 1.8 2005/07/25 15:39:47 ucacpco
 *          Exp $
 *  
 */
public class F_TestThresholdSelector extends JPanel {

    private JFrame frame;

    private JPanel thresholdSelectionPanel;

    private JPanel sliceSelectionPanel;

    private TopPanel topPanel;

    private JList<Integer> sliceList;

    private boolean autoUpdate = false;

    // set to true if any of the thresholds change
    private boolean updateNeeded = false;

    // controls threshold selection
    private FTestThresholdSelectionPanel thresholdPanel;

    // controls plane and zoom of view
    //   private ViewSelectionPanel viewPanel;
    // need to add this later

    // the image that is being displayed
    private BufferedImage slice;

    // used to set slice orientation
    private static final int XY = 100;

    private static final int XZ = 200;

    private static final int YZ = 300;

    private int sliceOrientation = -1;

    private double f1 = 0.0;

    private double f2 = 0.0;

    private double f3 = 0.0;

    private static final double DEFAULT_F1 = 1.0E-20;

    private static final double DEFAULT_F2 = 1.0E-6;

    private static final double DEFAULT_F3 = 1.0E-6;

    private double BACKGROUNDTHRESHOLD = 0.0;

    private double CSFTHRESHOLD = -1.0;

    private double[][][][] fStats;

    private int[][][] classes; // the actual data that is displayed

    private int xSize;

    private int ySize;

    private int zSize;

    private int order;

    // position of marker that controls slice positions
    private int xySliceMarker = 0;

    private int xzSliceMarker = 0;

    private int yzSliceMarker = 0;

    // used in working out which position a mouse click is in
    private double xScale = 0.0;

    private double yScale = 0.0;

    private double[] fractions = null; // fractions[o/2] contains fraction of

    // non-background voxels classified as
    // order o

    // listens for mouse clicks in the image
    private final MouseAdapter mouseClickListener = new MouseAdapter() {
        public void mouseClicked(MouseEvent e) {
            setSliceMarkers(e);
        }
    };

    /**
     * @param data
     *            from loadData.
     * @param ord
     *            max order.
     * @param backgroundThresh
     *            background threshold
     * @param csfThresh
     *            csf threshold.
     * @param f1
     *            f-test threshold (order 0/2).
     * @param f2
     *            f-test threshold (order 2/4).
     * @param f3
     *            f-test threshold (order 4/6).
     */
    public F_TestThresholdSelector(double[][][][] data, int ord, double backgroundThresh,
            double csfThresh) {
        this(data, CL_Initializer.maxOrder, CL_Initializer.BACKGROUNDTHRESHOLD,
                CL_Initializer.CSFTHRESHOLD, DEFAULT_F1, DEFAULT_F2, DEFAULT_F3);
    }

    /**
     * @param data
     *            from loadData.
     * @param ord
     *            max order.
     * @param backgroundThresh
     *            background threshold
     * @param csfThresh
     *            csf threshold.
     * @param f1
     *            f-test threshold (order 0/2).
     * @param f2
     *            f-test threshold (order 2/4).
     * @param f3
     *            f-test threshold (order 4/6).
     */
    public F_TestThresholdSelector(double[][][][] data, int ord, double backgroundThresh,
            double csfThresh, double f1, double f2, double f3) {

        this.f1 = f1;
        this.f2 = f2;
        this.f3 = f3;

        BACKGROUNDTHRESHOLD = backgroundThresh;
        CSFTHRESHOLD = csfThresh;

        fStats = data;
        xSize = fStats[0].length;
        ySize = fStats[0][0].length;
        zSize = fStats[0][0][0].length;
        order = ord;

        classes = new int[xSize][ySize][zSize];

        frame = new JFrame("F Test Threshold Selection");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        sliceOrientation = XY;
        initSliceSelectionPanel();

        addMouseListener(mouseClickListener);

        thresholdSelectionPanel = new FTestThresholdSelectionPanel(this);
        topPanel = new TopPanel(this);

        update();

        //Add stuff to the frame
        frame.getContentPane().add(this, BorderLayout.CENTER);
        frame.getContentPane().add(sliceSelectionPanel, BorderLayout.EAST);
        frame.getContentPane().add(thresholdSelectionPanel, BorderLayout.SOUTH);
        frame.getContentPane().add(topPanel, BorderLayout.NORTH);

        //Display the window.
        frame.pack();
        frame.setSize(600, 600);
        frame.setVisible(true);

    }

    public static void main(String[] args) {

	Locale.setDefault(Locale.UK);

        // The data file is always the first argument.
        //String inputFile = args[0];

        //CL_Initializer cl = new CL_Initializer(args);
        CL_Initializer.inputDataType = "double";
        CL_Initializer.CL_init(args);
        
        CL_Initializer.checkParsing(args);

	int xS = CL_Initializer.dataDims[0];
	int yS = CL_Initializer.dataDims[1];
	int zS = CL_Initializer.dataDims[2];

        // Test that all the required information was provided on the
        // command line.
        if (xS <= 0 || yS <= 0 || zS <= 0) {
            throw new RuntimeException(
                    "F_TestThresholdSelector -inputfile <data file> -datadims [x y z] [-order <max order>] [-bgthresh <voxels with q0 below this set to background] [-csfthresh <voxels with q0 above this set to order 0, default infinity>]");
        }

        int components = (CL_Initializer.maxOrder / 2)
                * (CL_Initializer.maxOrder / 2 + 1) / 2 + 2;

        double[][][][] data = loadData(CL_Initializer.inputFile, xS, yS, zS, components,
                CL_Initializer.inputDataType);

        if (CL_Initializer.f1 > 0.0) {
            double f1 = CL_Initializer.f1;
            double f2 = CL_Initializer.f2;
            double f3 = CL_Initializer.f3;

            F_TestThresholdSelector f = new F_TestThresholdSelector(data,
                    CL_Initializer.maxOrder, CL_Initializer.BACKGROUNDTHRESHOLD,
                    CL_Initializer.CSFTHRESHOLD, f1, f2, f3);
        }
        else {
            F_TestThresholdSelector f = new F_TestThresholdSelector(data,
                    CL_Initializer.maxOrder, CL_Initializer.BACKGROUNDTHRESHOLD,
                    CL_Initializer.CSFTHRESHOLD);
        }

    }

    /**
     * Loads the data into a 4D array. Requires doubles.
     * 
     * @param filename
     *            The name of the datafile.
     * 
     * @param xSize
     *            The number of voxels in the x-dimension of the data.
     * 
     * @param ySize
     *            The number of voxels in the y-dimension of the data.
     * 
     * @param zSize
     *            The number of voxels in the z-dimension of the data.
     * 
     * @param components
     *            The number of values in each voxel.
     * 
     * @return The array of data.
     */
    public static double[][][][] loadData(String filename, int xSize, int ySize,
            int zSize, int components, String datatype) {

        double[][][][] data = new double[components][xSize][ySize][zSize];

        DataSource source = ExternalDataSource.getDataSource(filename, components, datatype);

        try {
            for (int z = 0; z < zSize; z++) {
                for (int y = 0; y < ySize; y++) {
                    for (int x = 0; x < xSize; x++) {

                        double[] vox = source.nextVoxel();

                        for (int c = 0; c < components; c++) {
                            data[c][x][y][z] = vox[c];
                        }
                    }
                }
            }
        }
        catch (Exception e) {
            System.err.println(e);
            System.exit(1);
        }

        return data;
    }

    /**
     * Tell this object that the thresholds have changed and the image needs to
     * be recalculated.
     *  
     */
    private void setUpdateNeeded() {
        updateNeeded = true;

        if (autoUpdate) {
            update();
            //sets updateNeeded = false;
        }
    }

    /**
     * Returns status of auto update
     */
    public boolean getAutoUpdate() {
        return autoUpdate;
    }

    /**
     * Toggles automatic updating.
     */
    public void toggleAutoUpdate() {
        autoUpdate = !autoUpdate;
    }

    /**
     * Updates data and then display. Called when update button is pressed, or,
     * if auto updating is enabled, whenever the thresholds change.
     */
    public void update() {

        // trust if user says update is needed

        doModelSelection();
        setSlice();
        topPanel.setFractions(fractions);
        updateNeeded = false;
    }

    /**
     * Sets the thresholds for the F-test. Only supports orders up to 6 at
     * present.
     *  
     */
    public void setFTestThresholds(double f1, double f2, double f3) {
        this.f1 = f1;
        this.f2 = f2;
        this.f3 = f3;
        setUpdateNeeded();
    }

    /**
     * Sets the slice orientation to XY, XZ, or YZ, then updates the slice list
     * and display.
     *  
     */
    public void setSliceOrientation(int orientation) {
        sliceOrientation = orientation;
        setSliceListModel();
        setSlice();
    }

    /**
     * Called when user left-clicks in the image. Marks slices for easy
     * switching between planes.
     *  
     */
    public void setSliceMarkers(MouseEvent e) {

        if (e.getButton() == MouseEvent.BUTTON3) {
            xySliceMarker = 0;
            xzSliceMarker = 0;
            yzSliceMarker = 0;
            repaint();
            return;
        }

        int x = (int) Math.round(e.getX() / xScale);
        int y = (int) Math.round(e.getY() / yScale);

        switch (sliceOrientation) {
        case XY:
            xySliceMarker = sliceList.getSelectedIndex();
            xzSliceMarker = y < ySize - 1 ? ySize - y - 1 : 0;
            yzSliceMarker = x < xSize - 1 ? x : 0;
            break;
        case XZ:
            xzSliceMarker = sliceList.getSelectedIndex();
            xySliceMarker = y < zSize - 1 ? zSize - y - 1 : 0;
            yzSliceMarker = x < xSize - 1 ? x : 0;
            break;
        case YZ:
            yzSliceMarker = sliceList.getSelectedIndex();
            xySliceMarker = zSize - y - 1;
            xzSliceMarker = x < ySize - 1 ? x : 0;
            break;

        default:
            throw new java.lang.IllegalStateException("Unrecognised plane index "
                    + sliceOrientation);

        }

        repaint();

    }

    /**
     * Initialise the slice selection panel.
     */
    private void initSliceSelectionPanel() {
        sliceList = new JList<>();
        sliceSelectionPanel = new JPanel(new BorderLayout());

        setSliceListModel();
        sliceList.setBackground(Color.white);
        sliceList.setFixedCellWidth(22);
        sliceList.addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                if (sliceList.getSelectedIndex() > -1) {
                    // calls to setSliceListModel() will cause -1 to be selected
                    // this traps that situation and allows setSliceListModel()
                    // to correct it
                    setSlice();
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
        int noSlices = 0;
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

        switch (sliceOrientation) {
        case XY:
            sliceList.setSelectedIndex(xySliceMarker);
            break;
        case XZ:
            sliceList.setSelectedIndex(xzSliceMarker);
            break;
        case YZ:
            sliceList.setSelectedIndex(yzSliceMarker);
            break;

        default:
            throw new java.lang.IllegalStateException("Unrecognised plane index "
                    + sliceOrientation);

        }

    }

    /**
     * Called if slice is moved, the plane view is switched, or the data
     * changes.
     */
    public void setSlice() {
        int sl = sliceList.getSelectedIndex();

        switch (sliceOrientation) {
        case XY:
            slice = new BufferedImage(xSize, ySize, BufferedImage.TYPE_INT_RGB);

            //Initialise random image array.
            for (int i = 0; i < xSize; i++) {
                for (int j = 0; j < ySize; j++) {
                    int r = (int) (40 * (classes[i][j][sl] + 1));
                    int g = r;
                    int b = r;
                    int a = 128;

                    int index = b + 256 * (g + 256 * r);
                    slice.setRGB(i, ySize - 1 - j, index); // java draws top
                    // down
                }
            }

            break;
        case XZ:
            slice = new BufferedImage(xSize, zSize, BufferedImage.TYPE_INT_RGB);

            for (int i = 0; i < xSize; i++) {
                for (int j = 0; j < zSize; j++) {
                    int r = (int) (40 * (classes[i][sl][j] + 1));
                    int g = r;
                    int b = r;
                    int a = 128;

                    int index = b + 256 * (g + 256 * r);
                    slice.setRGB(i, zSize - 1 - j, index);
                }
            }

            break;
        case YZ:
            slice = new BufferedImage(ySize, zSize, BufferedImage.TYPE_INT_RGB);

            for (int i = 0; i < ySize; i++) {
                for (int j = 0; j < zSize; j++) {
                    int r = (int) (40 * (classes[sl][i][j] + 1));
                    int g = r;
                    int b = r;
                    int a = 128;

                    int index = b + 256 * (g + 256 * r);
                    slice.setRGB(i, zSize - 1 - j, index);
                }
            }

            break;

        default:
            throw new java.lang.IllegalStateException("Unrecognised plane index "
                    + sliceOrientation);

        }

        repaint();

    }

    /**
     * Does the voxel classification. Called from update. Shouldn't be called
     * directly.
     */
    private void doModelSelection() {
        int[] classSizes = new int[order / 2 + 1];
        int totalForeground = 0;
        for (int z = 0; z < zSize; z++) {
            for (int y = 0; y < ySize; y++) {
                for (int x = 0; x < xSize; x++) {

                    double exitCode = fStats[0][x][y][z];
                    double q0 = Math.exp(fStats[1][x][y][z]);

                    // Apply the background and CSF thresholds first.
                    if (exitCode < 0.0 || q0 < BACKGROUNDTHRESHOLD) {
                        classes[x][y][z] = -1;
                    }
                    else if (CSFTHRESHOLD > 0.0 && q0 > CSFTHRESHOLD) {
                        classes[x][y][z] = 0;
                    }
                    else {
                        double[] p = new double[fStats.length - 2];
                        for (int c = 0; c < p.length; c++) {
                            p[c] = fStats[c + 2][x][y][z];
                        }
                        classes[x][y][z] = EvenSphHarmFitter.selectModel(p, order, f1,
                                f2, f3);
                    }

                    if (classes[x][y][z] >= 0) {
                        classSizes[classes[x][y][z] / 2] += 1;
                        totalForeground += 1;
                    }

                }
            }
        }

        fractions = new double[classSizes.length];
        for (int i = 0; i < classSizes.length; i++) {
            fractions[i] = (double) classSizes[i] / (double) totalForeground;
        }
    }

    /**
     * Draws the image (and crosshairs, if needed).
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

        xScale = panelSize.getWidth() / slice.getWidth();
        yScale = panelSize.getHeight() / slice.getHeight();

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

        g2.setColor(Color.red);
        // draw lines
        //	g2.scale(1.0, 1.0);
        switch (sliceOrientation) {
        case XY:
            if (xzSliceMarker > 0) {
                g2.drawLine(0, (int) (ySize - 1 - xzSliceMarker),
                        (int) slice.getWidth() - 1, (int) (ySize - 1 - xzSliceMarker));
            }
            if (yzSliceMarker > 0) {
                g2.drawLine((int) (yzSliceMarker), 0, (int) (yzSliceMarker), (int) (slice
                        .getHeight() - 1));
            }
            break;
        case XZ:

            if (xySliceMarker > 0) {
                g2.drawLine(0, (int) ((zSize - 1 - xySliceMarker)), (int) (slice
                        .getWidth() - 1), (int) ((zSize - 1 - xySliceMarker)));
            }
            if (yzSliceMarker > 0) {
                g2.drawLine((int) (yzSliceMarker), 0, (int) (yzSliceMarker), (int) (slice
                        .getHeight() - 1));
            }

            break;
        case YZ:

            if (xySliceMarker > 0) {
                g2.drawLine(0, (int) ((zSize - 1 - xySliceMarker)), (int) (slice
                        .getWidth() - 1), (int) ((zSize - 1 - xySliceMarker)));
            }
            if (xzSliceMarker > 0) {
                g2.drawLine((int) (xzSliceMarker), 0, (int) (xzSliceMarker), (int) (slice
                        .getHeight()) - 1);

            }

            break;

        default:
            throw new java.lang.IllegalStateException("Unrecognised plane index "
                    + sliceOrientation);

        }

    }

    /**
     * Writes the current orders to a file.
     *  
     */
    private void saveVCImage(String file) {
        try {
            FileOutputStream fout = new FileOutputStream(file);
            DataOutputStream dout = new DataOutputStream(new BufferedOutputStream(fout,
                    1024 * 1024 * 4));

            for (int k = 0; k < zSize; k++) {
                for (int j = 0; j < ySize; j++) {
                    for (int i = 0; i < xSize; i++) {
                        dout.writeInt(classes[i][j][k]);
                    }
                }
            }

            dout.close();

        }
        catch (IOException e) {
            System.err.println(e);
        }

    }

    /**
     * Contains the scrollbars and text boxes that control the F-test
     * thresholds. Also contains the update button.
     */
    protected class FTestThresholdSelectionPanel extends JPanel {

        private double f1 = 0.0;

        private double f2 = 0.0;

        private double f3 = 0.0;

        private JScrollBar f1Controller;

        private JScrollBar f2Controller;

        private JScrollBar f3Controller;

        private JTextField f1Box;

        private JTextField f2Box;

        private JTextField f3Box;

        private int SCROLLRANGE = 50;

        private JButton updateButton;

        private JCheckBox autoUpdateCheck;

        private F_TestThresholdSelector container;

        protected FTestThresholdSelectionPanel(F_TestThresholdSelector container) {
            this.container = container;
            setLayout(new GridLayout(2, 4));

            f1 = container.f1;
            f2 = container.f2;
            f3 = container.f3;

            initComponents();
        }

        private void initComponents() {

            f1Controller = new JScrollBar(JScrollBar.HORIZONTAL, SCROLLRANGE / 2, 1, 0,
                    SCROLLRANGE);
            f2Controller = new JScrollBar(JScrollBar.HORIZONTAL, SCROLLRANGE / 2, 1, 0,
                    SCROLLRANGE);
            f3Controller = new JScrollBar(JScrollBar.HORIZONTAL, SCROLLRANGE / 2, 1, 0,
                    SCROLLRANGE);

            f1Controller.setUnitIncrement(1);
            f1Controller.setBlockIncrement(5);

            f2Controller.setUnitIncrement(1);
            f2Controller.setBlockIncrement(5);

            f3Controller.setUnitIncrement(1);
            f3Controller.setBlockIncrement(5);

            // initialise f thresholds to closest power of 10
            f1Controller.setValue((int) Math.round(Math.log(f1) / Math.log(10))
                    + SCROLLRANGE);
            f2Controller.setValue((int) Math.round(Math.log(f2) / Math.log(10))
                    + SCROLLRANGE);
            f3Controller.setValue((int) Math.round(Math.log(f3) / Math.log(10))
                    + SCROLLRANGE);

            f1 = Math.pow(10.0, f1Controller.getValue() - 50.0);
            f2 = Math.pow(10.0, f2Controller.getValue() - 50.0);
            f3 = Math.pow(10.0, f3Controller.getValue() - 50.0);

            f1Controller.addAdjustmentListener(new AdjustmentListener() {
                public void adjustmentValueChanged(AdjustmentEvent e) {
                    int newVal = e.getValue();
                    double newThresh = Math.pow(10.0, (double) (newVal - SCROLLRANGE));
                    f1 = newThresh;
                    syncBoxesToScrollers();
                    container.setFTestThresholds(f1, f2, f3);
                }
            });

            f2Controller.addAdjustmentListener(new AdjustmentListener() {
                public void adjustmentValueChanged(AdjustmentEvent e) {
                    int newVal = e.getValue();
                    double newThresh = Math.pow(10.0, (double) (newVal - SCROLLRANGE));
                    f2 = newThresh;
                    syncBoxesToScrollers();
                    container.setFTestThresholds(f1, f2, f3);

                }
            });

            f3Controller.addAdjustmentListener(new AdjustmentListener() {
                public void adjustmentValueChanged(AdjustmentEvent e) {
                    int newVal = e.getValue();
                    double newThresh = Math.pow(10.0, (double) (newVal - SCROLLRANGE));
                    f3 = newThresh;
                    syncBoxesToScrollers();
                    container.setFTestThresholds(f1, f2, f3);
                }
            });

            // top row of grid
            add(f1Controller);
            add(f2Controller);
            add(f3Controller);

            autoUpdateCheck = new JCheckBox("Auto-Update", container.getAutoUpdate());
            autoUpdateCheck.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.toggleAutoUpdate();
                }
            });

            add(autoUpdateCheck);

            // now add text boxes
            f1Box = new JTextField();

            // calling f*Controller.setValue() sets the scroll bar to the
            // nearest power of 10 to what's in the box, then resets the box to
            // the nearest power of 10.

            f1Box.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    try {
                        double val = Double.parseDouble(f1Box.getText());
                        f1Controller.setValue((int) Math.round(Math.log(val)
                                / Math.log(10))
                                + SCROLLRANGE); // confines values to the scroll
                        // range
                    }
                    catch (NumberFormatException ex) {
                        System.err.println("NumberFormatException when changing f1");
                    }

                }
            });

            // add event listener to update things when focus lost
            f1Box.addFocusListener(new FocusAdapter() {

                public void focusLost(FocusEvent e) {
                    double val = Double.parseDouble(f1Box.getText());
                    f1Controller.setValue((int) Math.round(Math.log(val) / Math.log(10))
                            + SCROLLRANGE);
                }

            });

            f2Box = new JTextField();

            f2Box.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    try {
                        double val = Double.parseDouble(f2Box.getText());
                        f2Controller.setValue((int) Math.round(Math.log(val)
                                / Math.log(10))
                                + SCROLLRANGE); // confines values to the scroll
                        // range
                    }
                    catch (NumberFormatException ex) {
                        System.err.println("NumberFormatException when changing f2");
                    }

                }
            });

            f2Box.addFocusListener(new FocusAdapter() {

                public void focusLost(FocusEvent e) {
                    double val = Double.parseDouble(f2Box.getText());
                    f2Controller.setValue((int) Math.round(Math.log(val) / Math.log(10))
                            + SCROLLRANGE);
                }

            });

            f3Box = new JTextField();

            f3Box.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    try {

                        double val = Double.parseDouble(f3Box.getText());
                        f3Controller.setValue((int) Math.round(Math.log(val)
                                / Math.log(10))
                                + SCROLLRANGE); // confines values to the scroll
                        // range
                    }
                    catch (NumberFormatException ex) {
                        System.err.println("NumberFormatException when changing f3");
                    }

                }
            });

            f3Box.addFocusListener(new FocusAdapter() {

                public void focusLost(FocusEvent e) {
                    double val = Double.parseDouble(f3Box.getText());
                    f3Controller.setValue((int) Math.round(Math.log(val) / Math.log(10))
                            + SCROLLRANGE);
                }

            });

            add(f1Box);
            add(f2Box);
            add(f3Box);

            updateButton = new JButton("Update");
            updateButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {

                    try {
                        f1 = Double.parseDouble(f1Box.getText());
                        f2 = Double.parseDouble(f2Box.getText());
                        f3 = Double.parseDouble(f3Box.getText());
                    }
                    catch (NumberFormatException ex) {
                        System.err.println(e);
                    }
                    container.setFTestThresholds(f1, f2, f3);
                    container.update();
                }
            });

            add(updateButton);

            syncBoxesToScrollers();
        }

        /**
         * Make sure that what's in the text box agrees with the scrollers.
         * Changes the text boxes.
         *  
         */
        private void syncBoxesToScrollers() {
            // should use decimalFormat, but for now, just dump text
            DecimalFormat df = new DecimalFormat("0.00E00");

            f1Box.setText(df.format(f1));
            f2Box.setText(df.format(f2));
            f3Box.setText(df.format(f3));
        }

    }

    /**
     * Contains buttons to switch slice, save the voxel classification, and
     * exit. Displays the fraction of non-background voxels that are classified
     * in each order.
     */
    protected class TopPanel extends JPanel {

        private F_TestThresholdSelector container;

        private JButton axialButton;

        private JButton coronalButton;

        private JButton sagittalButton;

        private JButton exitButton;

        private JButton saveButton;

        private JLabel fractionLabel;

        public TopPanel(F_TestThresholdSelector container) {
            this.container = container;
            setLayout(new GridLayout(2, 1));
            initComponents();
        }

        private void initComponents() {

            JPanel upperPanel = new JPanel();
            upperPanel.setLayout(new BoxLayout(upperPanel, BoxLayout.Y_AXIS));
            JPanel lowerPanel = new JPanel(new GridLayout(1, 5));

            axialButton = new JButton("AXIAL");

            axialButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.setSliceOrientation(F_TestThresholdSelector.XY);
                }
            });

            lowerPanel.add(axialButton);

            coronalButton = new JButton("CORONAL");

            coronalButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.setSliceOrientation(F_TestThresholdSelector.XZ);
                }
            });

            lowerPanel.add(coronalButton);

            sagittalButton = new JButton("SAGITTAL");

            sagittalButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.setSliceOrientation(F_TestThresholdSelector.YZ);
                }
            });

            lowerPanel.add(sagittalButton);

            saveButton = new JButton("Save");

            saveButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    save();
                }
            });

            lowerPanel.add(saveButton);

            exitButton = new JButton("Exit");

            exitButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    System.exit(0);
                }
            });

            lowerPanel.add(exitButton);

            fractionLabel = new JLabel();
            upperPanel.add(fractionLabel);

            add(upperPanel);
            add(lowerPanel);

        }

        public void save() {
            JFileChooser chooser = new JFileChooser();

            int r = chooser.showSaveDialog(TopPanel.this);
            if (r == JFileChooser.APPROVE_OPTION) {

                try {
                    container.saveVCImage(chooser.getSelectedFile().getCanonicalPath());
                }
                catch (IOException e) {
                    System.err.println(e);
                }
            }

        }

        /**
         * Called by container.update()
         */
        public void setFractions(double[] fractions) {

            DecimalFormat df = new DecimalFormat("00.00");

            String s = "<html>";

            Font f = fractionLabel.getFont();

            fractionLabel.setFont(new Font(f.getName(), Font.PLAIN, 14));

            double total = 0.0;

            // does order 0 2 4
            int max = 3;

            if (order / 2 + 1 < 3) {
                max = order / 2 + 1;
            }

            for (int i = 0; i < max; i++) {

                s = s + "<FONT COLOR=\"FF6633\"><b>Order " + i * 2 + ":</b> </FONT>"
                        + df.format(fractions[i] * 100.0) + " % ";

                total += fractions[i];

                s = s + " &nbsp &nbsp &nbsp &nbsp ";

            }

            s = s + "<FONT COLOR=\"FF6633\"><b>Others:</b> </FONT>"
                    + df.format((1.0 - total) * 100.0) + " % ";

            s = s + "</html>";

            fractionLabel.setText(s);
        }

    }

}

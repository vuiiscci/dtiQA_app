package apps;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;

import javax.swing.event.*;
import java.text.*;

 

/**
 * Top panel for PD_Orientation viewer. 
 * Contains buttons to switch slice and a zoomed out version of the brain.
 * The right panel (to the right of the zoom control) is customizable.
 * @author Philip Cook
 * @version $Id$
 * 
 */
class NavPanel extends JPanel {

    private PD_OrientationViewer container;

    private JButton axialButton;

    private JButton coronalButton;

    private JButton sagittalButton;

    private JButton saveButton;

    private JTextArea textPane;

    private JSlider zoomSlider;

    private DecimalFormat df = new DecimalFormat("0.0");

    private ZoomPanel zoomPanel;

    private JPanel rightPanel;

    public NavPanel(PD_OrientationViewer container) {
        this.container = container;
        setLayout(new GridLayout(1, 3));
        initComponents();
    }

    public int[] getZoomedRegion() {
        return zoomPanel.getZoomedRegion();
    }

    public void resetZoomPosition() {
        zoomPanel.resetZoomPosition();
    }

    
    public void setText(String text) {
	textPane.setText(text);
    }

    public void setSlice() {
        zoomPanel.setSlice();
    }

    private void initComponents() {


	// use BoxLayout. Add things in order left-right
	// Each set of tools (buttons, Zoom panel, etc) has its own layout 

	setLayout(new BoxLayout(this, BoxLayout.X_AXIS));

        JPanel leftPanel = new JPanel();
        leftPanel.setLayout(new GridLayout(3,1));

	JPanel buttonPanel = new JPanel();

	buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.X_AXIS));
       
        axialButton = new JButton("AXIAL");

        axialButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.setSliceOrientation(PD_OrientationViewer.XY);
                }
            });

        buttonPanel.add(axialButton);

        coronalButton = new JButton("CORONAL");

        coronalButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.setSliceOrientation(PD_OrientationViewer.XZ);
                }
            });

        buttonPanel.add(coronalButton);

        sagittalButton = new JButton("SAGITTAL");

        sagittalButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.setSliceOrientation(PD_OrientationViewer.YZ);
                }
            });

        buttonPanel.add(sagittalButton);

        saveButton = new JButton("SAVE RGB");

        saveButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.writeRGB();
                }
            });

        buttonPanel.add(saveButton);

        saveButton = new JButton("SAVE SLICES");

        saveButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.writeSlices();
                }
            });

        buttonPanel.add(saveButton);


	textPane = new JTextArea();

	leftPanel.add(buttonPanel);
	
	JPanel displayControlsPanel = new JPanel();

	displayControlsPanel.setLayout(new BoxLayout(displayControlsPanel, BoxLayout.X_AXIS));
	
	JCheckBox vectorsBox = new JCheckBox("Show vectors", true);
	
	vectorsBox.addChangeListener(new ChangeListener() {

		public void stateChanged(ChangeEvent e) {
		    container.showVectors = ! container.showVectors;
		    container.repaint();
		}
	    });

	displayControlsPanel.add(vectorsBox);

	int largestDimension = Math.max(container.xSize, container.ySize);

	largestDimension = Math.max(largestDimension, container.zSize);

	int maxZoomSliderValue = largestDimension;

	if (maxZoomSliderValue % 2 == 1) {
	    maxZoomSliderValue -= 1;
	}
	if (maxZoomSliderValue < 0) {
	    maxZoomSliderValue = 0;
	}
	int initialZoom = maxZoomSliderValue - container.preferredZoomedExtent;

	if (initialZoom < 0) {
	    initialZoom = 0;
	}

	zoomSlider = new JSlider(JSlider.HORIZONTAL, 0, maxZoomSliderValue, 
				     initialZoom);

	zoomSlider.setMajorTickSpacing(2);
	zoomSlider.setMinorTickSpacing(2);
	zoomSlider.setSnapToTicks(true);

	zoomSlider.addChangeListener(new ChangeListener() {

 		public void stateChanged(ChangeEvent e) {
		    
		    int extent = zoomSlider.getMaximum() - zoomSlider.getValue() - (zoomSlider.getValue() % 2);

		    container.preferredZoomedExtent = extent >= 2 ? extent : 2;
			
		    
		    NavPanel.this.zoomPanel.recentreBox();
		    container.setZoomedSlice();
		    container.repaint();
 		}
 	    });	

	displayControlsPanel.add(new JLabel("    Zoom"));

	displayControlsPanel.add(zoomSlider);

	leftPanel.add(displayControlsPanel);

	leftPanel.add(textPane);
	
        add(leftPanel);
	

        // next panel contains BufferedImage
        zoomPanel = new ZoomPanel(container);

        add(zoomPanel);
			
	JPanel rightPanel = createRightPanel();
	//	rightPanel.setMinimumSize(new Dimension(0,0));
	add(rightPanel);
        
    }


    /**
     * Override this method to put custom components in the top-right panel.
     *
     */
    protected JPanel createRightPanel() {

	final JSlider gsGammaSlider;
	final JSlider rgbGammaSlider;
	final JLabel gsGammaLabel;
	final JLabel rgbGammaLabel;


        JPanel rightPanel = new JPanel();

        rightPanel.setLayout(new GridLayout(3,1));


	rgbGammaLabel = new JLabel(df.format(container.image.rgbGamma()) + " ");	
	gsGammaLabel = new JLabel(df.format(container.image.scalarGamma()) + " ");


	// range 0 to 20 for gamma 0.0 to 2.0
	// if you change this change the listeners below
	gsGammaSlider = 
	    new JSlider(JSlider.HORIZONTAL, 0, 20, (int)(10.0 * container.image.scalarGamma()));

	gsGammaSlider.setMajorTickSpacing(1);
	gsGammaSlider.setMinorTickSpacing(1);
	gsGammaSlider.setSnapToTicks(true);
	gsGammaSlider.setPreferredSize(new Dimension(100,20));
	
	gsGammaSlider.addChangeListener(new ChangeListener() {
		
 		public void stateChanged(ChangeEvent e) {
		    
		    container.image.setScalarGamma(gsGammaSlider.getValue() / 10.0);
		    gsGammaLabel.setText(df.format(container.image.scalarGamma()));
		    container.updateDisplay();
		    container.repaint();
 		}
 	    });	



	rgbGammaSlider = 
	    new JSlider(JSlider.HORIZONTAL, 0, 20, (int)(10.0 * container.image.rgbGamma()));

	rgbGammaSlider.setMajorTickSpacing(1);
	rgbGammaSlider.setMinorTickSpacing(1);
	rgbGammaSlider.setSnapToTicks(true);
	rgbGammaSlider.setPreferredSize(new Dimension(100,20));

	rgbGammaSlider.addChangeListener(new ChangeListener() {

		public void stateChanged(ChangeEvent e) {
		    
		    container.image.setRGB_Gamma(rgbGammaSlider.getValue() / 10.0);

		    rgbGammaLabel.setText(df.format(container.image.rgbGamma()) + " ");

		    JSlider source = (JSlider)e.getSource();

		    // uncomment to adjust RGB only once mouse is released
		    //		    if (!source.getValueIsAdjusting()) {
			container.updateDisplay();
			container.repaint();
			//		    }
 		}
		
		
 	    });	
	

	// row 1
	JPanel rightPanelTop = new JPanel();
	rightPanelTop.setLayout(new BoxLayout(rightPanelTop, BoxLayout.X_AXIS));


	rightPanelTop.add(new JLabel(" Grey gamma"));
     	rightPanelTop.add(gsGammaSlider);
     	rightPanelTop.add(gsGammaLabel);

	// row 2
	JPanel rightPanelBottom = new JPanel();
	rightPanelBottom.setLayout(new BoxLayout(rightPanelBottom, BoxLayout.X_AXIS));

	rightPanelBottom.add(new JLabel("  RGB gamma"));
	rightPanelBottom.add(rgbGammaSlider);
	rightPanelBottom.add(rgbGammaLabel);

	rightPanel.add(rightPanelTop);
	rightPanel.add(rightPanelBottom);

	return rightPanel;

    }




}


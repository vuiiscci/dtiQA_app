package apps;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;

import javax.swing.event.*;
import java.text.*;

import imaging.DW_Scheme;
 

/**
 * Top panel for PD_Orientation viewer. 
 * Contains buttons to switch slice and a zoomed out version of the brain.
 * The right panel (to the right of the zoom control) is customizable.
 * @author Philip Cook
 * @version $Id$
 * 
 */
class SchemePanel extends JPanel {

    private PD_OrientationViewer container;

    private JRadioButton flipxButton;

    private JRadioButton flipyButton;

    private JRadioButton flipzButton;

    private JRadioButton noFlipButton;

    private String[] comboOptions = { "X Y Z", "X Z Y", "Y X Z", "Y Z X", "Z X Y", "Z Y X"};

    private JComboBox<String> swapGradsCombo;

    private JButton saveSchemeButton;

    public SchemePanel(PD_OrientationViewer container) {
        this.container = container;
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
	setBorder(BorderFactory.createTitledBorder("Scheme file options"));

        initComponents();
    }

    private void initComponents() {

	noFlipButton = new JRadioButton("no flip");
	noFlipButton.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
		    container.setFlip(PD_OrientationViewer.NO_FLIP);
                }
            });
	noFlipButton.setSelected(true);
	add(noFlipButton);

	flipxButton = new JRadioButton("flip x");
	flipxButton.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
		    container.setFlip(PD_OrientationViewer.FLIP_X);
                }
            });
	add(flipxButton);

	flipyButton = new JRadioButton("flip y");
	flipyButton.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
		    container.setFlip(PD_OrientationViewer.FLIP_Y);
		}
            });
	add(flipyButton);

	flipzButton = new JRadioButton("flip z");
	flipzButton.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
		    container.setFlip(PD_OrientationViewer.FLIP_Z);
                }
            });
	add(flipzButton);

	
	ButtonGroup bGroup = new ButtonGroup();
	
	bGroup.add(noFlipButton);
	bGroup.add(flipxButton);
	bGroup.add(flipyButton);
	bGroup.add(flipzButton);

	swapGradsCombo = new JComboBox<>(comboOptions);
	swapGradsCombo.setSelectedIndex(0);
	swapGradsCombo.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
		    int [] combo = new int [3];
		    switch(swapGradsCombo.getSelectedIndex()){
		    
		    case 0:
			combo=DW_Scheme.gradXYZ;
			break;
		    case 1:
			combo=DW_Scheme.gradXZY;
			break;
		    case 2:
			combo=DW_Scheme.gradYXZ;
			break;
		    case 3:
			combo=DW_Scheme.gradYZX;
			break;
		    case 4:
			combo=DW_Scheme.gradZXY;
			break;
		    case 5:
			combo=DW_Scheme.gradZYX;
			break;
		    }
                    container.swapDirs(combo);
                }
            });
	add(swapGradsCombo);

	saveSchemeButton = new JButton("SAVE SCHEME");
	//saveSchemeButton.setToolTipText("Saves scheme file.  Note: does not save pds");
	saveSchemeButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    container.writeScheme();
                }
            });
	add(saveSchemeButton);
	
        
    }
    protected void updateFlipButtons(int currentFlip) {

	// undo whatever flip we had originally
	flipxButton.setSelected(false);
	flipyButton.setSelected(false);
	flipzButton.setSelected(false);

	switch(currentFlip) {
	case PD_OrientationViewer.NO_FLIP:
	    noFlipButton.setSelected(true);
	    break;
	case PD_OrientationViewer.FLIP_X:
	    flipxButton.setSelected(true);
	    break;
	case PD_OrientationViewer.FLIP_Y:
	    flipyButton.setSelected(true);
	    break;
	case PD_OrientationViewer.FLIP_Z:
	    flipzButton.setSelected(true);
	    break;
	}
	
	container.repaint();
    }

    protected void resetButtons() {

	// resets buttons so that "no flip" radiobutton set and gradient set as "X Y Z"
	noFlipButton.setSelected(true);
	swapGradsCombo.setSelectedIndex(0);
	container.repaint();
    }

}


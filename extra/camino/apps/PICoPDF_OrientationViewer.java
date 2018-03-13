package apps;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.event.*;

import data.*;
import numerics.*;
import tools.*;


/**
 *
 * Simple tool to visualize principal directions of PICo data.
 * 
 * @author Philip Cook
 * @version $Id: PICoPDF_OrientationViewer.java 290 2007-05-23 20:13:26Z ucacpco $
 *  
 */
public class PICoPDF_OrientationViewer extends PD_OrientationViewer {

   
    /**
     * Creates and displays the viewer.
     */
    public PICoPDF_OrientationViewer(RGB_ScalarImage im, double scalarThresh) {
        
        super(im, scalarThresh);
	
	display();
    }


    // this seemingly pointless class exists so that we can easily add specific behavior to the viewer
    // by overriding the getTopPanel() method in the superclass

}


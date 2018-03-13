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
 * Simple tool to visualize principal directions of ball and stick data.
 * 
 * @author Philip Cook
 * @version $Id$
 *  
 */
public class BallSticksOrientationViewer extends PD_OrientationViewer {

   
    /**
     * Creates and displays the viewer.
     */
    public BallSticksOrientationViewer(RGB_ScalarImage im, double scalarThresh) {
        
        super(im, scalarThresh);
	
	display();
    }


    // this seemingly pointless class exists so that we can easily add specific behavior to the viewer
    // by overriding the getTopPanel() method in the superclass

}


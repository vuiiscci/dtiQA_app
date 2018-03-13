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
 * Simple tool to visualize principal directions of spherical function data.
 * 
 * @author Philip Cook
 * @version $Id$
 *  
 */
public class SphFuncPeaksOrientationViewer extends PD_OrientationViewer {

   
    /**
     * Creates and displays the viewer.
     */
    public SphFuncPeaksOrientationViewer(RGB_ScalarImage im, double scalarThresh) {
        
        super(im, scalarThresh);
	
	display();
    }


    // this seemingly pointless class exists so that we can easily add specific behavior to the viewer
    // by overriding the getTopPanel() method in the superclass

}

